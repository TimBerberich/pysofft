#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <fftw3.h>

#include <limits.h>
#include <stddef.h>
#include <stdint.h>

#define CAPSULE_NAME "pysofft.fftw_buffer"

/* ---------- free memeory on destruction ---------- */

static void fftw_capsule_destructor(PyObject *capsule)
{
    void *ptr = PyCapsule_GetPointer(capsule, CAPSULE_NAME);
    if (ptr != NULL) {
        fftw_free(ptr);
    } else {
        /* Destructor must not leak an exception */
        PyErr_Clear();
    }
}

/* ---------- helpers ---------- */

static int
parse_shape(PyObject *obj, npy_intp **dims_out, int *nd_out)
{
    npy_intp *dims = NULL;
    int nd = 0;

    if (PyLong_Check(obj)) {
        Py_ssize_t n = PyLong_AsSsize_t(obj);
        if (n == -1 && PyErr_Occurred()) {
            return -1;
        }
        if (n < 0) {
            PyErr_SetString(PyExc_ValueError, "shape entries must be >= 0");
            return -1;
        }

        dims = PyMem_Malloc(sizeof(npy_intp));
        if (dims == NULL) {
            PyErr_NoMemory();
            return -1;
        }
        dims[0] = (npy_intp)n;
        nd = 1;
    } else {
        PyObject *seq = PySequence_Fast(obj, "shape must be an int or a sequence of ints");
        if (seq == NULL) {
            return -1;
        }

        nd = (int)PySequence_Fast_GET_SIZE(seq);

        if (nd > NPY_MAXDIMS) {
            Py_DECREF(seq);
            PyErr_Format(PyExc_ValueError, "too many dimensions (max %d)", NPY_MAXDIMS);
            return -1;
        }

        dims = PyMem_Malloc((size_t)nd * sizeof(npy_intp));
        if (dims == NULL) {
            Py_DECREF(seq);
            PyErr_NoMemory();
            return -1;
        }

        for (int i = 0; i < nd; ++i) {
            PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
            Py_ssize_t v = PyLong_AsSsize_t(item);
            if (v == -1 && PyErr_Occurred()) {
                Py_DECREF(seq);
                PyMem_Free(dims);
                return -1;
            }
            if (v < 0) {
                Py_DECREF(seq);
                PyMem_Free(dims);
                PyErr_SetString(PyExc_ValueError, "shape entries must be >= 0");
                return -1;
            }
            dims[i] = (npy_intp)v;
        }

        Py_DECREF(seq);
    }

    *dims_out = dims;
    *nd_out = nd;
    return 0;
}

static int
compute_nbytes(const npy_intp *dims, int nd, size_t itemsize, size_t *nbytes_out)
{
    size_t nitems = 1;

    if (nd == 0) {
        *nbytes_out = itemsize;
        return 0;
    }

    for (int i = 0; i < nd; ++i) {
        if (dims[i] == 0) {
            *nbytes_out = 1;  /* allocate a harmless non-NULL pointer */
            return 0;
        }
        if ((size_t)dims[i] > SIZE_MAX / nitems) {
            return -1;
        }
        nitems *= (size_t)dims[i];
    }

    if (itemsize > SIZE_MAX / nitems) {
        return -1;
    }

    *nbytes_out = nitems * itemsize;
    return 0;
}

static PyObject *
make_aligned_array(PyObject *shape_obj, int typenum, size_t itemsize)
{
    npy_intp *dims = NULL;
    int nd = 0;
    size_t nbytes = 0;
    void *data = NULL;
    PyObject *arr = NULL;
    PyObject *capsule = NULL;

    if (parse_shape(shape_obj, &dims, &nd) < 0) {
        return NULL;
    }

    if (compute_nbytes(dims, nd, itemsize, &nbytes) < 0) {
        PyMem_Free(dims);
        PyErr_SetString(PyExc_OverflowError, "array is too large");
        return NULL;
    }

    data = fftw_malloc(nbytes);
    if (data == NULL) {
        PyMem_Free(dims);
        PyErr_NoMemory();
        return NULL;
    }

    /* C-contiguous ndarray backed by external FFTW memory */
    arr = PyArray_SimpleNewFromData(nd, dims, typenum, data);
    PyMem_Free(dims);

    if (arr == NULL) {
        fftw_free(data);
        return NULL;
    }

    capsule = PyCapsule_New(data, CAPSULE_NAME, fftw_capsule_destructor);
    if (capsule == NULL) {
        Py_DECREF(arr);
        fftw_free(data);
        return NULL;
    }

    if (PyArray_SetBaseObject((PyArrayObject *)arr, capsule) < 0) {
        /* steals reference only on success */
        Py_DECREF(capsule);
        Py_DECREF(arr);
        fftw_free(data);
        return NULL;
    }

    return arr;
}

/* ---------- Python-callable functions ---------- */

static PyObject *
py_create_float64(PyObject *self, PyObject *arg)
{
    (void)self;
    return make_aligned_array(arg, NPY_FLOAT64, sizeof(double));
}

static PyObject *
py_create_complex128(PyObject *self, PyObject *arg)
{
    (void)self;
    return make_aligned_array(arg, NPY_COMPLEX128, sizeof(fftw_complex));
}

/*  Returns the memory adress/pointer to check for alignment from python
	e.g. via  a = create_float64((32,32))
	          ptr = get_address(a)
			  aligned_16 = ptr%16==0
			  aligned_32 = ptr%32==0
			  aligned_64 = ptr%64==0
 */
static PyObject *
py_get_address(PyObject *self, PyObject *arg)
{
    (void)self;
    PyObject *arr_obj = NULL;
    PyArrayObject *arr = NULL;
    uintptr_t addr;

    if (!PyArg_ParseTuple(arg, "O", &arr_obj)) {
        return NULL;
    }

    arr = (PyArrayObject *)PyArray_FROM_O(arr_obj);
    if (arr == NULL) {
        return NULL;
    }

    addr = (uintptr_t)PyArray_DATA(arr);
    Py_DECREF(arr);

    return PyLong_FromUnsignedLongLong((unsigned long long)addr);
}

static PyMethodDef module_methods[] = {
    {
        "create_float64",
        (PyCFunction)py_create_float64,
        METH_O,
        PyDoc_STR("empty_float64(shape) -> aligned C-contiguous float64 ndarray")
    },
    {
        "create_complex128",
        (PyCFunction)py_create_complex128,
        METH_O,
        PyDoc_STR("empty_complex128(shape) -> aligned C-contiguous complex128 ndarray")
    },
    {
        "get_address",
        (PyCFunction)py_get_address,
        METH_VARARGS,
        PyDoc_STR("address(arr) -> integer data pointer address")
    },
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_fftw_aligned_alloc",
    "FFTW-aligned NumPy array allocation helpers.",
    -1,
    module_methods
};

PyMODINIT_FUNC
PyInit__fftw_aligned_alloc(void)
{
    PyObject *m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }

    if (PyArray_ImportNumPyAPI() < 0) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
