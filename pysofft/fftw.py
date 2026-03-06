from pysofft import _soft
from pathlib import Path
import os
from collections import namedtuple
import multiprocessing

fftw_flags = {n:int(getattr(_soft.softclass,n)) for n in dir(_soft.softclass) if n.startswith("fftw_")}
tmp = namedtuple("FFTW_Flags",list(fftw_flags.keys()))

flags = tmp(*fftw_flags.values())
del(fftw_flags)
del(tmp)

_wisdom_file_path = Path('~/.config/pysofft/fftw_wisdom.dat').expanduser()

def set_wisdom_file_path(file_path):
    global _wisdom_file_path
    _wisdom_file_path = Path(file_path).expanduser()

def load_wisdom_from_file():
    if _wisdom_file_path.exists():
        _soft.softclass.load_fftw_wisdom(str(_wisdom_file_path))

def save_wisdom_to_file():
    # make sure files are only written to in the main process
    # to prevent write conflicts.
    if multiprocessing.parent_process() is None:
        if not _wisdom_file_path.parent.exists():
            os.makedirs(str(_wisdom_file_path.parent))
        _soft.softclass.save_fftw_wisdom(str(_wisdom_file_path))
    
