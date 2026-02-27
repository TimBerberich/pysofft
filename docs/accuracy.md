# Accurracy

Currently PySOFFT uses the [three-term recursion proposed by Kostelec](wigner.md#kostelec-recurrence) during transforms.
Its accuracy has been verified, until a bandwidth of `bw=600`, by checking whether the comuted small Wigner-matrices $d^l_{m,n}(\beta)$ are orthogonal.  

![Wigner-d orthogonality](images/accuracy_kostelec.png){width = 300}

/// warning | Issues for `bw>600`
For degrees/bandwidth higher that 600 and certain $\beta$ values this recurrence brakes down, e.g. for 
$\beta =\frac{\pi}{4}$ and $l>800$:  

![Wigner-d orthogonality bad](images/accuracy_kostelec_bad.png){width = 300}
///
