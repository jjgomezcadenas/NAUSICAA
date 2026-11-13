using LinearAlgebra
"""
    ion_emin_moments(energy_keV, w_eV, fano, covariance, K; arrangement)

Compute mean collected ions, Fano variance, summed linear-CSA noise and S/N.
The covariance is in reconstructed-ion equivalents squared. `independent_pads`
sums one central bucket from each of K independent pads. `consecutive_buckets`
sums K central consecutive buckets on one pad, including all covariance terms.
The full deposited charge is assumed collected and summed. There is no ion
attachment or gain. Return additive-noise S/N and signal/total RMS separately.
The central submatrix keeps the sum away from acquisition boundaries.
"""
function ion_emin_moments(energy_keV,w_eV,fano,covariance,K;arrangement)
    n=size(covariance,1)
    1<=K<=n && size(covariance,2)==n || throw(ArgumentError("Invalid bucket count or covariance shape"))
    energy_keV>=0 && w_eV>0 && fano>=0 || throw(ArgumentError("Invalid ionisation parameters"))
    mid=cld(n,2);first=div(n-K,2)+1
    variance=if arrangement=="independent_pads"
        K*covariance[mid,mid]
    elseif arrangement=="consecutive_buckets"
        sum(@view covariance[first:first+K-1,first:first+K-1])
    else
        throw(ArgumentError("Unknown arrangement"))
    end
    variance>0 || throw(ArgumentError("Nonpositive summed variance"))
    ions=1000energy_keV/w_eV
    (;ions,fano_variance=fano*ions,noise_variance=variance,noise_rms_ions=sqrt(variance),
      total_rms_ions=sqrt(variance+fano*ions),sn_noise=ions/sqrt(variance),sn_total=ions/sqrt(variance+fano*ions))
end
