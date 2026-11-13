"""
    emin_moments(; energy_keV, w_eV, fano, survival, yield_per_e, J,
                   acceptance, pde, cells, dark_per_cell, noise_pe)

Return mean counts and separate variance contributions for a deposit summed
across `cells` SiPM–time buckets. The number of primary electron–ion pairs has
variance F*N. Attachment is binomial. Conditional EL variance is J*Y per
surviving electron, as in production_readout.jl. Optical collection and PDE
are successive binomial operations. Independent ASIC noise has variance
cells*noise_pe²; stationary primary dark counts have variance equal to their
mean cells*dark_per_cell. Subtract the known mean dark signal.

PE denotes fully recovered avalanche equivalents. This linear calculation
assumes no clipping, crosstalk or afterpulsing. Returns S/N for ASIC noise alone,
ASIC plus dark noise, and all signal/noise fluctuations. These are distinct
quantities; none is a false-trigger probability. Inputs use keV/eV and counts.
"""
function emin_moments(;energy_keV,w_eV,fano,survival,yield_per_e,J,acceptance,pde,cells,dark_per_cell,noise_pe)
    n=1000energy_keV/w_eV;s=survival;Y=yield_per_e;a=acceptance;p=pde
    pe=n*s*Y*a*p
    vf=n*fano*(s*Y*a*p)^2
    va=n*s*(1-s)*(Y*a*p)^2
    vel=n*s*J*Y*(a*p)^2
    vo=n*s*Y*a*(1-a)*p^2
    vp=n*s*Y*a*p*(1-p)
    vd=cells*dark_per_cell;vn=cells*noise_pe^2
    total=vf+va+vel+vo+vp+vd+vn
    (;primary_electrons=n,arriving_electrons=n*s,emitted_vuv=n*s*Y,photons_at_sipm=n*s*Y*a,pe,
      fano_variance=vf,attachment_variance=va,el_variance=vel,optical_variance=vo,pde_variance=vp,
      dark_mean_pe=vd,dark_variance=vd,asic_variance=vn,noise_rms_pe=sqrt(vd+vn),total_rms_pe=sqrt(total),
      sn_asic=pe/sqrt(vn),sn_noise=pe/sqrt(vd+vn),sn_total=pe/sqrt(total))
end
