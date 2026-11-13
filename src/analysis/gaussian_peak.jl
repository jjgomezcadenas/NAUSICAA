"""Gaussian peak fits with a fixed fit interval and explicit truncation normalization."""
module GaussianPeak
using Statistics, LinearAlgebra
using Distributions: Normal, cdf, pdf
export fit_peak

"""
    fit_peak(values; window_sigma=2.5)

Fit an unbinned Gaussian to the peak within median(values) ± window_sigma ×
1.4826*MAD. Keep that interval fixed during optimization. The likelihood is
normalized to the probability inside the interval, avoiding the width bias of
an ordinary standard deviation computed on a clipped sample.

Fit mean and log(sigma) with Newton steps and a decreasing-likelihood line search.
Return fit parameters, observed-likelihood errors, fixed bounds, included count,
and the robust centre/scale used to define the interval. Fit-window variations
must call this function on the same full sample. No event is removed from the
saved input. Values must be finite and have a nonzero peak width.
"""
function fit_peak(values;window_sigma=2.5)
    x=Float64.(values);all(isfinite,x)||throw(ArgumentError("Nonfinite data"))
    centre=median(x);scale=1.482602218505602*median(abs.(x.-centre))
    scale>0&&window_sigma>0||throw(ArgumentError("Zero width or invalid window"))
    lo=centre-window_sigma*scale;hi=centre+window_sigma*scale
    z=(filter(v->lo<=v<=hi,x).-centre)./scale
    n=length(z);n>=20||throw(ArgumentError("Too few peak entries"))
    m=mean(z);m2=mean(abs2,z);normal=Normal()
    function objective(q)
        mu,logs=q;sigma=exp(logs)
        isfinite(sigma)&&1e-4<sigma<1e4||return Inf
        prob=cdf(normal,(window_sigma-mu)/sigma)-cdf(normal,(-window_sigma-mu)/sigma)
        prob>1e-12||return Inf
        logs+(m2-2mu*m+mu^2)/(2sigma^2)+log(prob)
    end
    function derivatives(q)
        h=1e-4;f=objective(q);g=zeros(2);H=zeros(2,2)
        for i in 1:2
            d=zeros(2);d[i]=h
            g[i]=(objective(q+d)-objective(q-d))/(2h)
            H[i,i]=(objective(q+d)-2*f+objective(q-d))/h^2
        end
        a=[h,h];b=[h,-h]
        H[1,2]=H[2,1]=(objective(q+a)+objective(q-a)-objective(q+b)-objective(q-b))/(4h^2)
        g,H
    end
    q=[m,log(std(z;corrected=false))];converged=false
    for _ in 1:100
        grad,H=derivatives(q)
        if norm(grad)<1e-7;converged=true;break;end
        λ=minimum(eigvals(Symmetric(H)));λ<=1e-6&&(H+=I*(1e-6-λ))
        step=H\grad;factor=1.;old=objective(q)
        while factor>1e-8&&objective(q-factor*step)>=old;factor/=2;end
        factor>1e-8||break
        q-=factor*step
    end
    grad,H=derivatives(q);converged|=norm(grad)<1e-6
    converged||error("Gaussian peak fit did not converge: $grad")
    covariance=inv(H*n);all(diag(covariance).>0)||error("Invalid fit covariance")
    mu=centre+q[1]*scale;sigma=exp(q[2])*scale
    (;mu,sigma,mu_error=scale*sqrt(covariance[1,1]),sigma_error=sigma*sqrt(covariance[2,2]),
      lower=lo,upper=hi,entries=n,total_entries=length(x),centre,scale,window_sigma,
      relative_sigma_percent=100sigma/mu)
end
end
