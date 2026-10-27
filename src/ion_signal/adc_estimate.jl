module ADCEstimate
using HDF5,TOML,SHA
sha(p)=bytes2hex(sha256(read(p)))
"""Independent uniform quantization of two single readings, with no averaging credit."""
quantized_difference_rms(span,bits,alpha)=span/2.0^bits*sqrt((1+alpha^2)/12)
function analyze(root)
 cp=joinpath(root,"config/adc_estimate.toml");c=TOML.parsefile(cp)
 dir=joinpath(root,c["output"]);mkpath(dir);source=joinpath(root,c["source"])
 h5open(source,"r") do src
  a=attributes(src[c["source_group"]]);gain=read(a["gain_V_per_ion"])
  noise=read(a["noise_rms_V"]);tau=read(a["tau_s"]);alpha=exp(-c["frame_s"]/tau)
  h5open(joinpath(dir,"analysis.h5"),"w") do out
   at=attributes(out)
   for (k,v) in ("source_sha256"=>sha(source),"config_sha256"=>sha(cp),"analysis_sha256"=>sha(@__FILE__),
        "gain_V_per_ion"=>gain,"CSA_noise_rms_V"=>noise,"CSA_ENC_e_rms"=>noise/gain,
        "baseline_V"=>read(a["baseline_V"]),"tau_s"=>tau,"alpha"=>alpha,"frame_s"=>c["frame_s"],
        "adc_noise_fraction_budget"=>c["adc_noise_fraction_budget"],"total_enc_increase_budget"=>c["total_enc_increase_budget"])
    at[k]=v
   end
   at["assumptions"]="Uniform mean ion arrivals; per-circuit 300 us averaging and calibrated frame difference from the rounded paper circuit. Populations scale linearly from10ms to20ms at fixed arrival rate. Independent uniform single-reading quantization errors, no averaging reduction. Source gain at100ions;4000ion extrapolation initially unvalidated."
   p=create_group(out,"populations")
   counts=c["ions_per_reference_frame"].*(c["frame_s"]/c["population_reference_frame_s"])
   p["labels"]=c["population_labels"];p["ions_per_frame"]=counts
   p["estimated_frame_signal_V"]=counts.*gain;p["analog_snr"]=counts.*gain./noise
   for span in c["full_scale_spans_V"]
    g=create_group(out,"span_$(Int(span))V");a=attributes(g)
    a["span_V"]=span
    bits=c["bits"];q=span./2.0.^bits
    sigma=[quantized_difference_rms(span,b,alpha) for b in bits]
    total=hypot.(noise,sigma)
    g["bits"]=bits;g["LSB_V"]=q;g["ADC_noise_rms_V"]=sigma;g["ADC_ENC_e_rms"]=sigma./gain
    g["total_ENC_e_rms"]=total./gain;g["total_ENC_relative_increase"]=total./noise.-1
    g["ADC_to_CSA_noise_ratio"]=sigma./noise
    g["minimum_signal_LSBs"]=first(counts)*gain./q
    g["passes_ADC_noise_budget"]=UInt8.(sigma.<=c["adc_noise_fraction_budget"]*noise)
    g["passes_total_ENC_budget"]=UInt8.(total.<=noise*(1+c["total_enc_increase_budget"]))
    required_q=c["adc_noise_fraction_budget"]*noise*sqrt(12/(1+alpha^2))
    a["effective_bits_for_ADC_noise_budget"]=log2(span/required_q)
    a["integer_bits_for_ADC_noise_budget"]=ceil(Int,read(a["effective_bits_for_ADC_noise_budget"]))
   end
  end
 end
end
end
