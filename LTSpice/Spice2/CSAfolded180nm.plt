[Noise Spectral Density - (V/Hz½ or A/Hz½)]
{
   Npanes: 6
   Active Pane: 2
   {
      traces: 1 {1048579,0,"i(inoise)"}
      X: ('G',0,0.1,0,1e+09)
      Y[0]: ('p',0,0,8e-12,8.8e-11)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Units: "A/Hz½" ('p',0,0,0,0,8e-12,8.8e-11)
      Log: 1 0 0
      GridStyle: 1
   },
   {
      traces: 1 {8388612,0,"gain"}
      X: ('G',0,0.1,0,1e+09)
      Y[0]: ('T',0,0,1e+13,1.5e+14)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Units: "ohm" ('T',0,0,0,0,1e+13,1.5e+14)
      Log: 1 0 0
      GridStyle: 1
   },
   {
      traces: 2 {524295,0,"v(M1.id)"} {268959752,0,"v(M1.rd)"}
      X: ('G',0,0.1,0,1e+09)
      Y[0]: ('n',0,0,6e-08,6.6e-07)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Units: "V/Hz½" ('n',0,0,0,0,6e-08,6.6e-07)
      Log: 1 0 0
      GridStyle: 1
      Text: "V/Hz½" 1 (12195.7046015944,5.24383561643836e-07) ;v(M1,id) channel thermal noise, v(M1,rd) parasitic drain resistance thermal noise
      Text: "V/Hz½" 1 (8531678.52417281,3.07397260273973e-07) ;Transistor M1
   },
   {
      traces: 1 {524293,0,"v(M1.1overf)"}
      X: ('G',0,0.1,0,1e+09)
      Y[0]: ('µ',0,0,5e-05,0.0005)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Units: "V/Hz½" ('µ',0,0,0,0,5e-05,0.0005)
      Log: 1 0 0
      GridStyle: 1
      Text: "V/Hz½" 1 (11721.0229753348,0.000397260273972603) ;OUTPUT NOISE CONTRIBUTION FROM TRANSISTOR M1´s FLICKER (1/f) NOISE
   },
   {
      traces: 1 {524294,0,"v(M1)"}
      X: ('G',0,0.1,0,1e+09)
      Y[0]: ('µ',0,0,5e-05,0.0005)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Units: "V/Hz½" ('µ',0,0,0,0,5e-05,0.0005)
      Log: 1 0 0
      GridStyle: 1
   },
   {
      traces: 1 {524290,0,"v(onoise)"}
      X: ('G',0,0.1,0,1e+09)
      Y[0]: ('m',0,0,0.002,0.02)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Units: "V/Hz½" ('m',0,0,0,0,0.002,0.02)
      Log: 1 0 0
      GridStyle: 1
   }
}
[Transient Analysis]
{
   Npanes: 2
   {
      traces: 2 {524292,0,"V(in)"} {34603011,1,"I(Is1)"}
      X: (' ',0,0,1,10)
      Y[0]: ('m',1,0.0161,0.0007,0.0252)
      Y[1]: ('f',0,0,1e-14,1.3e-13)
      Volts: ('m',0,0,0,0.0161,0.0007,0.0252)
      Amps: ('f',0,0,0,0,1e-14,1.3e-13)
      Log: 0 0 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   },
   {
      traces: 1 {524290,0,"V(out)"}
      X: (' ',0,0,1,10)
      Y[0]: (' ',2,0.09,0.09,1.08)
      Y[1]: ('_',0,1e+308,0,-1e+308)
      Volts: (' ',0,0,0,0.09,0.09,1.08)
      Log: 0 0 0
      GridStyle: 1
      PltMag: 1
      PltPhi: 1 0
   }
}
