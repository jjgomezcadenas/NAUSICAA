"""LTspice schematic connectivity and binary data for the supplied CSA.

Limited intentionally to the six standard symbols in Rogelio's schematic.
Pin positions/order transcribed from the official LTspice 26 symbol library.
Unknown symbols, disconnected pins and ambiguous raw formats fail explicitly.
"""
from pathlib import Path
import re
import numpy as np

PINS = {
    'pmos4': [(48,0),(0,80),(48,96),(48,48)],
    'nmos4': [(48,0),(0,80),(48,96),(48,48)],
    'voltage': [(0,16),(0,96)], 'current': [(0,0),(0,80)],
    'cap': [(16,0),(16,64)], 'bv': [(0,16),(0,96)]}


def transform_pin(pin, origin, orientation):
    x,y=pin
    angle=int(orientation[1:])
    for _ in range(angle//90): x,y=-y,x
    if orientation.startswith('M'): x=-x
    return (origin[0]+x, origin[1]+y)


def schematic_netlist(path):
    """Convert wires/flags and known symbols; keep all active directives."""
    wires=[]; flags=[]; symbols=[]; directives=[]
    for line in Path(path).read_text(encoding='latin-1').splitlines():
        p=line.split()
        if not p: continue
        if p[0]=='WIRE': wires.append((tuple(map(int,p[1:3])),tuple(map(int,p[3:5]))))
        elif p[0]=='FLAG': flags.append((tuple(map(int,p[1:3])),p[3]))
        elif p[0]=='SYMBOL':
            name=p[1]
            if name not in PINS: raise ValueError(f'Unsupported symbol: {name}')
            symbols.append(dict(kind=name,pins=[transform_pin(pin,tuple(map(int,p[2:4])),p[4]) for pin in PINS[name]]))
        elif p[0]=='SYMATTR': symbols[-1][p[1]]=' '.join(p[2:])
        elif p[0]=='TEXT' and ' !' in line: directives.extend(line.split(' !',1)[1].split('\\n'))
    points={p for wire in wires for p in wire}|{p for p,n in flags}|{p for s in symbols for p in s['pins']}
    parent={p:p for p in points}
    def find(p):
        if parent[p]!=p: parent[p]=find(parent[p])
        return parent[p]
    def join(a,b): parent[find(a)]=find(b)
    for a,b in wires:
        for p in points:
            if ((b[0]-a[0])*(p[1]-a[1])==(b[1]-a[1])*(p[0]-a[0])
                and min(a[0],b[0])<=p[0]<=max(a[0],b[0])
                and min(a[1],b[1])<=p[1]<=max(a[1],b[1])): join(a,p)
    named={}
    for p,n in flags:
        if n in named: join(p,named[n])
        named[n]=p
    names={}
    for p,n in flags:
        root=find(p)
        if root not in names or n=='0': names[root]=n
    for p in sorted(points): names.setdefault(find(p),f'N{len(names):03d}')
    elements=[]
    for s in symbols:
        for pin in s['pins']:
            if not any(pin==p for w in wires for p in w) and pin not in dict(flags):
                raise ValueError(f'Pin lacks wire endpoint/flag: {s["InstName"]} {pin}')
        elements.append(' '.join([s['InstName'],*[names[find(p)] for p in s['pins']],s['Value'],s.get('Value2','')]).strip())
    return elements,directives


def read_raw(path):
    """Read ordinary (not stepped/FastAccess) native LTspice binary raw data."""
    b=Path(path).read_bytes()
    encoding='utf-16le' if b[1:2]==b'\x00' else 'utf-8'
    marker='Binary:\n'.encode(encoding); offset=b.find(marker)
    if offset<0: raise ValueError('Expected binary LTspice raw file')
    header=b[:offset].decode(encoding)
    flags=re.search(r'^Flags: (.*)$',header,re.M)[1].split()
    if any(x in flags for x in ('stepped','FastAccess')): raise ValueError(f'Unsupported raw flags {flags}')
    n=int(re.search(r'No. Variables:\s+(\d+)',header)[1]); m=int(re.search(r'No. Points:\s+(\d+)',header)[1])
    names=[l.split()[1] for l in header.split('Variables:\n',1)[1].splitlines() if l.strip()]
    payload=b[offset+len(marker):]
    if 'complex' in flags:
        data=np.frombuffer(payload,dtype='<c16').reshape(m,n).copy()
    elif len(payload)==m*n*8:
        data=np.frombuffer(payload,dtype='<f8').reshape(m,n).copy()
    elif len(payload)==m*(8+4*(n-1)):
        rows=np.frombuffer(payload,dtype=[('x','<f8'),('ys','<f4',(n-1,))])
        data=np.column_stack([rows['x'],rows['ys']])
    else: raise ValueError(f'Unexpected raw length {len(payload)} for {m}x{n}')
    if names[0]=='time':
        # LTspice compression encodes a marker in the time sign bit.
        data[:,0]=np.abs(data[:,0])
        if np.any(np.diff(data[:,0])<0): raise ValueError('Nonmonotone transient time')
    if len(names)!=n or not np.all(np.isfinite(data)): raise ValueError('Invalid raw data')
    return dict(zip(names,data.T)),header


def sampled_charge(time_s, output_V, sample_times_s, cf_F, tau_s):
    """Paper's point-sample estimator; subtract the DC output before calling."""
    times=np.asarray(sample_times_s)
    if times[0]<time_s[0] or times[-1]>time_s[-1]: raise ValueError('Samples outside transient')
    dt=np.diff(times); alpha=np.exp(-dt/tau_s)
    v=np.interp(times,time_s,output_V)
    return -cf_F*(dt/tau_s)/(-np.expm1(-dt/tau_s))*(v[1:]-alpha*v[:-1])


def charge_noise(frequency_Hz, output_asd_V_sqrtHz, dt_s, cf_F, tau_s, order=8):
    """Integrate one-sided PSD with point-sample difference filter.

    Log-linear PSD interpolation; Gauss panels include every spectral knot and
    every half period of the oscillatory covariance filter. Increasing order
    tests integration separately from increasing the SPICE frequency density.
    """
    f=np.asarray(frequency_Hz); asd=np.asarray(output_asd_V_sqrtHz)
    if np.any(np.diff(f)<=0) or np.any(f<=0) or np.any(asd<0): raise ValueError('Invalid noise spectrum')
    if not np.any(asd): return 0.
    alpha=np.exp(-dt_s/tau_s); k=cf_F*(dt_s/tau_s)/(-np.expm1(-dt_s/tau_s))
    boundaries=np.unique(np.r_[f,np.arange(f[0],f[-1],.5/dt_s)])
    x,w=np.polynomial.legendre.leggauss(order)
    half=np.diff(boundaries)/2; mid=(boundaries[:-1]+boundaries[1:])/2
    nodes=mid[:,None]+half[:,None]*x
    spectrum=np.exp(np.interp(np.log(nodes),np.log(f),2*np.log(np.maximum(asd,1e-150))))
    weight=(1-alpha)**2+4*alpha*np.sin(np.pi*nodes*dt_s)**2
    variance=np.sum(half[:,None]*w*spectrum*weight)
    return float(k*np.sqrt(variance))
