#!/usr/bin/env python3
"""Monte Carlo distance of the K X-ray to the main track in gamma events.

  python scripts/production/xray_truth_distance.py EVENT_LIST OUTPUT.txt NEXUS_FILE...

EVENT_LIST holds the numbers of the accepted events (one per line). For each of
them found in the NEXUS files, when the leading electron of the generated gamma
is a photoelectron and the absorption emitted a K X-ray (a gamma of 29.2-34.2
keV produced by the generated gamma), the script writes one row: event number,
distance from the energy-weighted centre of the X-ray deposits to the nearest
hit of the main track (mm), distance to the creation point of the main
electron (mm), energy deposited by the X-ray and its descendants (keV), X-ray
energy (keV), Z of the deposit (cm) and energy of the main track (keV). The
main track is the leading electron and its descendants not passing through a
gamma. Events whose leading electron is not a photoelectron, or with no K
X-ray, give a row with NaN distances and the creator process code
(1 photoelectric without K X-ray, 2 Compton, 3 pair production).
"""
import sys
import h5py,numpy as np

def rows(path,accepted):
    f=h5py.File(path,'r');P=f['MC/particles'][:];H=f['MC/hits'][:];name=np.char.decode(P['particle_name']);proc=np.char.decode(P['creator_proc']);out=[]
    order=np.argsort(P['event_id'],kind='stable');P=P[order];name=name[order];proc=proc[order];ho=np.argsort(H['event_id'],kind='stable');H=H[ho]
    ids,start=np.unique(P['event_id'],return_index=True);end=np.r_[start[1:],len(P)];hs=np.searchsorted(H['event_id'],ids);he=np.searchsorted(H['event_id'],ids,side='right')
    for ev,a,b,c,d in zip(ids,start,end,hs,he):
        if int(ev) not in accepted: continue
        p=P[a:b];nm=name[a:b];pr=proc[a:b];h=H[c:d];primary=p['particle_id'][p['primary']==1][0]
        el=np.flatnonzero((p['mother_id']==primary)&(nm=='e-'));lead=el[np.argmax(p['kin_energy'][el])]
        xr=np.flatnonzero((p['mother_id']==primary)&(nm=='gamma')&(p['kin_energy']>0.0292)&(p['kin_energy']<0.0342))
        if pr[lead]!='phot' or len(xr)==0:
            out.append((ev,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,{'phot':1,'compt':2,'conv':3}.get(pr[lead],9)));continue
        mother=dict(zip(p['particle_id'].tolist(),p['mother_id'].tolist()));gamma=dict(zip(p['particle_id'].tolist(),(nm=='gamma').tolist()))
        lid=int(p['particle_id'][lead]);xid=int(p['particle_id'][xr[0]]);main=[];xh=[]
        for j,pid in enumerate(h['particle_id'].tolist()):
            chain=[];i=pid
            while i!=0: chain.append(i);i=mother[i]
            if xid in chain: xh.append(j)
            elif lid in chain and not any(gamma[k] for k in chain[:chain.index(lid)]): main.append(j)
        if not xh: out.append((ev,np.nan,np.nan,0.0,1000*p['kin_energy'][xr[0]],np.nan,np.nan,0));continue
        xyz=np.c_[h['x'],h['y'],h['z']].astype(float);e=h['energy'].astype(float)*1000;centre=np.average(xyz[xh],axis=0,weights=e[xh])
        origin=np.array([p['initial_x'][lead],p['initial_y'][lead],p['initial_z'][lead]],float)
        out.append((ev,np.linalg.norm(xyz[main]-centre,axis=1).min(),np.linalg.norm(centre-origin),e[xh].sum(),1000*p['kin_energy'][xr[0]],centre[2]/10,e[main].sum(),0))
    return out

if __name__=='__main__':
    if len(sys.argv)<4: raise SystemExit(__doc__)
    accepted=set(np.loadtxt(sys.argv[1],dtype=np.int64).ravel().tolist());table=[r for path in sys.argv[3:] for r in rows(path,accepted)]
    np.savetxt(sys.argv[2],np.array(table),fmt='%.4f',header='event min_distance_to_main_track_mm distance_to_start_mm xray_deposit_keV xray_energy_keV xray_z_cm main_track_keV code(0 photoelectric with K X-ray, 1 photoelectric without, 2 Compton, 3 pair)')
    print(len(table),'events written to',sys.argv[2])
