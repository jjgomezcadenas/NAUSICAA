"""Build browser plots from Julia's completed MC-energy snapshot.

Only presentation geometry, ranges and colour limits are calculated here.
Energy sums, projections and contributing-hit counts are read from HDF5.
"""
import h5py
import numpy as np
from plotly.colors import get_colorscale


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def snapshot_payload(path, display):
    """Read one complete viewer-v1 snapshot and return JSON-ready Plotly figures.

Projection arrays arrive from h5py as (vertical,horizontal). Occupied voxel
indices remain one-based in storage and are converted to edge offsets here.
Each voxel has eight independent vertices and twelve triangles, keeping its
colour constant across all faces. Every occupied cell is rendered. Each axis
uses its own padded data bounds; the 3D scene maps them to a cube. Coordinates
and hover bounds remain in mm while the displayed axis scales are independent.
"""
    with h5py.File(path) as f:
        if text(f.attrs['schema']) != 'nausicaa.event.viewer.v1' or text(f.attrs['status']) != 'complete':
            raise ValueError('The event snapshot is incomplete or has an unsupported schema')
        edges = {a: f[f'display/axes/{a}_edges_mm'][:] for a in 'xyz'}
        ranges = {}
        for a in 'xyz':
            hits = f[f'event/hits/{a}_mm'][:]
            lo, hi = float(hits.min()), float(hits.max())
            mid = (lo + hi) / 2
            half = max((hi-lo)*float(display['extent_scale'])/2, .325)
            ranges[a] = [min(mid-half, float(edges[a][0])), max(mid+half, float(edges[a][-1]))]
        projections = {a: f[f'display/projections/{a}_energy_MeV'][:] for a in ('xy','xz','yz')}
        max_energy = max(float(v.max()) for v in projections.values()) or 1.0
        colour = get_colorscale(display['colormap'])
        base = dict(paper_bgcolor='#ffffff', plot_bgcolor='#fafbfc',
                    font=dict(family='system-ui, sans-serif', size=12, color='#25364b'),
                    margin=dict(l=57,r=88,t=12,b=47), hoverlabel=dict(namelength=0))
        figures = {}
        for plane, energy in projections.items():
            a,b = plane
            xa,ya = edges[a],edges[b]
            if energy.shape != (len(ya)-1,len(xa)-1):
                raise ValueError(f'Projection axis mismatch: {plane}')
            counts = f[f'display/projections/{plane}_hit_count'][:]
            custom = np.empty((*energy.shape,5),dtype=float)
            custom[:,:,0] = xa[:-1][None,:]; custom[:,:,1] = xa[1:][None,:]
            custom[:,:,2] = ya[:-1,None]; custom[:,:,3] = ya[1:,None]
            custom[:,:,4] = counts
            z = energy.astype(object); z[counts==0] = None
            figures[plane] = dict(data=[dict(type='heatmap',x=xa.tolist(),y=ya.tolist(),z=z.tolist(),
                customdata=custom.tolist(),colorscale=colour,zmin=0,zmax=max_energy,
                colorbar=dict(title=dict(text='MeV / pixel',side='right'),thickness=10,len=.8),
                hoverongaps=False,zsmooth=False,
                hovertemplate=(f'{a}: [%{{customdata[0]:.3f}}, %{{customdata[1]:.3f}}] mm<br>'
                               f'{b}: [%{{customdata[2]:.3f}}, %{{customdata[3]:.3f}}] mm<br>'
                               'Energy: %{z:.6g} MeV<br>Hits: %{customdata[4]:.0f}<extra></extra>'))],
                layout={**base,'xaxis':dict(title=dict(text=f'{a} [mm]'),range=ranges[a],zeroline=False),
                        'yaxis':dict(title=dict(text=f'{b} [mm]'),range=ranges[b],zeroline=False)})
        occ=f['display/occupied']
        ix,iy,iz=(occ[k][:]-1 for k in ('ix','iy','iz'))
        e=occ['energy_MeV'][:]; n=occ['hit_count'][:]
        corners=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
        vertices=np.stack([edges[a][idx[:,None]+corners[None,:,j]] for j,(a,idx) in enumerate(zip('xyz',(ix,iy,iz)))],axis=2).reshape(-1,3)
        triangles=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
                            [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        faces=(triangles[None,:,:]+8*np.arange(len(e))[:,None,None]).reshape(-1,3)
        bounds=np.column_stack([edges['x'][ix],edges['x'][ix+1],edges['y'][iy],edges['y'][iy+1],
                                edges['z'][iz],edges['z'][iz+1],e,n])
        mesh=dict(type='mesh3d',x=vertices[:,0].tolist(),y=vertices[:,1].tolist(),z=vertices[:,2].tolist(),
                  i=faces[:,0].tolist(),j=faces[:,1].tolist(),k=faces[:,2].tolist(),
                  intensity=np.repeat(e,8).tolist(),customdata=np.repeat(bounds,8,axis=0).tolist(),
                  colorscale=colour,cmin=0,cmax=float(e.max()) or 1.0,
                  opacity=display['voxel_opacity'],flatshading=True,
                  lighting=dict(ambient=1,diffuse=0,specular=0,fresnel=0,roughness=1),
                  colorbar=dict(title=dict(text='MeV / voxel',side='right'),thickness=10,len=.8),
                  hovertemplate='x: [%{customdata[0]:.3f}, %{customdata[1]:.3f}] mm<br>'
                    'y: [%{customdata[2]:.3f}, %{customdata[3]:.3f}] mm<br>'
                    'z: [%{customdata[4]:.3f}, %{customdata[5]:.3f}] mm<br>'
                    'Energy: %{customdata[6]:.6g} MeV<br>Hits: %{customdata[7]:.0f}<extra></extra>')
        camera=dict(eye=dict(x=1.7,y=1.7,z=1.4),up=dict(x=0,y=0,z=1))
        figures['xyz']=dict(data=[mesh],layout={**base,'margin':dict(l=10,r=88,t=4,b=18),
            'scene':dict(**{a+'axis':dict(title=dict(text=a+' [mm]'),range=ranges[a],showbackground=True,
                backgroundcolor='#fafbfc',gridcolor='#dce3ec',zeroline=False) for a in 'xyz'},
                aspectmode='cube',camera=camera,dragmode='orbit')})
        widths={a:[float(np.diff(v).min()),float(np.diff(v).max())] for a,v in edges.items()}
        return dict(sample=text(f['context/sample'][()]), event_id=str(f['event/event_id'][()]),
                    position=int(f['event/position'][()]),event_count=int(f['event/count'][()]),
                    hit_count=int(f['event/summary/hit_count'][()]),energy_MeV=float(f['event/summary/total_energy_MeV'][()]),
                    frame=text(f['event/frame'][()]),nbins=int(f['configuration/display_nbins'][()]),
                    bin_counts={a:len(v)-1 for a,v in edges.items()},bin_widths_mm=widths,
                    occupied_cells=len(e),figures=figures)
