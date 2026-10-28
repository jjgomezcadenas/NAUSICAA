"""Build browser plots from Julia's completed MC-energy snapshot.

Only presentation geometry, ranges and colour limits are calculated here.
Energy sums, projections and contributing-hit counts are read from HDF5.
"""
import h5py
import numpy as np
from plotly.colors import get_colorscale


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


COMPONENT_COLOURS=['#d35400','#0072b2','#009e73','#cc79a7','#b59a00','#56b4e9','#8f5a32','#6b55b8']


def graph_payload(f):
    """Read Julia component records and build overlay geometry in each projection.

Stored graph-local indices are one-based. Edge segments are separated by None;
backbones follow the saved node order. Summaries and cumulative distances are
read directly from HDF5. Component colours label topology; heatmaps retain the
energy colour scale. Trace metadata identifies the component and display layer.
"""
    r=f['reconstruction']
    if text(r['status'][()])!='complete':return dict(status='pending')
    settings={k:(text(v[()]) if v.dtype.kind in 'OSU' else v[()].item())
              for k,v in f['configuration/mc_graph'].items()}
    summary={k:v[()].item() for k,v in r['summary'].items()}
    components=[];overlays={p:[] for p in ('xy','xz','yz','xyz')}
    for key in sorted(r['components'],key=int):
        g=r['components'][key];cid=int(key);colour=COMPONENT_COLOURS[(cid-1)%len(COMPONENT_COLOURS)]
        nodes={a:g['nodes/'+a+'_mm'][:].tolist() for a in 'xyz'}
        energy=g['nodes/energy_MeV'][:].tolist()
        src=g['edges/src'][:]-1;dst=g['edges/dst'][:]-1
        path=g['backbone/path_indices'][:]-1;ends=g['backbone/endpoint_indices'][:]-1
        stats={k:v[()].item() for k,v in g['summary'].items()}
        endpoints=[[nodes[a][i] for a in 'xyz'] for i in ends]
        components.append(dict(component_id=cid,colour=colour,endpoints_mm=endpoints,**stats))
        for plane,traces in overlays.items():
            kind='scatter3d' if plane=='xyz' else 'scatter'
            def trace(layer,coords,**style):
                # Plotly's 2D vertical coordinate is y, also for XZ and YZ.
                positions={target:coords[a] for target,a in zip('xyz',plane)}
                return dict(type=kind,**positions,showlegend=False,name=f'Component {cid} · {layer}',
                            meta=dict(component_id=cid,layer=layer),**style)
            edge_coords={a:[v for i,j in zip(src,dst) for v in (nodes[a][i],nodes[a][j],None)] for a in plane}
            traces.append(trace('edges',edge_coords,mode='lines',line=dict(color=colour,width=2),hoverinfo='skip'))
            custom=[[i+1,energy[i]] for i in range(len(energy))]
            traces.append(trace('nodes',nodes,mode='markers',marker=dict(color=colour,size=4),customdata=custom,
                hovertemplate=f'Component {cid}<br>Node %{{customdata[0]}}<br>Energy: %{{customdata[1]:.6g}} MeV<extra></extra>'))
            traces.append(trace('backbone',{a:[nodes[a][i] for i in path] for a in plane},mode='lines',
                line=dict(color=colour,width=6),customdata=g['backbone/s_mm'][:].tolist(),
                hovertemplate=f'Component {cid}<br>Backbone s: %{{customdata:.3f}} mm<extra></extra>'))
            traces.append(trace('endpoints',{a:[nodes[a][i] for i in ends] for a in plane},mode='markers',
                marker=dict(color=colour,size=8,symbol='diamond'),text=['Start','End'] if len(ends) else [],
                hovertemplate=f'Component {cid} · %{{text}}<br>'+'<br>'.join(f'{a}: %{{{target}:.3f}} mm' for target,a in zip('xyz',plane))+'<extra></extra>'))
    return dict(status='complete',settings=settings,components=components,overlays=overlays,
                build_request_id=text(r['build_request_id'][()]),**summary)


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
        transported=bool(f['context/active_transport'][()]) if 'context/active_transport' in f else False
        pads=bool(f['context/active_pads'][()]) if 'context/active_pads' in f else False
        prefix='pads/display' if pads else 'transport/display' if transported else 'display'
        axis_unit={a:'ms' if pads and a=='z' else 'mm' for a in 'xyz'}
        axis_name={a:'t' if pads and a=='z' else a for a in 'xyz'}
        pad_info=dict(status='pending')
        if 'pads' in f:
            pad_info={k:(v[()].tolist() if v.shape else v[()].item()) for k,v in f['pads'].items() if isinstance(v,h5py.Dataset) and k not in ('xy_ions','xy_charge_C','slice_ions','z_edges_mm','time_edges_ms')}
            pad_info['status']='complete'
        quantity='charge_fC' if pads else 'ions' if transported else 'energy_MeV'
        unit='fC' if pads else 'ions' if transported else 'MeV'
        noun='Charge' if pads else 'Ions' if transported else 'Energy'
        count_noun='Ions' if transported else 'Hits'
        transport=dict(status='pending')
        if 'transport' in f:
            transport={k:(text(v[()]) if v.dtype.kind in 'OSU' else v[()].item()) for k,v in f['transport'].items() if isinstance(v,h5py.Dataset)}
            transport['status']='complete'
        edges = {a: f[f'{prefix}/axes/{a}_edges_{axis_unit[a]}'][:] for a in 'xyz'}
        ranges = {}
        for a in 'xyz':
            hits = f[f'event/hits/{a}_mm'][:]
            lo, hi = float(hits.min()), float(hits.max())
            if transported:
                lo,hi=float(edges[a][0]),float(edges[a][-1])
            if pads:
                index=f[prefix+'/occupied/i'+a][:]
                if len(index):lo,hi=float(edges[a][int(index.min())-1]),float(edges[a][int(index.max())])
            mid = (lo + hi) / 2
            half = max((hi-lo)*float(display['extent_scale'])/2, .325)
            ranges[a] = [min(mid-half, float(edges[a][0])), max(mid+half, float(edges[a][-1]))]
        projections = {a: f[f'{prefix}/projections/{a}_{quantity}'][:] for a in ('xy','xz','yz')}
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
            counts = f[f'{prefix}/projections/{plane}_hit_count'][:]
            custom = np.empty((*energy.shape,5),dtype=float)
            custom[:,:,0] = xa[:-1][None,:]; custom[:,:,1] = xa[1:][None,:]
            custom[:,:,2] = ya[:-1,None]; custom[:,:,3] = ya[1:,None]
            custom[:,:,4] = counts
            z = energy.astype(object); z[counts==0] = None
            figures[plane] = dict(data=[dict(type='heatmap',x=xa.tolist(),y=ya.tolist(),z=z.tolist(),
                customdata=custom.tolist(),colorscale=colour,zmin=0,zmax=max_energy,
                colorbar=dict(title=dict(text=unit+' / pixel',side='right'),thickness=10,len=.8),
                hoverongaps=False,zsmooth=False,
                hovertemplate=(f'{axis_name[a]}: [%{{customdata[0]:.3f}}, %{{customdata[1]:.3f}}] {axis_unit[a]}<br>'
                               f'{axis_name[b]}: [%{{customdata[2]:.3f}}, %{{customdata[3]:.3f}}] {axis_unit[b]}<br>'
                               +noun+': %{z:.6g} '+unit+'<br>'+count_noun+': %{customdata[4]:.0f}<extra></extra>'))],
                layout={**base,'xaxis':dict(title=dict(text=f'{axis_name[a]} [{axis_unit[a]}]'),range=ranges[a],zeroline=False),
                        'yaxis':dict(title=dict(text=f'{axis_name[b]} [{axis_unit[b]}]'),range=ranges[b],zeroline=False)})
        occ=f[prefix+'/occupied']
        ix,iy,iz=(occ[k][:]-1 for k in ('ix','iy','iz'))
        e=occ[quantity][:]; n=occ['hit_count'][:]
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
                  colorscale=colour,cmin=0,cmax=(float(e.max()) if len(e) else 1.0) or 1.0,
                  opacity=display['voxel_opacity'],flatshading=True,
                  lighting=dict(ambient=1,diffuse=0,specular=0,fresnel=0,roughness=1),
                  colorbar=dict(title=dict(text=unit+' / cell',side='right'),thickness=10,len=.8),
                  hovertemplate='x: [%{customdata[0]:.3f}, %{customdata[1]:.3f}] mm<br>'
                    'y: [%{customdata[2]:.3f}, %{customdata[3]:.3f}] mm<br>'
                    +('t: [%{customdata[4]:.3f}, %{customdata[5]:.3f}] ms<br>' if pads else 'z: [%{customdata[4]:.3f}, %{customdata[5]:.3f}] mm<br>')
                    +noun+': %{customdata[6]:.6g} '+unit+'<br>'+count_noun+': %{customdata[7]:.0f}<extra></extra>')
        camera=dict(eye=dict(x=1.7,y=1.7,z=1.4),up=dict(x=0,y=0,z=1))
        figures['xyz']=dict(data=[mesh],layout={**base,'margin':dict(l=10,r=88,t=4,b=18),
            'scene':dict(**{a+'axis':dict(title=dict(text=axis_name[a]+' ['+axis_unit[a]+']'),range=ranges[a],showbackground=True,
                backgroundcolor='#fafbfc',gridcolor='#dce3ec',zeroline=False) for a in 'xyz'},
                aspectmode='cube',camera=camera,dragmode='orbit')})
        widths={a:[float(np.diff(v).min()),float(np.diff(v).max())] for a,v in edges.items()}
        return dict(sample=text(f['context/sample'][()]), event_id=str(f['event/event_id'][()]),
                    position=int(f['event/position'][()]),event_count=int(f['event/count'][()]),
                    hit_count=int(f['event/summary/hit_count'][()]),energy_MeV=float(f['event/summary/total_energy_MeV'][()]),
                    frame=text(f['event/frame'][()]),nbins=int(f['configuration/display_nbins'][()]),
                    bin_counts={a:len(v)-1 for a,v in edges.items()},bin_widths=widths,axis_units=axis_unit,
                    occupied_cells=len(e),figures=figures,graph=graph_payload(f),transport=transport,pads=pad_info)
