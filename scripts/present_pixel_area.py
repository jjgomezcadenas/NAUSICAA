#!/usr/bin/env python3
"""Present the Julia area calculation, separating calculated geometry and reserves."""
from pathlib import Path
import os
import h5py
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'notes/physics/work/matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def main():
    result=ROOT/'results/pixel_area';out=ROOT/'latex_generated/pixel_area';out.mkdir(parents=True,exist_ok=True)
    with h5py.File(result/'analysis.h5') as h:
        a=dict(h.attrs);nom=dict(h['nominal'].attrs);wide=dict(h['wide_margins'].attrs)
        device=h['devices'];names=device['names'].asstr()[:]
        devices=dict(zip(names,zip(device['W_m'][:],device['L_m'][:])))
        rows=[]
        for name,(w,l) in devices.items():
            label=name.replace('_',r'\_')
            rows.append(f'{label} & {w*1e6:g} & {l*1e6:g}'+r' \\')
        scenarios=[]
        for key,label in [('compact','5 µm'),('nominal','10 µm'),('wide_margins','20 µm')]:
            b=dict(h[key].attrs)
            scenarios.append(f"{b['device_margin_m']*1e6:g} & {b['csa_width_m']*1e6:g} & {b['csa_height_m']*1e6:g} & {b['csa_area_m2']*1e6:.4f} & {b['budget_used_m2']*1e6:.4f}"+r' \\')
        region=h['regions'];rnames=region['names'].asstr()[:]
        rects=list(zip(rnames,region['x_m'][:],region['y_m'][:],region['w_m'][:],region['h_m'][:]))
        areas=dict(zip(rnames,region['area_m2'][:]))
        with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
            fig,axes=plt.subplots(1,2,figsize=(9,5.3),layout='constrained')
            ax=axes[0];offset=(a['pitch_m']-a['envelope_side_m'])/2
            ax.add_patch(Rectangle((0,0),a['pitch_m']*1e3,a['pitch_m']*1e3,fill=False,lw=1.5))
            ax.add_patch(Rectangle((offset*1e3,offset*1e3),a['envelope_side_m']*1e3,a['envelope_side_m']*1e3,fill=False,ls='--',lw=1.2))
            labels={'Pad and clearance':'Pad + clearance','Shared services':'Shared-services allowance','Buffer and selection':'Buffer / selection reserve','Control and routing':'Control / routing reserve'}
            for name,x,y,w,hh in rects:
                color='#e7e7e7' if name=='Pad and clearance' else '#f5ead2'
                if name=='CSA strip':color='#e8edf2'
                ax.add_patch(Rectangle(((offset+x)*1e3,(offset+y)*1e3),w*1e3,hh*1e3,facecolor=color,edgecolor='0.5',lw=.6,hatch='//' if name in labels and name!='Pad and clearance' else None))
                if name=='CSA strip':continue
                if name=='Pad and clearance':continue
                ax.text((offset+x+w/2)*1e3,(offset+y+hh/2)*1e3,labels[name],ha='center',va='center',rotation=90 if hh>w else 0,fontsize=7)
            center=a['pitch_m']/2;pad=a['pad_width_m']
            ax.add_patch(Rectangle(((center-pad/2)*1e3,(center-pad/2)*1e3),pad*1e3,pad*1e3,facecolor='#777777',edgecolor='black'))
            ax.text(center*1e3,center*1e3,'Pad\n0.5 × 0.5 mm',ha='center',va='center',color='white',fontsize=8)
            sx=offset+a['envelope_side_m']-nom['csa_width_m'];sy=offset+(a['envelope_side_m']-nom['csa_height_m'])/2
            ax.add_patch(Rectangle((sx*1e3,sy*1e3),nom['csa_width_m']*1e3,nom['csa_height_m']*1e3,facecolor='#31688e',edgecolor='black'))
            ax.text((sx+nom['csa_width_m']/2)*1e3,(sy+nom['csa_height_m']/2)*1e3,'CSA',ha='center',va='center',color='white',rotation=90)
            ax.text(1,.25,'Dashed envelope: 1 mm²\nPixel cell: 4 mm²',ha='center',fontsize=9)
            ax.set(xlim=(-.05,2.05),ylim=(-.05,2.05),xlabel='x (mm)',ylabel='y (mm)',title='Pixel area allocation');ax.set_aspect('equal')
            ax=axes[1];g=h['nominal'];m=nom['device_margin_m']
            for name,x,y,w,hh in zip(g['names'].asstr()[:],g['x_m'][:],g['y_m'][:],g['w_m'][:],g['h_m'][:]):
                ax.add_patch(Rectangle((x*1e6,y*1e6),w*1e6,hh*1e6,facecolor='#d9e5ef',edgecolor='#31688e'))
                if name in devices:
                    dw,dl=devices[name]
                    ax.add_patch(Rectangle(((x+m)*1e6,(y+m)*1e6),min(dw,dl)*1e6,max(dw,dl)*1e6,facecolor='#d8842f'))
                ax.text((x+w/2)*1e6,(y+hh/2)*1e6,name,ha='center',va='center',rotation=90 if hh>w else 0,fontsize=8)
            ax.set(xlim=(0,nom['csa_width_m']*1e6),ylim=(0,nom['csa_height_m']*1e6),xlabel='x (µm)',ylabel='y (µm)',title='CSA packing; 10 µm device margins');ax.set_aspect('equal')
            fig.savefig(out/'pixel_area.pdf',metadata={'CreationDate':None,'ModDate':None})
            fig.savefig(out/'pixel_area.png',dpi=180);plt.close(fig)
        text=r'''\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,booktabs,graphicx,hyperref}
\title{Presupuesto de área para un píxel de NAUSICA}
\author{Registro de cálculo}\date{\today}
\begin{document}\maketitle
El pitch de 2 mm corresponde a una celda cuadrada de @CELL@ mm$^2$.
Se estudia un presupuesto más exigente de @ENVELOPE@ mm$^2$, con el pad de
0,5 mm de lado centrado y toda la electrónica situada fuera de su zona de
separación. Se adopta la arquitectura prevista en el manuscrito: lectura
multiplexada y ADC compartidos fuera del píxel. Los bloques locales de lectura
y control todavía se representan por reservas de área.

\section*{CSA calculado a partir del circuito}
Las dimensiones se leen directamente del netlist M3 x4, M4 x4. Se incluyen
las tres ramas paralelas de M2. M5 utiliza los parámetros efectivos de su
línea del netlist. La suma de áreas de canal es @CHANNELS@ $\mu$m$^2$.
\begin{center}\begin{tabular}{lrr}\toprule
Dispositivo & $W$ ($\mu$m) & $L$ ($\mu$m)\\\midrule
@DEVICES@
\bottomrule\end{tabular}\end{center}
Cada canal se conserva como un rectángulo, con posibilidad de giro de 90 grados.
Se añade un margen por lado para contactos, difusiones, pozos y conexiones.
Los rectángulos se colocan en filas dentro de una franja de 200 $\mu$m de ancho.
Se añade una reserva de 50 por 50 $\mu$m para el condensador y conexiones de
realimentación; es una asignación de espacio, cuya geometría deberá ajustarse
al valor de 3 fF. Entre rectángulos se deja un espacio igual al margen individual,
y el borde exterior del bloque tiene el doble de ese margen.

\begin{center}\begin{tabular}{rrrrr}\toprule
Margen & Ancho CSA & Alto CSA & Área CSA & Presupuesto total\\
($\mu$m) & ($\mu$m) & ($\mu$m) & (mm$^2$) & (mm$^2$)\\\midrule
@SCENARIOS@
\bottomrule\end{tabular}\end{center}
Los márgenes de 5, 10 y 20 $\mu$m son hipótesis de planificación. Como
referencia de escala, el PDK público GF180 especifica contactos de 0,22 $\mu$m,
separación de 0,25 $\mu$m entre contactos y separaciones de pozos del orden
de micras. La verificación de reglas se hará con el PDK del proceso utilizado.

\section*{Reservas para el resto del píxel}
\begin{center}\begin{tabular}{lr}\toprule
Bloque & Área (mm$^2$)\\\midrule
Pad y separación de 50 $\mu$m por lado & @PADKEEP@\\
CSA, margen individual de 10 $\mu$m & @CSA@\\
Reserva de buffer y selección & @READOUT@\\
Reserva de control y conexiones & @CONTROL@\\
Asignación para servicios compartidos & @SHARED@\\\midrule
Total nominal & @TOTAL@\\
Total con margen individual de 20 $\mu$m & @WIDETOTAL@\\
\bottomrule\end{tabular}\end{center}
Las tres últimas asignaciones son presupuestos para bloques pendientes de
diseño. El resultado demuestra espacio geométrico para el CSA simulado y
mantiene reservas explícitas para completar la lectura. La extracción de
parásitos permitirá comprobar el ruido del circuito después del layout.

\begin{figure}[htbp]\centering
\includegraphics[width=\linewidth]{pixel_area.pdf}
\caption{Izquierda: presupuesto de 1 mm$^2$ dentro de una celda de pitch 2 mm.
Las regiones rayadas representan reservas de diseño. La región de servicios
compartidos representa su asignación por píxel. Derecha: colocación geométrica
de los canales y sus envolventes con márgenes de 10 $\mu$m.}
\end{figure}

\section*{Comparación con Topmetal-S y servicios compartidos}
Yang et al. publican un chip completo de 2,632 por 1,794 mm, es decir,
@YANG@ mm$^2$. Incluye electrodo, CSA, buffer, DAC de polarización, control y
51 pads de conexión. Liang et al. publican 2,6 por 1,5 mm, @LIANG@ mm$^2$,
con ADC externo. Son referencias de chip completo.

Como prueba del presupuesto compartido, se toma una matriz ilustrativa de
@N@ píxeles y un área común igual al chip completo de Yang. Su coste
amortizado es @YANGPER@ mm$^2$ por píxel, frente a la reserva de @SHARED@
mm$^2$. Esta comparación fija una escala de área; el circuito común de NAUSICA
se dimensionará según su propia lectura.

La matriz ilustrativa mide @TILE@ por @TILE@ mm. Su reserva compartida total,
@SHAREDTOTAL@ mm$^2$, cabe geométricamente en dos bandas exteriores de @BAND@
mm de ancho, dentro de los márgenes de @EDGE@ mm que quedan junto a los
envolventes centrales. Esta comprobación mantiene las posiciones de los pads
en la retícula de 2 mm. El borde de dado, protección ESD y encapsulado se
concretarán en el diseño de la matriz.

\section*{Reproducción y comprobaciones}
\texttt{config/pixel\_area.toml} contiene las hipótesis y fuentes;
\texttt{src/ion\_signal/pixel\_area.jl} calcula áreas y posiciones y guarda
HDF5. Python genera la figura y este informe a partir de ese archivo.
Se verifican la multiplicidad de M2, las dimensiones efectivas de M5, la suma
de áreas, los límites de cada región, la ausencia de solapamientos y la ley
de escala: duplicar todas las longitudes cuadruplica las áreas.
El README de resultados contiene los comandos completos.

\begin{thebibliography}{9}
\bibitem{yang} Y. Yang et al., JINST 19 (2024) C03031,
\url{https://doi.org/10.1088/1748-0221/19/03/C03031}.
\bibitem{liang} T. Liang et al., JINST 19 (2024) C04004,
\url{https://doi.org/10.1088/1748-0221/19/04/C04004}.
\bibitem{contact} GlobalFoundries, GF180MCU PDK, Contact rules,
\url{https://gf180mcu-pdk.readthedocs.io/en/latest/physical_verification/design_manual/drm_07_13.html}.
\bibitem{well} GlobalFoundries, GF180MCU PDK, Nwell rules,
\url{https://gf180mcu-pdk.readthedocs.io/en/latest/physical_verification/design_manual/drm_07_05.html}.
\end{thebibliography}
\end{document}
'''
        fmt=lambda x:f'{x:.6g}'
        values={'CELL':fmt(a['cell_area_m2']*1e6),'ENVELOPE':fmt(a['envelope_area_m2']*1e6),'CHANNELS':fmt(a['gate_area_m2']*1e12),'DEVICES':'\n'.join(rows),'SCENARIOS':'\n'.join(scenarios),'PADKEEP':fmt(areas['Pad and clearance']*1e6),'CSA':fmt(nom['csa_area_m2']*1e6),'READOUT':fmt(areas['Buffer and selection']*1e6),'CONTROL':fmt(areas['Control and routing']*1e6),'SHARED':fmt(areas['Shared services']*1e6),'TOTAL':fmt(nom['budget_used_m2']*1e6),'WIDETOTAL':fmt(wide['budget_used_m2']*1e6),'YANG':fmt(a['yang_chip_area_m2']*1e6),'LIANG':fmt(a['liang_chip_area_m2']*1e6),'N':str(a['array_pixels']),'YANGPER':fmt(a['yang_area_per_array_pixel_m2']*1e6),'TILE':fmt(a['tile_side_m']*1e3),'SHAREDTOTAL':fmt(a['shared_reserve_total_m2']*1e6),'BAND':fmt(a['shared_edge_band_m']*1e3),'EDGE':fmt(a['edge_free_margin_m']*1e3)}
        for k,v in values.items():text=text.replace('@'+k+'@',v)
        (out/'report.tex').write_text(text)
        (result/'README.md').write_text(f'''# Pixel area budget

Reproduce from the repository root:

```sh
julia --project=. tests/test_pixel_area.jl
julia --project=. scripts/build_pixel_area.jl
python scripts/present_pixel_area.py
```

Configuration: config/pixel_area.toml. Device dimensions are parsed directly
from LTSpice/variants/m4_x4_m3_x4/noise.cir, including multiplicity. Pitch and
pad width come from config/detector.toml. Julia saves numerical results,
rectangle coordinates and source hashes to analysis.h5. Python presents them.

The square 2 mm pitch gives {values['CELL']} mm2 per pixel. The stricter
1 mm2 budget includes the pad and 50 um clearance, transistor envelopes,
a 50 x 50 um feedback/passive reserve, and explicit unimplemented allocations
for buffer/selection, control/routing and shared services. ADCs are shared
outside the pixel. No circuitry is placed under the pad keepout.

Channel area: {values['CHANNELS']} um2. CSA footprint with 10 um device
margins: {values['CSA']} mm2. Total nominal budget: {values['TOTAL']} mm2;
with 20 um device margins: {values['WIDETOTAL']} mm2. Readout allocations
are planning reserves, not measured circuit footprints or completed layouts.
Public GF180 rules supply a scale comparison, not DRC rules for the supplied
180 nm transistor card. Routing, extraction, die edges, ESD and packaging
will be evaluated with the actual process and array design.

Topmetal-S complete-chip areas are {values['YANG']} mm2 (Yang 2024) and
{values['LIANG']} mm2 (Liang 2024). A complete Yang-sized common block
amortized over {values['N']} pixels would be {values['YANGPER']} mm2/pixel.
This is a reference-area comparison; it does not define the NAUSICA common
circuit. Two {values['BAND']} mm edge bands geometrically accommodate the
{values['SHAREDTOTAL']} mm2 common-area allocation in a {values['TILE']} mm
square example array while preserving the 2 mm pad pitch.

Report and figures: latex_generated/pixel_area/. Source generated without
compiling or modifying the manuscript. Tests check multiplicity, physical
areas, containment, overlap rejection and dimensional scaling.
''')
if __name__=='__main__':main()
