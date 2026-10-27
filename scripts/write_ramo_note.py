#!/usr/bin/env python3
"""Generate a brief Spanish clarification referring to the NAUSICA paper."""
from pathlib import Path
import h5py
ROOT=Path(__file__).resolve().parents[1]
with h5py.File(ROOT/'results/ramo_check/values.h5') as h:
    v={k:h[k][()] for k in h if k not in ('source_paths','source_sha256')}
fmt=lambda x,n=3:f'{x:.{n}f}'.rstrip('0').rstrip('.')
macros={
 'PadSide':fmt(v['pad_side_m']*1e3),'PadHalf':fmt(v['pad_half_side_m']*1e3),
 'GapSize':fmt(v['gap_m']*1e3),'GapSpeed':fmt(v['gap_speed_m_s']),
 'GapField':fmt(v['gap_field_V_m']/1e5),'Mobility':fmt(v['mobility_m2_V_s']*1e6),
 'PeakCurrent':fmt(v['peak_current_A']*1e15),'GapTime':fmt(v['gap_transit_s']*1e3),
 'GapPhi':fmt(v['gap_weighting_potential'],7),'GapFraction':fmt(v['gap_charge_fraction'],7),
 'FinalTime':fmt(v['final_90percent_full_charge_s']*1e6,1),'FrameTime':fmt(v['frame_s']*1e3),
 'Plateau':fmt(v['uniform_plateau_per_ion_A']*1e15,5),'ExampleCharge':fmt(v['rogelio_charge_for_10ke_C']*1e15,2),
 'ExampleCurrent':fmt(v['rogelio_parallel_plate_current_A']*1e15,2),
 'DerivativeError':f"{v['derivative_relative_error']:.2e}",'ChargeError':f"{v['charge_integral_error']:.2e}"}
text=r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[spanish,es-nodecimaldot,es-noshorthands]{babel}
\usepackage{lmodern,microtype,amsmath,amssymb,siunitx}
\usepackage[margin=23mm]{geometry}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\setlength{\emergencystretch}{2em}
\newcommand{\dd}{\mathrm{d}}
"""+''.join('\\newcommand{\\'+key+'}{'+value+'}\n' for key,value in macros.items())+r"""
\title{Corriente inducida en un pad de NAUSICA}
\author{J. J. G\'omez-Cadenas\\{\small DIPC}}
\date{\today}
\begin{document}
\maketitle
Rogelio, la derivación de la corriente está en la sección
\emph{Signal formation and response of one pad} del artículo
\emph{NAUSICAA for ITACA}. El apéndice \emph{Analytical single-ion response}
recoge las condiciones de contorno, el dibujo del potencial de ponderación
(\emph{weighting potential}) y las comprobaciones del cálculo.

En la expresión del pico de corriente, en \emph{Signal size and duration},
$a$ es la \textbf{mitad del lado del pad}. El pad es cuadrado y su lado mide
$2a=\SI{\PadSide}{mm}$, por tanto
\[
 \boxed{a=\SI{\PadHalf}{mm}=\SI{2.5e-4}{m}.}
\]
El pico de \SI{\PeakCurrent}{fA} corresponde a \textbf{un ion de carga $+e$}
justo antes de llegar al centro del pad, con la velocidad de
\SI{\GapSpeed}{m/s} adoptada en el artículo.

\section*{Comparación con la cuenta de tu esquema}
La expresión $I=Qv/d$ corresponde a un campo de ponderación uniforme,
$W=1/d$, como el de dos placas planas paralelas separadas por $d$.
En tu cálculo asignas al pad la carga total de unos $10^5$ iones.
La geometría de condensador utilizada y la asignación de esa carga total
al pad no representan la señal de entrada de un canal de NAUSICA.
Por tanto, la respuesta del CSA a esa excitación no permite evaluar su
sensibilidad a los iones que llegan a un pad del detector.

Para calcular esa señal necesitamos el campo de ponderación del pad y
los tiempos de llegada de los iones que recoge ese canal. El campo de
ponderación aumenta cerca del pad, y cada ion produce una corriente que
varía durante su aproximación. La sección \emph{Response to uniform arrivals}
del artículo explica cómo construimos la entrada de SPICE sumando estas
corrientes para el caso de llegadas uniformes durante un frame de
\SI{\FrameTime}{ms}. El pico de \SI{\PeakCurrent}{fA} de un ion no se
mantiene durante todo el frame.

Por separado, hay una inconsistencia en la anotación de carga: los
\SI{16}{fC} utilizados corresponden aproximadamente a $10^5$ iones,
es decir, $100\,\mathrm{ke}$. La etiqueta «$10\,\mathrm{ke}$» corresponde
a $10^4$ cargas elementales (\SI{\ExampleCharge}{fC}) y debería corregirse
a «$100\,\mathrm{ke}$».
\end{document}
"""
(ROOT/'LTSpice/ramo_para_rogelio.tex').write_text(text)
out=ROOT/'latex_generated/ramo_check';out.mkdir(parents=True,exist_ok=True)
(out/'report.tex').write_text(text)
print(macros)
