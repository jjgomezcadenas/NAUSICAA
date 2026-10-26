# Respuesta y ruido del CSA para NAUSICA

26 de octubre de 2026

Me he apoyado en Codex, versión Astra, high, para estos cálculos. Aunque están repasados, pueden haber errores.

El objetivo es determinar cuántos iones deben llegar a un pad durante 10 ms
para medir su carga con una relación señal/ruido mayor que cinco. Para ello
hemos aplicado al circuito LTspice recibido una corriente calculada a partir
del movimiento de los iones. Calculamos el voltaje de salida, aplicamos el
filtrado elegido para la lectura e integramos el ruido que queda después de
ese mismo filtrado.

Con una capacitancia externa de 1.25 pF (mucho más pequeña de la que tú suponías, pero es la que se corresponde a los pads que tendrá el circuito) y promedios de 300 µs de lectura obtenemos un ruido
equivalente de 40.5 iones RMS y un umbral de 203 iones para S/N > 5. Con 5 pF,
los valores son 51.5 iones RMS y 258 iones. Estos resultados corresponden al
ruido total exportado por nuestra ejecución de LTspice.

## 1. Circuito simulado

Partimos del esquema `CSAfolded180nm.asc` y de los modelos de
`180nm_professional.inc`. Conservamos los transistores, sus dimensiones y sus
polarizaciones, así como la capacitancia de realimentación de 3 fF. En
particular, M5 conserva las dimensiones que tiene asignadas en el esquema.

Para estudiar la señal de iones cambiamos la fuente de corriente y el valor
de la capacitancia externa `Cd1`. El circuito calcula su propia ganancia,
ancho de banda y descarga a través del transistor de realimentación.

Conviene distinguir las capacitancias que aparecen en esta comparación:

| Caso | Capacitancia externa Cd1 | Procedencia |
|---|---:|---|
| Esquema original (Rogelio) | 250 fF | Valor del fichero recibido |
| Primera referencia de nuestro estudio | 5 pF | Carga de referencia de un diseño anterior de Topmetal-S |
| Estimación para el pad de NAUSICA | 1.25 pF | Hipótesis de reducción por área |

En términos prácticos, la diferencia entre 5 pF y 1.25 pF no es relevante, ya que 5' e^- ENS es muy bueno.


Usamos LTspice 17.2.4 nativo en macOS. Las simulaciones de señal y ENC se
hicieron a 20 °C. La versión 26 para macOS con Wine no ejecuta bien.

## 2. Corriente de entrada

La obtengo a partir de los iones que llegan al pad y está documentado en el paper que te mando más tarde:

$$
i_1(t)=e\frac{d\phi_w}{dt},
\qquad
Q_g=\int i_1(t)\,dt=e\,[\phi_w(\mathrm{final})-\phi_w(\mathrm{inicio})].
$$

Para $N$ llegadas distribuidas uniformemente durante un intervalo
$\Delta=10\,\mathrm{ms}$, usamos la corriente media

$$
i_N(t)=\frac{N}{\Delta}\int_0^\Delta i_1(t-t_a)\,dt_a,
$$

donde $t_a$ es el instante de llegada al pad. Esta expresión suma las corrientes
de iones que llegan en distintos momentos.

El programa exporta la corriente como una tabla de tiempos y valores para una
fuente PWL de LTspice.

## 3. Medida de la señal

Estudiamos primero una inyección breve y después la señal de iones de 10 ms.
Las dos medidas usan el mismo circuito, pero necesitan una lectura distinta.

**Inyección breve.** Aplicamos un pulso de corriente de 1 µs con flancos de
10 ns. Ajustamos su amplitud para que su integral sea la carga deseada. A la
salida aplicamos un filtro trapezoidal con tiempo de subida $T$ y meseta $T$.
Medimos la amplitud en el instante del máximo de la respuesta sin ruido. Esto
permite comparar el ruido de carga para una inyección rápida, como en una
caracterización habitual del CSA.

**Iones durante 10 ms.** La carga llega durante todo el intervalo y el voltaje
se va acumulando. Para medir ese cambio promediamos el voltaje en ventanas de
duración $T$ centradas en los extremos del intervalo:

$$
\overline V(t_n)=\frac{1}{T}
\int_{t_n-T/2}^{t_n+T/2} V(t)\,dt.
$$

Formamos entonces la diferencia corregida por descarga

$$
D_n=-\left[\overline V(t_n)-\alpha\overline V(t_{n-1})\right],
\qquad
\alpha=\exp(-\Delta/\tau_0),
\qquad \tau_0=0.3\,\mathrm{s}.
$$

El signo hace positiva la señal para nuestra polaridad de salida. El factor
$\alpha$ procede del estimador del CSA ideal: descuenta la señal anterior que
habría sobrevivido a la descarga. Aquí conservamos ese estimador para comparar
los casos. **Los 0.3 s son un parámetro de la medida, no una constante de
descarga impuesta al circuito LTspice.**

La conversión a carga se calibra con el propio circuito. Llamamos $g$ al valor
de $D_n$ producido por una corriente uniforme normalizada a un ion. Si la
respuesta es lineal, $N$ iones producen aproximadamente $Ng$, y la estimación
del número de iones es $\widehat N=D_n/g$. Así se incluye la reducción de
voltaje debida a la ganancia finita del amplificador.

Un tiempo $T=300\,\mathrm{µs}$ significa, por tanto, dos cosas distintas:
subida del filtro trapezoidal en la inyección breve, y duración de cada promedio
en la lectura de iones. Las ventanas centradas requieren esperar $T/2$ después
del extremo del frame. También requieren acceso al voltaje durante la ventana;
una sola lectura puntual cada 10 ms no basta para construirlas.

## 4. Del espectro de ruido al ENC

LTspice proporciona la densidad espectral de amplitud de ruido de salida,
que denotamos $a_V(f)$, en $\mathrm{V}/\sqrt{\mathrm{Hz}}$. Para obtener la
varianza de nuestra medida elevamos esa densidad al cuadrado y aplicamos el
mismo filtro que a la señal:

$$
\sigma_D^2=\int_{0.1\,\mathrm{Hz}}^{1\,\mathrm{MHz}}
a_V^2(f)\,|H(f)|^2\,df.
$$

En la lectura de frames, el módulo del filtro es

$$
|H(f)|=\left|\operatorname{sinc}(fT)
\left[1-\alpha e^{-2\pi i f\Delta}\right]\right|,
\qquad
\operatorname{sinc}(x)=\frac{\sin(\pi x)}{\pi x}.
$$

El primer factor corresponde al promedio de duración $T$; el segundo, a la
diferencia entre las dos lecturas. Esta expresión incluye la correlación del
ruido entre ventanas. No sumamos las dos varianzas como si fueran lecturas
independientes. Como $\alpha$ no es exactamente uno, el estimador tampoco
rechaza por completo una componente constante; el límite inferior de
integración forma parte de las hipótesis del cálculo.

El ruido equivalente referido a la entrada es

$$
\mathrm{ENC}=\frac{\sigma_D}{g}.
$$

Para la señal de iones, esta razón se expresa en iones de carga $+e$ RMS,
calibrados con la distribución uniforme de llegadas. Para la inyección breve
usamos la amplitud filtrada por carga elemental inyectada. El ENC depende, por
tanto, de cómo llega la carga y de cómo se mide la salida.

En el régimen lineal,

$$
\frac{S}{N}=\frac{Ng}{\sigma_D}=\frac{N}{\mathrm{ENC}}.
$$

Esto da una primera estimación del número de iones necesario para superar
S/N=5. Después simulamos directamente números enteros de iones a ambos lados
del umbral para comprobarlo. El ruido es el del punto de trabajo sin señal;
no hemos incluido ruido de ADC, cuantización o ruido externo de la placa.

## 5. Resultados y efecto del filtrado

| Cd1 externo | T | ENC inyección, e RMS | ENC iones, 10 ms, iones RMS | N mínimo para S/N > 5 |
|---:|---:|---:|---:|---:|
| 5 pF | 180 µs | 44.61 | 58.41 | 293 |
| 5 pF | 200 µs | 42.45 | 56.75 | 284 |
| 5 pF | 250 µs | 38.28 | 53.64 | 269 |
| 5 pF | 300 µs | 35.27 | 51.45 | 258 |
| 5 pF | 350 µs | 33.00 | 49.83 | 250 |
| 1.25 pF | 300 µs | 15.42 | 40.47 | 203 |

La primera estimación, con dos lecturas puntuales y 5 pF, daba unos 315 iones
RMS. Los promedios cortos reducen el ruido rápido que entraba en esas lecturas.
La mejora debe evaluarse después de aplicar el filtro tanto a la señal como
al ruido; eso es lo que hace la tabla.

Al pasar de 5 a 1.25 pF también mejora la respuesta a la carga. La reducción de
ENC no es proporcional a la reducción de capacitancia: permanecen las
capacitancias de los transistores y cambia la respuesta del circuito. Hemos
recalculado los transitorios y el espectro de ruido para la nueva carga externa.

Los resultados sitúan el ENC de la lectura de iones alrededor de 40–50 iones
RMS para los casos estudiados. No debe usarse el valor de 15.42 de la inyección
breve para calcular el umbral de una señal distribuida durante 10 ms. En este
último caso, la simulación da S/N=4.989 para 202 iones y S/N=5.014 para 203.

También comprobamos la respuesta a una población igual de iones desplazada
un frame. Para 1.25 pF y 300 µs, una población en el frame anterior contribuye
un 2.88% de la señal de la población del frame de interés; una en el siguiente,
un 1.01%. Son pruebas separadas de la respuesta a pequeñas señales. Recogen
el efecto de la inducción antes de la llegada, la memoria del circuito y el
promedio temporal; deberán tenerse en cuenta al reconstruir pistas.

## 6. Qué hemos comprobado sobre el ruido flicker

El PDF recibido, `Noise Analysis CSA 180 nm.pdf`, muestra el ruido total de
tensión de salida y el ruido equivalente de corriente de entrada. Hemos
repetido ese cálculo con el esquema original: `Cd1=250 fF`, temperatura por
defecto de 27 °C y barrido de 1 Hz a 1 GHz con 1000 puntos por década.
Esta reproducción es distinta del cálculo de ENC, que usa las condiciones
indicadas anteriormente.

| Frecuencia | Ruido de salida reproducido | Ruido equivalente de corriente de entrada |
|---:|---:|---:|
| 1 Hz | 2.167 mV/√Hz | 0.09287 fA/√Hz |
| 10 Hz | 0.2187 mV/√Hz | 0.09288 fA/√Hz |
| 100 Hz | 0.02207 mV/√Hz | 0.09373 fA/√Hz |

El resultado a 1 Hz es compatible con los aproximadamente 2.2 mV/√Hz que se
ven en el PDF. Es una comparación a la resolución de la figura: no disponemos
del fichero numérico con el que se generó.

En nuestra ejecución, las siete contribuciones que LTspice identifica como
`1overf` son cero. El fichero de modelos recibido no especifica `KF`,
`NOIA/NOIB/NOIC` ni `NOIMOD`. **No hemos desactivado el flicker ni eliminado
contribuciones del espectro: el ENC se calcula con el ruido total `V(onoise)`
que devuelve LTspice.**

Para comprobar la lectura de los datos hicimos una prueba aparte, añadiendo
`KF=1e-25` a ambos modelos en una copia. Las siete contribuciones `1overf`
pasaron a ser no nulas. Ese valor es arbitrario y sirve solo como diagnóstico;
no se utiliza en los resultados de señal y ENC. Además, en ambas ejecuciones
la suma de las varianzas por transistor reproduce la varianza total.

La caída del ruido de salida hacia frecuencias altas no identifica por sí sola
el flicker. En nuestra reproducción, la corriente de ruido equivalente es casi
constante entre 1 y 10 Hz, mientras que el ruido de salida cae aproximadamente
un factor diez. La respuesta del circuito transforma ese ruido de entrada y
produce la caída observada.

Por ello, debemos precisar nuestra afirmación anterior: **nuestra ejecución
de LTspice devuelve contribución flicker cero**. El PDF confirma que se realizó
un análisis de ruido, pero no separa los términos necesarios para resolver
esta cuestión. Conviene comparar con Rogelio la versión de LTspice, el netlist,
los modelos y las contribuciones `1overf` de su ejecución antes de modificar
el modelo o atribuir al circuito físico el ENC calculado.

## 7. Reproducción y material adjunto

Los resultados de señal y ENC están guardados en la rama `LTSpice`, commit
`93694c3`. Los de la comprobación de flicker se añaden en esta entrega. Los
siete ficheros originales de Rogelio permanecen intactos.

Los programas comprueban la integral de la corriente, la convergencia al
reducir el paso temporal y la convergencia de la integración del espectro.
Los tests de filtrado incluyen comparaciones independientes en los dominios
del tiempo y de la frecuencia. Los 52 tests numéricos del repositorio y las
comprobaciones de física pasaron en el estudio de capacitancia. La auditoría
adicional comprueba los términos flicker y la suma de varianzas.

Desde la raíz del repositorio, con Python >=3.11 y las dependencias de
`requirements.txt`:

```sh
PYTHONPATH=src python scripts/audit_ltspice_noise.py --simulator /ruta/al/ejecutable/LTspice
PYTHONPATH=src python scripts/build_shaping_study.py --simulator /ruta/al/ejecutable/LTspice --include-ion-frames
PYTHONPATH=src python scripts/build_capacitance_study.py --simulator /ruta/al/ejecutable/LTspice
make PYTHON=python test physics
```

| Directorio | Contenido |
|---|---|
| `results/ltspice/shaping/` | Barrido de tiempos, tablas, figuras, señales y comprobaciones |
| `results/ltspice/capacitance/pad_1.25pF/` | Comparación de capacitancia y umbral de 203 iones |
| `results/ltspice/noise_audit/` | Reproducción del ruido original y prueba independiente del término flicker |
| `LTSpice/` | Ficheros originales y notas de trabajo |

Los JSON contienen los valores completos, las condiciones y las huellas de
los archivos de entrada. El paquete de revisión incluye los espectros y los
archivos necesarios para repetir los cálculos; el paquete completo incluye
también los transitorios generados.
