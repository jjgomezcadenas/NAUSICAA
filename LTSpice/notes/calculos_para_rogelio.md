# CSA para NAUSICA: puntos pendientes

Los cálculos se han realizado con apoyo de Codex, versión Astra, high.

El cálculo aplica a tu circuito una corriente de iones distribuidos uniformemente durante 10 ms, con una capacitancia externa de 1.25 pF y una capacitancia de realimentación de 3 fF. LTspice calcula el voltaje y el espectro de ruido. Después se promedian los voltajes durante 300 µs en los extremos del frame y se calcula su diferencia, corregida por descarga. El mismo filtro se aplica al espectro de ruido. Se obtiene un ENC de unos 40 electrones RMS y un umbral de 203 iones para S/N > 5.

Para completar el estudio, te pido:

Confirmar los modelos de transistores y el espectro de ruido que deben usarse. La ejecución con los archivos que me mandaste da una contribución flicker nula. Necesito que me pases el espectro numérico de ruido de salida entre 0.1 Hz y 1 MHz para esta configuración, a 20 °C, junto con los modelos y la versión de LTspice que has utilizado. Con esos datos se puede comprobar el ENC.

El artículo describe el cálculo. Los scripts, los modelos y los resultados están en el repositorio.
