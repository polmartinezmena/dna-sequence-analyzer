"""Vista aprende."""

import streamlit as st


def render(record, item, stats, checks, results, params, settings, exports):
    settings.caption("Sin parámetros para esta vista.")
    st.markdown(
        """
### Contenido en GC

Es la fracción de las bases contadas que son G o C. Resume la composición de una
secuencia. No mide, por sí solo, lo estable que es la doble hélice.

Las letras ambiguas, como N, no entran en el porcentaje: no sabemos si son G o C.
Sí cuentan en la longitud y en la columna de ambiguas. Una secuencia hecha solo de
N tiene el GC indefinido, no 0 %. Cero significaría "todo es A o T".

### Apilamiento y puentes de hidrógeno

En el dúplex, cada base está emparejada con la de enfrente y apilada sobre el par
de al lado. La estabilidad depende sobre todo de ese **apilamiento**. Los pasos que
contienen G o C suelen ser más estables, y por eso un ADN rico en GC tiende a
separar sus hebras a más temperatura. El orden importa: `GCGCGCGC` y `GGGGCCCC`
tienen el mismo GC y no la misma Tm.

Los puentes de hidrógeno existen (tres entre G y C, dos entre A y T). Aportan menos
a la diferencia neta de lo que sugiere ese recuento, porque al abrir la hélice la
base sigue formando puentes con el agua. Por eso la Tm se calcula sumando cada
dinucleótido, con un modelo de vecinos más próximos, y no con una fórmula que solo
mire el porcentaje de GC.

### Islas CpG

CpG es el dinucleótido C seguida de G. En muchos genomas de vertebrados ese par
es más raro de lo que predice el azar. El cociente observado/esperado compara los
CG que hay con los que tocarían si C y G se colocaran al azar. Un cribado clásico
busca tramos de más de 200 bp, con GC por encima del 50 % y observado/esperado por
encima de 0,6. La vista Resumen muestra los números. No etiqueta islas sola.

### Sesgo GC y replicación bacteriana

El sesgo (G − C) / (G + C) dice si una ventana tiene más G que C. En muchas
bacterias la hebra líder de la replicación es más rica en G. La suma acumulada
baja a lo largo de un replichoro y sube a lo largo del otro. El mínimo suele caer
cerca del origen (oriC) y el máximo cerca del término (ter). Es una pista para
mirar el cromosoma, no una prueba, y no tiene sentido en un cebador de 20 bases.

### Sesgo técnico en PCR e Illumina

Un molde muy rico en GC cuesta más de abrir en una PCR y a menudo amplifica peor.
Uno muy rico en AT se abre con facilidad y puede amplificar con otra eficiencia.
En Illumina la cobertura suele bajar en los dos extremos. El aviso `low-gc` (por
debajo del 25 %) y `high-gc` (por encima del 65 %) es un corte práctico de este
programa, no un umbral del fabricante. Sirve para mirar otra vez esa secuencia
antes de diseñar cebadores o de fiarse de una cobertura.

### ARN

Si hay uracilo y no hay timina, la secuencia se marca como ARN. La U cuenta como
una base que no es GC, igual que la T. Los tallos del ARN también se estabilizan
por apilamiento, y además el ARN forma pares G·U, llamados wobble, que no son el
par Watson-Crick pero aparecen en muchos tallos. Esta aplicación no pliega el ARN.
        """
    )
