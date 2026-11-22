# dna-sequence-analyzer

Herramienta de línea de comandos para analizar secuencias de ADN.
Esta versión, la **0.6.0**, lee un archivo FASTA y, para cada secuencia,
dice cómo se llama, cuántos símbolos tiene, cuál es su contenido en GC
y cuántas bases ambiguas lleva. También parte la secuencia en ventanas,
calcula el sesgo GC y estima la temperatura de fusión de un oligo corto.

El programa está pensado para instalarse con `pip` y usarse en cualquier
ordenador. Las etapas siguientes (BLAST y publicación
en PyPI) todavía no están hechas: la hoja de ruta está más abajo.

Autor: Pol Martinez Mena, estudiante de Biotecnología en la UAB.

## Qué hace

### Marcos abiertos de lectura (ORFs)

```powershell
uv run dna-analyzer orfs data/examples.fasta --min-aa 30 --table 1 --format csv
```

Busca candidatos completos en los seis marcos: un inicio y la primera parada
inequívoca en el mismo marco. Sigue las tablas genéticas de Biopython y los
criterios de búsqueda de [NCBI ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/).
Con `--alternative-starts` habilitas los inicios de la tabla. Con
`--exclude-nested` conservas el inicio más largo entre los que comparten marco y parada.
La pestaña ORFs muestra dónde se sitúan y permite descargar CSV, FASTA de
proteínas o GFF3. En la CLI usa `--format fasta` o `--format gff`.
Las coordenadas exportadas empiezan en 1 y ambos extremos son inclusivos.
El intervalo incluye la parada y la proteína no la incluye.

**Ejemplo reproducible.** En la app, pulsa **Cargar ejemplo**, abre **ORFs**
y selecciona `NM_000207.3`. Con la tabla 1, un mínimo de 30 aa e inicios
solo ATG, el marco `+3` muestra un candidato de 110 aa entre las posiciones
60 y 392. Es un intervalo de 333 nucleótidos: 330 codifican los 110
aminoácidos y los 3 restantes forman el codón de parada.

Para ver ese mismo resultado desde PowerShell:

```powershell
uv run dna-analyzer orfs data/examples.fasta --min-aa 30 --table 1 --format csv |
    ConvertFrom-Csv |
    Where-Object { $_.id -eq 'NM_000207.3' -and $_.frame -eq '+3' } |
    Format-Table id, frame, start_1based, end_inclusive, length_aa
```

Verás también un candidato de 106 aa que empieza en 72 y termina en 392.
Al desactivar **Incluir ORFs anidados** (o usar `--exclude-nested` en la CLI)
queda el de 110 aa. Esa diferencia muestra el efecto del filtro de inicios.

### Entrada FASTA

Un archivo FASTA guarda una o varias secuencias. Cada una empieza por
una línea de cabecera que comienza con `>`:

```text
>NM_000207.3 Homo sapiens insulin (INS), transcript variant 1, mRNA
AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGATCACT...
```

El identificador es la primera palabra después de `>`. En el ejemplo,
`NM_000207.3`. El resto de la línea es la descripción. La secuencia
puede ocupar varias líneas: esos saltos de línea no forman parte de la
molécula, y el programa los junta.

El comando `dna-analyzer gc` recorre todos los registros del archivo e
imprime una tabla con tres columnas:

| Columna | Significado |
| --- | --- |
| `id` | Nombre de la secuencia (la primera palabra de la cabecera). |
| `length` | Longitud: cuántos símbolos tiene, contando también `N`, huecos y cualquier otra letra. |
| `gc_percent` | Contenido en GC, en tanto por ciento, con dos decimales. `NA` si no se puede calcular. |
| `ambiguous` | Cuántas letras ambiguas hay (`N`, `R`, `Y`, `K`, `M`, `B`, `D`, `H`, `V`, `X`). Siguen fuera del porcentaje. |
| `type` | `dna` si hay T y no hay U, `rna` si hay U y no hay T, `mixed` si hay las dos, `unspecified` si no hay ni T ni U. |
| `flag` | `low-gc` por debajo del 25 % o `high-gc` por encima del 65 %. Vacío en el medio. |

## La biología del contenido en GC

El contenido en GC es la fracción de las bases contadas que son G o C.
Resume la composición. No mide, por sí solo, lo estable que es la doble
hélice.

En el dúplex, cada base está emparejada con la de enfrente y, además,
apilada sobre el par de al lado. La estabilidad depende sobre todo de
ese **apilamiento** entre pares adyacentes (en inglés, *base stacking*).
Los pasos que contienen G o C suelen apilarse con más fuerza, y por eso
un ADN rico en GC tiende a separar sus hebras a más temperatura. "Suele"
no significa "siempre": el orden importa. `GCGCGCGC` y `GGGGCCCC` tienen
las dos un 100 % de GC y no se funden a la misma temperatura, porque los
contactos entre vecinos no son los mismos.

Los puentes de hidrógeno también existen: tres entre G y C, dos entre A
y T. Aportan menos a la estabilidad **neta** de lo que sugiere ese
recuento. Una base que ya no está emparejada con su complementaria sigue
formando puentes de hidrógeno con el agua. Al abrir la hélice se pierden
los puentes entre bases y se ganan otros con el disolvente. La diferencia
que queda es pequeña comparada con el apilamiento.

Por eso la temperatura de fusión (Tm) no sale de una fórmula que solo
mire el porcentaje de GC. Se calcula con un modelo de **vecinos más
próximos** (*nearest-neighbor*): se suma la contribución medida de cada
dinucleótido, AA, AT, GC, CG y el resto, y se corrige por la sal y por
la concentración del oligo. El comando `dna-analyzer tm` hace ese
cálculo. Sirve para oligos y cebadores cortos, no para un ARNm de
cientos de bases.

El porcentaje de GC sigue siendo un buen resumen de composición:

- **Comparar organismos.** Hay genomas muy ricos en GC y genomas muy
  ricos en AT. *Streptomyces* ronda el 70 % de GC. *Plasmodium
  falciparum*, el parásito de la malaria, ronda el 20 %. El genoma
  humano anda cerca del 41 %.
- **Caracterizar una región.** Un trozo con un GC muy distinto del de
  sus vecinos puede ser una pista de transferencia horizontal, de una
  isla CpG o, simplemente, de un gen con otro uso de codones. El número
  no demuestra la causa: invita a mirar más de cerca.
- **No sustituye a la Tm de un cebador.** Dos oligos con el mismo GC
  pueden fundir a temperaturas distintas. Para eso está `dna-analyzer tm`.

El contenido en GC de **un transcrito no es el GC del genoma entero**.
Los dos ejemplos de `data/examples.fasta` lo enseñan: el ARNm de la
insulina humana sale al 63,87 %, bastante por encima del genoma humano,
y el ARNm de la GFP de la medusa *Aequorea victoria* sale al 36,23 %.
Son la composición de esos mensajeros, no una ficha del organismo.

## Cómo se calcula

La fórmula es:

```text
GC% = 100 * (G + C + S) / (A + C + G + T + U + S + W)
```

No todas las letras de un FASTA son una base decidida. La regla, la
misma que usa Biopython en `gc_fraction(..., ambiguous="remove")`, es:

| Letras | Qué significan | Qué hace el programa |
| --- | --- | --- |
| `A`, `C`, `G`, `T` | Bases de ADN decididas. | Entran en la fórmula. G y C suman al numerador. |
| `U` | Uracilo, el equivalente del T en el ARN. | Cuenta como T: entra en el total y no suma como GC. Así un FASTA de ARN también da un resultado sensato. |
| `S` | La base es G o C, pero no se supo cuál. El apodo antiguo de esta letra es "strong". | Cuenta como GC porque la base es G o C. Ese apodo es histórico: no describe de dónde sale la estabilidad del dúplex. |
| `W` | La base es A o T, pero no se supo cuál. El apodo antiguo es "weak". | Entra en el total y no cuenta como GC, porque es A o T. |
| `N`, `R`, `Y`, `K`, `M`, `B`, `D`, `H`, `V`, `X` | Podría ser un par GC o un par AT. `N` es "cualquier base". `X` aparece a veces con el mismo sentido. | No entran ni arriba ni abajo. Contarlas abajo haría bajar el porcentaje sin motivo. Contarlas como G o C sería inventar el dato. |
| `-`, `.`, espacios y el resto | Huecos de un alineamiento u otros símbolos. | No entran en el porcentaje. Sí cuentan en la longitud. |
| minúsculas | La misma base. A menudo marcan regiones repetidas. | Se tratan igual que las mayúsculas: `atgc` y `ATGC` dan el mismo GC. |

Si después de aplicar la regla no queda ninguna letra contable, el GC
**no es 0 %**. Cero significaría "todo es A o T". El programa devuelve
un valor vacío y la tabla muestra `NA`. Biopython, en ese mismo caso,
devuelve 0; aquí se separan las dos situaciones a propósito.

Cuatro secuencias cortas, para ver la regla en acción:

```text
ATGC      longitud 4    GC% = 100 * 2 / 4 = 50
ATGCNNNN  longitud 8    GC% = 100 * 2 / 4 = 50    (las N no entran)
SSWW      longitud 4    GC% = 100 * 2 / 4 = 50    (S es GC, W es AT)
NNNN      longitud 4    GC% = NA
AAAA      longitud 4    GC% = 0
```

## Instalación

Hace falta Python 3.10 o posterior. Desde la carpeta del repositorio:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

En Windows, el entorno se activa con `.venv\Scripts\activate`.

`pip install -e .` instala el programa en modo editable: si cambias el
código, el comando `dna-analyzer` usa ese cambio sin reinstalar. El
extra `[dev]` añade pytest para lanzar las pruebas. El extra `[plot]`
añade matplotlib, solo si quieres los PNG de `--plot`. Sin esos extras
el comando también funciona.

Comprueba la instalación:

```bash
dna-analyzer --version
```

```text
dna-analyzer 0.6.0
```

## Uso

El archivo de ejemplo, `data/examples.fasta`, contiene dos secuencias
públicas copiadas de NCBI Nucleotide, sin modificar las bases:

- `NM_000207.3`, *Homo sapiens* insulin (INS), transcript variant 1, mRNA.
  465 bases. RefSeq.
- `L29345.1`, *Aequorea victoria* green-fluorescent protein (GFP) mRNA,
  complete cds. 922 bases. GenBank.

```bash
dna-analyzer gc data/examples.fasta
```

```text
id           length  gc_percent  ambiguous  type  flag
-----------  ------  ----------  ---------  ----  ----
NM_000207.3     465       63.87          0   dna
L29345.1        922       36.23          0   dna
```

La misma información en CSV, por ejemplo para abrirla en una hoja de
cálculo:

```bash
dna-analyzer gc data/examples.fasta --csv
```

```text
id,length,gc_percent,ambiguous,type,flag
NM_000207.3,465,63.87,0,dna,
L29345.1,922,36.23,0,dna,
```

Si una secuencia no tiene ninguna base contable, la tabla escribe `NA`
y, por la salida de error, una nota. El CSV deja esa celda vacía para
que un programa no confunda el texto `NA` con un número:

```text
id           length  gc_percent  ambiguous         type  flag
-----------  ------  ----------  ---------  -----------  ----
ejemplo           8       50.00          4          dna
desconocido       4          NA          4  unspecified
Note: GC% is NA when a record has no bases that can be counted (A, C, G, T, U, S or W).
```

Los mensajes del programa están en inglés para que cualquiera pueda
usar la herramienta. Esta documentación está en español.

### Temperatura de fusión

`dna-analyzer tm` llama a `Bio.SeqUtils.MeltingTemp.Tm_NN`. La tabla
termodinámica es DNA_NN3 (Allawi y SantaLucia, 1997), el defecto de
Biopython: un modelo de ADN contra ADN. Una `U` se lee como `T`; no se
usan las tablas de ARN.

Estas dos secuencias tienen el mismo GC y distinta Tm, porque el orden
cambia el apilamiento:

```bash
dna-analyzer tm --sequence GCGCGCGC
dna-analyzer tm --sequence GGGGCCCC
```

En un archivo, con esos nombres:

```text
>alterno
GCGCGCGC
>bloque
GGGGCCCC
```

```bash
dna-analyzer tm oligos.fasta
```

```text
id       length  tm_celsius
-------  ------  ----------
alterno       8       36.46
bloque        8       24.87
```

Un oligo de tamaño más parecido al de un cebador, el ejemplo que usa
la propia documentación de Biopython
(`CGTTCCAAAGATGTGGGCATGAGCTTAC`, 28 bases):

```bash
dna-analyzer tm --sequence CGTTCCAAAGATGTGGGCATGAGCTTAC
```

```text
id        length  tm_celsius
--------  ------  ----------
sequence      28       60.32
```

La temperatura sale en grados Celsius. Las condiciones por defecto son
las de `Tm_NN`. Las concentraciones de iones van en milimolar (mM) y
las de las hebras en nanomolar (nM):

| Opción | Qué es | Defecto |
| --- | --- | --- |
| `--na` | Na+ | 50 mM |
| `--k` | K+ | 0 |
| `--tris` | Tris | 0 |
| `--mg` | Mg2+ | 0 |
| `--dntps` | dNTPs | 0 |
| `--dnac1` | Hebra más concentrada, normalmente el cebador | 25 nM |
| `--dnac2` | Hebra menos concentrada, normalmente la plantilla | 25 nM |
| `--saltcorr` | Corrección de sal. 5 ajusta la entropía según SantaLucia 1998. 0 la apaga. Vale 0-7, las fórmulas de `salt_correction` en Biopython. | 5 |
| `--self-complementary` | El oligo se empareja consigo mismo. Entonces `--dnac2` no se usa. | desactivado |
| `--csv` | La misma tabla en CSV. | desactivado |

En Biopython, 25 nM y 25 nM corresponde a lo que Primer3Plus llama una
concentración total de oligo de 50 nM.

Cada nucleótido lleva un fosfato, y a pH fisiológico ese fosfato está
cargado negativo. Las dos hebras se repelen. Esa repulsión favorece
que el dúplex se abra, así que baja la Tm. Los cationes (Na+, K+ y,
con más efecto por ion, Mg2+) se acumulan alrededor de los fosfatos y
apantallan la repulsión: las hebras se estorban menos y el dúplex
aguanta hasta una temperatura más alta. Los dNTPs también están
cargados negativo y secuestran Mg2+, así que bajan el magnesio libre
que estaba haciendo de pantalla. Por eso `--mg` y `--dntps` son
opciones distintas. Subir la concentración del oligo (`--dnac1` y
`--dnac2`) también sube un poco la Tm, porque a más hebra en el tubo
más favorecido está el dúplex.

El mismo oligo de 28 bases, `CGTTCCAAAGATGTGGGCATGAGCTTAC`, lo enseña:

| Condiciones | Tm |
| --- | --- |
| Defecto: 50 mM Na+, sin magnesio | 60,32 °C |
| `--na 200` | 67,34 °C |
| `--mg 2` | 67,82 °C |
| `--mg 2 --dntps 1.5` | 65,31 °C |

El magnesio sube la Tm. Parte de ese efecto se pierde cuando los dNTPs
se quedan con parte del Mg2+.

El cálculo es para oligos y cebadores, del orden de 60 nucleótidos o
menos. `data/examples.fasta` son ARNm largos. El comando igual produce
un número, y avisa de que no hay que leerlo como temperatura de PCR:

```text
id           length  tm_celsius
-----------  ------  ----------
NM_000207.3     465       88.30
L29345.1        922       76.31
Note: Tm_NN is meant for short DNA oligos and primers (about 60 nt or fewer), not for long molecules such as mRNA. A longer record is still calculated, but do not use that number as a PCR melting temperature.
```

88 °C no es la temperatura a la que pondrías una PCR de la insulina.
La fórmula se puede evaluar con una secuencia larga; el modelo no
describe un mensajero de ese tamaño.

Si una secuencia tiene una `N`, un hueco o menos de dos bases, esa fila
sale `NA`. Las letras ambiguas no se borran antes de calcular: borrarlas
sería puntuar otro oligo. En CSV la celda queda vacía.

También desde Python:

```python
from dna_sequence_analyzer import melting_temperature

print(round(melting_temperature("GCGCGCGC"), 2))
print(round(melting_temperature("GGGGCCCC"), 2))
```

```text
36.46
24.87
```

También se puede usar como biblioteca de Python:

```python
from dna_sequence_analyzer import analyze_fasta, gc_percent

print(gc_percent("ATGC"))

for sequence in analyze_fasta("data/examples.fasta"):
    print(sequence.id, sequence.length, round(sequence.gc_percent, 2))
```

```text
50.0
NM_000207.3 465 63.87
L29345.1 922 36.23
```

`analyze_fasta` devuelve una lista de `SequenceStats`. Además del
nombre, la longitud y el porcentaje, cada elemento guarda cuántas
bases contaron como GC (`gc_bases`), cuántas entraron en el total
(`unambiguous_bases`), cuántas letras ambiguas se dejaron fuera
(`ambiguous_bases`), si la secuencia parece ADN o ARN (`nucleic_acid`)
y si el GC es extremo (`gc_flag`). El campo `name` es el mismo que `id`.

### Ventanas de GC y CpG

El GC de todo el registro es una media. Una isla CpG o un tramo rico
en AT pueden quedar escondidos dentro. `gc-window` corta cada registro
en ventanas de longitud `--window` que avanzan `--step` bases. Las
posiciones son desde 1 e incluyen los dos extremos. Una ventana que no
cabe entera al final no se imprime.

```bash
dna-analyzer gc-window data/examples.fasta --window 200 --step 200
```

```text
id           start  end  gc_percent  cpg_oe  ambiguous
-----------  -----  ---  ----------  ------  ---------
NM_000207.3      1  200       63.50  0.2999          0
NM_000207.3    201  400       64.00  0.2447          0
L29345.1         1  200       39.50  0.3953          0
L29345.1       201  400       38.50  0.5405          0
L29345.1       401  600       33.50  0.3571          0
L29345.1       601  800       40.50  0.2551          0
```

`cpg_oe` es el cociente observado/esperado del dinucleótido CG:

```text
O/E = (número de CG * longitud de la ventana) / (número de C * número de G)
```

Solo cuenta una C seguida de una G. Si en la ventana no hay C o no hay
G, el esperado es cero y la celda sale `NA`. Un cribado clásico de
islas CpG (Gardiner-Garden y Frommer, 1987) busca tramos de más de
200 bp con GC por encima del 50 % y O/E por encima de 0,6. Este
comando no decide si una ventana es una isla: deja los dos números
para que los mires. En el ARNm de la insulina el GC de esas ventanas
es alto, pero el O/E está cerca de 0,3, muy por debajo de 0,6.

`--step` menor que `--window` hace que las ventanas se solapen.
`--csv` escribe la misma tabla. `--plot figura.png` guarda un PNG del
GC a lo largo de la secuencia. El dibujo necesita matplotlib, que no
viene en la instalación básica:

```bash
pip install -e ".[plot]"
dna-analyzer gc-window data/examples.fasta --window 100 --step 100 --plot gc.png
```

### Sesgo GC

El sesgo GC de una ventana es:

```text
sesgo = (G - C) / (G + C)
```

Vale +1 si solo hay G, -1 si solo hay C, y 0 si hay las mismas. La S
no entra: es G o C, pero el signo depende de cuál, y eso no se sabe.
El sesgo acumulado es la suma de esos valores desde el principio del
registro. Cada registro vuelve a empezar de cero.

En muchos cromosomas bacterianos la hebra líder de la replicación es
más rica en G que en C. La curva acumulada baja a lo largo de un
replichoro y sube a lo largo del otro. El mínimo suele caer cerca del
origen de replicación (oriC) y el máximo cerca del término (ter). Es
una pista, no una demostración, y tiene sentido en un cromosoma
bacteriano, no en un cebador de 20 bases.

Un ejemplo corto, solo para ver la aritmética. Ocho G seguidas de
ocho C, en ventanas de 4:

```text
>chrom
GGGGGGGGCCCCCCCC
```

```bash
dna-analyzer skew cromosoma.fasta --window 4 --step 4
```

```text
id     start  end  gc_skew  cumulative_skew
-----  -----  ---  -------  ---------------
chrom      1    4   1.0000           1.0000
chrom      5    8   1.0000           2.0000
chrom      9   12  -1.0000           1.0000
chrom     13   16  -1.0000           0.0000
```

`--plot sesgo.png` dibuja el sesgo de cada ventana y, debajo, la suma
acumulada. También pide el extra `plot`.

### ARN, apilamiento y pares G·U

Si la secuencia lleva uracilo (U) y no lleva timina, la columna `type`
dice `rna`. La U cuenta como base conocida que no es GC, igual que la
T. Si aparecen las dos letras, el tipo es `mixed`. El comando `tm`
lee la U como T para poder usar la tabla de ADN contra ADN, y avisa de
que ese número no es una temperatura de plegamiento.

La estructura secundaria del ARN, los tallos y los bucles, también
depende del apilamiento: un tallo rico en pares G-C suele ser más
estable. Además el ARN forma pares G·U, llamados wobble. No son el
par Watson-Crick G-C, pero son frecuentes y lo bastante estables para
aparecer en un tallo. Este programa no pliega la secuencia y no puntúa
esos pares. Solo los nombra para que no se confunda el GC% con un
predictor de estructura.

### Sesgo técnico del GC

El GC extremo no es solo una curiosidad de composición. En una PCR, un
molde muy rico en GC cuesta más de desnaturalizar y a menudo amplifica
peor. Un molde muy rico en AT se abre con facilidad y puede amplificar
con otra eficiencia. En la secuenciación Illumina la cobertura suele
bajar en las regiones muy ricas en GC y en las muy ricas en AT: la
librería no las representa igual que el resto del genoma.

`dna-analyzer gc` marca `low-gc` por debajo del 25 % y `high-gc` por
encima del 65 %. Esos cortes son un aviso práctico de este programa,
no un umbral del fabricante del secuenciador. Un registro dentro de la
banda también puede dar problemas, y uno fuera de la banda no está
condenado. La marca sirve para mirar otra vez esa secuencia antes de
diseñar cebadores o de fiarse de una cobertura.

### Cuando algo sale mal

Ruta que no existe (el programa termina con código 1):

```bash
dna-analyzer gc data/no-existe.fasta
```

```text
Error: FASTA file not found: data/no-existe.fasta
```

Un archivo vacío da `Error: FASTA file is empty: ...`. Un archivo que
no empieza por `>` da `Error: Could not parse ... as FASTA`. Una carpeta
da `Error: Path is not a file: ...`.

## App web

Carga una muestra FASTA, pega bases o abre el ejemplo desde la bienvenida. El panel lateral reúne la entrada, los parámetros y la descarga del informe HTML. El control de calidad señala registros vacíos, caracteres inválidos, bases ambiguas y GC extremo. Cada informe incluye los resultados de las pestañas, el identificador de ejecución, la fecha UTC, la huella de la muestra, los parámetros y la versión. Cada tabla de resultados se puede descargar en CSV.

La misma herramienta se puede abrir en el navegador, sin escribir comandos
de análisis. La interfaz está en español y llama a las funciones del
paquete: no tiene un cálculo aparte.

En tu ordenador, desde la carpeta del repositorio:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[web]"
dna-analyzer web
```

Se abre en `http://localhost:8501`. El extra `[web]` instala Streamlit, Plotly y Jinja2. Si ya tenías el entorno de las pruebas, basta con
`pip install -e ".[web]"`.

También puedes lanzarla con el archivo que usa Streamlit Cloud:

```bash
streamlit run app/streamlit_app.py
```

Las seis pestañas son **Resumen GC**, **GC por ventanas**, **GC skew**, **Tm de primers**, **Traducción** y **Aprende**. La Tm utiliza su propio oligo y sus concentraciones, independientemente de la muestra FASTA.

### Publicarla gratis en Streamlit Community Cloud

Hace falta una cuenta de GitHub y que el repositorio sea visible para esa
cuenta. El archivo `requirements.txt` de la raíz lista lo que Cloud instala.
El archivo de la app es `app/streamlit_app.py`.

1. Entra en [share.streamlit.io](https://share.streamlit.io) e inicia sesión con GitHub.
2. Pulsa **Create app**.
3. Elige el repositorio `polmartinezmena/dna-sequence-analyzer`.
4. Elige la rama que quieras publicar. Para probar esta versión, `ui-lab-redesign`. Cuando esté en `main`, usa `main`.
5. En **Main file path** escribe `app/streamlit_app.py`.
6. Pulsa **Deploy**. El primer arranque tarda un poco: Cloud instala Biopython, Streamlit y Plotly.
7. Cuando termine, te da una dirección pública. Quien la abra ve la app, no la terminal.

Si Cloud pide un requirements file y no lo detecta, el nombre tiene que ser
exactamente `requirements.txt` en la raíz del repositorio. El tema de color
está en `.streamlit/config.toml`.

## Cómo está organizado el código

```text
pyproject.toml                  nombre del paquete, dependencias y el comando
src/dna_sequence_analyzer/
    __init__.py                 lo que se importa desde fuera
    gc.py                       lectura del FASTA, GC, ambiguas y tipo ADN/ARN
    windows.py                  ventanas de GC, CpG O/E y sesgo GC
    tm.py                       Tm por vecinos más próximos, para oligos
    plots.py                    PNG opcionales; matplotlib no es obligatorio
    cli.py                      las tablas y los comandos, incluido web
    webapp.py                   la interfaz de Streamlit
    web_input.py                pegar texto o leer un archivo subido
    __main__.py                 permite python -m dna_sequence_analyzer
app/streamlit_app.py            archivo que abre Streamlit Cloud
data/examples.fasta             las dos secuencias públicas de ejemplo
tests/                          pruebas automáticas
.github/workflows/tests.yml     ejecuta las pruebas en GitHub
```

El recorrido de una llamada es corto:

1. `cli.py` recibe la ruta del archivo.
2. `analyze_fasta`, en `gc.py`, comprueba que la ruta exista y que sea
   un archivo con contenido. Biopython (`Bio.SeqIO`) separa los
   registros FASTA y junta las líneas de cada secuencia.
3. Para cada secuencia, un bucle mira letra a letra y aplica la tabla
   de más arriba. El porcentaje es `100 * gc_bases / unambiguous_bases`.
4. `cli.py` imprime la tabla, o el CSV si se pasó `--csv`.
5. `dna-analyzer tm` usa el mismo lector y llama a `Tm_NN`. Na+, Mg2+
   y los dNTPs entran en la corrección de sal de Biopython.
6. `gc-window` y `skew` recorren la misma secuencia en bloques. El
   sesgo acumulado se reinicia en cada registro.

Los comentarios del código están en inglés y explican la biología de
cada decisión, no solo la sintaxis.

## Pruebas

Con el extra `[dev]` instalado:

```bash
pytest
```

Las pruebas cubren una secuencia normal, varias secuencias en el mismo
archivo, mayúsculas y minúsculas, las letras `N`, `S` y `W`, el uracilo,
los huecos, un registro vacío, un archivo vacío, una ruta que no existe,
una carpeta y un archivo que no es FASTA. También comprueban el formato
de la tabla y del CSV, y que los dos ejemplos de NCBI dan 297 G/C de
465 y 334 G/C de 922.

Para la Tm comprueban que el número coincide con `Tm_NN`, que dos oligos
con el mismo GC% no dan la misma temperatura, que subir Na+ o Mg2+ sube
la Tm, que los dNTPs le quitan parte de ese efecto al magnesio, que una
`N` no se borra en silencio, y que una secuencia larga se calcula y queda
marcada como demasiado larga para un cebador.

Para las ventanas comprueban el O/E de CpG, el sesgo, que la suma
acumulada vuelve a cero al cambiar de registro, y que una ventana que
no cabe no se inventa. Si matplotlib está instalado, también se
comprueba que `--plot` escribe un PNG.

En GitHub, el workflow `Tests` ejecuta lo mismo con Python 3.10, 3.12
y 3.13 en cada pull request.

## Transcripción y traducción en seis marcos

Al transcribir la hebra codificante de ADN escrita de 5′ a 3′, cada T se
sustituye por U. Por ejemplo, `ATGAAATAG` da el ARN `AUGAAAUAG`. El
código genético asigna un aminoácido o una señal de parada a cada grupo
de tres bases (codón). En la tabla estándar, esos tres codones dan
metionina, lisina y una parada, que se muestran como `MK*`.

Hay tres formas de agrupar una hebra en tripletes según se empiece en
la primera, segunda o tercera base. Son los marcos +1, +2 y +3.
La hebra opuesta es antiparalela. Para leerla también de 5′ a 3′ se
invierte el orden y se toma la base complementaria (A con T, C con G).
Esa **complementaria inversa** admite otros tres marcos, −1, −2 y −3.

```console
dna-analyzer translate data/examples.fasta
dna-analyzer translate data/examples.fasta --frame -2 --table 2
```

La tabla 1 es la estándar y la 2 corresponde a las mitocondrias de
vertebrados. Cada fila cuenta los aminoácidos (incluida X), las paradas
y la longitud del tramo más largo entre paradas o extremos. Una X indica
que el codón ambiguo admite varias interpretaciones, como ocurre con NNN.
Se omiten las bases finales que no forman un triplete completo. En la
pestaña **Traducción** puedes consultar el ARN y ver cada parada resaltada
en los seis marcos de la secuencia elegida.

## Hoja de ruta

1. **Lectura de FASTA, GC, ventanas, sesgo GC y Tm de oligos.** Hecho.
   El GC es la 0.1.0, la Tm la 0.1.1, y las ventanas, el sesgo, la
   columna de ambiguas y el aviso de GC extremo son la 0.2.0.
   La app web, para usar el programa sin terminal, es la 0.3.0.
2. Transcripción a ARN y traducción a proteína en los seis marcos de lectura (disponible desde la 0.4.0).
3. Búsqueda de ORFs (open reading frames, marcos de lectura abiertos).
4. BLAST remoto contra NCBI.
5. Línea de comandos pulida, README de publicación e instalación desde PyPI.

BLAST y la publicación en PyPI siguen pendientes. La Tm no
los sustituye: es un cálculo para oligos cortos, apoyado en el mismo
lector de FASTA.

## Qué he aprendido

Esta sección es tuya, Pol. Sustituye este párrafo por lo que de verdad
te haya quedado claro al construir la etapa 1. Si necesitas un sitio
por donde empezar, puedes contar:

- cómo está organizado un archivo FASTA y qué hace Biopython al leerlo;
- por qué G y C se cuentan juntas, y por qué una `N` no es lo mismo que
  un 0 %;
- por qué el porcentaje de GC no predice la Tm, y qué aporta el
  apilamiento frente a los puentes de hidrógeno;
- por qué el sodio y el magnesio suben la Tm (la repulsión de los
  fosfatos), y para qué sirve el sesgo GC en una bacteria;
- qué tiene que haber en `pyproject.toml` para que `pip install -e .`
  cree el comando `dna-analyzer`;
- qué prueba falló la primera vez y qué te enseñó el fallo.

Cuando lo escribas, borra esta lista.

## Licencia

MIT. Ver [LICENSE](LICENSE).
