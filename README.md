# dna-sequence-analyzer

Herramienta de línea de comandos para analizar secuencias de ADN.
Esta versión, la **0.1.1**, es la etapa 1 del proyecto: lee un archivo
FASTA y, para cada secuencia, dice cómo se llama, cuántos símbolos
tiene y cuál es su contenido en GC. También estima la temperatura de
fusión de un oligo corto.

El programa está pensado para instalarse con `pip` y usarse en cualquier
ordenador. Las etapas siguientes (traducción, ORFs, BLAST y publicación
en PyPI) todavía no están hechas: la hoja de ruta está más abajo.

Autor: Pol Martinez Mena, estudiante de Biotecnología en la UAB.

## Qué hace

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
extra `[dev]` añade pytest para lanzar las pruebas. Sin `[dev]` el
comando también funciona; simplemente no tendrás pytest.

Comprueba la instalación:

```bash
dna-analyzer --version
```

```text
dna-analyzer 0.1.1
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
id           length  gc_percent
-----------  ------  ----------
NM_000207.3     465       63.87
L29345.1        922       36.23
```

La misma información en CSV, por ejemplo para abrirla en una hoja de
cálculo:

```bash
dna-analyzer gc data/examples.fasta --csv
```

```text
id,length,gc_percent
NM_000207.3,465,63.87
L29345.1,922,36.23
```

Si una secuencia no tiene ninguna base contable, la tabla escribe `NA`
y, por la salida de error, una nota. El CSV deja esa celda vacía para
que un programa no confunda el texto `NA` con un número:

```text
id           length  gc_percent
-----------  ------  ----------
ejemplo           8       50.00
desconocido       4          NA
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
concentración total de oligo de 50 nM. Para subir solo el sodio:

```bash
dna-analyzer tm --sequence GCGCGCGC --na 200
```

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
(`unambiguous_bases`) y cuántas letras ambiguas se dejaron fuera
(`ambiguous_bases`). El campo `name` es el mismo que `id`.

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

## Cómo está organizado el código

```text
pyproject.toml                  nombre del paquete, dependencias y el comando
src/dna_sequence_analyzer/
    __init__.py                 lo que se importa desde fuera
    gc.py                       lectura del FASTA y cálculo del GC
    tm.py                       Tm por vecinos más próximos, para oligos
    cli.py                      las tablas y los comandos gc y tm
    __main__.py                 permite python -m dna_sequence_analyzer
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
5. `dna-analyzer tm` usa el mismo lector y, para cada secuencia válida,
   llama a `Tm_NN` con las condiciones que hayas pasado.

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
con el mismo GC% no dan la misma temperatura, que una `N` no se borra en
silencio, y que una secuencia larga se calcula y queda marcada como
demasiado larga para un cebador.

En GitHub, el workflow `Tests` ejecuta lo mismo con Python 3.10, 3.12
y 3.13 en cada pull request.

## Hoja de ruta

1. **Lectura de FASTA, longitud, contenido en GC y Tm de oligos.** Hecho.
   El GC es la 0.1.0. La Tm por vecinos más próximos es la 0.1.1.
2. Transcripción a ARN y traducción a proteína en los seis marcos de lectura.
3. Búsqueda de ORFs (open reading frames, marcos de lectura abiertos).
4. BLAST remoto contra NCBI.
5. Línea de comandos pulida, README de publicación e instalación desde PyPI.

La transcripción, los ORFs, BLAST y PyPI siguen sin hacerse. La Tm no
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
- qué tiene que haber en `pyproject.toml` para que `pip install -e .`
  cree el comando `dna-analyzer`;
- qué prueba falló la primera vez y qué te enseñó el fallo.

Cuando lo escribas, borra esta lista.

## Licencia

MIT. Ver [LICENSE](LICENSE).
