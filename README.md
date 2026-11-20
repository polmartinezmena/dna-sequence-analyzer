# dna-sequence-analyzer

Herramienta de línea de comandos para analizar secuencias de ADN.
Esta versión, la **0.1.0**, es la etapa 1 del proyecto: lee un archivo
FASTA y, para cada secuencia, dice cómo se llama, cuántos símbolos
tiene y cuál es su contenido en GC.

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

En el ADN, guanina (G) se empareja con citosina (C), y adenina (A) con
timina (T). El par G-C se mantiene con **tres puentes de hidrógeno**. El
par A-T se mantiene con **dos**. Por eso un ADN rico en G y C suele
necesitar más energía para separar las dos hebras: su temperatura de
fusión (Tm) tiende a ser más alta.

Eso no lo explica todo. El apilamiento entre bases vecinas también
cambia la estabilidad, así que el porcentaje de GC es un resumen útil,
no una descripción completa de la molécula. Aun así es uno de los
primeros números que se miran de una secuencia, porque sirve para
varias cosas:

- **Comparar organismos.** Hay genomas muy ricos en GC y genomas muy
  ricos en AT. *Streptomyces* ronda el 70 % de GC. *Plasmodium
  falciparum*, el parásito de la malaria, ronda el 20 %. El genoma
  humano anda cerca del 41 %.
- **Caracterizar una región.** Un trozo con un GC muy distinto del de
  sus vecinos puede ser una pista de transferencia horizontal, de una
  isla CpG o, simplemente, de un gen con otro uso de codones. El número
  no demuestra la causa: invita a mirar más de cerca.
- **Diseñar un experimento.** La Tm de un cebador de PCR depende en
  parte de cuántas G y C lleva. Dos cebadores con un GC muy distinto
  pueden no funcionar bien en la misma reacción.

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
| `S` | "Strong": la base es G o C, pero no se supo cuál. | Cuenta como GC. Las dos opciones tienen tres puentes de hidrógeno. |
| `W` | "Weak": la base es A o T, pero no se supo cuál. | Entra en el total y no cuenta como GC. |
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
dna-analyzer 0.1.0
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
    cli.py                      la tabla y el comando dna-analyzer gc
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

En GitHub, el workflow `Tests` ejecuta lo mismo con Python 3.10, 3.12
y 3.13 en cada pull request.

## Hoja de ruta

1. **Lectura de FASTA, longitud y contenido en GC.** Hecho en la versión 0.1.0.
2. Transcripción a ARN y traducción a proteína en los seis marcos de lectura.
3. Búsqueda de ORFs (open reading frames, marcos de lectura abiertos).
4. BLAST remoto contra NCBI.
5. Línea de comandos pulida, README de publicación e instalación desde PyPI.

Cada etapa se apoya en la anterior. Esta versión se queda en la primera
a propósito: el paquete, las pruebas y el comando tienen que estar
firmes antes de añadir biología más complicada.

## Qué he aprendido

Esta sección es tuya, Pol. Sustituye este párrafo por lo que de verdad
te haya quedado claro al construir la etapa 1. Si necesitas un sitio
por donde empezar, puedes contar:

- cómo está organizado un archivo FASTA y qué hace Biopython al leerlo;
- por qué G y C se cuentan juntas, y por qué una `N` no es lo mismo que
  un 0 %;
- qué tiene que haber en `pyproject.toml` para que `pip install -e .`
  cree el comando `dna-analyzer`;
- qué prueba falló la primera vez y qué te enseñó el fallo.

Cuando lo escribas, borra esta lista.

## Licencia

MIT. Ver [LICENSE](LICENSE).
