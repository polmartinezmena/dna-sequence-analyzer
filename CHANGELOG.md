# Historial de cambios

## Próximos cambios

- Rediseño de la web estática con menú lateral plegable, métricas del registro y control de calidad junto al visor.
- Mapa de ORFs por hebra y perfil de GC conectados al selector de región. El visor permite copiar, complementar o traducir la región activa.
- La web ya no muestra una Tm para secuencias de más de 60 pb: indica "No aplica", igual que el límite orientativo de `dna-analyzer tm`.
- El README empieza por la web pública, resume qué hace la herramienta y explica cómo instalarla desde GitHub mientras no esté en PyPI.

## 1.0.0 (sin publicar)

- Los ORFs empiezan por ATG por defecto y omiten los inicios anidados salvo que se soliciten. La CLI y la app permiten activar los inicios alternativos y muestran el codón, la hebra y la posición de inicio.
- Interfaz clara de escritorio con siete vistas, importación en diálogo, registro compartido y exportación desde la cabecera.
- Visor con casillas por base, regla, perfil GC centrado, tramos en seis marcos y anotaciones de ORFs.
- El GC global pondera por bases definidas y las tablas muestran números con formato español.
- La lectura de FASTA omite los comentarios iniciales. La app los indica en el control de calidad y la CLI muestra los errores de entrada en español.
- Búsqueda de ORFs en ambas hebras con coordenadas inclusivas desde 1, inicios alternativos, tramos parciales y exportación FASTA, CSV y GFF3.
- BLAST remoto con MegaBLAST por defecto para nucleótidos, caché local y peticiones serializadas. La app limita las consultas por tamaño y a una cada cinco minutos por sesión.
- Informe HTML con metadatos reproducibles, CI para Python 3.10 a 3.13 y publicación en PyPI mediante etiquetas y Trusted Publishing.
