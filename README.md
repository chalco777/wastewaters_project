
## Carpeta `pathogen_abundance`

El propósito es obtener patrones temporales de abundancia de los patógenos de interés.

El script `temporal_figures.Rmd` genera las figuras de abundancia temporal de los patógenos de interés.

La data de la cual hace uso el scriptse encuentra en la carpeta `data` e incluye:

- `abundance_data_export.csv` contiene los resultados del análisis de abundancia realizado con `Sylph` (sin aplicar `ConQR`) para todos los secuenciamientos, exceptuando `seq21-11-24` que se encuentra en `final_data`.
- `final_data` contiene el resultado del conteo de reads de `Sylph` para el secuenciamiento `seq21-11-24`.
- Los archivos `countreads` contiene el conteo total de reads de cada muestra de todos los secuenciamientos.
