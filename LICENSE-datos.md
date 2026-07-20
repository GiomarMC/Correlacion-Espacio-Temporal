# Licencia y atribución de los datos

Este repositorio redistribuye datos procesados derivados de fuentes de terceros. El código propio
del proyecto puede usarse libremente con fines académicos; los datos se rigen por las licencias de
sus fuentes originales.

## OpenDengue (casos de dengue)

- Licencia: **Creative Commons CC BY-SA**.
- Los archivos de `datos_dengue/` derivan del *Spatial extract* v1.3 de OpenDengue mediante el ETL de
  `src/20_preparar_dengue.py`.
- **Cita obligatoria:**
  Clarke, J., Lim, A., Gupte, P., Pigott, D. M., van Panhuis, W. G., & Brady, O. J. (2024).
  *A global dataset of publicly available dengue case count data.* Scientific Data, 11, 296.
  <https://doi.org/10.1038/s41597-024-03120-7>
- Sitio del proyecto: <https://opendengue.org>

## GADM (geometría de provincias)

- Base de datos de áreas administrativas globales, versión 4.1 (<https://gadm.org>).
- Uso académico y no comercial conforme a los términos de GADM. No redistribuir con fines comerciales.

## Nota

`datos_dengue/provincias_dengue.geojson` contiene geometría derivada de GADM; los conteos de casos
provienen de OpenDengue. Cualquier reutilización debe conservar estas atribuciones.
