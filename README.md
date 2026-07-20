# Correlación Cruzada Espacial con Retardo Temporal — Difusión del dengue en el Perú

Trabajo de Investigación Formativa (TIF) — Análisis Exploratorio de Datos Espaciales (AEDE),
Escuela Profesional de Ciencia de la Computación, Universidad Nacional de San Agustín de Arequipa.

Este repositorio contiene el **código, los datos procesados y los resultados** del análisis de la
difusión espacio-temporal del dengue entre las provincias del Perú mediante el **Índice Cruzado de
Moran con retardo temporal**.

## Idea

El Análisis Exploratorio de Datos Espaciales tradicional estudia la autocorrelación espacial en un
único instante. Muchos fenómenos, sin embargo, se **propagan** entre regiones vecinas con un desfase
temporal. El dengue es un caso claro: los brotes surgen en focos y se difunden a las provincias
vecinas en las semanas siguientes. Se extiende el Índice de Moran a una versión **cruzada con retardo**:

```
I(k) = promedio_t [ z_t · (W · z_{t-k}) / N ]
```

donde `z_t` es la incidencia estandarizada entre provincias en la semana `t`, `W` la matriz de
contigüidad Queen (fila-estandarizada) y `k` el retardo en semanas.

## Resultados principales

- **I(0) = 0.187** (p < 0.001, z ≈ 113): la incidencia de dengue se **agrupa** en provincias contiguas.
- `I(k)` se mantiene **positivo y significativo hasta 26 semanas**, decayendo de forma gradual →
  **difusión espacio-temporal**.
- **Focos de difusión** (LISA): eje **Tumbes–Piura** y **selva central** (epicentros reales del dengue peruano).
- **Validación:** el mismo método aplicado a indicadores socioeconómicos da heterogeneidad (no difusión),
  lo que muestra que el método **discrimina** difusión de heterogeneidad.

> Las tablas y figuras (curva I(k), mapas de focos, secuencia del brote, etc.) se generan al ejecutar
> los scripts; se guardan en `resultados_dengue/` (no versionado, ver más abajo).

## Instalación

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Ejecución

```bash
# 1) Descarga OpenDengue + GADM (si faltan) y reconstruye el panel + la matriz W.
python src/20_preparar_dengue.py

# 2) Índice Cruzado de Moran con retardo (0–26 semanas) + LISA local. Genera moran_temporal_dengue.csv.
python src/21_moran_cruzado_dengue.py

# 3) Figuras de la difusión espacio-temporal (secuencia de mapas, correlación foco–vecinas, scatter).
python src/22_visualizacion_difusion.py
```

El repositorio contiene **solo el código y la documentación** (lo esencial para replicar). Los datos y
los resultados **no se versionan** porque se regeneran por completo con los scripts: ejecuta el paso 1
(descarga OpenDengue + GADM y reconstruye `datos_dengue/`) y luego los pasos 2 y 3
(reconstruyen `resultados_dengue/`). Requiere conexión a internet en el paso 1.

## Estructura

```
src/
  depuru.py                    # utilidades: normalización de nombres, adyacencia de referencia
  aede.py                      # núcleo estadístico (Moran, estandarización, permutación)
  20_preparar_dengue.py        # ETL: descarga OpenDengue + GADM -> panel + W
  21_moran_cruzado_dengue.py   # Índice Cruzado de Moran con retardo + LISA
  22_visualizacion_difusion.py # visualizaciones de difusión espacio-temporal
DATASET.md                     # diccionario de datos y proceso de normalización
LICENSE-datos.md               # licencias de OpenDengue y GADM
requirements.txt
# (datos_dengue/ y resultados_dengue/ se generan al correr los scripts; no versionados)
```

## Datos y licencia

- Casos de dengue: **OpenDengue** v1.3 (CC BY-SA). Cita obligatoria:
  Clarke, J., Lim, A., Gupte, P., et al. (2024). *A global dataset of publicly available dengue case
  count data.* Scientific Data, 11, 296. <https://doi.org/10.1038/s41597-024-03120-7>
- Geometría: **GADM** v4.1, nivel 2 (<https://gadm.org>).

Ver `DATASET.md` para el detalle del origen, tipo de dato y normalización, y `LICENSE-datos.md`.
