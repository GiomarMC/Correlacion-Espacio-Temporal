# Correlación Cruzada Espacial con Retardo Temporal — Difusión del dengue en el Perú

Trabajo de Investigación Formativa (TIF) — Análisis Exploratorio de Datos Espaciales (AEDE),
Escuela Profesional de Ciencia de la Computación, Universidad Nacional de San Agustín de Arequipa.

Este repositorio contiene el **código y la documentación** del análisis de la
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

## Requisitos

- **Python 3.11 o superior.**
- **Conexión a internet** en los pasos 1 y 4 (descargan OpenDengue y GADM).
- **~2 GB de RAM libres** y **~500 MB de disco**: el paso 1 descarga un extracto
  comprimido de 55 MB que se expande a 427 MB, y lo procesa en memoria.

### Instalación

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows:  .venv\Scripts\activate
pip install -r requirements.txt
```

Todos los comandos siguientes asumen el entorno **activado**; si no lo está, usa
`python3` en lugar de `python` (en Linux y macOS `python` a secas puede no existir).

## Ejecución

> **El orden importa.** Casi ningún archivo de datos se versiona: cada script
> consume lo que genera el anterior. Si te saltas un paso, el script te dirá cuál
> falta y cómo generarlo.

```bash
# 1) Descarga OpenDengue + GADM (si faltan) y reconstruye el panel + la matriz W.
python src/20_preparar_dengue.py

# 2) Índice Cruzado de Moran con retardo (0–26 semanas) + LISA local. Genera moran_temporal_dengue.csv.
python src/21_moran_cruzado_dengue.py

# 3) Figuras de la difusión espacio-temporal (secuencia de mapas, correlación foco–vecinas, scatter).
python src/22_visualizacion_difusion.py
```

### Caso de contraste (validación del método)

Para reproducir la validación metodológica —el mismo método aplicado a indicadores
socioeconómicos, que **no** se difunden— se ejecuta además:

```bash
# 4) Matriz W de los 24 departamentos (descarga la geometría GADM nivel 1 si falta).
python src/01_construir_matriz_W.py

# 5) Moran global, LISA y Moran cruzado 2007→2017 sobre los indicadores del INEI.
python src/02_correlacion_cruzada_moran.py

# 6) Figura comparativa dengue (difusión) vs socioeconómico (heterogeneidad).
python src/23_contraste.py
```

> Los pasos 4 y 5 corresponden a los archivos `01_…` y `02_…`: la numeración de los
> archivos agrupa por caso de estudio (`0x` socioeconómico, `2x` dengue), no por orden
> de ejecución. **El paso 4 debe ejecutarse antes que el paso 5.**

El repositorio contiene **solo lo esencial para replicar**: el código, la documentación y un único
archivo de datos (`datos_consolidados/analisis_tecnologia_educacion.csv`, los indicadores del INEI del
caso de contraste, que no son descargables por script). Todo lo demás **no se versiona** porque se
regenera con los scripts: los pasos 1 y 4 descargan OpenDengue y GADM, y los pasos 2, 3, 5 y 6
reconstruyen `resultados_dengue/` y `resultados/`. Requiere conexión a internet en los pasos 1 y 4.

### Si falla la descarga de la geometría

`01_construir_matriz_W.py` tiene un respaldo: si no puede descargar la geometría oficial,
construye la matriz W con una lista de fronteras escrita a mano. Es útil para no bloquear
el pipeline, pero **da 53 aristas en vez de 52 y los índices de Moran no coinciden con los
publicados**. Cuando eso ocurre el script lo advierte de forma destacada al terminar; en ese
caso, restablece la conexión y vuelve a ejecutarlo antes de comparar resultados.

## Estructura

```
src/
  depuru.py                    # utilidades: normalización de nombres, adyacencia de referencia
  aede.py                      # núcleo estadístico (Moran, estandarización, permutación)
  # --- caso principal: dengue ---
  20_preparar_dengue.py        # ETL: descarga OpenDengue + GADM -> panel + W
  21_moran_cruzado_dengue.py   # Índice Cruzado de Moran con retardo + LISA
  22_visualizacion_difusion.py # visualizaciones de difusión espacio-temporal
  # --- caso de contraste: socioeconómico ---
  01_construir_matriz_W.py     # matriz W Queen de los 24 departamentos
  02_correlacion_cruzada_moran.py  # Moran global, LISA y cruzado 2007->2017
  23_contraste.py              # figura comparativa difusión vs heterogeneidad
datos_consolidados/
  analisis_tecnologia_educacion.csv  # indicadores INEI 2007/2017 (único dato versionado)
DATASET.md                     # diccionario de datos y proceso de normalización
LICENSE-datos.md               # licencias de OpenDengue y GADM
requirements.txt
# (datos_dengue/, resultados_dengue/ y resultados/ se generan al correr los scripts)
```

## Datos y licencia

- Casos de dengue: **OpenDengue** v1.3 (CC BY-SA). Cita obligatoria:
  Clarke, J., Lim, A., Gupte, P., et al. (2024). *A global dataset of publicly available dengue case
  count data.* Scientific Data, 11, 296. <https://doi.org/10.1038/s41597-024-03120-7>
- Geometría: **GADM** v4.1, nivel 2 (<https://gadm.org>).

Ver `DATASET.md` para el detalle del origen, tipo de dato y normalización, y `LICENSE-datos.md`.
