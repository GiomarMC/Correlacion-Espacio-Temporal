# Documentación del conjunto de datos

**Proyecto:** Correlación Cruzada Espacial con Retardo Temporal — Difusión espacio-temporal del dengue
en el Perú (TIF, curso AEDE, UNSA).

Este documento describe el origen, el tipo de dato y el proceso de normalización del conjunto de datos
usado en el análisis.

---

## 1. Origen y naturaleza de los datos

| Componente | Fuente | Licencia |
|---|---|---|
| Casos de dengue (provincia × semana) | **OpenDengue**, *Spatial extract* v1.3 | CC BY-SA |
| Geometría de provincias del Perú (Admin2) | **GADM** v4.1, nivel 2 | Uso académico (ver GADM) |

- **Fuente primaria (casos):** OpenDengue — base global curada de conteos de dengue de acceso público,
  construida a partir de ministerios de salud, literatura revisada y reportes de vigilancia
  (Clarke et al., 2024, *Scientific Data* 11:296, <https://doi.org/10.1038/s41597-024-03120-7>).
- **Unidad de análisis:** **provincia-semana** (nivel agregado / ecológico). No son personas.
- **Cobertura:** 113 provincias del Perú con notificación de dengue; 1252 semanas (2000-01-02 a 2023-12-24).
- **Total de casos en el panel:** 754 732.

> Los datos crudos **no se versionan** en este repositorio (el extracto de OpenDengue pesa ~54 MB).
> El script `src/20_preparar_dengue.py` los **descarga automáticamente** si faltan, junto con la
> geometría de GADM, y regenera los archivos procesados de `datos_dengue/`.

---

## 2. Diccionario del extracto crudo de OpenDengue (`Spatial_extract_V1_3.csv`)

Campos relevantes (una fila por unidad espacial × ventana temporal):

| Campo | Tipo | Descripción |
|---|---|---|
| `adm_0_name` | texto | País (se filtra `PERU`). |
| `adm_1_name` | texto | Departamento (Admin1). |
| `adm_2_name` | texto | **Provincia (Admin2)** — unidad espacial usada. |
| `calendar_start_date` | fecha | Inicio de la ventana temporal (se usa como fecha de la semana). |
| `calendar_end_date` | fecha | Fin de la ventana temporal. |
| `Year` | entero | Año. |
| `dengue_total` | numérico | **Número de casos de dengue** (variable de interés). |
| `S_res` | texto | Resolución espacial (`Admin0`/`Admin2`…); se filtra `Admin2`. |
| `T_res` | texto | Resolución temporal (`Week`/`Year`); se filtra `Week`. |
| `UUID` | texto | Identificador de la fuente/serie original. |

---

## 3. Proceso de normalización (ETL)

Implementado en `src/20_preparar_dengue.py`. Pasos:

1. **Filtrado.** Se conservan solo las filas con `adm_0_name = PERU`, `S_res = Admin2` y `T_res = Week`.
   → 30 991 filas, 115 provincias.
2. **Normalización de nombres de provincia.** Cada `adm_2_name` se normaliza (minúsculas, sin tildes ni
   espacios) para emparejarlo con la geometría GADM (`NAME_2`). De 115 provincias:
   - **111 casan directamente.**
   - **2 requieren alias:** la provincia de Huánuco figura en GADM como `Huenuco`; `Marañón` llega con
     un carácter corrupto (cirílico `Ѐ`) que se repara a `Ñ`.
   - **2 se descartan** por no existir en GADM v4.1 (provincias de creación reciente): **Datem del
     Marañón** y **Putumayo**.
   → quedan **113 provincias**.
3. **Agregación.** Se agrupa por (provincia, semana) sumando `dengue_total`; se ordena por fecha.
4. **Relleno de huecos.** Las semanas sin notificación se completan con **0 casos** (ausencia de
   reporte = sin casos notificados).
5. **Construcción de la matriz de pesos `W`.** Sobre la geometría de las 113 provincias se calcula la
   contigüidad de primer orden tipo **Queen** (comparten frontera o vértice) con `libpysal`. La matriz
   se **estandariza por filas** (Σ_j w_ij = 1). Las provincias que quedaran aisladas se conectan a su
   vecina más cercana (evita filas nulas). Resultado: 226 aristas, 4.0 vecinos en promedio.

> **Nota sobre estandarización estadística (en el análisis, no en el ETL):** al calcular el Índice de
> Moran (`src/21`, `src/22`), cada semana se transforma con `log(1 + casos)` y luego se estandariza
> **entre provincias** (media 0, desv 1 por semana). Esta *estandarización espacial por periodo*
> implementa el Moran correcto (desviaciones respecto de la media de cada mapa) y elimina el
> componente estacional común, aislando el patrón espacial de difusión.

---

## 4. Archivos procesados (`datos_dengue/`)

| Archivo | Contenido | Dimensión |
|---|---|---|
| `panel_dengue_semanal.csv` | Panel de casos de dengue. Índice = fecha de semana; columnas = provincias. | 1252 × 113 |
| `provincias_dengue.geojson` | Geometría (polígonos) de las 113 provincias, con `provincia` y `departamento`. | 113 |
| `matriz_W_dengue.csv` | Matriz de pesos espaciales W (Queen, fila-estandarizada), etiquetada por provincia. | 113 × 113 |

## 5. Resultados (`resultados_dengue/`)

| Archivo | Contenido |
|---|---|
| `moran_temporal_dengue.csv` | Índice Cruzado de Moran I(k) por retardo (0–26 semanas): I, nº de pares, z, p. |
| `figuras/serie_nacional.png` | Serie temporal nacional de casos (2000–2023) con picos 2017/2023. |
| `figuras/incidencia_media.png` | Coropleta de incidencia media por provincia (con etiquetas). |
| `figuras/I_vs_lag_dengue.png` | Curva I(k) vs retardo temporal (figura principal). |
| `figuras/lisa_dengue.png` | Mapa de focos de difusión (LISA cruzado local, retardo 1 semana). |
| `figuras/scatter_moran_dengue.png` | Diagrama de dispersión de Moran contemporáneo. |
| `figuras/difusion_secuencia.png` | Secuencia semanal de mapas del brote de 2023 (difusión). |
| `figuras/mapa_onset_2023.png` | Frente de difusión: provincias coloreadas por semana de pico en 2023. |
| `figuras/correlacion_foco_vecinas.png` | Serie foco (Piura) vs vecinas + función de correlación cruzada. |

---

## 6. Cómo regenerar todo desde cero

```bash
python src/20_preparar_dengue.py     # descarga OpenDengue + GADM y reconstruye datos_dengue/
python src/21_moran_cruzado_dengue.py # Índice Cruzado de Moran con retardo + LISA
python src/22_visualizacion_difusion.py # figuras de difusión espacio-temporal
```
