#!/usr/bin/env python3
"""
20_preparar_dengue.py
=====================
Prepara el panel espacio-temporal de DENGUE en el Perú a nivel de PROVINCIA
(Admin2) y semana, a partir de la base curada OpenDengue (v1.3), y construye la
matriz de pesos espaciales W (contigüidad Queen entre provincias).

El dengue es un fenómeno que se PROPAGA entre provincias vecinas con retardo
(vector Aedes aegypti + movilidad humana), por lo que es idóneo para demostrar la
correlación cruzada espacial con retardo temporal.

Fuente: OpenDengue, Spatial_extract_V1_3 (https://opendengue.org).
Geometría: GADM 4.1 nivel 2 (provincias del Perú).

Entradas (ya descargadas en datos_dengue/):
  - spatial_extract_V1_3.csv
  - gadm_per_2.json

Salidas (en datos_dengue/):
  - panel_dengue_semanal.csv     (index=semana, columnas=provincias, casos)
  - provincias_dengue.geojson    (geometría de las provincias del estudio)
  - matriz_W_dengue.csv          (Queen fila-estandarizada entre provincias)
"""
from __future__ import annotations

import os
import sys
import unicodedata

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import depuru as d  # noqa: E402  (reutilizamos normaliza())

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(RAIZ, "datos_dengue")

# Alias para nombres de provincia que no casan directamente con GADM.
# clave = nombre normalizado de OpenDengue -> valor = NAME_2 de GADM.
ALIAS = {
    "huanuco": "Huenuco",     # GADM escribe la provincia de Huánuco como "Huenuco"
    "maranon": "Marañón",     # en OpenDengue viene con un carácter corrupto
}


def limpia_provincia(s: str) -> str:
    """Normaliza un nombre de provincia y repara el carácter cirílico corrupto
    (Ѐ, U+0400) que OpenDengue trae en 'MARAÑON'."""
    if not isinstance(s, str):
        return ""
    s = s.replace("Ѐ", "N").replace("Ѐ", "N")  # Marañón corrupto
    return d.normaliza(s)


def cargar_gadm():
    import geopandas as gpd
    g = gpd.read_file(os.path.join(DATOS, "gadm_per_2.json"))[["NAME_1", "NAME_2", "geometry"]]
    g = g.rename(columns={"NAME_1": "departamento", "NAME_2": "provincia"})
    g["_key"] = g["provincia"].map(d.normaliza)
    return g


# URLs de las fuentes públicas (para reproducibilidad sin subir datos crudos al repo).
URL_OPENDENGUE = ("https://github.com/OpenDengue/master-repo/raw/main/"
                  "data/releases/V1.3/Spatial_extract_V1_3.zip")
URL_GADM = "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_PER_2.json"


def asegurar_datos():
    """Descarga (si faltan) el extracto espacial de OpenDengue y la geometría GADM
    nivel 2. Así el repositorio no necesita versionar los datos crudos (54 MB)."""
    import io
    import zipfile

    import requests

    os.makedirs(DATOS, exist_ok=True)
    csv = os.path.join(DATOS, "spatial_extract_V1_3.csv")
    gjson = os.path.join(DATOS, "gadm_per_2.json")

    if not os.path.exists(csv):
        print(f"[descarga] OpenDengue Spatial_extract V1.3…\n           {URL_OPENDENGUE}")
        r = requests.get(URL_OPENDENGUE, headers={"User-Agent": "Mozilla/5.0"}, timeout=180)
        r.raise_for_status()
        z = zipfile.ZipFile(io.BytesIO(r.content))
        nombre = [n for n in z.namelist() if n.lower().endswith(".csv")][0]
        with z.open(nombre) as f, open(csv, "wb") as out:
            out.write(f.read())
        print("[descarga] OpenDengue OK.")
    if not os.path.exists(gjson):
        print(f"[descarga] Geometría GADM nivel 2…\n           {URL_GADM}")
        r = requests.get(URL_GADM, timeout=180)
        r.raise_for_status()
        with open(gjson, "wb") as out:
            out.write(r.content)
        print("[descarga] GADM OK.")


def main():
    import geopandas as gpd

    print("=" * 72)
    print(" PREPARACIÓN DEL PANEL ESPACIO-TEMPORAL DE DENGUE (provincia × semana) ")
    print("=" * 72)

    asegurar_datos()

    # ---- 1) OpenDengue: Peru, Admin2, semanal ----
    df = pd.read_csv(os.path.join(DATOS, "spatial_extract_V1_3.csv"), low_memory=False)
    pe = df[(df.adm_0_name.str.upper() == "PERU") & (df.S_res == "Admin2") & (df.T_res == "Week")].copy()
    pe["semana"] = pd.to_datetime(pe["calendar_start_date"])
    pe["_key"] = pe["adm_2_name"].map(limpia_provincia)
    print(f"[1] OpenDengue Perú Admin2×Week: {len(pe)} filas, "
          f"{pe['_key'].nunique()} provincias, {pe.semana.dt.year.min()}–{pe.semana.dt.year.max()}")

    # ---- 2) Geometría GADM nivel 2 y emparejamiento de nombres ----
    g = cargar_gadm()
    gkeys = set(g["_key"])

    def a_gadm_key(k):
        if k in gkeys:
            return k
        if k in ALIAS:
            return d.normaliza(ALIAS[k])
        return None

    pe["_gkey"] = pe["_key"].map(a_gadm_key)
    sin_match = sorted(pe.loc[pe["_gkey"].isna(), "adm_2_name"].unique())
    print(f"[2] Provincias sin geometría (se descartan): {sin_match}")
    pe = pe.dropna(subset=["_gkey"])

    provincias = sorted(pe["_gkey"].unique())
    g = g[g["_key"].isin(provincias)].drop_duplicates("_key").set_index("_key").loc[provincias].reset_index()
    print(f"    provincias con geometría: {len(provincias)}")

    # ---- 3) Panel provincia × semana (casos de dengue) ----
    panel = (pe.groupby(["semana", "_gkey"])["dengue_total"].sum()
             .unstack("_gkey").reindex(columns=provincias).sort_index())
    panel.index.name = "semana"
    # Rellenar huecos internos con 0 (semana sin reporte = sin casos notificados)
    panel = panel.asfreq("W-SUN") if panel.index.inferred_freq else panel
    panel = panel.fillna(0.0)
    print(f"[3] Panel: {panel.shape[0]} semanas × {panel.shape[1]} provincias "
          f"({panel.index.min().date()}..{panel.index.max().date()})")
    print(f"    casos totales de dengue en el panel: {int(panel.to_numpy().sum()):,}")

    # ---- 4) Matriz W (Queen entre provincias del estudio) ----
    gg = gpd.GeoDataFrame(g[["provincia", "departamento", "geometry"]], geometry="geometry", crs=g.crs)
    A = queen(gg)
    # Conectar posibles provincias aisladas a su vecina más cercana (KNN-1)
    A = conecta_aislados(A, gg)
    W = A / A.sum(axis=1, keepdims=True)

    # ---- Guardado ----
    prov_names = g["provincia"].tolist()
    panel.columns = prov_names
    panel.to_csv(os.path.join(DATOS, "panel_dengue_semanal.csv"))
    gg.assign(provincia=prov_names).to_file(os.path.join(DATOS, "provincias_dengue.geojson"), driver="GeoJSON")
    pd.DataFrame(W, index=prov_names, columns=prov_names).to_csv(os.path.join(DATOS, "matriz_W_dengue.csv"))

    print(f"[4] W Queen construida ({A.sum() // 2:.0f} aristas, "
          f"vecinos: min={int(A.sum(1).min())}, medio={A.sum(1).mean():.1f}, max={int(A.sum(1).max())})")
    print("Guardado en datos_dengue/: panel_dengue_semanal.csv, provincias_dengue.geojson, matriz_W_dengue.csv")
    print("=" * 72)


def queen(gdf) -> np.ndarray:
    from libpysal.weights import Queen
    w = Queen.from_dataframe(gdf, use_index=False)
    A = (np.asarray(w.full()[0]) > 0).astype(float)
    np.fill_diagonal(A, 0)
    return A


def conecta_aislados(A: np.ndarray, gdf) -> np.ndarray:
    aislados = np.nonzero(A.sum(1) == 0)[0]
    if len(aislados) == 0:
        return A
    cent = gdf.geometry.to_crs("EPSG:32718").centroid
    pts = np.c_[cent.x.values, cent.y.values]
    for i in aislados:
        dist = np.linalg.norm(pts - pts[i], axis=1)
        dist[i] = np.inf
        j = int(np.argmin(dist))
        A[i, j] = A[j, i] = 1.0
    print(f"    (conectadas {len(aislados)} provincias aisladas a su vecina más cercana)")
    return A


if __name__ == "__main__":
    main()
