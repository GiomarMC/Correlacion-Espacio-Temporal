#!/usr/bin/env python3
"""
01_construir_matriz_W.py
========================
Construye la MATRIZ DE PESOS ESPACIALES W (contigüidad tipo *Queen*) para los 24
departamentos del Perú, a partir de la geometría OFICIAL de límites departamentales.

Flujo:
  1. Descarga (y cachea) la geometría oficial de departamentos (GADM nivel 1).
  2. Wrangling: excluye Callao, fusiona (dissolve) Lima en un único registro, y
     enlaza cada polígono con los 24 departamentos canónicos del estudio.
  3. Deriva la contigüidad *Queen* (comparten frontera o vértice) -> adyacencia A.
  4. Valida A contra la adyacencia de fronteras verificada a mano (depuru.py).
  5. Estandariza W por filas ( Σ_j w_ij = 1 ) y persiste ambos archivos.

Si la descarga falla, usa la adyacencia verificada a mano como respaldo, de modo
que el pipeline nunca queda bloqueado (con aviso explícito en consola).

Salidas (en datos_consolidados/):
  - departamentos.geojson              (geometría cacheada, 24 polígonos)
  - matriz_adyacencia_queen.csv        (binaria 24x24, etiquetada)
  - matriz_W_estandarizada.csv         (fila-estandarizada 24x24, etiquetada)
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import depuru as d  # noqa: E402

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(RAIZ, "datos_consolidados")
GEOJSON_CACHE = os.path.join(DATOS, "departamentos.geojson")

# Fuentes candidatas de geometría oficial (se prueban en orden).
FUENTES_GEOMETRIA = [
    "https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_PER_1.json",
    "https://raw.githubusercontent.com/juaneladio/peru-geojson/master/peru_departamental_simple.geojson",
]


# ---------------------------------------------------------------------------
def descargar_geometria() -> "gpd.GeoDataFrame | None":  # noqa: F821
    """Descarga la geometría oficial (con caché local). Devuelve GeoDataFrame en
    EPSG:4326 o None si no fue posible obtenerla."""
    import geopandas as gpd

    if os.path.exists(GEOJSON_CACHE):
        print(f"[geo] Usando geometría cacheada: {GEOJSON_CACHE}")
        return gpd.read_file(GEOJSON_CACHE)

    import requests

    for url in FUENTES_GEOMETRIA:
        try:
            print(f"[geo] Descargando geometría oficial:\n      {url}")
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            tmp = os.path.join(DATOS, "_descarga_tmp.json")
            with open(tmp, "wb") as f:
                f.write(r.content)
            gdf = gpd.read_file(tmp)
            os.remove(tmp)
            if gdf.crs is None:
                gdf.set_crs("EPSG:4326", inplace=True)
            else:
                gdf = gdf.to_crs("EPSG:4326")
            print(f"[geo] OK — {len(gdf)} unidades crudas descargadas.")
            return gdf
        except Exception as e:  # noqa: BLE001
            print(f"[geo] Falló esa fuente ({type(e).__name__}: {e}). Probando otra…")
    return None


def _columna_nombre(gdf) -> str:
    """Detecta la columna con el nombre del departamento en el GeoDataFrame."""
    for c in ["NAME_1", "NOMBDEP", "nombre", "name", "DEPARTAMEN", "departamento"]:
        if c in gdf.columns:
            return c
    # heurística: primera columna de texto
    for c in gdf.columns:
        if gdf[c].dtype == object and c.lower() != "geometry":
            return c
    raise ValueError(f"No se encontró columna de nombre. Columnas: {list(gdf.columns)}")


def preparar_24_departamentos(gdf):
    """Excluye Callao, fusiona Lima y deja exactamente 24 polígonos en orden
    canónico. Devuelve GeoDataFrame indexado por nombre canónico."""
    import geopandas as gpd

    col = _columna_nombre(gdf)
    # Trabaja sólo con el nombre y la geometría para evitar colisiones de columnas
    # cuando el archivo de entrada ya fue preparado antes (idempotencia con la caché).
    canon = gdf[col].map(d.a_canonico)  # None para Callao / no-departamentos
    base = gpd.GeoDataFrame({"_canon": canon.values}, geometry=gdf.geometry.values, crs=gdf.crs)

    n_desc = int(base["_canon"].isna().sum())
    if n_desc:
        descartados = sorted(set(gdf.loc[canon.isna(), col]))
        print(f"[geo] Descartados (no pertenecen a los 24): {descartados}")
    base = base.dropna(subset=["_canon"])

    # Fusiona geometrías que caen en el mismo departamento canónico (p.ej. Lima).
    base = base.dissolve(by="_canon", as_index=False)

    faltan = set(d.DEPARTAMENTOS) - set(base["_canon"])
    if faltan:
        raise ValueError(f"Faltan departamentos tras el wrangling: {sorted(faltan)}")
    if len(base) != d.N:
        raise ValueError(f"Se esperaban {d.N} polígonos, hay {len(base)}.")

    base = base.set_index("_canon").loc[d.DEPARTAMENTOS].reset_index()
    out = gpd.GeoDataFrame(
        {
            "departamento": base["_canon"].values,
            "ubigeo": [d.UBIGEO[x] for x in base["_canon"].values],
        },
        geometry=base.geometry.values,
        crs=base.crs,
    )
    return out


def queen_desde_geometria(gdf) -> np.ndarray:
    """Matriz de adyacencia binaria 24x24 por contigüidad Queen."""
    try:
        from libpysal.weights import Queen

        w = Queen.from_dataframe(gdf, use_index=False)
        A = w.full()[0]
        A = (np.asarray(A) > 0).astype(int)
        print("[queen] Contigüidad calculada con libpysal.weights.Queen.")
    except Exception as e:  # noqa: BLE001
        print(f"[queen] libpysal no disponible ({e}); usando shapely .touches / intersección.")
        geoms = list(gdf.geometry)
        A = np.zeros((d.N, d.N), dtype=int)
        for i in range(d.N):
            for j in range(i + 1, d.N):
                gi, gj = geoms[i], geoms[j]
                # Queen: comparten cualquier punto de frontera (borde o vértice)
                if gi.touches(gj) or (gi.intersection(gj).length > 0) or gi.intersects(gj):
                    inter = gi.intersection(gj)
                    if not inter.is_empty:
                        A[i, j] = A[j, i] = 1
    np.fill_diagonal(A, 0)
    return A


def matriz_desde_fallback() -> np.ndarray:
    print("[fallback] Usando la adyacencia de fronteras VERIFICADA A MANO (depuru.py).")
    return d.adyacencia_a_matriz(d.ADYACENCIA_VERIFICADA)


def comparar_con_verificada(A: np.ndarray) -> None:
    """Reporta diferencias entre la contigüidad derivada y la verificada a mano."""
    Aref = d.adyacencia_a_matriz(d.ADYACENCIA_VERIFICADA)
    dif = np.argwhere(np.triu(A != Aref, 1))
    if len(dif) == 0:
        print("[valida] La contigüidad Queen COINCIDE 100% con la verificada a mano. ✔")
        return
    print(f"[valida] {len(dif)} pares difieren de la adyacencia verificada a mano:")
    for i, j in dif:
        a, b = d.DEPARTAMENTOS[i], d.DEPARTAMENTOS[j]
        en_geo = "vecinos" if A[i, j] else "NO vecinos"
        en_ref = "vecinos" if Aref[i, j] else "NO vecinos"
        print(f"         {a:<14} ↔ {b:<14} | geometría: {en_geo:<10} | verificada: {en_ref}")


def estandarizar_por_filas(A: np.ndarray) -> np.ndarray:
    grados = A.sum(axis=1, keepdims=True)
    grados[grados == 0] = 1  # evita división por cero (no debería haber aislados)
    return A / grados


# ---------------------------------------------------------------------------
def main() -> None:
    os.makedirs(DATOS, exist_ok=True)
    print("=" * 70)
    print(" CONSTRUCCIÓN DE LA MATRIZ DE PESOS ESPACIALES W (Queen, 24 deptos) ")
    print("=" * 70)

    gdf = descargar_geometria()
    if gdf is not None:
        try:
            gdf = preparar_24_departamentos(gdf)
            gdf.to_file(GEOJSON_CACHE, driver="GeoJSON")
            print(f"[geo] Geometría de 24 departamentos guardada en {GEOJSON_CACHE}")
            A = queen_desde_geometria(gdf)
            comparar_con_verificada(A)
        except Exception as e:  # noqa: BLE001
            print(f"[geo] Error procesando la geometría ({type(e).__name__}: {e}).")
            A = matriz_desde_fallback()
    else:
        print("[geo] No se pudo descargar la geometría oficial.")
        A = matriz_desde_fallback()

    # Chequeos de integridad
    assert np.array_equal(A, A.T), "La matriz de adyacencia NO es simétrica."
    assert A.diagonal().sum() == 0, "La diagonal debe ser 0 (una región no es su propia vecina)."
    aislados = [d.DEPARTAMENTOS[i] for i in range(d.N) if A[i].sum() == 0]
    if aislados:
        print(f"[WARN] Departamentos sin vecinos: {aislados}")

    W = estandarizar_por_filas(A)

    df_A = pd.DataFrame(A, index=d.DEPARTAMENTOS, columns=d.DEPARTAMENTOS)
    df_W = pd.DataFrame(W, index=d.DEPARTAMENTOS, columns=d.DEPARTAMENTOS)
    df_A.to_csv(os.path.join(DATOS, "matriz_adyacencia_queen.csv"))
    df_W.to_csv(os.path.join(DATOS, "matriz_W_estandarizada.csv"))

    print("-" * 70)
    print(f"Vecinos por departamento:\n{df_A.sum(1).to_string()}")
    print("-" * 70)
    print(f"Total de conexiones (aristas): {int(A.sum() // 2)}")
    print(f"Suma de filas de W (debe ser 1): min={W.sum(1).min():.4f}, max={W.sum(1).max():.4f}")
    print("Guardado:")
    print(f"  - {os.path.join(DATOS, 'matriz_adyacencia_queen.csv')}")
    print(f"  - {os.path.join(DATOS, 'matriz_W_estandarizada.csv')}")
    print("=" * 70)


if __name__ == "__main__":
    main()
