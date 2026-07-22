#!/usr/bin/env python3
"""
02_correlacion_cruzada_moran.py
===============================
Análisis Exploratorio de Datos Espaciales (AEDE) con RETARDO TEMPORAL para el TIF.

A partir de:
  - datos_consolidados/analisis_tecnologia_educacion.csv   (indicadores 2007/2017)
  - datos_consolidados/matriz_W_estandarizada.csv          (W Queen fila-estandarizada)

calcula:
  1. Índice de Moran UNIVARIADO (estático) de cada indicador, como caracterización.
  2. ÍNDICE CRUZADO DE MORAN CON RETARDO TEMPORAL para 4 emparejamientos
     pasado(2007) -> presente(2017):
        (a) Tecnología07 -> Tecnología17     (persistencia digital rezagada)
        (b) Analfabetismo07 -> Analfabetismo17  (persistencia educativa rezagada)
        (c) Tecnología07 -> Analfabetismo17  (cruce tecnología->educación)
        (d) Analfabetismo07 -> Tecnología17  (cruce educación->tecnología)
     Estadístico:  I = (z_presente · W·z_pasado) / n
  3. Significancia por PERMUTACIÓN (999, aleatorización condicional) -> pseudo p-valor
     y z-score respecto de la distribución nula.
  4. LISA CRUZADO local por departamento (I_i = z_presente_i · (W·z_pasado)_i) con
     clasificación Alto-Alto / Bajo-Bajo / Alto-Bajo / Bajo-Alto y p-valor local.

Salidas en resultados/ :
  - moran_global.csv               (univariados + 4 cruzados)
  - lisa_<par>.csv                 (uno por emparejamiento)
  - resumen_capitulo_iv.csv        (tabla lista para el informe)
  - figuras/scatter_<par>.png      (diagrama de dispersión de Moran)
  - figuras/lisa_<par>.png         (mapa de clusters LISA, si hay geometría)
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import depuru as d  # noqa: E402
from aede import (  # noqa: E402
    N_PERM, SEMILLA, moran_cruzado, moran_univariado, permutacion_global,
    pseudo_p, zscore,
)

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(RAIZ, "datos_consolidados")
RES = os.path.join(RAIZ, "resultados")
FIG = os.path.join(RES, "figuras")

# Variables del estudio (columnas del CSV).
V_TEC07 = "TEC07_sin_TIC_%"      # % hogares sin servicio de info/comunicación (proxy 2007)
V_TEC17 = "TEC17_sin_pc_%"       # % hogares sin PC/laptop (2017)
V_EDU07 = "EDU07_analfab_%"      # tasa de analfabetismo 2007
V_EDU17 = "EDU17_analfab_%"      # tasa de analfabetismo 2017

# Emparejamientos pasado -> presente: (clave, etiqueta, var_pasado_2007, var_presente_2017)
PARES = [
    ("tec_tec", "Tecnología'07 → Tecnología'17", V_TEC07, V_TEC17),
    ("edu_edu", "Analfabetismo'07 → Analfabetismo'17", V_EDU07, V_EDU17),
    ("tec_edu", "Tecnología'07 → Analfabetismo'17", V_TEC07, V_EDU17),
    ("edu_tec", "Analfabetismo'07 → Tecnología'17", V_EDU07, V_TEC17),
]


# ---------------------------------------------------------------------------
# Utilidades estadísticas (zscore, moran_*, pseudo_p, permutacion_global
# provienen de aede.py). Aquí solo el LISA cruzado local.
# ---------------------------------------------------------------------------
def lisa_cruzado(z_pres, z_pas, W, seed=SEMILLA, M=N_PERM):
    """LISA cruzado local con permutación condicional por ubicación.
    Devuelve (I_local, lag, p_local, cuadrante[str])."""
    n = len(z_pres)
    lag = W @ z_pas                      # rezago espacial del pasado
    I_local = z_pres * lag               # contribución local
    rng = np.random.default_rng(seed)

    p_local = np.empty(n)
    for i in range(n):
        vecinos = np.nonzero(W[i] > 0)[0]
        pesos = W[i, vecinos]
        # valores del pasado disponibles para reasignar (todos menos i)
        resto = np.delete(z_pas, i)
        sims = np.empty(M)
        ki = len(vecinos)
        for k in range(M):
            muestra = rng.choice(resto, size=ki, replace=False)
            sims[k] = z_pres[i] * float(pesos @ muestra)
        obs = I_local[i]
        larger = int(np.sum(sims >= obs))
        if (M - larger) < larger:
            larger = M - larger
        p_local[i] = (larger + 1.0) / (M + 1.0)

    # Cuadrante: eje x = presente en i, eje y = rezago del pasado en vecinos
    cuad = []
    for i in range(n):
        hp = z_pres[i] > 0
        hl = lag[i] > 0
        cuad.append(
            "Alto-Alto" if (hp and hl) else
            "Bajo-Bajo" if (not hp and not hl) else
            "Alto-Bajo" if (hp and not hl) else
            "Bajo-Alto"
        )
    return I_local, lag, p_local, cuad


# ---------------------------------------------------------------------------
# Carga de datos
# ---------------------------------------------------------------------------
def cargar():
    csv = os.path.join(DATOS, "analisis_tecnologia_educacion.csv")
    df = pd.read_csv(csv)
    df["__c"] = df["Departamento"].map(d.a_canonico)
    faltan = set(d.DEPARTAMENTOS) - set(df["__c"])
    if faltan:
        raise ValueError(f"Faltan departamentos en el CSV: {sorted(faltan)}")
    df = df.set_index("__c").loc[d.DEPARTAMENTOS]

    Wdf = pd.read_csv(os.path.join(DATOS, "matriz_W_estandarizada.csv"), index_col=0)
    Wdf = Wdf.loc[d.DEPARTAMENTOS, d.DEPARTAMENTOS]
    W = Wdf.to_numpy(dtype=float)
    # Verificación de completitud (diseño censal: N=24 sin faltantes)
    for v in (V_TEC07, V_TEC17, V_EDU07, V_EDU17):
        n_na = int(df[v].isna().sum())
        if n_na:
            raise ValueError(f"La variable {v} tiene {n_na} faltantes; se esperaba censo completo.")
    return df, W


# ---------------------------------------------------------------------------
# Figuras
# ---------------------------------------------------------------------------
COLOR_CUAD = {
    "Alto-Alto": "#c0392b", "Bajo-Bajo": "#2c6fbb",
    "Alto-Bajo": "#e59866", "Bajo-Alto": "#7fb3d5", "No signif.": "#d5d8dc",
}


def fig_scatter(clave, etiqueta, z_pres, lag, I, p, cuad, sig):
    fig, ax = plt.subplots(figsize=(7, 6))
    for i, dep in enumerate(d.DEPARTAMENTOS):
        c = COLOR_CUAD[cuad[i]] if sig[i] else COLOR_CUAD["No signif."]
        ax.scatter(z_pres[i], lag[i], color=c, s=60, edgecolor="k", linewidth=0.4, zorder=3)
        ax.annotate(dep, (z_pres[i], lag[i]), fontsize=6.5, xytext=(3, 3),
                    textcoords="offset points")
    ax.axhline(0, color="grey", lw=0.8); ax.axvline(0, color="grey", lw=0.8)
    xs = np.linspace(z_pres.min() - 0.3, z_pres.max() + 0.3, 50)
    ax.plot(xs, I * xs, color="black", lw=1.6, label=f"pendiente = I = {I:.3f}")
    ax.set_xlabel("Presente (z) — valor 2017 en la región")
    ax.set_ylabel("Rezago espacial del pasado — W·z(2007) de las vecinas")
    ax.set_title(f"Moran cruzado con retardo temporal\n{etiqueta}   (I={I:.3f}, p={p:.3f})")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"scatter_{clave}.png"), dpi=140)
    plt.close(fig)


def fig_lisa_mapa(clave, etiqueta, cuad, sig, gdf):
    if gdf is None:
        return
    cat = [cuad[i] if sig[i] else "No signif." for i in range(d.N)]
    g = gdf.copy()
    g["cluster"] = pd.Categorical(cat, categories=list(COLOR_CUAD.keys()))
    fig, ax = plt.subplots(figsize=(7, 8))
    for categoria, color in COLOR_CUAD.items():
        sub = g[g["cluster"] == categoria]
        if len(sub):
            sub.plot(ax=ax, color=color, edgecolor="white", linewidth=0.5,
                     label=f"{categoria} ({len(sub)})")
    ax.set_title(f"Clusters LISA cruzados (retardo 2007→2017)\n{etiqueta}")
    ax.axis("off")
    ax.legend(loc="lower left", fontsize=8, title="Cuadrante (p<0.05)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f"lisa_{clave}.png"), dpi=140)
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    os.makedirs(RES, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    df, W = cargar()
    n = d.N

    try:
        import geopandas as gpd
        gpath = os.path.join(DATOS, "departamentos.geojson")
        gdf = gpd.read_file(gpath).set_index("departamento").loc[d.DEPARTAMENTOS].reset_index() \
            if os.path.exists(gpath) else None
    except Exception:  # noqa: BLE001
        gdf = None

    # z-scores de todas las variables
    Z = {v: zscore(df[v].to_numpy()) for v in (V_TEC07, V_TEC17, V_EDU07, V_EDU17)}
    EI = -1.0 / (n - 1)  # esperanza bajo H0

    print("=" * 74)
    print(" ÍNDICE DE MORAN — AUTOCORRELACIÓN Y CORRELACIÓN CRUZADA CON RETARDO ")
    print(f" N = {n} departamentos | {N_PERM} permutaciones | E[I]=-1/(N-1)={EI:.4f}")
    print("=" * 74)

    filas = []

    # ---- 1) Moran univariado (caracterización estática) ----
    print("\n[1] Moran univariado (estático):")
    etiquetas_uni = {V_TEC07: "Tecnología 2007", V_TEC17: "Tecnología 2017",
                     V_EDU07: "Analfabetismo 2007", V_EDU17: "Analfabetismo 2017"}
    for v in (V_TEC07, V_TEC17, V_EDU07, V_EDU17):
        z = Z[v]
        I = moran_univariado(z, W)
        perms = permutacion_global(z, z, W, cruzado=False)
        p, mean, sd, zsc = pseudo_p(I, perms)
        filas.append(dict(tipo="Univariado", indicador=etiquetas_uni[v], I=I, E_I=EI,
                          media_perm=mean, sd_perm=sd, z_score=zsc, p_sim=p))
        print(f"    {etiquetas_uni[v]:<22} I={I:+.4f}  z={zsc:+.2f}  p={p:.3f}")

    # ---- 2) Moran cruzado con retardo temporal + 4) LISA ----
    print("\n[2] Índice cruzado de Moran con retardo temporal (pasado→presente):")
    for clave, etiqueta, vp, vq in PARES:
        z_pas, z_pres = Z[vp], Z[vq]
        I = moran_cruzado(z_pres, z_pas, W)
        perms = permutacion_global(z_pres, z_pas, W, cruzado=True)
        p, mean, sd, zsc = pseudo_p(I, perms)
        filas.append(dict(tipo="Cruzado", indicador=etiqueta, I=I, E_I=EI,
                          media_perm=mean, sd_perm=sd, z_score=zsc, p_sim=p))
        print(f"    {etiqueta:<38} I={I:+.4f}  z={zsc:+.2f}  p={p:.3f}")

        # LISA cruzado local
        I_loc, lag, p_loc, cuad = lisa_cruzado(z_pres, z_pas, W)
        sig = p_loc < 0.05
        lisa_df = pd.DataFrame({
            "departamento": d.DEPARTAMENTOS,
            "ubigeo": [d.UBIGEO[x] for x in d.DEPARTAMENTOS],
            "z_presente": z_pres, "rezago_pasado": lag,
            "I_local": I_loc, "p_local": p_loc,
            "cuadrante": cuad, "significativo_5pct": sig,
        })
        lisa_df.to_csv(os.path.join(RES, f"lisa_{clave}.csv"), index=False)

        fig_scatter(clave, etiqueta, z_pres, lag, I, p, cuad, sig)
        fig_lisa_mapa(clave, etiqueta, cuad, sig, gdf)

        hotspots = lisa_df.loc[sig, ["departamento", "cuadrante", "p_local"]]
        if len(hotspots):
            print("        clusters significativos (p<0.05): " +
                  ", ".join(f"{r.departamento}[{r.cuadrante}]" for r in hotspots.itertuples()))
        else:
            print("        sin clusters locales significativos al 5%.")

    # ---- Verificación cruzada con esda (si está instalado) ----
    verificar_con_esda(Z, W)

    # ---- Guardado de tablas ----
    glob = pd.DataFrame(filas)
    glob.to_csv(os.path.join(RES, "moran_global.csv"), index=False)

    resumen = glob.copy()
    resumen["significancia"] = np.where(resumen["p_sim"] < 0.01, "** (p<0.01)",
                               np.where(resumen["p_sim"] < 0.05, "* (p<0.05)",
                               np.where(resumen["p_sim"] < 0.10, ". (p<0.10)", "n.s.")))
    resumen["patron"] = np.where(resumen["I"] > resumen["E_I"], "Agrupamiento (positivo)",
                                 "Dispersión (negativo)")
    resumen[["tipo", "indicador", "I", "E_I", "z_score", "p_sim", "significancia", "patron"]] \
        .to_csv(os.path.join(RES, "resumen_capitulo_iv.csv"), index=False)

    print("\n" + "=" * 74)
    print("Resultados guardados en resultados/:")
    print("  - moran_global.csv, resumen_capitulo_iv.csv")
    print("  - lisa_<par>.csv (4)")
    print("  - figuras/scatter_<par>.png (4)" +
          ("  y  figuras/lisa_<par>.png (4)" if gdf is not None else ""))
    print("=" * 74)


def verificar_con_esda(Z, W):
    """Comprobación cruzada del Moran univariado contra esda (tolerancia 1e-6)."""
    try:
        from esda.moran import Moran
        from libpysal.weights import full2W
    except Exception:  # noqa: BLE001
        print("\n[check] esda/libpysal no disponibles: se omite la verificación cruzada.")
        return
    w = full2W(W)
    w.transform = "r"
    ok = True
    for v, z in Z.items():
        mine = moran_univariado(z, W)
        theirs = Moran(z, w, permutations=0).I
        if abs(mine - theirs) > 1e-6:
            ok = False
            print(f"[check] DISCREPANCIA en {v}: numpy={mine:.6f} vs esda={theirs:.6f}")
    print(f"\n[check] Verificación cruzada con esda: {'coincide ✔' if ok else 'REVISAR �’'}")


if __name__ == "__main__":
    main()
