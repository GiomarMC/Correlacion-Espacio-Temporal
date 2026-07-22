#!/usr/bin/env python3
"""
21_moran_cruzado_dengue.py
==========================
CORRELACIÓN CRUZADA ESPACIAL CON RETARDO TEMPORAL aplicada al DENGUE en el Perú
(provincias × semanas, 2000–2023). Demuestra la difusión espacial del dengue entre
provincias vecinas con desfase temporal.

Para cada retardo k (en semanas):

        I(k) = promedio_t [ z_t · ( W · z_{t-k} ) / N ]

donde z_t estandariza la incidencia ENTRE provincias dentro de cada semana t
(log(1+casos) y luego z espacial). Esta estandarización por semana es el Moran
correcto (deviaciones respecto de la media de cada mapa) y, además, elimina el
componente estacional común (todas las provincias suben en verano), aislando el
patrón espacial de dónde están los brotes respecto del promedio nacional de esa
semana. I(k) mide si un brote en las vecinas hace k semanas anticipa un brote aquí.

Significancia: permutación espacial independiente por semana (999), nulo en E[I].

Salidas (en resultados_dengue/):
  - moran_temporal_dengue.csv
  - figuras/I_vs_lag_dengue.png      (curva I(k) — figura principal)
  - figuras/lisa_dengue.png          (mapa de clusters LISA promedio)
  - figuras/incidencia_media.png     (coropleta de incidencia media)
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import depuru as d  # noqa: E402  (reutilizamos exigir())


def etiquetar_provincias(ax, g, valores, n=8):
    """Rotula en el centroide las n provincias de mayor 'valores', con halo blanco."""
    orden = np.argsort(valores)[::-1][:n]
    cent = g.geometry.representative_point()
    for i in orden:
        ax.annotate(g["provincia"].iloc[i], (cent.x.iloc[i], cent.y.iloc[i]),
                    fontsize=7, ha="center", va="center", color="black",
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATOS = os.path.join(RAIZ, "datos_dengue")
RES = os.path.join(RAIZ, "resultados_dengue")
FIG = os.path.join(RES, "figuras")

SEMILLA = 42
N_PERM = 999
LAGS = list(range(0, 27))     # retardos de 0 a 26 semanas (medio año)
EI = None                     # se fija en main según N


def z_espacial_por_semana(M: np.ndarray) -> np.ndarray:
    """Estandariza cada FILA (semana) entre provincias, sobre log(1+casos).
    Devuelve la matriz z y un booleano de semanas válidas (varianza > 0)."""
    L = np.log1p(M)
    mu = L.mean(axis=1, keepdims=True)
    sd = L.std(axis=1, ddof=0, keepdims=True)
    valida = (sd[:, 0] > 1e-9)
    sd[sd < 1e-9] = 1.0
    return (L - mu) / sd, valida


def pares_validos(fechas, valida, k):
    """Índices (i,j) de pares (semana t, semana t-k) ambas válidas y presentes."""
    idx = {f: n for n, f in enumerate(fechas)}
    ii, jj = [], []
    for t in fechas:
        t0 = t - pd.Timedelta(weeks=k)
        i, j = idx.get(t), idx.get(t0)
        if i is not None and j is not None and valida[i] and valida[j]:
            ii.append(i); jj.append(j)
    return np.array(ii), np.array(jj)


def I_de_lag(Z, W, ii, jj):
    n = W.shape[0]
    Zi, WZj = Z[ii], (Z[jj] @ W.T)
    return float(np.mean(np.sum(Zi * WZj, axis=1) / n))


def permutar(Z, W, ii, jj, obs, seed, M=N_PERM):
    n = W.shape[0]
    Zi, Zj, Wt = Z[ii], Z[jj], W.T
    P = len(ii)
    rng = np.random.default_rng(seed)
    perms = np.empty(M)
    for r in range(M):
        orden = np.argsort(rng.random((P, n)), axis=1)
        WZ = np.take_along_axis(Zj, orden, axis=1) @ Wt
        perms[r] = np.mean(np.sum(Zi * WZ, axis=1) / n)
    mean, sd = perms.mean(), perms.std(ddof=0)
    larger = int(np.sum(perms >= obs))
    if (M - larger) < larger:
        larger = M - larger
    return (larger + 1) / (M + 1), (obs - mean) / sd if sd > 0 else np.nan


def main():
    global EI
    d.salida_utf8()
    d.exigir({
        os.path.join(DATOS, "panel_dengue_semanal.csv"): "python src/20_preparar_dengue.py",
        os.path.join(DATOS, "matriz_W_dengue.csv"): "python src/20_preparar_dengue.py",
    })
    os.makedirs(FIG, exist_ok=True)
    panel = pd.read_csv(os.path.join(DATOS, "panel_dengue_semanal.csv"),
                        index_col=0, parse_dates=True)
    W = pd.read_csv(os.path.join(DATOS, "matriz_W_dengue.csv"), index_col=0)
    prov = list(panel.columns)
    W = W.loc[prov, prov].to_numpy(float)
    n = len(prov)
    EI = -1.0 / (n - 1)

    fechas = panel.index
    Z, valida = z_espacial_por_semana(panel.to_numpy(float))

    print("=" * 74)
    print(" CORRELACIÓN CRUZADA ESPACIAL CON RETARDO TEMPORAL — DENGUE PERÚ ")
    print(f" {n} provincias | {len(fechas)} semanas ({fechas.min().date()}..{fechas.max().date()})")
    print(f" {int(valida.sum())} semanas con brotes (varianza espacial>0) | E[I]={EI:.4f}")
    print("=" * 74)
    print(f"{'lag(sem)':>8} {'I(k)':>9} {'pares':>7} {'z':>8} {'p_sim':>7}")

    filas = []
    for k in LAGS:
        ii, jj = pares_validos(fechas, valida, k)
        I = I_de_lag(Z, W, ii, jj)
        p, z = permutar(Z, W, ii, jj, I, seed=SEMILLA + k)
        filas.append(dict(lag_semanas=k, I=I, n_pares=len(ii), z=z, p_sim=p))
        marca = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"{k:>8} {I:>9.4f} {len(ii):>7} {z:>8.2f} {p:>7.3f} {marca}")

    tab = pd.DataFrame(filas)
    tab.to_csv(os.path.join(RES, "moran_temporal_dengue.csv"), index=False)

    # ---- Figura principal: I(k) vs retardo ----
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(tab.lag_semanas, tab.I, "o-", color="#117a65", lw=2, zorder=3)
    sig = tab[tab.p_sim < 0.05]
    ax.scatter(sig.lag_semanas, sig.I, s=90, facecolor="none", edgecolor="black",
               linewidth=1.3, label="significativo (p<0.05)", zorder=4)
    ax.axhline(0, color="grey", lw=0.8)
    ax.axhline(EI, color="grey", ls="--", lw=1, label="E[I] bajo H0")
    ax.set_xlabel("Retardo temporal k (semanas)")
    ax.set_ylabel("Índice Cruzado de Moran  I(k)")
    ax.set_title("Correlación cruzada espacial con retardo temporal — Dengue, provincias del Perú\n"
                 "(difusión: la incidencia presente se asocia con la de las vecinas semanas atrás)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "I_vs_lag_dengue.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)

    # ---- Mapas (incidencia media y LISA local promedio a lag 1) ----
    try:
        import geopandas as gpd
        g = gpd.read_file(os.path.join(DATOS, "provincias_dengue.geojson"))
        g = g.set_index("provincia").loc[prov].reset_index()

        inc = np.log1p(panel.to_numpy(float)).mean(axis=0)
        g["inc_media"] = inc
        fig, ax = plt.subplots(figsize=(7.5, 8.5))
        g.plot(column="inc_media", cmap="YlOrRd", legend=True, edgecolor="white",
               linewidth=0.2, ax=ax,
               legend_kwds={"label": "incidencia media  log(1+casos)", "shrink": 0.6})
        etiquetar_provincias(ax, g, inc, n=8)
        ax.set_title("Incidencia media de dengue por provincia, 2000–2023\n"
                     "(más oscuro = mayor carga histórica de dengue)", fontsize=12)
        ax.axis("off")
        fig.savefig(os.path.join(FIG, "incidencia_media.png"), dpi=140, bbox_inches="tight")
        plt.close(fig)

        # LISA local promedio a lag 1 semana: I_i = prom_t z_i(t)·(W z(t-1))_i
        ii, jj = pares_validos(fechas, valida, 1)
        Ii = (Z[ii] * (Z[jj] @ W.T)).mean(axis=0)
        g["lisa"] = Ii
        fig, ax = plt.subplots(figsize=(7.5, 8.5))
        g.plot(column="lisa", cmap="RdBu_r", legend=True, edgecolor="white",
               linewidth=0.2, ax=ax, vmin=-np.abs(Ii).max(), vmax=np.abs(Ii).max(),
               legend_kwds={"label": "contribución local a I (retardo 1 semana)", "shrink": 0.6})
        etiquetar_provincias(ax, g, Ii, n=7)
        ax.set_title("Focos de difusión del dengue — LISA cruzado local (retardo de 1 semana)\n"
                     "en rojo, provincias con alta incidencia rodeadas de vecinas también altas "
                     "la semana previa", fontsize=11)
        ax.axis("off")
        fig.savefig(os.path.join(FIG, "lisa_dengue.png"), dpi=140, bbox_inches="tight")
        plt.close(fig)
        top = pd.Series(Ii, index=prov).sort_values(ascending=False).head(8)
        print("-" * 74)
        print("Provincias-foco de difusión (mayor contribución local a lag 1 semana):")
        print(top.round(3).to_string())
    except Exception as e:  # noqa: BLE001
        print(f"[mapas] omitidos ({type(e).__name__}: {e})")

    print("-" * 74)
    I0 = tab.loc[tab.lag_semanas == 0, "I"].iloc[0]
    print(f"Moran contemporáneo I(0) = {I0:.3f}. Curva I(k) guardada en resultados_dengue/.")
    print("=" * 74)


if __name__ == "__main__":
    main()
