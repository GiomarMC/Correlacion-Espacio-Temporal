#!/usr/bin/env python3
"""
22_visualizacion_difusion.py
============================
Visualizaciones de la DIFUSIÓN ESPACIO-TEMPORAL del dengue en el Perú, para incluir
en el informe. Genera tres figuras (matplotlib, PNG):

  1. difusion_secuencia.png       — secuencia de mapas semanales de un brote real
                                     (temporada 2023) mostrando cómo se propaga a las
                                     provincias vecinas: el MAPA + la difusión.
  2. correlacion_foco_vecinas.png — cómo se correlaciona UNA provincia-foco (Piura)
                                     con sus vecinas a lo largo del tiempo: serie
                                     lead-lag + función de correlación cruzada (TLCC).
  3. scatter_moran_dengue.png     — dispersión de Moran contemporáneo (pendiente=I(0)).

Entradas: datos_dengue/{panel_dengue_semanal.csv, provincias_dengue.geojson,
matriz_W_dengue.csv}. Salidas: resultados_dengue/figuras/.
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
from aede import z_espacial_por_periodo  # noqa: E402


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

FOCO = "Piura"
TEMPORADA = 2023          # temporada del brote a visualizar (la mayor de la serie)
N_MAPAS = 6               # nº de mapas en la secuencia


def cargar():
    d.exigir({
        os.path.join(DATOS, "panel_dengue_semanal.csv"): "python src/20_preparar_dengue.py",
        os.path.join(DATOS, "matriz_W_dengue.csv"): "python src/20_preparar_dengue.py",
        os.path.join(DATOS, "provincias_dengue.geojson"): "python src/20_preparar_dengue.py",
    })
    panel = pd.read_csv(os.path.join(DATOS, "panel_dengue_semanal.csv"),
                        index_col=0, parse_dates=True)
    W = pd.read_csv(os.path.join(DATOS, "matriz_W_dengue.csv"), index_col=0)
    import geopandas as gpd
    g = gpd.read_file(os.path.join(DATOS, "provincias_dengue.geojson"))
    prov = list(panel.columns)
    g = g.set_index("provincia").loc[prov].reset_index()
    W = W.loc[prov, prov]
    return panel, W, g, prov


# ---------------------------------------------------------------------------
def fig_secuencia(panel, g, prov):
    """Secuencia de mapas semanales del brote 2023 (escala de color común)."""
    año = panel[panel.index.year == TEMPORADA]
    tot = año.sum(axis=1)
    pico = tot.idxmax()
    # 6 semanas desde ~14 semanas antes del pico hasta el pico (fase de crecimiento)
    semanas = pd.date_range(end=pico, periods=N_MAPAS, freq="3W")
    semanas = [panel.index[np.argmin(np.abs(panel.index - s))] for s in semanas]

    valores = np.log1p(panel.loc[semanas].to_numpy(float))
    vmax = np.percentile(valores[valores > 0], 99) if (valores > 0).any() else 1
    fig, axes = plt.subplots(2, 3, figsize=(13, 9))
    for ax, sem in zip(axes.ravel(), semanas):
        gg = g.copy()
        gg["v"] = np.log1p(panel.loc[sem].to_numpy(float))
        gg.plot(column="v", cmap="YlOrRd", vmin=0, vmax=vmax, edgecolor="0.8",
                linewidth=0.15, ax=ax)
        casos = int(panel.loc[sem].sum())
        ax.set_title(f"{sem.date()}  —  {casos:,} casos", fontsize=10)
        ax.axis("off")
    # barra de color compartida
    sm = plt.cm.ScalarMappable(cmap="YlOrRd", norm=plt.Normalize(vmin=0, vmax=vmax))
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("log(1 + casos de dengue) por provincia")
    fig.suptitle(f"Difusión espacio-temporal del dengue — brote {TEMPORADA} (Perú)\n"
                 "el brote aparece en el norte y se propaga a las provincias vecinas",
                 fontsize=13)
    fig.savefig(os.path.join(FIG, "difusion_secuencia.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    return semanas


# ---------------------------------------------------------------------------
def fig_foco_vecinas(panel, W, prov):
    """Serie lead-lag foco vs vecinas + función de correlación cruzada con retardo."""
    vecinas = [c for c in prov if W.loc[FOCO, c] > 0]
    foco = panel[FOCO]
    vec_media = panel[vecinas].mean(axis=1)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # (izq) temporada del brote: foco vs media de vecinas
    m = (panel.index.year >= TEMPORADA - 1) & (panel.index.year <= TEMPORADA)
    ax1.plot(panel.index[m], foco[m], color="#c0392b", lw=2, label=f"{FOCO} (foco)")
    ax1.plot(panel.index[m], vec_media[m], color="#2c6fbb", lw=2,
             label="media de provincias vecinas")
    ax1.set_title(f"Incidencia de dengue: {FOCO} vs sus vecinas ({TEMPORADA-1}–{TEMPORADA})")
    ax1.set_ylabel("casos por semana"); ax1.legend(); ax1.tick_params(axis="x", rotation=30)

    # (der) TLCC: corr( vecinas_{t-k} , foco_t ) sobre log-incidencia estandarizada
    zf = (np.log1p(foco) - np.log1p(foco).mean()) / np.log1p(foco).std()
    zv = (np.log1p(vec_media) - np.log1p(vec_media).mean()) / np.log1p(vec_media).std()
    lags = range(-8, 9)
    cc = [zf.corr(zv.shift(k)) for k in lags]   # k>0: vecinas anteceden al foco
    kbest = list(lags)[int(np.argmax(cc))]
    ax2.axvline(0, color="grey", lw=0.8)
    ax2.bar(list(lags), cc, color=["#117a65" if k == kbest else "#a9c7bd" for k in lags])
    ax2.set_title("Correlación cruzada con retardo (TLCC)\n"
                  f"máximo en k={kbest} sem  →  " +
                  ("las vecinas anteceden al foco" if kbest > 0 else
                   "el foco antecede a las vecinas" if kbest < 0 else "sincronía"))
    ax2.set_xlabel("retardo k (semanas)"); ax2.set_ylabel(f"corr( {FOCO}$_t$ , vecinas$_{{t-k}}$ )")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "correlacion_foco_vecinas.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    return vecinas, kbest


# ---------------------------------------------------------------------------
def fig_scatter_moran(panel, W, prov):
    """Dispersión de Moran contemporáneo agrupada sobre todas las semanas."""
    Z, valida = z_espacial_por_periodo(panel.to_numpy(float))
    Z = Z[valida]
    Wm = W.to_numpy(float)
    lag = Z @ Wm.T
    x, y = Z.ravel(), lag.ravel()
    I0 = float(np.mean(np.sum(Z * lag, axis=1) / Wm.shape[0]))

    fig, ax = plt.subplots(figsize=(6.5, 6))
    hb = ax.hexbin(x, y, gridsize=60, cmap="Blues", bins="log", mincnt=1)
    xs = np.linspace(np.percentile(x, 0.5), np.percentile(x, 99.5), 50)
    ax.plot(xs, I0 * xs, color="#c0392b", lw=2, label=f"pendiente = I(0) = {I0:.3f}")
    ax.axhline(0, color="grey", lw=0.7); ax.axvline(0, color="grey", lw=0.7)
    ax.set_xlabel("incidencia estandarizada de la provincia  z (semana t)")
    ax.set_ylabel("incidencia rezagada de las vecinas  W·z (semana t)")
    ax.set_title("Diagrama de dispersión de Moran (contemporáneo) — Dengue Perú")
    fig.colorbar(hb, ax=ax, label="nº de observaciones (log)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "scatter_moran_dengue.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)
    return I0


# ---------------------------------------------------------------------------
def fig_serie_nacional(panel):
    """Serie temporal semanal de casos de dengue a nivel nacional, con los
    principales picos epidémicos anotados."""
    tot = panel.sum(axis=1)
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.fill_between(tot.index, tot.values, color="#117a65", alpha=0.35)
    ax.plot(tot.index, tot.values, color="#0e6252", lw=1)
    # anotar los picos de las mayores temporadas
    for anio, texto in [(2017, "2017\n(El Niño costero)"), (2023, "2023\n(máximo histórico)")]:
        s = tot[tot.index.year == anio]
        if len(s):
            pk = s.idxmax()
            ax.annotate(texto, (pk, s.max()), xytext=(0, 18), textcoords="offset points",
                        ha="center", fontsize=9, fontweight="bold",
                        arrowprops=dict(arrowstyle="->", color="grey"))
    ax.set_xlabel("Semana"); ax.set_ylabel("Casos de dengue por semana (Perú)")
    ax.set_title("Incidencia semanal de dengue en el Perú, 2000–2023\n"
                 "las epidemias llegan en olas, con máximos en 2017 y 2023", fontsize=12)
    ax.margins(x=0.01)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "serie_nacional.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_onset_2023(panel, g, prov):
    """Mapa del 'frente de difusión': cada provincia con brote relevante en 2023
    se colorea por la SEMANA en que alcanzó su pico; las de escaso brote, en gris.
    Muestra qué provincias se encienden primero y cómo avanza la ola."""
    import matplotlib.dates as mdates

    año = panel[panel.index.year == 2023]
    total_prov = año.sum(axis=0)
    umbral = 50  # casos mínimos en 2023 para considerar 'brote relevante'
    semana_pico = año.idxmax(axis=0)          # semana del máximo por provincia
    con_brote = total_prov >= umbral

    gg = g.copy()
    gg["pico"] = [mdates.date2num(semana_pico[p]) if con_brote[p] else np.nan for p in prov]

    fig, ax = plt.subplots(figsize=(7.5, 8.5))
    # provincias sin brote relevante: gris
    gg[gg["pico"].isna()].plot(ax=ax, color="#e5e7e9", edgecolor="white", linewidth=0.2)
    # provincias con brote: coloreadas por semana de pico
    sub = gg[~gg["pico"].isna()].reset_index(drop=True)
    sub.plot(column="pico", cmap="viridis", ax=ax, edgecolor="white", linewidth=0.2)
    # rotula las provincias con mayor brote en 2023 (las más reconocibles)
    sub_total = total_prov.loc[sub["provincia"]].to_numpy()
    etiquetar_provincias(ax, sub, sub_total, n=6)
    # barra de color con fechas
    vmin, vmax = sub["pico"].min(), sub["pico"].max()
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    cb = fig.colorbar(sm, ax=ax, shrink=0.6)
    loc = mdates.AutoDateLocator(); cb.ax.yaxis.set_major_locator(loc)
    cb.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
    cb.set_label("semana del pico de dengue en 2023")
    ax.set_title("Frente de difusión del brote de dengue de 2023\n"
                 "cada provincia coloreada por la semana de su pico "
                 "(morado = temprano · amarillo = tardío); gris = sin brote relevante", fontsize=11)
    ax.axis("off")
    fig.savefig(os.path.join(FIG, "mapa_onset_2023.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    d.salida_utf8()
    os.makedirs(FIG, exist_ok=True)
    panel, W, g, prov = cargar()
    print("Generando visualizaciones de difusión espacio-temporal…")
    sem = fig_secuencia(panel, g, prov)
    print(f"  [1] difusion_secuencia.png  (semanas: {[s.date().isoformat() for s in sem]})")
    vecinas, kbest = fig_foco_vecinas(panel, W, prov)
    print(f"  [2] correlacion_foco_vecinas.png  (foco={FOCO}, vecinas={vecinas}, TLCC máx k={kbest})")
    I0 = fig_scatter_moran(panel, W, prov)
    print(f"  [3] scatter_moran_dengue.png  (I(0)={I0:.3f})")
    fig_serie_nacional(panel)
    print("  [4] serie_nacional.png")
    fig_onset_2023(panel, g, prov)
    print("  [5] mapa_onset_2023.png")
    print(f"Figuras guardadas en {FIG}")


if __name__ == "__main__":
    main()
