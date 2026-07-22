#!/usr/bin/env python3
"""
23_contraste.py
===============
Figura de VALIDACIÓN METODOLÓGICA: contrasta el Índice Cruzado de Moran del dengue
(fenómeno que se difunde) frente al de los indicadores socioeconómicos (que no se
difunden), demostrando que el método discrimina difusión de heterogeneidad.

Requiere haber ejecutado antes 21_moran_cruzado_dengue.py y 02_correlacion_cruzada_moran.py.

Lee:  resultados_dengue/moran_temporal_dengue.csv  (dengue)
      resultados/resumen_capitulo_iv.csv           (socioeconómico, filas 'Cruzado')
Salida: resultados/figuras/contraste_metodo.png
        (y también informe/figuras/ si ese directorio existe)
"""
from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import depuru as d  # noqa: E402  (reutilizamos exigir())

RAIZ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_RES = os.path.join(RAIZ, "resultados", "figuras")
FIG_INFORME = os.path.join(RAIZ, "informe", "figuras")


def destinos() -> list[str]:
    """Directorios donde guardar la figura: siempre resultados/figuras/ y, si existe
    (copia local del informe), también informe/figuras/."""
    salidas = [FIG_RES]
    if os.path.isdir(FIG_INFORME):
        salidas.append(FIG_INFORME)
    return salidas


def main():
    d.salida_utf8()
    d.exigir({
        os.path.join(RAIZ, "resultados_dengue", "moran_temporal_dengue.csv"):
            "python src/21_moran_cruzado_dengue.py",
        os.path.join(RAIZ, "resultados", "resumen_capitulo_iv.csv"):
            "python src/02_correlacion_cruzada_moran.py",
    })
    os.makedirs(FIG_RES, exist_ok=True)

    den = pd.read_csv(os.path.join(RAIZ, "resultados_dengue", "moran_temporal_dengue.csv"))
    i0 = float(den.loc[den.lag_semanas == 0, "I"].iloc[0])
    i4 = float(den.loc[den.lag_semanas == 4, "I"].iloc[0])

    soc = pd.read_csv(os.path.join(RAIZ, "resultados", "resumen_capitulo_iv.csv"))
    soc = soc[soc.tipo == "Cruzado"]

    etiquetas, valores, colores, sig = [], [], [], []
    # dengue (difusión)
    for lab, val in [("Dengue\n(contemporáneo)", i0), ("Dengue\n(retardo 4 sem)", i4)]:
        etiquetas.append(lab); valores.append(val); colores.append("#117a65"); sig.append("***")
    # socioeconómico (heterogeneidad)
    mapa = {"Tecnología'07 → Tecnología'17": "Brecha digital\n(2007→2017)",
            "Analfabetismo'07 → Analfabetismo'17": "Analfabetismo\n(2007→2017)",
            "Tecnología'07 → Analfabetismo'17": "Digital→Educación",
            "Analfabetismo'07 → Tecnología'17": "Educación→Digital"}
    for _, r in soc.iterrows():
        etiquetas.append(mapa.get(r["indicador"], r["indicador"]))
        valores.append(float(r["I"])); colores.append("#95a5a6"); sig.append("n.s.")

    ei = -1.0 / (113 - 1)

    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = range(len(valores))
    barras = ax.bar(x, valores, color=colores, edgecolor="black", linewidth=0.4)
    for xi, v, s in zip(x, valores, sig):
        ax.text(xi, v + (0.008 if v >= 0 else -0.008), s, ha="center",
                va="bottom" if v >= 0 else "top", fontsize=9, fontweight="bold")
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(ei, color="grey", ls="--", lw=1)
    ax.text(len(valores) - 0.4, ei, "  E[I] bajo H0", va="center", color="grey", fontsize=8)
    ax.set_ylim(top=max(valores) + 0.075)   # espacio para los rótulos de sección
    ax.axvline(1.5, color="grey", lw=0.8, ls=":")
    ax.text(0.167, 0.965, "DIFUSIÓN (dengue)", ha="center", va="top", transform=ax.transAxes,
            fontsize=10, color="#117a65", fontweight="bold")
    ax.text(0.667, 0.965, "HETEROGENEIDAD (socioeconómico)", ha="center",
            va="top", transform=ax.transAxes, fontsize=10, color="#7f8c8d", fontweight="bold")
    ax.set_xticks(list(x)); ax.set_xticklabels(etiquetas, fontsize=8)
    ax.set_ylabel("Índice Cruzado de Moran  I")
    ax.set_title("El método discrimina difusión de heterogeneidad\n"
                 "el Índice Cruzado es positivo y significativo solo cuando el fenómeno "
                 "realmente se propaga entre vecinos", fontsize=12)
    fig.tight_layout()
    for carpeta in destinos():
        ruta = os.path.join(carpeta, "contraste_metodo.png")
        fig.savefig(ruta, dpi=140, bbox_inches="tight")
        print(f"Guardado: {ruta}")
    plt.close(fig)
    print(f"  dengue I(0)={i0:.3f}, I(4)={i4:.3f}; socioeconómico cruzados: "
          f"{[round(v,3) for v in valores[2:]]}")


if __name__ == "__main__":
    main()
