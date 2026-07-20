"""
aede.py — Núcleo estadístico del Análisis Exploratorio de Datos Espaciales con
retardo temporal. Funciones compartidas por los scripts 02 (análisis principal)
y 03 (robustez ante distintas definiciones de vecindad W).
"""
from __future__ import annotations

import numpy as np

SEMILLA = 42
N_PERM = 999


def zscore(x) -> np.ndarray:
    """Estandariza a media 0 y desviación 1 (ddof=0 -> zᵀz = n)."""
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / x.std(ddof=0)


def fila_estandariza(W: np.ndarray) -> np.ndarray:
    """Estandarización por filas: Σ_j w_ij = 1 (filas nulas se dejan en 0)."""
    W = np.asarray(W, dtype=float).copy()
    s = W.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return W / s


def z_espacial_por_periodo(M: np.ndarray, log1p: bool = True):
    """Estandariza cada FILA (periodo: día/semana) ENTRE unidades espaciales.
    Con log1p=True aplica log(1+x) antes (útil para conteos como casos de dengue).
    Devuelve (Z, valida) donde 'valida' marca los periodos con varianza espacial>0.
    Es la estandarización correcta del Moran (desviaciones respecto de la media de
    cada mapa) y elimina el componente común del periodo (p.ej. estacionalidad)."""
    L = np.log1p(M) if log1p else np.asarray(M, dtype=float)
    mu = L.mean(axis=1, keepdims=True)
    sd = L.std(axis=1, ddof=0, keepdims=True)
    valida = (sd[:, 0] > 1e-9)
    sd[sd < 1e-9] = 1.0
    return (L - mu) / sd, valida


def moran_univariado(z: np.ndarray, W: np.ndarray) -> float:
    """I = (zᵀ W z)/(zᵀ z)."""
    return float(z @ (W @ z) / (z @ z))


def moran_cruzado(z_pres: np.ndarray, z_pas: np.ndarray, W: np.ndarray) -> float:
    """Índice cruzado con retardo: I = z_presente · (W · z_pasado) / n."""
    return float(z_pres @ (W @ z_pas) / len(z_pres))


def pseudo_p(obs: float, perms: np.ndarray):
    """Pseudo p-valor (estilo esda, plegado a la cola de la observación) +
    media, desviación y z-score de la distribución nula."""
    M = len(perms)
    mean, sd = float(perms.mean()), float(perms.std(ddof=0))
    larger = int(np.sum(perms >= obs))
    if (M - larger) < larger:
        larger = M - larger
    p = (larger + 1.0) / (M + 1.0)
    z = (obs - mean) / sd if sd > 0 else np.nan
    return p, mean, sd, z


def permutacion_global(z_pres, z_pas, W, cruzado=True, seed=SEMILLA, M=N_PERM):
    """Distribución nula por permutación del vector que se rezaga (el 'pasado')."""
    rng = np.random.default_rng(seed)
    perms = np.empty(M)
    for k in range(M):
        zp = rng.permutation(z_pas)
        perms[k] = moran_cruzado(z_pres, zp, W) if cruzado else moran_univariado(zp, W)
    return perms
