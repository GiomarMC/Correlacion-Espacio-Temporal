"""
depuru.py — Utilidades compartidas del TIF (Correlación Cruzada Espacial con Retardo Temporal).

Contiene la definición canónica de los 24 departamentos del Perú usados como unidades
espaciales del estudio (Ubigeo 01–25, excluyendo Callao=07; Lima como registro único),
funciones de normalización de nombres y la lista de adyacencia de fronteras verificada a
mano, que se usa como validación y como respaldo (fallback) si falla la descarga de la
geometría oficial.
"""
from __future__ import annotations

import unicodedata

# ---------------------------------------------------------------------------
# 24 departamentos en orden canónico (por Ubigeo). Nombres tal como aparecen en
# datos_consolidados/analisis_tecnologia_educacion.csv (con "Madre De Dios").
# ---------------------------------------------------------------------------
DEPARTAMENTOS = [
    "Amazonas", "Ancash", "Apurimac", "Arequipa", "Ayacucho", "Cajamarca",
    "Cusco", "Huancavelica", "Huanuco", "Ica", "Junin", "La Libertad",
    "Lambayeque", "Lima", "Loreto", "Madre De Dios", "Moquegua", "Pasco",
    "Piura", "Puno", "San Martin", "Tacna", "Tumbes", "Ucayali",
]

# Ubigeo departamental oficial del INEI (Callao=07 excluido intencionalmente).
UBIGEO = {
    "Amazonas": "01", "Ancash": "02", "Apurimac": "03", "Arequipa": "04",
    "Ayacucho": "05", "Cajamarca": "06", "Cusco": "08", "Huancavelica": "09",
    "Huanuco": "10", "Ica": "11", "Junin": "12", "La Libertad": "13",
    "Lambayeque": "14", "Lima": "15", "Loreto": "16", "Madre De Dios": "17",
    "Moquegua": "18", "Pasco": "19", "Piura": "20", "Puno": "21",
    "San Martin": "22", "Tacna": "23", "Tumbes": "24", "Ucayali": "25",
}

N = len(DEPARTAMENTOS)  # 24


def normaliza(nombre: str) -> str:
    """Clave canónica de un nombre de departamento: sin tildes/ñ, minúsculas, sin
    espacios extra. Unifica variantes de Lima (Lima Province/Provincias -> lima)."""
    if nombre is None:
        return ""
    s = unicodedata.normalize("NFKD", str(nombre))
    s = "".join(c for c in s if not unicodedata.combining(c))  # quita tildes
    s = s.lower().strip()
    s = " ".join(s.split())
    # Equivalencias de escritura frecuentes
    s = s.replace("provincia constitucional del ", "")
    # Fusión de Lima Metropolitana / Lima Province(s) en un único "lima"
    if s.startswith("lima"):
        s = "lima"
    # Quita TODO espacio: unifica "La Libertad" con "LaLibertad" (nomenclatura GADM)
    s = s.replace(" ", "")
    return s


# Mapa: clave normalizada -> nombre canónico del estudio.
_CANON = {normaliza(d): d for d in DEPARTAMENTOS}


def a_canonico(nombre: str) -> str | None:
    """Devuelve el nombre canónico del estudio para cualquier variante, o None si
    no pertenece a los 24 departamentos (p.ej. 'Callao')."""
    return _CANON.get(normaliza(nombre))


# ---------------------------------------------------------------------------
# Lista de adyacencia de fronteras (contigüidad de primer orden), verificada a
# mano a partir de los límites políticos departamentales del Perú. Es simétrica.
# Se usa para: (a) validar la matriz derivada de la geometría oficial, y
# (b) respaldo si la descarga de geometría falla.
# ---------------------------------------------------------------------------
ADYACENCIA_VERIFICADA: dict[str, list[str]] = {
    "Amazonas": ["Cajamarca", "La Libertad", "San Martin", "Loreto"],
    "Ancash": ["La Libertad", "Huanuco", "Pasco", "Lima"],
    "Apurimac": ["Cusco", "Ayacucho", "Arequipa"],
    "Arequipa": ["Ica", "Ayacucho", "Apurimac", "Cusco", "Puno", "Moquegua"],
    "Ayacucho": ["Junin", "Huancavelica", "Ica", "Arequipa", "Apurimac", "Cusco"],
    "Cajamarca": ["Piura", "Lambayeque", "La Libertad", "Amazonas"],
    "Cusco": ["Junin", "Ucayali", "Madre De Dios", "Puno", "Arequipa", "Apurimac", "Ayacucho"],
    "Huancavelica": ["Junin", "Lima", "Ica", "Ayacucho"],
    "Huanuco": ["La Libertad", "San Martin", "Ucayali", "Pasco", "Ancash", "Lima"],
    "Ica": ["Lima", "Huancavelica", "Ayacucho", "Arequipa"],
    "Junin": ["Pasco", "Ucayali", "Cusco", "Ayacucho", "Huancavelica", "Lima"],
    "La Libertad": ["Lambayeque", "Cajamarca", "Amazonas", "San Martin", "Huanuco", "Ancash"],
    "Lambayeque": ["Piura", "Cajamarca", "La Libertad"],
    "Lima": ["Ancash", "Huanuco", "Pasco", "Junin", "Huancavelica", "Ica"],
    "Loreto": ["Amazonas", "San Martin", "Ucayali"],
    "Madre De Dios": ["Ucayali", "Cusco", "Puno"],
    "Moquegua": ["Arequipa", "Puno", "Tacna"],
    "Pasco": ["Huanuco", "Ucayali", "Junin", "Lima", "Ancash"],
    "Piura": ["Tumbes", "Cajamarca", "Lambayeque"],
    "Puno": ["Madre De Dios", "Cusco", "Arequipa", "Moquegua", "Tacna"],
    "San Martin": ["Loreto", "Amazonas", "La Libertad", "Huanuco", "Ucayali"],
    "Tacna": ["Moquegua", "Puno"],
    "Tumbes": ["Piura"],
    "Ucayali": ["Loreto", "Huanuco", "Pasco", "Junin", "Cusco", "Madre De Dios", "San Martin"],
}


def adyacencia_a_matriz(ady: dict[str, list[str]]):
    """Convierte una lista de adyacencia en matriz binaria 24x24 (numpy), en el
    orden canónico de DEPARTAMENTOS. Verifica simetría."""
    import numpy as np

    idx = {d: i for i, d in enumerate(DEPARTAMENTOS)}
    A = np.zeros((N, N), dtype=int)
    for dep, vecinos in ady.items():
        i = idx[dep]
        for v in vecinos:
            A[i, idx[v]] = 1
    if not np.array_equal(A, A.T):
        raise ValueError("La lista de adyacencia verificada NO es simétrica.")
    return A
