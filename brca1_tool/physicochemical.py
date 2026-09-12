"""
physicochemical.py
-------------------
Caracterización de cambios de aminoácido mediante propiedades fisicoquímicas
básicas: carga, polaridad, hidrofobicidad (escala Kyte-Doolittle) y tamaño.

Objetivo específico 4.
"""

from dataclasses import dataclass

# Escala de hidrofobicidad de Kyte & Doolittle (1982)
HIDROFOBICIDAD_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Carga a pH fisiológico aproximado
CARGA = {
    "D": "negativa", "E": "negativa",
    "K": "positiva", "R": "positiva", "H": "positiva (débil)",
}
for _aa in "ACFGILMNPQSTVWY":
    CARGA[_aa] = "neutra"

# Polaridad
POLARIDAD = {
    "A": "apolar", "V": "apolar", "L": "apolar", "I": "apolar", "P": "apolar",
    "F": "apolar", "M": "apolar", "W": "apolar", "G": "apolar", "C": "apolar",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar", "Y": "polar",
    "D": "polar-cargada", "E": "polar-cargada",
    "K": "polar-cargada", "R": "polar-cargada", "H": "polar-cargada",
}

# Tamaño relativo de la cadena lateral (clasificación aproximada)
TAMANO = {
    "G": "pequeño", "A": "pequeño", "S": "pequeño", "C": "pequeño",
    "D": "pequeño", "N": "pequeño", "P": "pequeño", "T": "pequeño",
    "E": "mediano", "Q": "mediano", "H": "mediano", "V": "mediano",
    "M": "mediano", "I": "mediano", "L": "mediano", "K": "mediano",
    "R": "grande", "F": "grande", "Y": "grande", "W": "grande",
}


@dataclass
class PropiedadesAA:
    aminoacido: str
    carga: str
    polaridad: str
    hidrofobicidad: float
    tamano: str


def propiedades(aa: str) -> PropiedadesAA:
    aa = aa.upper()
    if aa not in HIDROFOBICIDAD_KD:
        # Aminoácido no estándar (p.ej. 'X') -> propiedades desconocidas
        return PropiedadesAA(aa, "desconocida", "desconocida", 0.0, "desconocido")
    return PropiedadesAA(
        aminoacido=aa,
        carga=CARGA[aa],
        polaridad=POLARIDAD[aa],
        hidrofobicidad=HIDROFOBICIDAD_KD[aa],
        tamano=TAMANO[aa],
    )


@dataclass
class CambioFisicoquimico:
    aa_referencia: str
    aa_variante: str
    prop_referencia: PropiedadesAA
    prop_variante: PropiedadesAA
    delta_hidrofobicidad: float
    cambia_carga: bool
    cambia_polaridad: bool
    cambia_tamano: bool
    clasificacion: str  # "conservativo" | "no conservativo"


def caracterizar_sustitucion(aa_ref: str, aa_var: str) -> CambioFisicoquimico:
    """
    Compara las propiedades fisicoquímicas entre el aminoácido de referencia
    y el aminoácido variante, y determina si el cambio es conservativo
    (propiedades similares) o no conservativo (propiedades muy distintas).
    """
    p_ref = propiedades(aa_ref)
    p_var = propiedades(aa_var)

    delta_h = round(p_var.hidrofobicidad - p_ref.hidrofobicidad, 2)
    cambia_carga = p_ref.carga != p_var.carga
    cambia_polaridad = p_ref.polaridad != p_var.polaridad
    cambia_tamano = p_ref.tamano != p_var.tamano

    # Regla simple de conservación: se considera "no conservativo" si cambia
    # la carga, o si cambia la polaridad y además el salto de hidrofobicidad
    # es grande (>|4|, escala KD va de -4.5 a 4.5).
    if cambia_carga or (cambia_polaridad and abs(delta_h) > 4):
        clasificacion = "no conservativo"
    elif abs(delta_h) > 6:
        clasificacion = "no conservativo"
    else:
        clasificacion = "conservativo"

    return CambioFisicoquimico(
        aa_referencia=aa_ref,
        aa_variante=aa_var,
        prop_referencia=p_ref,
        prop_variante=p_var,
        delta_hidrofobicidad=delta_h,
        cambia_carga=cambia_carga,
        cambia_polaridad=cambia_polaridad,
        cambia_tamano=cambia_tamano,
        clasificacion=clasificacion,
    )
