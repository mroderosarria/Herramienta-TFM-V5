"""
conservation.py
----------------
Análisis de la conservación de las posiciones afectadas por las variantes,
comparando el residuo de referencia con el residuo presente en la misma
posición de cada secuencia homóloga (obtenida mediante BLAST).

Objetivo específico 6.

Estrategia: cada secuencia homóloga se alinea individualmente contra la
referencia (alineamiento global). Para cada variante de tipo "sustitución"
o "delección" (posiciones concretas de la referencia), se localiza la
columna del alineamiento que corresponde a esa posición de referencia y se
extrae el residuo del homólogo en esa columna. Con ello se calcula qué
porcentaje de homólogos conserva el aminoácido de referencia en esa
posición (score de conservación).
"""

from dataclasses import dataclass
from typing import List, Dict
from collections import Counter

from .variant_detection import _crear_alineador


@dataclass
class ConservacionPosicion:
    posicion_ref: str
    aa_referencia: str
    residuos_homologos: Dict[str, str]     # id_homologo -> residuo observado
    frecuencias: Dict[str, int]            # residuo -> nº de homólogos
    pct_conservacion: float                # % de homólogos con el mismo aa que la referencia
    residuo_consenso: str


def _mapa_posicion_referencia(referencia: str, homologo: str):
    """
    Alinea `referencia` con `homologo` y devuelve un diccionario que mapea
    posición 1-based de la referencia -> residuo del homólogo en esa columna
    (o '-' si en el homólogo hay una delección en esa posición).
    """
    aligner = _crear_alineador()
    alineamiento = aligner.align(referencia, homologo)[0]
    aln_ref, aln_hom = str(alineamiento[0]), str(alineamiento[1])

    mapa = {}
    pos_ref = 0
    for r, h in zip(aln_ref, aln_hom):
        if r != "-":
            pos_ref += 1
            mapa[pos_ref] = h
    return mapa


def analizar_conservacion(referencia: str, variantes, homologos: List) -> List[ConservacionPosicion]:
    """
    referencia: secuencia de referencia (str)
    variantes: lista de objetos Variante (variant_detection.Variante) de tipo
               'sustitucion' o 'delecion' con posición numérica simple.
    homologos: lista de objetos HitHomologo (blast_module.HitHomologo)

    Devuelve una lista de ConservacionPosicion, una por cada posición de
    referencia afectada por al menos una variante puntual.
    """
    # Precalcular el mapa posición->residuo para cada homólogo (una sola vez)
    mapas = {}
    for hom in homologos:
        try:
            mapas[hom.id_hit] = _mapa_posicion_referencia(referencia, hom.secuencia)
        except Exception:
            mapas[hom.id_hit] = {}

    resultados = []
    posiciones_vistas = set()

    for var in variantes:
        if var.tipo not in ("sustitucion", "delecion"):
            continue
        # Para delecciones de varios residuos, evaluar solo la primera posición
        try:
            pos = int(var.posicion_ref.split("-")[0])
        except ValueError:
            continue
        if pos in posiciones_vistas:
            continue
        posiciones_vistas.add(pos)

        aa_ref = referencia[pos - 1] if 0 < pos <= len(referencia) else "?"

        residuos_por_homologo = {}
        for hom in homologos:
            residuo = mapas.get(hom.id_hit, {}).get(pos, "?")
            residuos_por_homologo[hom.id_hit] = residuo

        conteo = Counter(residuos_por_homologo.values())
        total = sum(conteo.values()) or 1
        conservados = conteo.get(aa_ref, 0)
        pct = round(100 * conservados / total, 1)
        consenso = conteo.most_common(1)[0][0] if conteo else "?"

        resultados.append(ConservacionPosicion(
            posicion_ref=str(pos),
            aa_referencia=aa_ref,
            residuos_homologos=residuos_por_homologo,
            frecuencias=dict(conteo),
            pct_conservacion=pct,
            residuo_consenso=consenso,
        ))

    return resultados
