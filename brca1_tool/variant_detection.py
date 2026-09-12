"""
variant_detection.py
---------------------
Alineamiento de la secuencia de referencia frente a cada secuencia variante
y detección/clasificación de las diferencias (sustitución, inserción,
deleción), con su posición respecto a la referencia.

Objetivos específicos 2 y 3.
"""

from dataclasses import dataclass, field
from typing import List
from Bio import Align
from Bio.Align import substitution_matrices


def _crear_alineador() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.end_insertion_score = 0.0
    aligner.end_deletion_score = 0.0
    return aligner


@dataclass
class Variante:
    tipo: str                 # "sustitucion" | "insercion" | "delecion"
    posicion_ref: str         # posición (1-based) en la referencia, o rango/"entre X-Y"
    aa_referencia: str        # '-' si es inserción
    aa_variante: str          # '-' si es delección
    longitud: int = 1


@dataclass
class ResultadoComparacion:
    id_variante: str
    alineamiento_ref: str
    alineamiento_var: str
    identidad_pct: float
    variantes: List[Variante] = field(default_factory=list)


def alinear(referencia: str, variante: str):
    aligner = _crear_alineador()
    alineamientos = aligner.align(referencia, variante)
    return alineamientos[0]  # mejor alineamiento


def _pct_identidad(aln_ref: str, aln_var: str) -> float:
    columnas = [(a, b) for a, b in zip(aln_ref, aln_var) if a != "-" and b != "-"]
    if not columnas:
        return 0.0
    iguales = sum(1 for a, b in columnas if a == b)
    return round(100 * iguales / len(columnas), 2)


def detectar_variantes(referencia: str, variante_seq: str, id_variante: str) -> ResultadoComparacion:
    """
    Alinea `referencia` con `variante_seq` y extrae la lista de eventos
    (sustitución, inserción, delección) con su posición respecto a la
    referencia (numeración 1-based, contando solo residuos de la referencia).
    """
    alineamiento = alinear(referencia, variante_seq)
    aln_ref, aln_var = str(alineamiento[0]), str(alineamiento[1])

    eventos: List[Variante] = []
    pos_ref = 0  # posición 1-based en la secuencia de referencia (sin gaps)

    i = 0
    n = len(aln_ref)
    while i < n:
        r, v = aln_ref[i], aln_var[i]

        if r == "-":
            # Inserción: acumular todas las columnas de inserción consecutivas.
            # pos_ref todavía apunta al último residuo de referencia ya emitido.
            ins_seq = ""
            while i < n and aln_ref[i] == "-":
                ins_seq += aln_var[i]
                i += 1
            eventos.append(Variante(
                tipo="insercion",
                posicion_ref=f"entre {pos_ref} y {pos_ref + 1}",
                aa_referencia="-",
                aa_variante=ins_seq,
                longitud=len(ins_seq),
            ))
            continue

        elif v == "-":
            # Delección: acumular columnas de delección consecutivas
            del_seq = ""
            pos_inicio = pos_ref + 1
            while i < n and aln_var[i] == "-" and aln_ref[i] != "-":
                pos_ref += 1
                del_seq += aln_ref[i]
                i += 1
            pos_fin = pos_ref
            rango = f"{pos_inicio}" if pos_inicio == pos_fin else f"{pos_inicio}-{pos_fin}"
            eventos.append(Variante(
                tipo="delecion",
                posicion_ref=rango,
                aa_referencia=del_seq,
                aa_variante="-",
                longitud=len(del_seq),
            ))
            continue

        else:
            pos_ref += 1
            if r != v:
                eventos.append(Variante(
                    tipo="sustitucion",
                    posicion_ref=str(pos_ref),
                    aa_referencia=r,
                    aa_variante=v,
                    longitud=1,
                ))

        i += 1

    return ResultadoComparacion(
        id_variante=id_variante,
        alineamiento_ref=aln_ref,
        alineamiento_var=aln_var,
        identidad_pct=_pct_identidad(aln_ref, aln_var),
        variantes=eventos,
    )


def comparar_todas(referencia: str, variantes: list) -> List[ResultadoComparacion]:
    """
    `variantes` es una lista de objetos SecuenciaProteica (ver fasta_io.py).
    """
    return [detectar_variantes(referencia, v.secuencia, v.id) for v in variantes]
