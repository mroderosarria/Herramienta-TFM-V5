"""
fasta_io.py
-----------
Lectura y validación de secuencias proteicas en formato FASTA.

Objetivo específico 1: módulo de lectura y validación de secuencias FASTA,
permitiendo introducir una secuencia de referencia de BRCA1 y una o varias
secuencias variantes.
"""

from dataclasses import dataclass
from typing import List
from Bio import SeqIO

AMINOACIDOS_VALIDOS = set("ACDEFGHIKLMNPQRSTVWY")  # 20 aa estándar
CARACTERES_TOLERADOS = AMINOACIDOS_VALIDOS | {"X", "*"}  # X: desconocido, *: stop


@dataclass
class SecuenciaProteica:
    id: str
    descripcion: str
    secuencia: str

    def __len__(self):
        return len(self.secuencia)


class FastaValidationError(Exception):
    """Error lanzado cuando una secuencia FASTA no es válida."""
    pass


def _validar_secuencia(seq: str, origen: str) -> str:
    seq = seq.strip().upper().replace("\n", "")
    if not seq:
        raise FastaValidationError(f"La secuencia '{origen}' está vacía.")

    caracteres_invalidos = set(seq) - CARACTERES_TOLERADOS
    if caracteres_invalidos:
        raise FastaValidationError(
            f"La secuencia '{origen}' contiene caracteres no válidos para una "
            f"proteína: {sorted(caracteres_invalidos)}. Se permiten los 20 "
            f"aminoácidos estándar, 'X' (residuo indeterminado) y '*' (stop)."
        )
    return seq


def leer_fasta(ruta: str) -> List[SecuenciaProteica]:
    """
    Lee un archivo FASTA (puede contener una o varias secuencias) y devuelve
    una lista de objetos SecuenciaProteica ya validados.
    """
    registros = list(SeqIO.parse(ruta, "fasta"))
    if not registros:
        raise FastaValidationError(
            f"No se ha encontrado ninguna secuencia en el archivo '{ruta}'. "
            "Comprueba que el archivo está en formato FASTA (encabezado '>')."
        )

    secuencias = []
    for r in registros:
        seq_validada = _validar_secuencia(str(r.seq), origen=r.id)
        secuencias.append(
            SecuenciaProteica(id=r.id, descripcion=r.description, secuencia=seq_validada)
        )
    return secuencias


def leer_referencia(ruta: str) -> SecuenciaProteica:
    """
    Lee un archivo FASTA que debe contener EXACTAMENTE una secuencia de
    referencia (p. ej. BRCA1 canónica).
    """
    secuencias = leer_fasta(ruta)
    if len(secuencias) != 1:
        raise FastaValidationError(
            f"El archivo de referencia '{ruta}' debe contener una única secuencia, "
            f"pero contiene {len(secuencias)}."
        )
    return secuencias[0]


def leer_variantes(rutas: List[str]) -> List[SecuenciaProteica]:
    """
    Lee una lista de rutas de archivos FASTA de variantes. Cada archivo puede
    contener una o varias secuencias variantes.
    """
    variantes = []
    for ruta in rutas:
        variantes.extend(leer_fasta(ruta))
    if not variantes:
        raise FastaValidationError("No se ha proporcionado ninguna secuencia variante.")
    return variantes
