"""
blast_module.py
----------------
Búsqueda de secuencias proteicas homólogas a BRCA1 mediante BLASTp.

Objetivo específico 5.

Se ofrecen dos modos, seleccionables por el usuario:

1. "local" (recomendado quando sea posible): ejecuta el binario `blastp`
   de BLAST+ contra una base de datos ya descargada en la propia máquina
   (p. ej. con `update_blastdb.pl swissprot`). No depende de la red ni de
   la carga del servidor de NCBI: tarda típicamente menos de un segundo
   por consulta. Requiere tener BLAST+ instalado (`blastp` en el PATH) y
   la base de datos ya descargada localmente.

2. "remote": lanza la búsqueda contra el servidor público de NCBI mediante
   Bio.Blast.NCBIWWW.qblast. No requiere instalar nada, pero depende de la
   conexión a internet y de la carga/límites de NCBI, y puede tardar desde
   segundos hasta muchos minutos (NCBI limita los recursos de cómputo por
   dirección IP, por lo que redes compartidas -p. ej. universidades o
   proxies de salida compartidos en la nube- pueden sufrir tiempos de
   espera mucho más largos).
"""

import shutil
import subprocess
import tempfile
import os
from dataclasses import dataclass
from typing import List


@dataclass
class HitHomologo:
    id_hit: str
    descripcion: str
    identidad_pct: float
    longitud_alineamiento: int
    score: float
    evalue: float
    secuencia: str  # secuencia del homólogo alineada con la query (para análisis de conservación)


def blastp_disponible() -> bool:
    """Comprueba si el ejecutable blastp está disponible en el PATH."""
    return shutil.which("blastp") is not None


def buscar_homologos_local(secuencia: str, db_path: str, max_hits: int = 10) -> List[HitHomologo]:
    """
    Ejecuta blastp local contra una base de datos ya descargada/construida
    (p. ej. con `update_blastdb.pl swissprot` o `makeblastdb`).

    Parámetros
    ----------
    secuencia : str
        Secuencia proteica problema.
    db_path : str
        Ruta base de la base de datos BLAST local (sin extensión), p. ej.
        "/home/usuario/blastdb/swissprot".
    max_hits : int
        Número máximo de homólogos a recuperar.
    """
    if not blastp_disponible():
        raise RuntimeError(
            "El ejecutable 'blastp' no está disponible en el PATH. "
            "Instala BLAST+ (p. ej. 'conda install -c bioconda blast')."
        )
    if not os.path.exists(db_path + ".pin") and not os.path.exists(db_path + ".pdb"):
        raise RuntimeError(
            f"No se encuentra una base de datos BLAST válida en '{db_path}'. "
            "Comprueba la ruta (sin extensión) y que la base de datos esté descargada."
        )

    with tempfile.TemporaryDirectory() as tmp:
        query_path = os.path.join(tmp, "query.fasta")
        with open(query_path, "w") as f:
            f.write(">query\n" + secuencia + "\n")

        cmd = [
            "blastp", "-query", query_path, "-db", db_path,
            "-max_target_seqs", str(max_hits),
            "-outfmt", "6 sseqid stitle pident length score evalue sseq",
        ]
        resultado = subprocess.run(cmd, capture_output=True, text=True)
        if resultado.returncode != 0:
            raise RuntimeError(f"Error ejecutando blastp: {resultado.stderr}")

        hits = []
        for linea in resultado.stdout.strip().splitlines():
            campos = linea.split("\t")
            if len(campos) < 7:
                continue
            sseqid, stitle, pident, length, score, evalue, sseq = campos[:7]
            hits.append(HitHomologo(
                id_hit=sseqid,
                descripcion=stitle if stitle else sseqid,
                identidad_pct=float(pident),
                longitud_alineamiento=int(length),
                score=float(score),
                evalue=float(evalue),
                secuencia=sseq.replace("-", ""),
            ))
        return hits


def buscar_homologos_remoto(
    secuencia: str,
    hitlist_size: int = 10,
    database: str = "nr",
    entrez_query: str = "",
) -> List[HitHomologo]:
    """
    Lanza una búsqueda BLASTp remota contra NCBI y devuelve los homólogos
    encontrados. Requiere conexión a internet; ver notas de rendimiento en
    el docstring del módulo.
    """
    from Bio.Blast import NCBIWWW, NCBIXML

    kwargs = {"hitlist_size": hitlist_size}
    if entrez_query:
        kwargs["entrez_query"] = entrez_query

    result_handle = NCBIWWW.qblast("blastp", database, secuencia, **kwargs)
    registros = NCBIXML.read(result_handle)

    hits = []
    for alignment in registros.alignments:
        hsp = alignment.hsps[0]  # mejor HSP de cada hit
        identidad_pct = round(100 * hsp.identities / hsp.align_length, 2)
        hits.append(HitHomologo(
            id_hit=alignment.hit_id,
            descripcion=alignment.hit_def,
            identidad_pct=identidad_pct,
            longitud_alineamiento=hsp.align_length,
            score=hsp.score,
            evalue=hsp.expect,
            secuencia=hsp.sbjct.replace("-", ""),
        ))
    return hits


def buscar_homologos(secuencia: str, modo: str = "local", **kwargs) -> List[HitHomologo]:
    """
    Interfaz única de búsqueda de homólogos.

    modo:
        "local"  -> kwargs: db_path (obligatorio), max_hits (opcional)
        "remote" -> kwargs: hitlist_size, database, entrez_query (opcionales)
    """
    if modo == "local":
        return buscar_homologos_local(secuencia, **kwargs)
    elif modo == "remote":
        return buscar_homologos_remoto(secuencia, **kwargs)
    else:
        raise ValueError(f"Modo de búsqueda desconocido: '{modo}'. Usa 'local' o 'remote'.")
