#!/usr/bin/env python3
"""Concatenate per-locus FASTA alignments into a partitioned supermatrix.

Usage: concat_aln.py OUT_FASTA OUT_PARTITION locus1.aln locus2.aln ...
Taxa missing a locus are padded with gaps for that locus's aligned length, so the
5-locus MLST supermatrix tolerates the partial locus recovery typical of low-coverage
metagenomic Wolbachia. Partition file is RAxML/IQ-TREE style (one charset per locus).
"""
import sys
from pathlib import Path


def read_fasta(p: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main() -> None:
    out_fa, out_part = sys.argv[1], sys.argv[2]
    aln_files = [Path(x) for x in sys.argv[3:]]
    per_locus = []          # list of (locus_name, {taxon: seq}, aln_len)
    taxa: set[str] = set()
    for f in aln_files:
        d = read_fasta(f)
        if not d:
            continue
        aln_len = len(next(iter(d.values())))
        # sanity: all same length (true for an alignment)
        for t, s in d.items():
            if len(s) != aln_len:
                raise SystemExit(f"{f}: ragged alignment ({t} len {len(s)} != {aln_len})")
        per_locus.append((f.stem.replace(".aln", ""), d, aln_len))
        taxa.update(d.keys())

    taxa_sorted = sorted(taxa)
    # Build concatenated sequence per taxon
    cat: dict[str, list[str]] = {t: [] for t in taxa_sorted}
    parts = []
    pos = 1
    for locus, d, aln_len in per_locus:
        for t in taxa_sorted:
            cat[t].append(d.get(t, "-" * aln_len))
        parts.append(f"DNA, {locus} = {pos}-{pos + aln_len - 1}")
        pos += aln_len

    with open(out_fa, "w") as fh:
        for t in taxa_sorted:
            fh.write(f">{t}\n{''.join(cat[t])}\n")
    Path(out_part).write_text("\n".join(parts) + "\n")
    total = pos - 1
    sys.stderr.write(
        f"[concat] {len(taxa_sorted)} taxa x {len(per_locus)} loci = {total} bp\n"
    )


if __name__ == "__main__":
    main()
