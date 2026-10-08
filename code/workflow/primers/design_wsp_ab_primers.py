"""Design wCchiA / wCchiB strain-specific primers on the wsp locus.

Rationale for the acceptance rules (they follow the house convention established in
SINGLE_WASP_AB_PRIMERS.md, applied here to the two locally assembled genomes rather
than to external supergroup references):

  * Taq extension needs a paired 3' terminus, so a primer only discriminates if a
    DIAGNOSTIC position (A != B) sits in its 3'-terminal window. Internal mismatches
    alone do not block amplification.
  * Every candidate is additionally run through a 3'-anchored in-silico PCR against
    BOTH complete genomes: it must give exactly one product on its own genome and
    zero on the other. This catches off-target sites the local alignment cannot see.
  * Product length is the gel band length (5' end of F to 5' end of R), never the
    3'-to-3' distance.

Outputs (pd/out/): wsp_ab_primers.csv, wsp_ab_amplicons.fa, wsp_ab_design.json
"""
import json, re, itertools
import numpy as np
import primer3
from Bio.Align import PairwiseAligner
import insilico_pcr as ip

FLANK = 150
LEN_MIN, LEN_MAX = 18, 26
TM_MIN, TM_MAX = 56.0, 63.0
GC_MIN, GC_MAX = 30.0, 65.0
HAIRPIN_DG = -3000.0        # cal/mol, primer3 returns cal/mol
HOMODIMER_DG = -6000.0
HETERODIMER_DG = -6000.0
DIAG_3P_WINDOW = 3          # >=1 diagnostic base within the last 3 nt
DIAG_MIN_TOTAL = 3          # >=3 diagnostic mismatches across the footprint

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def load_fa(path):
    d, n = {}, None
    for line in open(path):
        if line.startswith(">"):
            n = line[1:].split()[0]
            d[n] = []
        else:
            d[n].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def align_alleles(tA, tB):
    """Global alignment -> (diagnostic mask on A, A-index -> B-index map)."""
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                         open_gap_score=-8, extend_gap_score=-1)
    a = al.align(tA, tB)[0]
    X, Y = str(a[0]), str(a[1])
    diag = np.zeros(len(tA), bool)
    amap = {}
    ia = ib = -1
    for x, y in zip(X, Y):
        if x != "-":
            ia += 1
        if y != "-":
            ib += 1
        if x != "-":
            if x != y:
                diag[ia] = True
            amap[ia] = ib if y != "-" else None
    return diag, amap, X, Y


def _bad_run(s, n=5):
    return re.search(r"(.)\1{%d,}" % (n - 1), s) is not None


def thermo_ok(seq):
    tm = primer3.calc_tm(seq)
    gc = 100.0 * (seq.count("G") + seq.count("C")) / len(seq)
    if not (TM_MIN <= tm <= TM_MAX and GC_MIN <= gc <= GC_MAX):
        return None
    if _bad_run(seq):
        return None
    if primer3.calc_hairpin(seq).dg < HAIRPIN_DG:
        return None
    if primer3.calc_homodimer(seq).dg < HOMODIMER_DG:
        return None
    return dict(tm=round(tm, 1), gc=round(gc, 1))


def candidates(template, diag, orient):
    """Enumerate primers on `template` (A-strand coordinates).

    orient 'F': primer sequence == template[i:j], 3' end at j-1.
    orient 'R': primer sequence == rc(template[i:j]), 3' end at i.
    """
    out = []
    n = len(template)
    for L in range(LEN_MIN, LEN_MAX + 1):
        for i in range(0, n - L):
            j = i + L
            foot = template[i:j]
            if "N" in foot:
                continue
            seq = foot if orient == "F" else rc(foot)
            # diagnostic positions inside the footprint, and inside the 3' window
            d_all = int(diag[i:j].sum())
            if orient == "F":
                d_3p = int(diag[j - DIAG_3P_WINDOW:j].sum())
                d_last = bool(diag[j - 1])
                p3 = j - 1
            else:
                d_3p = int(diag[i:i + DIAG_3P_WINDOW].sum())
                d_last = bool(diag[i])
                p3 = i
            t = thermo_ok(seq)
            if t is None:
                continue
            out.append(dict(seq=seq, i=i, j=j, orient=orient, p3=p3,
                            diag_total=d_all, diag_3p=d_3p, diag_last=d_last, **t))
    return out


def specific(c):
    return c["diag_3p"] >= 1 and c["diag_total"] >= DIAG_MIN_TOTAL


def conserved(c):
    return c["diag_total"] == 0


def pair_ok(f, r):
    if abs(f["tm"] - r["tm"]) > 2.5:
        return False
    return primer3.calc_heterodimer(f["seq"], r["seq"]).dg >= HETERODIMER_DG


def score(f, r, size, target_size):
    s = 0.0
    s += 3.0 * (f["diag_last"] + r["diag_last"])
    s += 1.0 * (min(f["diag_3p"], 3) + min(r["diag_3p"], 3))
    s += 0.4 * (min(f["diag_total"], 8) + min(r["diag_total"], 8))
    s -= 0.5 * (abs(f["tm"] - 60) + abs(r["tm"] - 60))
    s -= 1.5 * abs(f["tm"] - r["tm"])
    s -= 0.004 * abs(size - target_size)
    return s


def build_pairs(cands_f, cands_r, lo, hi, target, need, genome_own, genome_other,
                topn=8):
    """Return validated pairs, best first. `need` = 'specific' | 'conserved'."""
    filt = specific if need == "specific" else conserved
    F = [c for c in cands_f if filt(c)]
    R = [c for c in cands_r if filt(c)]
    F.sort(key=lambda c: (-c["diag_last"], -c["diag_3p"], abs(c["tm"] - 60)))
    R.sort(key=lambda c: (-c["diag_last"], -c["diag_3p"], abs(c["tm"] - 60)))
    F, R = F[:400], R[:400]
    cand_pairs = []
    for f in F:
        for r in R:
            size = r["j"] - f["i"]          # 5'(F) .. 5'(R) inclusive length
            if not (lo <= size <= hi):
                continue
            if not pair_ok(f, r):
                continue
            cand_pairs.append((score(f, r, size, target), f, r, size))
    cand_pairs.sort(key=lambda x: -x[0])
    out = []
    for sc, f, r, size in cand_pairs:
        own = ip.pcr_multicontig(genome_own, f["seq"], r["seq"])
        oth = ip.pcr_multicontig(genome_other, f["seq"], r["seq"])
        if len(own) != 1 or len(oth) != 0:
            continue
        out.append(dict(score=round(sc, 2), size=own[0]["size"], F=f, R=r,
                        own_products=len(own), other_products=len(oth)))
        if len(out) >= topn:
            break
    return out
