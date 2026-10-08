"""Design Wolbachia supergroup A- and B-specific PCR primers for single-wasp typing.

Strategy
--------
Allele-specific PCR: the primer 3' end must sit on a supergroup-diagnostic position,
because a 3'-terminal mismatch is what actually blocks extension by Taq. A single
internal mismatch is NOT reliable discrimination, so every candidate is additionally
required to carry >=3 diagnostic mismatches across its footprint.

"Diagnostic" is defined against TWO references per supergroup
(A: wMel NC_002978.6 + wRi NC_012416.1; B: wPip NC_010981.1 + wTpre NZ_CM003641.1):
a column counts only if both A refs agree, both B refs agree, and A != B. That makes a
site supergroup-diagnostic rather than strain-specific — the property needed for primers
that must work on an uncharacterised local strain.

Loci are restricted to genes that are single-copy in all four genomes
(gatB, ftsZ, fbpA, groEL, 16S). wsp is deliberately excluded: it has 8-14 paralogous
copies in these genomes and is recombination-prone, so it cannot support a
supergroup-specific assay even though the classic supergroup primers target it.
"""
import primer3

MIN_LEN, MAX_LEN = 18, 25
TM_LO, TM_HI = 57.0, 64.0
GC_LO, GC_HI = 35.0, 65.0
MIN_DIAG = 3          # diagnostic mismatches anywhere in the footprint
DIAG_3P_WINDOW = 3    # ...of which at least one must be in the last N bases
AMP_LO, AMP_HI = 120, 500
HAIRPIN_DG = -3000.0  # cal/mol, primer3 returns cal/mol
HOMODIMER_DG = -6000.0
HETERODIMER_DG = -6000.0

COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def ref_template(al, ref_id, mask, same_sg_refs):
    """Ungapped template taken from a REAL reference, with per-base diagnostic flags.

    Designing on a cross-strain "consensus" is wrong: at columns where two same-supergroup
    references disagree the consensus is not a subsequence of either genome, so a primer
    spanning such a column matches neither. Instead the template IS one reference genome's
    sequence; the alignment only supplies (a) which bases are supergroup-diagnostic and
    (b) which bases are conserved across the other reference of the same supergroup.

    Returns (seq, diag_flags, conserved_flags, columns). `conserved_flags[i]` is False where
    any same-supergroup reference differs, and candidate primers must avoid those bases so
    the assay does not become strain-specific.
    """
    seq, diag, cons, cols = [], [], [], []
    for i, ch in enumerate(al[ref_id]):
        if ch == "-":
            continue
        seq.append(ch)
        diag.append(mask[i])
        cons.append(all(al[r][i] == ch for r in same_sg_refs))
        cols.append(i)
    return "".join(seq), diag, cons, cols


def gc(s):
    return 100.0 * (s.count("G") + s.count("C")) / len(s)


def has_run(s, n=5):
    run = 1
    for a, b in zip(s, s[1:]):
        run = run + 1 if a == b else 1
        if run >= n:
            return True
    return False


def candidates(seq, flags, strand, cons=None):
    """Yield primer candidates on one strand of a reference template.

    strand '+' -> forward primer, sequence as-is, 3' end at the right edge.
    strand '-' -> reverse primer, reverse complement, 3' end at the LEFT edge of the window.
    `cons` (optional) is the within-supergroup conservation flag per base; any footprint
    containing a non-conserved base is rejected, so the primer matches every reference of
    its own supergroup and the assay generalises to an uncharacterised local strain.
    """
    out = []
    n = len(seq)
    for L in range(MIN_LEN, MAX_LEN + 1):
        for i in range(n - L + 1):
            win_flags = flags[i:i + L]
            ndiag = sum(win_flags)
            if ndiag < MIN_DIAG:
                continue
            if cons is not None and not all(cons[i:i + L]):
                continue
            if strand == "+":
                if not any(win_flags[-DIAG_3P_WINDOW:]):
                    continue
                pseq = seq[i:i + L]
                p3 = i + L - 1          # genomic index of the 3' base
            else:
                if not any(win_flags[:DIAG_3P_WINDOW]):
                    continue
                pseq = rc(seq[i:i + L])
                p3 = i
            if "N" in pseq or has_run(pseq):
                continue
            g = gc(pseq)
            if not (GC_LO <= g <= GC_HI):
                continue
            tm = primer3.calc_tm(pseq)
            if not (TM_LO <= tm <= TM_HI):
                continue
            if primer3.calc_hairpin(pseq).dg < HAIRPIN_DG:
                continue
            if primer3.calc_homodimer(pseq).dg < HOMODIMER_DG:
                continue
            out.append({"seq": pseq, "len": L, "tm": round(tm, 1), "gc": round(g, 1),
                        "start": i, "end": i + L - 1, "p3_idx": p3,
                        "ndiag": ndiag, "diag_3p": sum(win_flags[-DIAG_3P_WINDOW:]) if strand == "+"
                        else sum(win_flags[:DIAG_3P_WINDOW]), "strand": strand})
    return out


def pair_up(fwd, rev, locus, sg, max_pairs=400):
    """Pair forward and reverse candidates into amplicons."""
    pairs = []
    rev_sorted = sorted(rev, key=lambda r: r["p3_idx"])
    for f in fwd:
        for r in rev_sorted:
            amp = r["p3_idx"] - f["p3_idx"] + 1
            if not (AMP_LO <= amp <= AMP_HI):
                continue
            if abs(f["tm"] - r["tm"]) > 2.0:
                continue
            if primer3.calc_heterodimer(f["seq"], r["seq"]).dg < HETERODIMER_DG:
                continue
            pairs.append({"locus": locus, "supergroup": sg, "amplicon_bp": amp,
                          "F": f["seq"], "F_tm": f["tm"], "F_gc": f["gc"], "F_ndiag": f["ndiag"],
                          "F_diag3p": f["diag_3p"],
                          "R": r["seq"], "R_tm": r["tm"], "R_gc": r["gc"], "R_ndiag": r["ndiag"],
                          "R_diag3p": r["diag_3p"],
                          "min_ndiag": min(f["ndiag"], r["ndiag"]),
                          "f_start": f["start"], "r_end": r["end"]})
            if len(pairs) >= max_pairs:
                return pairs
    return pairs
