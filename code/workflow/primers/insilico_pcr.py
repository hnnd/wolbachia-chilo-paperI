"""3'-anchored in-silico PCR for allele-specific primer validation.

Why not blastn: with `blastn-short -word_size 7` a 20-mer returns ~10^5 HSPs per genome,
almost all of them spurious partial alignments, and an HSP that omits the primer's 3'
terminus is meaningless for PCR. Extension by Taq requires the 3' end to be paired, so
the assay is modelled directly: locate an exact seed at the 3' terminus, then count
mismatches back along the footprint.

Extension model (standard allele-specific PCR heuristic):
  amplifies  - 3'-terminal SEED_3P bases match exactly AND total footprint mismatches <= MAX_MM
  blocked    - any mismatch inside the 3'-terminal BLOCK_3P bases, or > MAX_MM mismatches
Mismatches 5+ bases from the 3' end are largely tolerated, which is why total mismatch
count alone does not decide specificity.
"""
SEED_3P = 8
BLOCK_3P = 3
MAX_MM = 3
AMP_MIN, AMP_MAX = 50, 1000

COMP = str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB")


def rc(s):
    return s.translate(COMP)[::-1]


def _seed_positions(hay, seed):
    out, i = [], hay.find(seed)
    while i != -1:
        out.append(i)
        i = hay.find(seed, i + 1)
    return out


def _mm(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def find_sites(template, primer):
    """Return binding sites of `primer` on `template`.

    Each site: {p3, orient, mm, mm_3p, extends}
      orient '+' : primer anneals to the minus strand and extends rightward (forward primer);
                   p3 is the template index of the primer's 3' base.
      orient '-' : primer anneals to the plus strand and extends leftward (reverse primer);
                   p3 is likewise the template index of the primer's 3' base.
    """
    L = len(primer)
    sites = []
    # forward orientation: primer sequence itself occurs on the plus strand
    seed = primer[-SEED_3P:]
    for pos in _seed_positions(template, seed):
        start = pos - (L - SEED_3P)
        if start < 0:
            continue
        foot = template[start:start + L]
        if len(foot) < L:
            continue
        mm = _mm(foot, primer)
        mm3 = _mm(foot[-BLOCK_3P:], primer[-BLOCK_3P:])
        sites.append({"p3": start + L - 1, "orient": "+", "mm": mm, "mm_3p": mm3,
                      "extends": mm3 == 0 and mm <= MAX_MM})
    # reverse orientation: reverse complement of the primer occurs on the plus strand,
    # and the primer's 3' base maps to the LEFT edge of that match
    r = rc(primer)
    seed_r = r[:SEED_3P]
    for pos in _seed_positions(template, seed_r):
        foot = template[pos:pos + L]
        if len(foot) < L:
            continue
        mm = _mm(foot, r)
        mm3 = _mm(foot[:BLOCK_3P], r[:BLOCK_3P])
        sites.append({"p3": pos, "orient": "-", "mm": mm, "mm_3p": mm3,
                      "extends": mm3 == 0 and mm <= MAX_MM})
    return sites


def pcr(template, fwd, rev):
    """Return list of predicted products: {size, f_mm, r_mm, f_p3, r_p3, strand}.

    BOTH amplicon orientations are tested. Primers are designed on gene-oriented CDS
    sequence, so for a gene on the template's minus strand the "forward" primer is found
    in '-' orientation and the "reverse" primer in '+'. Testing only the plus-strand case
    silently reports "no product" for every minus-strand locus -- which would make a
    non-specific pair look specific on the off-target genome.
    """
    lf, lr = len(fwd), len(rev)
    sites_f = [s for s in find_sites(template, fwd) if s["extends"]]
    sites_r = [s for s in find_sites(template, rev) if s["extends"]]
    prods = []
    for f in sites_f:
        for r in sites_r:
            # `size` is the PRODUCT length -- 5' end of the forward primer to 5' end of the
            # reverse primer, i.e. what runs on a gel. The 3'-to-3' distance is shorter by
            # (lf + lr - 2) and is NOT the band size.
            if f["orient"] == "+" and r["orient"] == "-":
                lo, hi, strand = f["p3"] - lf + 1, r["p3"] + lr - 1, "+"
            elif f["orient"] == "-" and r["orient"] == "+":
                lo, hi, strand = r["p3"] - lr + 1, f["p3"] + lf - 1, "-"
            else:
                continue
            if lo < 0 or hi >= len(template):
                continue
            size = hi - lo + 1
            if AMP_MIN <= size <= AMP_MAX:
                prods.append({"size": size, "lo": lo, "hi": hi, "strand": strand,
                              "f_mm": f["mm"], "r_mm": r["mm"],
                              "f_p3": f["p3"], "r_p3": r["p3"]})
    return prods


def pcr_multicontig(contigs, fwd, rev):
    """contigs: {name: seq}. Products cannot span contigs."""
    out = []
    for name, seq in contigs.items():
        for p in pcr(seq, fwd, rev):
            p = dict(p); p["contig"] = name
            out.append(p)
    return out


def worst_block(template, primer):
    """Best (most permissive) site for `primer` anywhere on `template`, either orientation.

    Used to show WHY a primer fails on the other supergroup: returns the site with the
    fewest total mismatches among those with a paired 3' terminus, else the fewest
    mismatches overall.
    """
    sites = find_sites(template, primer)
    if not sites:
        return None
    ext = [s for s in sites if s["mm_3p"] == 0]
    pool = ext if ext else sites
    return min(pool, key=lambda s: s["mm"])
