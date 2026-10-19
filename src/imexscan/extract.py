"""Extract filamentous prophages for one genome: cluster marker genes,
split clusters at dif-like sites, pick attL/attR boundaries, write sequences."""

import csv
import glob
import os
from collections import defaultdict

from .dif import open_maybe_gz
from .markers import MARKER_KIND, INOVIRUS_SIGNAL, gene_name

GAP = 10000
MAX_DIF_DIST = 20000
FLANK = 5000

SUMMARY_HEADER = ["strain", "prophage", "type", "contig", "start", "end", "length",
                  "boundary", "quality", "markers",
                  "attL_S", "attL_E", "attL_type", "attL_seq",
                  "attR_S", "attR_E", "attR_type", "attR_seq",
                  "gene_S", "gene_S_start", "gene_S_end",
                  "gene_E", "gene_E_start", "gene_E_end"]


def read_fasta(path):
    seqs, cur, curid = {}, [], None
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if curid is not None:
                    seqs[curid] = "".join(cur)
                curid = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.strip())
    if curid is not None:
        seqs[curid] = "".join(cur)
    return seqs


def load_dif(dif_dir, strain):
    path = os.path.join(dif_dir, f"{strain}-find_Vibrio_dif_motif-output.csv")
    dif = defaultdict(list)
    if not os.path.exists(path):
        return dif
    with open(path) as fh:
        for row in csv.DictReader(fh):
            c = row["Sequence_ID"].split()[0]
            s, e, t = int(row["Start_Position"]), int(row["End_Position"]), row["Type_of_DIF-like_att"]
            dif[c].append((s, e, t, row["Match_Sequence"]))
    for c in dif:
        dif[c].sort()
    return dif


def load_contig_map(fna):
    cmap, i = {}, 0
    with open_maybe_gz(fna) as fh:
        for line in fh:
            if line.startswith(">"):
                i += 1
                hdr = line[1:].split()[0]
                cmap[str(i)] = hdr
                cmap[f"seq{i}"] = hdr
                cmap[f"sequence{i}"] = hdr
    return cmap


def load_proteins(fna, prodigal_dir, stem):
    gff_path = os.path.join(prodigal_dir, f"{stem}.gff")
    cmap = load_contig_map(fna)
    prot = {}
    if os.path.exists(gff_path):
        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue
                contig = cmap.get(parts[0], parts[0])
                start, end = int(parts[3]), int(parts[4])
                strand = parts[6]
                attrs = dict(p.split("=", 1) for p in parts[8].split(";") if "=" in p)
                pid = attrs.get("ID", "")
                if pid:
                    prot[pid] = (contig, start, end, strand)
    return prot


def load_markers(dom_dir, strain):
    marks = defaultdict(list)
    for f in glob.glob(os.path.join(dom_dir, f"{strain}.*.domtbl")):
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.split()
                if len(cols) < 4:
                    continue
                marks[cols[0]].append(cols[3])
    return marks


def classify_cluster(marker_names, min_hallmark):
    kinds = {MARKER_KIND.get(m, "unknown") for m in marker_names}
    n_inov = len(kinds & (INOVIRUS_SIGNAL | {"shared"}))
    return "inovirus" if n_inov >= min_hallmark else None


def find_bounds(contig, core_start, core_end, dif):
    difs = dif.get(contig, [])
    left = [d for d in difs if d[0] <= core_start and d[2] == "full-length"] or \
           [d for d in difs if d[0] <= core_start]
    right = [d for d in difs if d[1] >= core_end and d[2] == "full-length"] or \
            [d for d in difs if d[1] >= core_end]
    attL = max(left, key=lambda d: d[0]) if left else None
    attR = min(right, key=lambda d: d[1]) if right else None
    l = attL[0] if attL else None
    r = attR[1] if attR else None
    if l is not None and core_start - l > MAX_DIF_DIST:
        l, attL = None, None
    if r is not None and r - core_end > MAX_DIF_DIST:
        r, attR = None, None
    if l is not None and r is not None:
        return l, r, "both", attL, attR
    if l is not None:
        return l, core_end + FLANK, "left", attL, None
    if r is not None:
        return core_start - FLANK, r, "right", None, attR
    return core_start - FLANK, core_end + FLANK, "none", None, None


def merge_dif_sites(difs, window=100):
    sites = sorted(d for d in difs if d[2] == "full-length")
    for d in sorted(d for d in difs if d[2] != "full-length"):
        if not any(abs(s[0] - d[0]) <= window for s in sites):
            sites.append(d)
    return sorted(sites)


def split_by_difs(cl, sites):
    if len(sites) < 2:
        return [(cl, None, None)]
    pieces = []
    for i in range(len(sites) - 1):
        lo, hi = sites[i][0], sites[i + 1][1]
        sub = [g for g in cl if g[1] >= lo and g[0] <= hi]
        if sub:
            pieces.append((sub, sites[i], sites[i + 1]))
    pre = [g for g in cl if g[1] < sites[0][0]]
    if pre:
        pieces.append((pre, None, sites[0]))
    post = [g for g in cl if g[0] > sites[-1][1]]
    if post:
        pieces.append((post, sites[-1], None))
    pieces.sort(key=lambda p: min(g[0] for g in p[0]))
    return pieces if pieces else [(cl, None, None)]


def quality(length, boundary):
    if length < 5000:
        return "fragment"
    if boundary in ("left", "right", "none"):
        return "incomplete"
    return "complete"


def extract_genome(stem, fna, prodigal_dir, dom_dir, dif_dir, seq_dir, min_hallmark):
    seqs = read_fasta(fna)
    dif = load_dif(dif_dir, stem)
    prot = load_proteins(fna, prodigal_dir, stem)
    marks = load_markers(dom_dir, stem)

    by_contig = defaultdict(list)
    for pid, mlist in marks.items():
        if pid in prot:
            c, s, e, st = prot[pid]
            by_contig[c].append((s, e, pid, mlist))

    prophages = []
    for c, genes in by_contig.items():
        genes.sort()
        clusters = [[genes[0]]]
        for g_ in genes[1:]:
            if g_[0] - clusters[-1][-1][1] <= GAP:
                clusters[-1].append(g_)
            else:
                clusters.append([g_])
        for cl in clusters:
            lo0, hi0 = cl[0][0] - 1000, cl[-1][1] + 1000
            nearby = [d for d in dif.get(c, []) if d[0] <= hi0 and d[1] >= lo0]
            sites = merge_dif_sites(nearby)
            pieces = split_by_difs(cl, sites)
            for sub, attL, attR in pieces:
                all_names = [m for g_ in sub for m in g_[3]]
                ptype = classify_cluster(all_names, min_hallmark)
                if ptype is None:
                    continue
                if attL is not None and attR is not None:
                    start, end, bc = attL[0], attR[1], "both"
                else:
                    cs = min(g_[0] for g_ in sub)
                    ce = max(g_[1] for g_ in sub)
                    start, end, bc, attL, attR = find_bounds(c, cs, ce, dif)
                start = max(1, start)
                end = min(len(seqs.get(c, "")), end)
                if end <= start:
                    continue
                prophages.append((c, ptype, start, end, bc, sub, attL, attR))

    prophages.sort(key=lambda x: (x[0], x[2]))
    rows = []
    cnt_fp = 0
    for c, ptype, start, end, bc, cl, attL, attR in prophages:
        length = end - start + 1
        cnt_fp += 1
        suffix = "abcdefghijklmnopqrstuvwxyz"[cnt_fp - 1] if cnt_fp > 1 else ""
        name = f"FPphi{stem}{suffix}"
        q = quality(length, bc)
        seq = seqs[c][start - 1:end]
        with open(os.path.join(seq_dir, f"{name}.fasta"), "w") as fh:
            fh.write(f">{name} {c}:{start}-{end} len={length} type={ptype} boundary={bc} quality={q}\n")
            fh.write(seq + "\n")
        msum = ",".join(sorted({gene_name(m) for g_ in cl for m in g_[3]}))
        cl_sorted = sorted(cl, key=lambda g: g[0])
        gene_names_all = lambda g: {gene_name(m) for m in g[3]}
        rstR_cands = [g for g in cl if "RstR" in gene_names_all(g)]
        rstA_cands = [g for g in cl if "RstA" in gene_names_all(g)]
        if rstR_cands:
            rstR_gene = rstR_cands[0]
        elif rstA_cands:
            rstR_gene = rstA_cands[0]
        else:
            rstR_gene = cl_sorted[0]
        xaft_gene = next((g for g in cl if "XafT" in gene_names_all(g)), None)
        zot_gene = next((g for g in cl if "Zot" in gene_names_all(g)), None)
        end_gene = xaft_gene if xaft_gene else (zot_gene if zot_gene else rstR_gene)
        gS, gE = rstR_gene, end_gene
        geneS = "/".join(sorted(gene_names_all(gS)))
        geneE = "/".join(sorted(gene_names_all(gE)))
        aL = attL if attL else ("N/A", "N/A", "N/A", "N/A")
        aR = attR if attR else ("N/A", "N/A", "N/A", "N/A")
        rows.append([stem, name, ptype, c, start, end, length, bc, q, msum,
                     aL[0], aL[1], aL[2], aL[3], aR[0], aR[1], aR[2], aR[3],
                     geneS, gS[0], gS[1], geneE, gE[0], gE[1]])
    return rows
