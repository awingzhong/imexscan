"""Table-level steps: merge per-genome summaries, expand genes, extract pHVR."""

import csv
import os
from collections import defaultdict

from .dif import open_maybe_gz
from .extract import SUMMARY_HEADER
from .markers import gene_name

GENE_HEADER = ["Bac_Contig", "Bac_Species", "Prophage_Name", "Prophage_Type",
               "Gene_Start", "Gene_End", "Strand",
               "Hallmarks", "Prokka_Gene_Accession", "Gene_Annotation", "Protein_Sequences"]

PHVR_HEADER = ["Bac_Contig", "Bac_Species", "Prophage_Name", "Gene_Start", "Gene_End",
               "Strand", "Prokka_Gene_Accession", "Gene_Annotation", "Protein_Sequences"]


def merge_summary(rows, table_dir):
    out = os.path.join(table_dir, "FPphi_summary.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(SUMMARY_HEADER) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"FPphi_summary.tsv: {len(rows)} prophages")


def _contig_map(fna):
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


def _load_prodigal(strain, fna, prodigal_dir):
    gff_path = os.path.join(prodigal_dir, f"{strain}.gff")
    faa_path = os.path.join(prodigal_dir, f"{strain}.faa")
    cmap = _contig_map(fna) if fna and os.path.exists(fna) else {}
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
                    prot[pid] = [contig, start, end, strand, ""]
    if os.path.exists(faa_path):
        curid, cur = None, []
        with open(faa_path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if curid and curid in prot:
                        prot[curid][4] = "".join(cur)
                    curid = line[1:].split()[0]
                    cur = []
                else:
                    cur.append(line.strip())
            if curid and curid in prot:
                prot[curid][4] = "".join(cur)
    return prot


def _load_markers(strain, marker_index):
    marks = defaultdict(list)
    for f in marker_index.get(strain, []):
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.split()
                if len(cols) < 4:
                    continue
                marks[cols[0]].append(cols[3])
    return marks


def _species_of(fna):
    if not fna or not os.path.exists(fna):
        return ""
    with open_maybe_gz(fna) as fh:
        for line in fh:
            if line.startswith(">"):
                toks = line[1:].split()
                return " ".join(toks[1:3])
    return ""


def expand_genes(table_dir, prodigal_dir, dom_dir, manifest):
    summary = os.path.join(table_dir, "FPphi_summary.tsv")
    out_tsv = os.path.join(table_dir, "FPphi_gene_annotation.tsv")

    marker_index = defaultdict(list)
    for f in os.listdir(dom_dir):
        if f.endswith(".domtbl"):
            marker_index[f.rsplit(".", 2)[0]].append(os.path.join(dom_dir, f))

    groups = defaultdict(list)
    with open(summary) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            groups[row["strain"]].append(row)

    out_rows, gene_map = [], {}
    for strain, rows in groups.items():
        fna = manifest.get(strain, "")
        prot = _load_prodigal(strain, fna, prodigal_dir)
        marks = _load_markers(strain, marker_index)
        species = _species_of(fna)

        for row in rows:
            name, ptype = row["prophage"], row["type"]
            contig, start, end = row["contig"], int(row["start"]), int(row["end"])
            genes = [(pid, v) for pid, v in prot.items()
                     if v[0] == contig and v[1] <= end and v[2] >= start]
            genes.sort(key=lambda pv: pv[1][1])
            for pid, (c, gs, ge, strand, seq) in genes:
                hs = sorted({gene_name(m) for m in marks.get(pid, [])})
                row_data = [contig, species, name, ptype,
                            gs, ge, "+" if strand in ("1", "+") else "-",
                            "/".join(hs), pid, "", seq]
                if pid in gene_map:
                    out_rows[gene_map[pid]] = row_data
                else:
                    gene_map[pid] = len(out_rows)
                    out_rows.append(row_data)

    with open(out_tsv, "w") as fh:
        fh.write("\t".join(GENE_HEADER) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"FPphi_gene_annotation.tsv: {len(out_rows)} genes")


def _backfill_summary(summary, has_hvr):
    with open(summary) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        if "has_pHVR" not in fieldnames:
            fieldnames.append("has_pHVR")
        rows = []
        for r in reader:
            r["has_pHVR"] = "yes" if r["prophage"] in has_hvr else "no"
            rows.append(r)
    with open(summary, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def extract_phvr(table_dir):
    summary = os.path.join(table_dir, "FPphi_summary.tsv")
    gene_ann = os.path.join(table_dir, "FPphi_gene_annotation.tsv")
    out_tsv = os.path.join(table_dir, "FPphi_pHVR_annotation.tsv")
    out_faa = os.path.join(table_dir, "FPphi_pHVR.faa")

    by_pro = defaultdict(list)
    with open(gene_ann) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = dict(r)
            g["Gene_Start"] = int(g["Gene_Start"])
            g["Gene_End"] = int(g["Gene_End"])
            g["_hs"] = set(g["Hallmarks"].split("/")) if g["Hallmarks"] else set()
            by_pro[g["Prophage_Name"]].append(g)
    for pro in by_pro:
        by_pro[pro].sort(key=lambda g: g["Gene_Start"])

    out_rows, faa_rows, has_hvr = [], [], set()
    for pro, genes in by_pro.items():
        zots = [g for g in genes if "Zot" in g["_hs"]]
        xafts = [g for g in genes if "XafT" in g["_hs"]]
        if not zots or not xafts:
            continue
        zot = zots[0]
        xaft = xafts[0]
        lo = min(zot["Gene_Start"], xaft["Gene_Start"])
        hi = max(zot["Gene_Start"], xaft["Gene_Start"])
        exclude = {zot["Prokka_Gene_Accession"], xaft["Prokka_Gene_Accession"]}
        hvrs = [g for g in genes
                if lo < g["Gene_Start"] < hi and g["Prokka_Gene_Accession"] not in exclude]
        if not hvrs:
            continue
        has_hvr.add(pro)
        for g in hvrs:
            out_rows.append([g["Bac_Contig"], g["Bac_Species"], g["Prophage_Name"],
                             g["Gene_Start"], g["Gene_End"], g["Strand"],
                             g["Prokka_Gene_Accession"], g["Gene_Annotation"],
                             g["Protein_Sequences"]])
            faa_rows.append((g["Prokka_Gene_Accession"], g["Protein_Sequences"]))

    with open(out_tsv, "w") as fh:
        fh.write("\t".join(PHVR_HEADER) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(out_faa, "w") as fh:
        for acc, seq in faa_rows:
            fh.write(f">{acc}\n{seq}\n")
    _backfill_summary(summary, has_hvr)

    print(f"FPphi_pHVR_annotation.tsv: {len(out_rows)} genes from {len(has_hvr)} prophages")
    print(f"FPphi_pHVR.faa: {len(faa_rows)} proteins")
