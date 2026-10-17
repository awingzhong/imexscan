"""End-to-end orchestration: prodigal -> dif scan + hmmsearch -> extraction -> tables."""

import glob
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

from . import dif as dif_mod
from . import extract as extract_mod
from . import tables as tables_mod

GENOME_EXTS = (".fa", ".fasta", ".fna")
DEFAULT_HMM_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "Vibrio_Inovirus")
PRODIGAL_ATTEMPTS = 3


def find_genomes(genome_dir):
    found = []
    for f in sorted(os.listdir(genome_dir)):
        base = f[:-3] if f.endswith(".gz") else f
        for ext in GENOME_EXTS:
            if base.endswith(ext):
                found.append((base[: -len(ext)], os.path.join(genome_dir, f)))
                break
    seen = {}
    for stem, path in found:
        if stem in seen:
            raise SystemExit(f"Error: duplicated genome name '{stem}': {seen[stem]} and {path}")
        seen[stem] = path
    return found


def _plain_fasta(path, workdir):
    if not path.endswith(".gz"):
        return path
    out = os.path.join(workdir, os.path.basename(path)[:-3])
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def _prefix_faa(raw, out, stem):
    with open(raw) as fi, open(out, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                i = line.rfind("ID=")
                if i != -1:
                    pid = line[i + 3:].split(";", 1)[0]
                    line = re.sub(r"^>\S+", f">{stem}_{pid}", line, count=1)
                    line = re.sub(r"ID=[^;]+", f"ID={stem}_{pid}", line, count=1)
            fo.write(line)


def _prefix_gff(raw, out, stem):
    with open(raw) as fi, open(out, "w") as fo:
        for line in fi:
            if not line.startswith("#"):
                line = line.replace("ID=", f"ID={stem}_", 1)
            fo.write(line)


def run_prodigal(binary, fna, stem, prodigal_dir):
    faa = os.path.join(prodigal_dir, f"{stem}.faa")
    gff = os.path.join(prodigal_dir, f"{stem}.gff")
    if os.path.exists(faa) and os.path.exists(gff):
        return False
    with tempfile.TemporaryDirectory() as tmp:
        plain = _plain_fasta(fna, tmp)
        raw_faa, raw_gff = os.path.join(tmp, "raw.faa"), os.path.join(tmp, "raw.gff")
        for attempt in range(1, PRODIGAL_ATTEMPTS + 1):
            for p in (raw_faa, raw_gff):
                if os.path.exists(p):
                    os.remove(p)
            proc = subprocess.run(
                [binary, "-i", plain, "-a", raw_faa, "-f", "gff", "-o", raw_gff, "-p", "meta", "-q"],
                capture_output=True, text=True)
            if proc.returncode == 0:
                break
            # prodigal on macOS/arm64 occasionally segfaults after writing truncated output
            print(f"[warn] prodigal exited {proc.returncode} for {stem}, "
                  f"retrying ({attempt}/{PRODIGAL_ATTEMPTS})", file=sys.stderr, flush=True)
        else:
            raise RuntimeError(f"prodigal failed (exit {proc.returncode}): "
                               f"{proc.stderr.strip()[-300:]}")
        _prefix_faa(raw_faa, faa, stem)
        _prefix_gff(raw_gff, gff, stem)
    return True


def run_dif(fna, stem, dif_dir):
    out = os.path.join(dif_dir, f"{stem}-find_Vibrio_dif_motif-output.csv")
    if os.path.exists(out):
        return False
    dif_mod.write_csv(dif_mod.scan(fna), out)
    return True


def run_hmmsearch(binary, faa, stem, dom_dir, hmm_files, evalue):
    for hmm in hmm_files:
        name = os.path.splitext(os.path.basename(hmm))[0]
        out = os.path.join(dom_dir, f"{stem}.{name}.domtbl")
        if os.path.exists(out):
            continue
        proc = subprocess.run([binary, "--domtblout", out, "-E", str(evalue), hmm, faa],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            if os.path.exists(out):
                os.remove(out)
            print(f"[warn] hmmsearch failed for {stem} / {name}",
                  file=sys.stderr, flush=True)


def _process_genome(job):
    stem, fna = job["genome"]
    prodigal_dir, dom_dir, dif_dir, seq_dir = job["dirs"]
    try:
        run_prodigal(job["prodigal"], fna, stem, prodigal_dir)
        if not job["skip_dif"]:
            run_dif(fna, stem, dif_dir)
        faa = os.path.join(prodigal_dir, f"{stem}.faa")
        if os.path.exists(faa):
            run_hmmsearch(job["hmmsearch"], faa, stem, dom_dir,
                          job["hmm_files"], job["evalue"])
        rows = extract_mod.extract_genome(stem, fna, prodigal_dir, dom_dir,
                                          dif_dir, seq_dir, job["min_hallmark"])
        return stem, rows, None
    except Exception as exc:
        return stem, [], f"{type(exc).__name__}: {exc}"


def run(input_dir, out_dir, threads=4, hmm_dir=None, min_hallmark=4, evalue="1e-5",
        skip_dif=False, prodigal="prodigal", hmmsearch="hmmsearch"):
    genomes = find_genomes(input_dir)
    if not genomes:
        raise SystemExit(f"Error: no genome files (*.fa/*.fasta/*.fna[.gz]) in {input_dir}")

    hmm_dir = hmm_dir or DEFAULT_HMM_DIR
    hmm_files = sorted(glob.glob(os.path.join(hmm_dir, "*.hmm")))
    if not hmm_files:
        raise SystemExit(f"Error: no .hmm files in {hmm_dir}")

    prodigal_dir = os.path.join(out_dir, "prodigal")
    dom_dir = os.path.join(out_dir, "dom")
    dif_dir = os.path.join(out_dir, "dif")
    seq_dir = os.path.join(out_dir, "results", "FPphi")
    table_dir = os.path.join(out_dir, "tables")
    for d in (prodigal_dir, dom_dir, dif_dir, seq_dir, table_dir):
        os.makedirs(d, exist_ok=True)

    jobs = [{"genome": g, "dirs": (prodigal_dir, dom_dir, dif_dir, seq_dir),
             "hmm_files": hmm_files, "evalue": evalue, "min_hallmark": min_hallmark,
             "skip_dif": skip_dif, "prodigal": prodigal, "hmmsearch": hmmsearch}
            for g in genomes]

    print(f"InoMiner: {len(genomes)} genomes, {threads} parallel workers, "
          f"{len(hmm_files)} marker HMMs")
    by_stem, n_fail = {}, 0
    with ProcessPoolExecutor(max_workers=threads) as pool:
        futures = [pool.submit(_process_genome, j) for j in jobs]
        for i, fut in enumerate(as_completed(futures), 1):
            stem, rows, err = fut.result()
            by_stem[stem] = rows
            if err:
                n_fail += 1
                print(f"[{i}/{len(jobs)}] {stem}: FAILED ({err})",
                      file=sys.stderr, flush=True)
            else:
                print(f"[{i}/{len(jobs)}] {stem}: {len(rows)} prophages", flush=True)

    rows = [r for stem, _ in genomes for r in by_stem.get(stem, [])]

    print("Merging tables ...")
    tables_mod.merge_summary(rows, table_dir)
    tables_mod.expand_genes(table_dir, prodigal_dir, dom_dir, dict(genomes))
    tables_mod.extract_phvr(table_dir)

    if n_fail:
        print(f"WARNING: {n_fail} genome(s) failed; see messages above", file=sys.stderr)
    print(f"Done. {len(rows)} filamentous prophages -> {out_dir}")


def _split_faa(faa, chunk_dir, chunk_size):
    os.makedirs(chunk_dir, exist_ok=True)
    chunks, n, fh = [], 0, None
    try:
        with open(faa) as src:
            for line in src:
                if line.startswith(">"):
                    n += 1
                    if (n - 1) % chunk_size == 0:
                        if fh:
                            fh.close()
                        path = os.path.join(chunk_dir,
                                            f"chunk_{(n - 1) // chunk_size + 1:04d}.faa")
                        chunks.append(path)
                        fh = open(path, "w")
                if fh:
                    fh.write(line)
    finally:
        if fh:
            fh.close()
    return chunks


def _merge_defense_tables(out_dir):
    merged = []
    for kind in ("hmmer", "genes", "systems"):
        parts = sorted(glob.glob(os.path.join(out_dir, "**",
                                              f"*_defense_finder_{kind}.tsv"),
                                 recursive=True))
        if not parts:
            continue
        dest = os.path.join(out_dir, f"defense_finder_{kind}.tsv")
        with open(dest, "w") as fo:
            with open(parts[0]) as fi:
                fo.write(fi.readline())
            for p in parts:
                with open(p) as fi:
                    fi.readline()
                    shutil.copyfileobj(fi, fo)
        merged.append(os.path.basename(dest))
    return merged


def run_defense(input_dir, out_dir, threads=8, chunk_size=5000):
    faa = os.path.join(input_dir, "tables", "FPphi_pHVR.faa")
    if not os.path.exists(faa):
        raise SystemExit(f"Error: {faa} not found (run `inomine run` first)")
    exe = shutil.which("defense-finder")
    if not exe:
        raise SystemExit("Error: defense-finder not found in PATH.\n"
                         "Install it with: conda install -c bioconda defensefinder")

    chunks = _split_faa(faa, os.path.join(out_dir, "chunks"), chunk_size)
    print(f"DefenseFinder: {len(chunks)} chunk(s) from {faa}")
    for chunk in chunks:
        name = os.path.splitext(os.path.basename(chunk))[0]
        cdir = os.path.join(out_dir, name)
        if os.path.exists(cdir):
            continue
        print(f"  {name} ...", flush=True)
        proc = subprocess.run([exe, "run", "-w", str(threads), "-a", "-o", cdir, chunk],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise SystemExit(f"defense-finder failed on {name}: "
                             f"{proc.stderr.strip()[-300:]}")

    merged = _merge_defense_tables(out_dir)
    if merged:
        print(f"Done. Merged tables in {out_dir}: {', '.join(merged)}")
    else:
        print(f"Done. No result tables found under {out_dir}")
