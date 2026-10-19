# IMEXscan

Mining filamentous (Inovirus) prophages and their pHVR regions from bacterial
genomes.

> Status: v0.1.0 released. Bioconda packaging is in progress; until then use
> the development installation below.

## Workflow

Given a directory of bacterial genome assemblies, IMEXscan runs:

1. **Gene calling** — Prodigal (metagenomic mode) with strain-prefixed gene IDs.
2. **Marker search** — hmmsearch with 16 curated inovirus marker HMMs
   (Zot, RstA/RstB, Ace, G3P/G8P/G2P/G5P, CtxA/CtxB, RstR, XafT, ...) plus a
   scan for the *Vibrio* dif-like *att* motif.
3. **Prophage extraction** — marker genes are clustered along contigs, clusters
   are split at dif-like sites, and clusters carrying ≥ 4 inovirus hallmark
   signals are reported as filamentous prophages. Boundaries are inferred from
   the flanking dif-like *attL/attR* sites when available.
4. **Tables** — per-prophage summary, per-gene annotation, and the pHVR (genes
   between *zot* and *xaft*) with protein sequences.
5. **Defense annotation** (optional) — DefenseFinder on pHVR proteins.

## Requirements

- Python ≥ 3.9
- [`prodigal`](https://github.com/hyattpd/Prodigal) and
  [`hmmsearch`](http://hmmer.org/) in `PATH`
- Optional: [DefenseFinder](https://github.com/mdmparis/defense-finder) for
  `imexscan defense`

## Installation

Development install:

```bash
git clone https://github.com/awingzhong/imexscan.git
cd imexscan
conda create -n imexscan -c conda-forge -c bioconda python=3.11 prodigal hmmer
conda activate imexscan
pip install -e .
```

After the Bioconda release:

```bash
conda create -n imexscan -c conda-forge -c bioconda imexscan
```

## Usage

```bash
imexscan run -i genomes/ -o results/ -t 16
```

Input genomes may be `.fa`, `.fasta`, `.fna`, optionally gzipped. Processing is
parallel across genomes; existing intermediate files are reused on re-runs
(resume-friendly).

Optional downstream defense-system annotation of pHVR proteins:

```bash
imexscan defense -i results/ -o defense_out/
```

Scan only the dif-like *att* motif of one genome:

```bash
imexscan dif -i genome.fna -o genome-dif.csv
```

### Options for `imexscan run`

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | required | genome directory |
| `-o, --output` | required | output directory |
| `-t, --threads` | 4 | genomes processed in parallel |
| `--hmm` | bundled | marker HMM directory |
| `--min-hallmark` | 4 | minimum inovirus hallmark signals |
| `--evalue` | 1e-5 | hmmsearch E-value cutoff |
| `--skip-dif` | off | skip dif scan (boundaries default to core ± 5 kb) |
| `--prodigal` | `prodigal` | prodigal executable |
| `--hmmsearch` | `hmmsearch` | hmmsearch executable |

## Output

| File | Content |
|------|---------|
| `results/FPphi/*.fasta` | extracted filamentous prophage sequences |
| `tables/FPphi_summary.tsv` | one row per prophage (coordinates, *attL/attR*, quality, hallmarks) |
| `tables/FPphi_gene_annotation.tsv` | all genes inside each prophage |
| `tables/FPphi_pHVR_annotation.tsv` | genes between *zot* and *xaft* (pHVR) |
| `tables/FPphi_pHVR.faa` | pHVR protein sequences |
| `prodigal/`, `dom/`, `dif/` | intermediates (reused on re-runs) |

## Bundled marker HMMs

The 16 marker profiles in `src/imexscan/data/Vibrio_Inovirus/` comprise 14
profiles derived from [InterPro](https://www.ebi.ac.uk/interpro/)/Pfam
(CC0 1.0) and 2 self-built profiles. InterPro/Pfam downloadable data are
redistributed under the
[CC0 1.0 Universal](https://interpro-documentation.readthedocs.io/en/latest/license.html)
dedication.

## Citation

Citation information will be added upon publication. The dif-like *att* motif
pattern builds on the *Vibrio* dif motif analysis in this project.

## License

MIT
