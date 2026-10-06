# curated_gene_models

Hand-curated venom gene models from published work, for benchmarking gene callers.
Most of these models are not in GenBank or RefSeq.

## Content (v1)

| dataset | genes | complete / partial / pseudogene | source |
|---|---|---|---|
| `3ftx` | 163 | 126 / 30 / 7 | Koludarov et al. 2023, Nat Commun, doi:10.1038/s41467-023-40550-0 |
| `sp` | 120 | 115 / 3 / 2 | Barua et al. 2021, BMC Biol, doi:10.1186/s12915-021-01191-1 |
| `pla2` | 429 | 427 / 1 / 1 | Koludarov et al. 2019, bioRxiv, doi:10.1101/583344 |
| `bee` | 154 | 140 / 12 / 2 | Koludarov et al. 2023, BMC Biol, doi:10.1186/s12915-023-01656-5 |
| `ant` | 1,787 | 1,737 / 50 / 0 | Weitz et al. 2026, bioRxiv 2026.02.12.705515 (pre-publication) |

The full dataset is on Zenodo (DOI: TBA). This repository holds everything except the FASTA files.

Each `data/<dataset>/` contains:

- `*.fasta.gz` (Zenodo only): genome windows (gene cluster ± 50 kb), named `<sequence>:<start>-<end>`
- `*.gff3`: gene / mRNA / exon / CDS, with `tier=` and genome coordinates in the attributes
- `*.genes.tsv`: one row per gene: genome accession, coordinates, tier, protein
- `*.excluded.tsv`: models left out, with the reason
- `pla2/*.non_ncbi_sequences.tsv`: species and source of the PLA2 sequences that have no NCBI accession

`MD5SUMS` lists checksums of the uncompressed content.

## Notes

- Models are CDS only (no UTRs). Score with CDS as the exon feature and ignore predictions outside
  annotated genes, e.g. `gene-benchmark run --gt-exon-feature-type CDS -Q`.
- Every protein was checked against the protein sequences published with its paper.
- `ant` is a snapshot of an unpublished curation (export of 2026-09-30). It will change in later versions.
- Bee `Xylocopa violacea` windows come from the assembly in Zenodo doi:10.5281/zenodo.8052397 (Additional
  file 41), not from NCBI.
- `scripts/` holds the build scripts, for provenance. They need inputs that are not in this repository.

## Citation and licence

Cite the version DOI on Zenodo, plus the source paper of each dataset you use. CC BY 4.0.
