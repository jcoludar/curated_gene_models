"""Build the ant benchmark dataset (locus-window FASTA + GFF3 + genes TSV) from the Ant_Venoms release.

Usage:
  python _crescendo_build_ant_20261006.py --release DIR --genome-roots DIR [DIR ...] --dataset ID
      --paper DOI --out-dir DIR [--flank 50000] [--merge-gap 100000]

Inputs (read-only)
  release/gff/ALL_SPECIES.gff3            functional gene -> mRNA -> CDS (genome-true coordinates)
                                          header lines '# assembly: <species> <ACC> (<name>)'
  release/EXPORT_MANIFEST.tsv             gm, species, family, is_venom, structural_status, roundtrip
  release/proteins/ant_venom_proteins.faa reference proteins (checked against our re-translation)
  genome-roots                            searched recursively for '<ACC>_*genomic.fna' (or '<ACC>_*.fna')
Tier mapping from structural_status:
  complete-functional -> complete ; partial-functional / assembly-fragmented / fragment /
  flagged-not-counted -> partial ; anything pseudogene-like -> pseudogene (only functional GFF genes are read).
Windows: per (assembly, scaffold) cluster genes (gap <= merge-gap), window = cluster +- flank, clipped.
FASTA seqid = '<scaffold>:<start>-<end>'. Uses the .fai next to the genome when present, else
streams the FASTA once and keeps only the scaffolds needed. Round trip: re-splice from the written
window, translate, compare with the release protein (EQUAL / EQUAL_TRIMMED_STOP / MISMATCH / NO_REF).
Writes only into --out-dir (refuses if non-empty).
"""
import argparse
import collections
import csv
import glob
import os
import re
import sys

from Bio.Seq import Seq

COMP = str.maketrans("ACGTNRYKMSWBDHVacgtnrykmswbdhv", "TGCANYRMKSWVHDBtgcanyrmkswvhdb")
TIER = {"complete-functional": "complete", "partial-functional": "partial", "assembly-fragmented": "partial",
        "fragment": "partial", "flagged-not-counted": "partial"}


def rc(s):
    return s.translate(COMP)[::-1]


def find_genome(roots, acc):
    for root in roots:
        hits = glob.glob(os.path.join(root, "**", f"{acc}_*genomic.fna"), recursive=True)
        hits += glob.glob(os.path.join(root, "**", f"{acc}_*.fna"), recursive=True)
        hits = [h for h in hits if os.path.getsize(h) > 0]
        if hits:
            return sorted(hits, key=len)[0]
    return None


def fetch_scaffolds(path, wanted):
    """Return {scaffold: sequence} for the wanted scaffolds (uppercased)."""
    out = {}
    fai = path + ".fai"
    if os.path.exists(fai):
        idx = {}
        for line in open(fai, encoding="utf-8"):
            f = line.split("\t")
            idx[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        with open(path, "rb") as fh:
            for name in wanted:
                if name not in idx:
                    continue
                length, offset, lb, lw = idx[name]
                nlines = (length + lb - 1) // lb
                fh.seek(offset)
                raw = fh.read(nlines * lw)
                seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:length]
                out[name] = seq.decode("ascii").upper()
        return out
    cur, buf = None, []
    with open(path, encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in wanted:
                    out[cur] = "".join(buf).upper()
                cur, buf = line[1:].split()[0], []
                if len(out) == len(wanted):
                    break
            elif cur in wanted:
                buf.append(line.strip())
        if cur in wanted and cur not in out:
            out[cur] = "".join(buf).upper()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", required=True)
    ap.add_argument("--genome-roots", nargs="+", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--paper", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--merge-gap", type=int, default=100000)
    a = ap.parse_args()
    if os.path.isdir(a.out_dir) and os.listdir(a.out_dir):
        sys.exit(f"refusing: {a.out_dir} is not empty")
    os.makedirs(a.out_dir, exist_ok=True)

    man = {r["gm"]: r for r in csv.DictReader(open(os.path.join(a.release, "EXPORT_MANIFEST.tsv"),
                                                    encoding="utf-8"), delimiter="\t")}
    ref = {}
    cur = None
    for line in open(os.path.join(a.release, "proteins", "ant_venom_proteins.faa"), encoding="utf-8"):
        if line.startswith(">"):
            m = re.search(r"(GM_\d+)", line)
            cur = m.group(1) if m else None
            ref[cur] = []
        elif cur:
            ref[cur].append(line.strip())
    ref = {k: "".join(v).upper() for k, v in ref.items()}

    species_acc, genes, mrna_gene, cds = {}, {}, {}, collections.defaultdict(list)
    phase_of = {}
    for line in open(os.path.join(a.release, "gff", "ALL_SPECIES.gff3"), encoding="utf-8"):
        m = re.match(r"# assembly: (.+?) (GC[AF]_\d+\.\d+) \((.+)\)", line)
        if m:
            species_acc[m.group(1)] = (m.group(2), m.group(3))
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            genes[at["ID"]] = {"seqid": f[0], "strand": f[6], "species": at.get("species", ""),
                               "family": at.get("family", ""), "status": at.get("status", "")}
        elif f[2] == "mRNA":
            mrna_gene[at["ID"]] = at["Parent"]
        elif f[2] == "CDS":
            gid_ = mrna_gene[at["Parent"]]
            cds[gid_].append((int(f[3]), int(f[4])))
            phase_of[(gid_, int(f[3]), int(f[4]))] = int(f[7]) if f[7] in "012" else 0

    by_asm = collections.defaultdict(lambda: collections.defaultdict(list))
    missing_asm = collections.Counter()
    for gid, g in genes.items():
        acc = species_acc.get(g["species"], (None, None))[0]
        if not acc:
            missing_asm[g["species"]] += 1
            continue
        g["acc"], g["asm_name"] = species_acc[g["species"]]
        g["cds"] = sorted(cds[gid])
        by_asm[acc][g["seqid"]].append(gid)

    fasta_out = open(os.path.join(a.out_dir, f"{a.dataset}.fasta"), "w", encoding="utf-8")
    gff_lines = ["##gff-version 3"]
    rows, problems = [], []
    for acc, scaffolds in sorted(by_asm.items()):
        path = find_genome(a.genome_roots, acc)
        if not path:
            problems.append(f"NO_GENOME {acc} genes={sum(len(v) for v in scaffolds.values())}")
            continue
        seqs = fetch_scaffolds(path, set(scaffolds))
        for scaf, gids in sorted(scaffolds.items()):
            if scaf not in seqs:
                problems.append(f"NO_SCAFFOLD {acc} {scaf} genes={len(gids)}")
                continue
            pseq = seqs[scaf]
            gids.sort(key=lambda x: genes[x]["cds"][0][0])
            clusters, cl = [], [gids[0]]
            for x in gids[1:]:
                if genes[x]["cds"][0][0] - max(genes[y]["cds"][-1][1] for y in cl) <= a.merge_gap:
                    cl.append(x)
                else:
                    clusters.append(cl)
                    cl = [x]
            clusters.append(cl)
            for cl in clusters:
                ws = max(1, min(genes[x]["cds"][0][0] for x in cl) - a.flank)
                we = min(len(pseq), max(genes[x]["cds"][-1][1] for x in cl) + a.flank)
                sid = f"{scaf}:{ws}-{we}"
                wseq = pseq[ws - 1:we]
                fasta_out.write(f">{sid} {acc}\n")
                for i in range(0, len(wseq), 80):
                    fasta_out.write(wseq[i:i + 80] + "\n")
                gff_lines.append(f"##sequence-region {sid} 1 {len(wseq)}")
                for gid in cl:
                    g, mf = genes[gid], man.get(gid, {})
                    off = ws - 1
                    c = [(s - off, e - off) for s, e in g["cds"]]
                    tier = TIER.get(g["status"], "pseudogene")
                    attrs = (f"Name={gid};species={g['species'].replace(' ', '_')};family={g['family']};"
                             f"is_venom={mf.get('is_venom', '')};structural_status={g['status']};tier={tier};"
                             f"genome={acc};scaffold={scaf};genome_start={g['cds'][0][0]};"
                             f"genome_end={g['cds'][-1][1]};source={a.paper}")
                    gff_lines.append(f"{sid}\tcurated\tgene\t{c[0][0]}\t{c[-1][1]}\t.\t{g['strand']}\t.\tID=gene-{gid};{attrs}")
                    gff_lines.append(f"{sid}\tcurated\tmRNA\t{c[0][0]}\t{c[-1][1]}\t.\t{g['strand']}\t.\tID=mrna-{gid};Parent=gene-{gid};Name={gid};tier={tier}")
                    # phases are carried from the release GFF (a 5'-partial gene starts off-frame)
                    phase = {(s - off, e - off): phase_of[(gid, s, e)] for s, e in g["cds"]}
                    first = c[0] if g["strand"] == "+" else c[-1]
                    p0 = phase[first]
                    for k, (s, e) in enumerate(c, 1):
                        gff_lines.append(f"{sid}\tcurated\texon\t{s}\t{e}\t.\t{g['strand']}\t.\tID=exon-{gid}-{k};Parent=mrna-{gid}")
                    for s, e in c:
                        gff_lines.append(f"{sid}\tcurated\tCDS\t{s}\t{e}\t.\t{g['strand']}\t{phase[(s, e)]}\tID=cds-{gid};Parent=mrna-{gid}")
                    nt = "".join(wseq[s - 1:e] for s, e in c)
                    if g["strand"] == "-":
                        nt = rc(nt)
                    core = nt[p0:]
                    prot = str(Seq(core[: len(core) - len(core) % 3]).translate())
                    body = prot[:-1] if prot.endswith("*") else prot
                    r = ref.get(gid, "")
                    rt = ("NO_REF" if not r else "EQUAL" if body == r else
                          "EQUAL_TRIMMED_STOP" if body.rstrip("*") == r.rstrip("*") else "MISMATCH")
                    rows.append([a.dataset, f"gene-{gid}", gid, g["species"], g["family"], mf.get("is_venom", ""),
                                 g["status"], tier, acc, g["asm_name"], scaf, str(g["cds"][0][0]),
                                 str(g["cds"][-1][1]), g["strand"], str(len(c)), str(len(nt)), str(len(body)),
                                 rt, sid, body])
    fasta_out.close()
    with open(os.path.join(a.out_dir, f"{a.dataset}.gff3"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(gff_lines) + "\n")
    with open(os.path.join(a.out_dir, f"{a.dataset}.genes.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(["dataset", "gene_id", "name", "species", "family", "is_venom", "structural_status",
                            "tier", "assembly", "assembly_name", "scaffold", "genome_start", "genome_end",
                            "strand", "n_cds", "cds_nt", "protein_aa", "vs_release_protein", "window",
                            "protein"]) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"genes_in_gff={len(genes)} written={len(rows)} tiers={dict(collections.Counter(r[7] for r in rows))}")
    print(f"vs_release_protein={dict(collections.Counter(r[17] for r in rows))}")
    print(f"assemblies={len(by_asm)} missing_assembly_species={dict(missing_asm)}")
    for p in problems:
        print(p)


if __name__ == "__main__":
    main()
