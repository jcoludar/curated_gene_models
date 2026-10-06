"""Build a benchmark dataset (locus-window FASTA + GFF3 + genes TSV) from Geneious exports.

Usage:
  python _crescendo_build_dataset_20261006.py --exports DIR --locate TSV --records DIR
      --published NAMES_TXT --norm {3ftx,sp} --dataset ID --paper DOI --out-dir DIR
      [--flank 50000] [--merge-gap 100000] [--max-shift 5000]

Inputs
  exports    Geneious <doc>.gff + <doc>.fasta (USER export)
  locate     output of _crescendo_locate_docs_20261006.py (doc -> parent, strand, start)
  records    NCBI parent FASTAs (<parent>.fasta)
  published  one published gene name per line (SD1 identifiers / SP final-tree names)

Per curated CDS group (Name without NCBI suffix, shared ID):
  1. normalise the name; keep it only if it is in the published list (else: excluded, reason)
  2. place EVERY CDS interval on the parent by exact match of the interval's sequence,
     searched within +-max-shift of the position projected from the doc offset
     (absorbs small indels between the Geneious copy and the current NCBI version);
     an interval that is not found exactly and uniquely -> gene excluded (reason recorded)
  3. drop duplicate models (identical genomic CDS sets) keeping the first published name
  4. tier: complete (M..stop, no internal stop, len%3==0) / partial (no internal stop, but
     missing M or stop) / pseudogene (psi-named, internal stop, or len%3!=0)
Then genes are clustered per parent (gap <= merge-gap), each cluster becomes a window
[min-flank, max+flank] clipped to the record; FASTA seqid = '<acc.ver>:<start>-<end>'.
GFF3 (window coordinates): gene / mRNA / exon / CDS (exon == CDS, no UTRs known).
Round trip: every CDS is re-spliced FROM THE WRITTEN FASTA+GFF3 and re-translated; it must
equal the translation from the Geneious doc. Writes only into --out-dir (refuses if the
out-dir already has files).
"""
import argparse
import collections
import os
import re
import sys
import urllib.parse

from Bio.Seq import Seq

NCBI_SUFFIXES = (" CDS", " mRNA", " gene")
COMP = str.maketrans("ACGTNRYKMSWBDHV", "TGCANYRMKSWVHDB")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta_one(path):
    hdr, buf = "", []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr:
                    break
                hdr = line[1:].strip()
            else:
                buf.append(line.strip())
    return hdr, "".join(buf).upper()


def find_all(hay, needle):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def norm_name(name, mode):
    if mode == "sp":
        return name.strip().split(" ")[0].replace("_", "")
    return name.strip().replace(" ", "_")


def translate(nt):
    core = nt[: len(nt) - len(nt) % 3]
    return str(Seq(core).translate())


def tier_of(name, nt, prot):
    body = prot[:-1] if prot.endswith("*") else prot
    if name.lower().startswith("psi") or " psi" in name.lower() or "*" in body or len(nt) % 3:
        return "pseudogene"
    if body.startswith("M") and prot.endswith("*"):
        return "complete"
    return "partial"


def load_groups(gff_path):
    groups = collections.defaultdict(list)
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            name = urllib.parse.unquote(attrs.get("Name", "?"))
            if name.endswith(NCBI_SUFFIXES) or name == "?":
                continue
            groups[(attrs.get("ID", name), name)].append((int(f[3]), int(f[4]), f[6]))
    # Geneious exports can repeat a CDS line verbatim (PLA2: every interval twice) -> de-duplicate
    return {k: sorted(set(v)) for k, v in groups.items()}


def main():
    ap = argparse.ArgumentParser()
    for a_ in ("exports", "locate", "records", "norm", "dataset", "paper", "out_dir"):
        ap.add_argument("--" + a_.replace("_", "-"), required=True)
    ap.add_argument("--published", default=None, help="published gene names, one per line")
    ap.add_argument("--published-fasta", default=None,
                    help="published proteins; a model is published if its protein equals / contains / is "
                         "contained (>=90%% length) in one of them (gate by SEQUENCE, survives renames)")
    ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--merge-gap", type=int, default=100000)
    ap.add_argument("--max-shift", type=int, default=5000)
    ap.add_argument("--min-contain-frac", type=float, default=0.9,
                    help="published-vs-model containment must cover >= this fraction of the longer one "
                         "(lower it when the published set is MATURE/trimmed proteins, e.g. PLA2 SM9)")
    ap.add_argument("--gate-aln", action="store_true",
                    help="fallback gate: local alignment with 100%% identity over >=95%% of a published "
                         "protein (>=60 aa); for published sets that are trimmed alignment rows")
    ap.add_argument("--extra", help="USER-approved names NOT in the published table (flagged in_published_table=no)")
    ap.add_argument("--alias", help="TSV old_name<TAB>published_name (same model, renamed at publication)")
    a = ap.parse_args()
    extra = {l.strip() for l in open(a.extra, encoding="utf-8") if l.strip()} if a.extra else set()
    alias = {}
    if a.alias:
        for l in open(a.alias, encoding="utf-8"):
            if l.strip():
                k, v = l.rstrip("\n").split("\t")
                alias[k] = v
    if os.path.isdir(a.out_dir) and os.listdir(a.out_dir):
        sys.exit(f"refusing: {a.out_dir} is not empty")
    os.makedirs(a.out_dir, exist_ok=True)
    published = ({l.strip() for l in open(a.published, encoding="utf-8") if l.strip()}
                 if a.published else set())
    pub_seqs = {}
    if a.published_fasta:
        cur = None
        for l in open(a.published_fasta, encoding="utf-8", errors="replace"):
            if l.startswith(">"):
                cur = l[1:].strip().split()[0]
                pub_seqs[cur] = []
            elif cur:
                pub_seqs[cur].append(l.strip().upper().replace("-", ""))
        pub_seqs = {k: "".join(v).rstrip("*") for k, v in pub_seqs.items() if v}

    kmer_index = collections.defaultdict(set)
    if a.gate_aln:
        for k_, p_ in pub_seqs.items():
            for i in range(len(p_) - 9):
                kmer_index[p_[i:i + 10]].add(k_)
        from Bio import Align
        from Bio.Align import substitution_matrices
        aligner = Align.PairwiseAligner(mode="local", open_gap_score=-10, extend_gap_score=-0.5)
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aln_stats = collections.Counter()

    def aln_match(body):
        """100 % identity over >= 95 % of a published protein (>= 60 aa), gaps allowed but counted."""
        cands = set()
        for i in range(len(body) - 9):
            cands |= kmer_index.get(body[i:i + 10], set())
        best = ("", 0, 0)
        for k_ in cands:
            p_ = pub_seqs[k_]
            al_ = aligner.align(body.replace("*", "X"), p_)[0]
            ident = sum(1 for (a1, a2), (b1, b2) in zip(*al_.aligned) for j in range(a2 - a1)
                        if body[a1 + j] == p_[b1 + j])
            alen = sum(a2 - a1 for a1, a2 in al_.aligned[0])
            cov = sum(b2 - b1 for b1, b2 in al_.aligned[1]) / len(p_)
            if ident == alen and alen >= 60 and cov >= 0.95 and alen > best[1]:
                best = (k_, alen, len(al_.aligned[0]))
        if best[0]:
            aln_stats["gapped" if best[2] > 1 else "ungapped"] += 1
        return best[0]

    def seq_match(body):
        if not body:
            return ""
        for k, p in pub_seqs.items():
            if body == p:
                return k
        for k, p in pub_seqs.items():
            f = a.min_contain_frac
            if (len(p) >= 30 and p in body and len(p) >= f * len(body)) or \
               (len(body) >= 30 and body in p and len(body) >= f * len(p)):
                return k
        return aln_match(body) if a.gate_aln else ""

    loc = {}
    with open(a.locate, encoding="utf-8") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            loc[r["doc"]] = r

    genes, excluded = [], []
    parent_seqs = {}
    for doc, r in sorted(loc.items()):
        if r["status"] == "NO_PARENT_RECORD" or not r["parent"]:
            excluded.append([doc, "*", f"doc not placed ({r['status']})"])
            continue
        # EXACT -> projected local search; otherwise (assembly version differs) -> global search
        global_mode = r["status"] != "EXACT"
        doc_fa = os.path.join(a.exports, doc + ".fasta")
        if not os.path.exists(doc_fa) and doc == r["parent"]:
            doc_fa = os.path.join(a.records, doc + ".fasta")  # doc IS the record (identity locate)
        _, dseq = read_fasta_one(doc_fa)
        if r["parent"] not in parent_seqs:
            phdr, pseq = read_fasta_one(os.path.join(a.records, r["parent"] + ".fasta"))
            parent_seqs[r["parent"]] = (phdr.split()[0], phdr, pseq, rc(pseq))
        pacc, phdr, pseq, pseq_rc = parent_seqs[r["parent"]]
        dstrand = r["strand"] or "+"
        dstart = int(r["start_1based"]) if r["start_1based"] else 1
        dlen = len(dseq)
        for (gid, name), ivs in load_groups(os.path.join(a.exports, doc + ".gff")).items():
            nname = norm_name(name, a.norm)
            nname = alias.get(nname, nname)
            in_pub = nname in published
            pub_id = nname if in_pub else ""
            if pub_seqs and not in_pub:
                ivs_ = sorted(ivs)
                nt_ = "".join(dseq[s - 1:e] for s, e, _ in ivs_)
                if ivs_ and ivs_[0][2] == "-":
                    nt_ = rc(nt_)
                pr_ = translate(nt_)
                pub_id = seq_match(pr_[:-1] if pr_.endswith("*") else pr_)
                in_pub = bool(pub_id)
            if not in_pub and nname not in extra:
                excluded.append([doc, name, "not in published list"])
                continue
            ivs.sort()
            strands = {s for _, _, s in ivs}
            if len(strands) != 1:
                excluded.append([doc, name, "mixed strand"])
                continue
            gstrand_doc = strands.pop()
            # transcript-order nt from the doc
            seg_doc = [dseq[s - 1:e] for s, e, _ in ivs]
            nt_doc = "".join(seg_doc)
            if gstrand_doc == "-":
                nt_doc = rc(nt_doc)
            placed, bad = [], ""
            if global_mode:
                orient = set()
                for (s, e, _), seg in zip(ivs, seg_doc):
                    hit = None
                    for flank in (0, 30, 100, 300):
                        q = dseq[max(0, s - 1 - flank):min(dlen, e + flank)]
                        lead = (s - 1) - max(0, s - 1 - flank)
                        fw = find_all(pseq, q)
                        rv = find_all(pseq_rc, q)
                        if len(fw) + len(rv) == 1:
                            if fw:
                                hit = ("+", fw[0] + lead + 1)
                            else:
                                # position on the rc string -> forward coordinate of the segment
                                rc_start = rv[0] + lead
                                hit = ("-", len(pseq) - (rc_start + (e - s + 1)) + 1)
                            break
                        if len(fw) + len(rv) == 0:
                            break
                    if hit is None:
                        bad = f"interval {s}-{e}: no unique exact genome-wide hit"
                        break
                    orient.add(hit[0])
                    placed.append((hit[1], hit[1] + (e - s)))
                if not bad and len(orient) != 1:
                    bad = "intervals map to both parent strands"
                if not bad and max(p[1] for p in placed) - min(p[0] for p in placed) > 300000:
                    bad = "placed intervals span >300 kb"
                if bad:
                    excluded.append([doc, name, bad])
                    continue
                dstrand = orient.pop()
                gstrand = gstrand_doc if dstrand == "+" else ("-" if gstrand_doc == "+" else "+")
                placed.sort()
                prot = translate(nt_doc)
                genes.append({
                    "name": name.strip(), "norm": nname, "doc": doc, "parent": r["parent"], "acc": pacc,
                    "species": phdr.split(" ", 1)[1] if " " in phdr else "",
                    "strand": gstrand, "cds": placed, "nt": nt_doc, "prot": prot,
                    "tier": tier_of(name, nt_doc, prot), "in_pub": in_pub, "pub_id": pub_id, "placement": "global",
                })
                continue
            for (s, e, _), seg in zip(ivs, seg_doc):
                found = None
                for flank in (0, 30, 100, 300):
                    xs, xe = max(1, s - flank), min(dlen, e + flank)
                    ext = dseq[xs - 1:xe]
                    if dstrand == "+":
                        exp, q, lead = dstart + xs - 1, ext, s - xs
                    else:
                        exp, q, lead = dstart + (dlen - xe), rc(ext), xe - e
                    lo = max(0, exp - 1 - a.max_shift)
                    hi = min(len(pseq), exp - 1 + len(q) + a.max_shift)
                    hits = [lo + i + 1 for i in find_all(pseq[lo:hi], q)]
                    if len(hits) == 1:
                        found = hits[0] + lead
                        break
                    if not hits:
                        break
                if found is None:
                    bad = f"interval {s}-{e}: no unique exact hit near {exp}"
                    break
                placed.append((found, found + (e - s)))
            if bad:
                excluded.append([doc, name, bad])
                continue
            gstrand = gstrand_doc if dstrand == "+" else ("-" if gstrand_doc == "+" else "+")
            placed.sort()
            prot = translate(nt_doc)
            genes.append({
                "name": name.strip(), "norm": nname, "doc": doc, "parent": r["parent"], "acc": pacc,
                "species": phdr.split(" ", 1)[1] if " " in phdr else "",
                "strand": gstrand, "cds": placed, "nt": nt_doc, "prot": prot,
                "tier": tier_of(name, nt_doc, prot), "in_pub": in_pub, "pub_id": pub_id, "placement": "local",
            })

    # dedup identical genomic models
    seen, uniq = {}, []
    for g in genes:
        key = (g["acc"], g["strand"], tuple(g["cds"]))
        if key in seen:
            excluded.append([g["doc"], g["name"], f"duplicate of {seen[key]}"])
            continue
        seen[key] = g["name"]
        uniq.append(g)
    genes = uniq

    # windows
    by_parent = collections.defaultdict(list)
    for g in genes:
        by_parent[g["parent"]].append(g)
    fasta_out = open(os.path.join(a.out_dir, f"{a.dataset}.fasta"), "w", encoding="utf-8")
    gff_lines = ["##gff-version 3"]
    windows = {}
    for parent, gs in sorted(by_parent.items()):
        pacc, _, pseq, _ = parent_seqs[parent]
        gs.sort(key=lambda g: g["cds"][0][0])
        clusters, cur = [], [gs[0]]
        for g in gs[1:]:
            if g["cds"][0][0] - max(x["cds"][-1][1] for x in cur) <= a.merge_gap:
                cur.append(g)
            else:
                clusters.append(cur)
                cur = [g]
        clusters.append(cur)
        for cl in clusters:
            ws = max(1, min(g["cds"][0][0] for g in cl) - a.flank)
            we = min(len(pseq), max(g["cds"][-1][1] for g in cl) + a.flank)
            sid = f"{pacc}:{ws}-{we}"
            seq = pseq[ws - 1:we]
            windows[sid] = seq
            fasta_out.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                fasta_out.write(seq[i:i + 80] + "\n")
            gff_lines.append(f"##sequence-region {sid} 1 {len(seq)}")
            for g in cl:
                g["window"], g["wstart"] = sid, ws
    fasta_out.close()

    rows = []
    used_ids = collections.Counter()
    for g in sorted(genes, key=lambda g: (g["window"], g["cds"][0][0])):
        off = g["wstart"] - 1
        cds = [(s - off, e - off) for s, e in g["cds"]]
        base = g["norm"] if (g["in_pub"] and g["norm"] == g["pub_id"]) or not g["pub_id"] else g["pub_id"]
        gid = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{a.dataset}:{base}")
        used_ids[gid] += 1
        if used_ids[gid] > 1:  # names are not unique (generic Geneious names / shared mature peptides)
            gid = f"{gid}.{used_ids[gid]}"
        attrs = (f"Name={g['norm']};original_name={urllib.parse.quote(g['name'])};tier={g['tier']};"
                 f"genome={g['acc']};genome_start={g['cds'][0][0]};genome_end={g['cds'][-1][1]};"
                 f"in_published_table={'yes' if g['in_pub'] else 'no'};"
                 f"published_id={urllib.parse.quote(g['pub_id']) or '.'};"
                 f"placement={g['placement']};source={urllib.parse.quote(a.paper)}")
        gs_, ge_ = cds[0][0], cds[-1][1]
        gff_lines.append(f"{g['window']}\tcurated\tgene\t{gs_}\t{ge_}\t.\t{g['strand']}\t.\tID=gene-{gid};{attrs}")
        gff_lines.append(f"{g['window']}\tcurated\tmRNA\t{gs_}\t{ge_}\t.\t{g['strand']}\t.\tID=mrna-{gid};Parent=gene-{gid};Name={g['norm']};tier={g['tier']}")
        order = cds if g["strand"] == "+" else cds[::-1]
        done = 0
        phases = {}
        for s, e in order:
            phases[(s, e)] = (3 - done % 3) % 3
            done += e - s + 1
        for k, (s, e) in enumerate(cds, 1):
            gff_lines.append(f"{g['window']}\tcurated\texon\t{s}\t{e}\t.\t{g['strand']}\t.\tID=exon-{gid}-{k};Parent=mrna-{gid}")
        for s, e in cds:
            gff_lines.append(f"{g['window']}\tcurated\tCDS\t{s}\t{e}\t.\t{g['strand']}\t{phases[(s, e)]}\tID=cds-{gid};Parent=mrna-{gid}")
        # round trip from the written window sequence
        wseq = windows[g["window"]]
        nt = "".join(wseq[s - 1:e] for s, e in cds)
        if g["strand"] == "-":
            nt = rc(nt)
        rt = "OK" if translate(nt) == g["prot"] else "MISMATCH"
        body = g["prot"][:-1] if g["prot"].endswith("*") else g["prot"]
        rows.append([a.dataset, gid, g["norm"], g["name"], g["species"], g["tier"], g["acc"],
                     str(g["cds"][0][0]), str(g["cds"][-1][1]), g["strand"], str(len(cds)),
                     str(len(g["nt"])), str(len(body)), rt, g["window"], g["placement"],
                     "yes" if g["in_pub"] else "no", g["pub_id"], body])
    with open(os.path.join(a.out_dir, f"{a.dataset}.gff3"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(gff_lines) + "\n")
    with open(os.path.join(a.out_dir, f"{a.dataset}.genes.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(["dataset", "gene_id", "name", "original_name", "genome_record_title", "tier",
                            "genome_acc", "genome_start", "genome_end", "strand", "n_cds", "cds_nt",
                            "protein_aa", "roundtrip", "window", "placement", "in_published_table",
                            "published_id", "protein"]) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(os.path.join(a.out_dir, f"{a.dataset}.excluded.tsv"), "w", encoding="utf-8") as fh:
        fh.write("doc\tname\treason\n")
        for r in excluded:
            fh.write("\t".join(r) + "\n")
    tiers = collections.Counter(r[5] for r in rows)
    reasons = collections.Counter(r[2].split(":")[0].split(" (")[0] for r in excluded)
    print(f"genes={len(rows)} tiers={dict(tiers)} windows={len(windows)} "
          f"window_bp={sum(len(s) for s in windows.values())} roundtrip_ok={sum(r[13] == 'OK' for r in rows)}")
    print(f"excluded={len(excluded)} reasons={dict(reasons)}")
    if a.gate_aln:
        print(f"alignment-gate matches={dict(aln_stats)}")


if __name__ == "__main__":
    main()
