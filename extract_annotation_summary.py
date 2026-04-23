#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import re

PROJECT = Path("/home/benla/life748report")
ANNOT_DIR = PROJECT / "annotations"
OUT_CSV = PROJECT / "annotation_summary_real.csv"


def fasta_stats(path: Path):
    lengths = []
    current = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
            else:
                current += len(line)
        if current > 0:
            lengths.append(current)
    return len(lengths), sum(lengths)


def parse_runtime_from_log(log_file: Path):
    try:
        lines = log_file.read_text(errors="ignore").splitlines()
        times = []
        for line in lines:
            m = re.search(r"\[(\d{2}:\d{2}:\d{2})\]", line)
            if m:
                times.append(m.group(1))
        if len(times) >= 2:
            from datetime import datetime
            t0 = datetime.strptime(times[0], "%H:%M:%S")
            t1 = datetime.strptime(times[-1], "%H:%M:%S")
            runtime_min = (t1 - t0).total_seconds() / 60.0
            if runtime_min >= 0:
                return runtime_min
    except Exception:
        pass
    return None


def find_first_existing(directory: Path, patterns):
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def parse_prokka(prokka_dir: Path):
    tsv_file = find_first_existing(prokka_dir, ["*.tsv"])
    fna_file = find_first_existing(prokka_dir, ["*.fna", "*.fsa", "*.fa", "*.fasta"])
    log_file = find_first_existing(prokka_dir, ["*.log"])

    if tsv_file is None:
        raise FileNotFoundError(f"No Prokka TSV file found in {prokka_dir}")

    tsv = pd.read_csv(tsv_file, sep="\t")
    tsv.columns = [str(c).strip() for c in tsv.columns]

    possible_ftype = [c for c in tsv.columns if c.lower() in ["ftype", "type", "feature", "feature_type"]]
    if not possible_ftype:
        print("Prokka columns found:", list(tsv.columns))
        raise ValueError(f"Could not find feature-type column in {tsv_file}")
    fcol = possible_ftype[0]

    possible_product = [c for c in tsv.columns if c.lower() == "product"]
    pcol = possible_product[0] if possible_product else None

    feature_series = tsv[fcol].astype(str).str.strip().str.lower()

    n_cds = int((feature_series == "cds").sum())
    n_rrna = int((feature_series == "rrna").sum())
    n_trna = int((feature_series == "trna").sum())

    frac_function = None
    if pcol is not None:
        cds = tsv[feature_series == "cds"].copy()
        if len(cds) > 0:
            cds["has_function"] = (
                cds[pcol].notna()
                & (cds[pcol].astype(str).str.strip().str.lower() != "hypothetical protein")
            )
            frac_function = float(cds["has_function"].mean())

    runtime_min = parse_runtime_from_log(log_file) if log_file else None

    total_bp = None
    n_contigs = None
    if fna_file:
        n_contigs, total_bp = fasta_stats(fna_file)

    return {
        "assembler": "Flye",
        "annotator": "Prokka",
        "total_bp": total_bp,
        "n_contigs": n_contigs,
        "n_CDS": n_cds,
        "n_rRNA": n_rrna,
        "n_tRNA": n_trna,
        "frac_function": frac_function,
        "runtime_min": runtime_min,
    }


def choose_bakta_main_tsv(bakta_dir: Path):
    candidates = sorted(bakta_dir.glob("*.tsv"))
    filtered = []
    for f in candidates:
        name = f.name.lower()
        if "hypothet" in name:
            continue
        if "inference" in name:
            continue
        if "plot" in name:
            continue
        filtered.append(f)

    if not filtered:
        raise FileNotFoundError(
            f"No suitable Bakta main TSV found in {bakta_dir}. Found: {[f.name for f in candidates]}"
        )

    filtered = sorted(filtered, key=lambda x: (len(x.name), x.name))
    return filtered[0]


def parse_bakta(bakta_dir: Path):
    tsv_file = choose_bakta_main_tsv(bakta_dir)
    fasta_file = find_first_existing(bakta_dir, ["*.fna", "*.fasta", "*.fa"])
    log_file = find_first_existing(bakta_dir, ["*.log"])

    # Bakta TSV in your folder has no header row
    tsv = pd.read_csv(
        tsv_file,
        sep="\t",
        comment="#",
        header=None,
        names=["contig", "type", "start", "end", "strand", "locus", "gene", "product", "dbxref"],
    )

    feature_series = tsv["type"].astype(str).str.strip().str.lower()

    n_cds = int((feature_series == "cds").sum())
    n_rrna = int((feature_series == "rrna").sum())
    n_trna = int((feature_series == "trna").sum())

    frac_function = None
    cds = tsv[feature_series == "cds"].copy()
    if len(cds) > 0:
        cds["has_function"] = (
            cds["product"].notna()
            & (cds["product"].astype(str).str.strip().str.lower() != "hypothetical protein")
        )
        frac_function = float(cds["has_function"].mean())

    runtime_min = parse_runtime_from_log(log_file) if log_file else None

    total_bp = None
    n_contigs = None
    if fasta_file:
        n_contigs, total_bp = fasta_stats(fasta_file)

    return {
        "assembler": "Flye",
        "annotator": "Bakta",
        "total_bp": total_bp,
        "n_contigs": n_contigs,
        "n_CDS": n_cds,
        "n_rRNA": n_rrna,
        "n_tRNA": n_trna,
        "frac_function": frac_function,
        "runtime_min": runtime_min,
    }


def main():
    rows = []

    prokka_dir = ANNOT_DIR / "prokka_flye"
    bakta_dir = ANNOT_DIR / "bakta_flye"

    rows.append(parse_prokka(prokka_dir))
    rows.append(parse_bakta(bakta_dir))

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)

    print("\nAnnotation summary:")
    print(df.to_string(index=False))
    print(f"\nSaved to: {OUT_CSV}")


if __name__ == "__main__":
    main()
