#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


PROJECT = Path("/home/benla/life748report")
QUAST_TSV = PROJECT / "quast" / "GN9_compare" / "report.tsv"
FIG_DIR = PROJECT / "figures"

# Edit these if your file names differ
ASSEMBLY_FASTAS = {
    "Flye": PROJECT / "assemblies" / "flye" / "assembly.fasta",
    "SPAdes": PROJECT / "assemblies" / "spades" / "GN9_spades_fast" / "contigs.fasta",
}

# Fill in from your logs / report text
RUNTIME_MIN = {
    "Flye": 20,   # replace with your actual Flye runtime
    "SPAdes": 4,  # from your successful run
}

PEAK_MEM_GB = {
    "Flye": 4.0,   # replace with your actual/estimated Flye peak RAM
    "SPAdes": 2.3, # from log max ~2316M
}


def read_quast_report(tsv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    df = df.rename(columns={df.columns[0]: "metric"})
    df = df.set_index("metric")
    return df


def get_metric(df: pd.DataFrame, candidates: list[str]) -> pd.Series:
    for name in candidates:
        if name in df.index:
            return pd.to_numeric(df.loc[name], errors="coerce")
    raise KeyError(f"Could not find any of these metrics in QUAST report: {candidates}")


def read_fasta_lengths(path: Path) -> list[int]:
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
    return sorted(lengths, reverse=True)


def plot_assembly_metrics(df: pd.DataFrame, outpath: Path) -> None:
    total_len = get_metric(df, ["Total length", "Total length (>= 0 bp)"])
    n_contigs = get_metric(df, ["# contigs", "# contigs (>= 0 bp)"])
    n50 = get_metric(df, ["N50"])
    largest = get_metric(df, ["Largest contig"])
    gc = get_metric(df, ["GC (%)"])

    metrics = {
        "Total length (bp)": total_len,
        "Contigs": n_contigs,
        "N50 (bp)": n50,
        "Largest contig (bp)": largest,
        "GC %": gc,
    }

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()

    for ax, (title, series) in zip(axes, metrics.items()):
        series.plot(kind="bar", ax=ax)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=0)

    axes[-1].axis("off")
    fig.suptitle("Assembly comparison: Flye vs SPAdes", fontsize=14)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_quast_quality(df: pd.DataFrame, outpath: Path) -> None:
    possible = {
        "Genome fraction (%)": ["Genome fraction (%)"],
        "Misassemblies": ["# misassemblies"],
        "Mismatches per 100 kbp": ["# mismatches per 100 kbp"],
        "Indels per 100 kbp": ["# indels per 100 kbp"],
        "Duplication ratio": ["Duplication ratio"],
    }

    found = {}
    for pretty, candidates in possible.items():
        try:
            found[pretty] = get_metric(df, candidates)
        except KeyError:
            pass

    if not found:
        print("No correctness/completeness metrics found in QUAST report; skipping quality plot.")
        return

    fig, axes = plt.subplots(1, len(found), figsize=(4 * len(found), 4))
    if len(found) == 1:
        axes = [axes]

    for ax, (title, series) in zip(axes, found.items()):
        series.plot(kind="bar", ax=ax)
        ax.set_title(title)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=0)

    fig.suptitle("Assembly completeness / correctness metrics", fontsize=14)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_runtime_memory(outpath: Path) -> None:
    df = pd.DataFrame({
        "Assembler": list(RUNTIME_MIN.keys()),
        "Runtime_min": list(RUNTIME_MIN.values()),
        "Peak_RAM_GB": [PEAK_MEM_GB[k] for k in RUNTIME_MIN.keys()],
    })

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    df.plot(x="Assembler", y="Runtime_min", kind="bar", ax=axes[0], legend=False)
    axes[0].set_title("Runtime")
    axes[0].set_ylabel("Minutes")
    axes[0].tick_params(axis="x", rotation=0)

    df.plot(x="Assembler", y="Peak_RAM_GB", kind="bar", ax=axes[1], legend=False)
    axes[1].set_title("Peak memory")
    axes[1].set_ylabel("GB")
    axes[1].tick_params(axis="x", rotation=0)

    fig.suptitle("Computational efficiency", fontsize=14)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cumulative_contig_lengths(outpath: Path) -> None:
    plt.figure(figsize=(8, 5))

    for label, fasta in ASSEMBLY_FASTAS.items():
        lengths = read_fasta_lengths(fasta)
        cumulative = pd.Series(lengths).cumsum()
        ranks = range(1, len(lengths) + 1)
        plt.plot(ranks, cumulative, marker="o", label=label)

    plt.xlabel("Contig rank")
    plt.ylabel("Cumulative assembly length (bp)")
    plt.title("Cumulative contig length")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    FIG_DIR.mkdir(exist_ok=True)

    df = read_quast_report(QUAST_TSV)

    plot_assembly_metrics(df, FIG_DIR / "assembly_metrics_panel.png")
    plot_quast_quality(df, FIG_DIR / "assembly_quality_panel.png")
    plot_runtime_memory(FIG_DIR / "runtime_memory_panel.png")
    plot_cumulative_contig_lengths(FIG_DIR / "cumulative_contig_length.png")

    print("Saved assembly figures to:", FIG_DIR)


if __name__ == "__main__":
    main()
