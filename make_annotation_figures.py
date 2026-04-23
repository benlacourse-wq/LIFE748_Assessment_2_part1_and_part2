#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

PROJECT = Path("/home/benla/life748report")
FIG_DIR = PROJECT / "figures"
CSV = PROJECT / "annotation_summary_real.csv"


def plot_feature_counts(df: pd.DataFrame, outpath: Path) -> None:
    plot_df = df.copy()
    plot_df["label"] = plot_df["assembler"] + "_" + plot_df["annotator"]

    cols = ["n_CDS", "n_rRNA", "n_tRNA"]
    long_df = plot_df.melt(
        id_vars="label",
        value_vars=cols,
        var_name="feature",
        value_name="count",
    )

    pivot = long_df.pivot(index="label", columns="feature", values="count")
    pivot.plot(kind="bar", figsize=(10, 6))
    plt.title("Annotation feature counts")
    plt.xlabel("Pipeline")
    plt.ylabel("Count")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_functional_coverage(df: pd.DataFrame, outpath: Path) -> None:
    plot_df = df.copy()
    plot_df["label"] = plot_df["assembler"] + "_" + plot_df["annotator"]

    plt.figure(figsize=(8, 5))
    plt.bar(plot_df["label"], plot_df["frac_function"])
    plt.ylabel("Fraction of CDS with functional annotation")
    plt.xlabel("Pipeline")
    plt.title("Functional annotation coverage")
    plt.xticks(rotation=0)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_runtime(df: pd.DataFrame, outpath: Path) -> None:
    plot_df = df.copy()
    plot_df["label"] = plot_df["assembler"] + "_" + plot_df["annotator"]

    plt.figure(figsize=(8, 5))
    plt.bar(plot_df["label"], plot_df["runtime_min"])
    plt.ylabel("Runtime (min)")
    plt.xlabel("Pipeline")
    plt.title("Annotation runtime")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_gene_density(df: pd.DataFrame, outpath: Path) -> None:
    plot_df = df.copy()
    plot_df["label"] = plot_df["assembler"] + "_" + plot_df["annotator"]
    plot_df["cds_per_mb"] = plot_df["n_CDS"] / (plot_df["total_bp"] / 1_000_000)

    plt.figure(figsize=(8, 5))
    plt.bar(plot_df["label"], plot_df["cds_per_mb"])
    plt.ylabel("CDS per Mb")
    plt.xlabel("Pipeline")
    plt.title("Gene density")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_fragmentation_vs_annotation(df: pd.DataFrame, outpath: Path) -> None:
    plt.figure(figsize=(7, 5))
    for _, row in df.iterrows():
        label = f"{row['assembler']}_{row['annotator']}"
        plt.scatter(row["n_contigs"], row["frac_function"], s=80)
        plt.text(row["n_contigs"], row["frac_function"], label, fontsize=9)

    plt.xlabel("Number of contigs")
    plt.ylabel("Fraction functionally annotated")
    plt.title("Assembly fragmentation vs annotation yield")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    FIG_DIR.mkdir(exist_ok=True)
    df = pd.read_csv(CSV)

    plot_feature_counts(df, FIG_DIR / "annotation_feature_counts.png")
    plot_functional_coverage(df, FIG_DIR / "annotation_functional_coverage.png")
    plot_runtime(df, FIG_DIR / "annotation_runtime.png")
    plot_gene_density(df, FIG_DIR / "annotation_gene_density.png")
    plot_fragmentation_vs_annotation(df, FIG_DIR / "fragmentation_vs_annotation.png")

    print("Saved annotation figures to:", FIG_DIR)


if __name__ == "__main__":
    main()
