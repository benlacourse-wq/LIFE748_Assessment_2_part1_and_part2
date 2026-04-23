import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
import subprocess
from pathlib import Path

# -----------------------------
# helpers
# -----------------------------
def find_one(pattern):
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No file matched: {pattern}")
    return matches[0]

def count_gff_feature(gff_file, feature):
    count = 0
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2] == feature:
                count += 1
    return count

# make figures folder
Path("figures").mkdir(exist_ok=True)

# -----------------------------
# 1. QUAST assembly metrics
# -----------------------------
quast = pd.read_csv("quast/report.tsv", sep="\t")

# first column is metric name; second column is your assembly
metric_col = quast.columns[0]
value_col = quast.columns[1]

flye = quast[[metric_col, value_col]].copy()
flye.columns = ["Metric", "Value"]

wanted = ["# contigs", "Total length", "Largest contig", "N50", "GC (%)"]
metrics = flye[flye["Metric"].isin(wanted)].copy()

# convert numeric safely
metrics["Value"] = pd.to_numeric(metrics["Value"], errors="coerce")
metrics = metrics.dropna(subset=["Value"])

# nicer order
order = ["Total length", "N50", "Largest contig", "# contigs", "GC (%)"]
metrics["Metric"] = pd.Categorical(metrics["Metric"], categories=order, ordered=True)
metrics = metrics.sort_values("Metric")

plt.figure(figsize=(9, 5))
plt.bar(metrics["Metric"], metrics["Value"])
plt.xticks(rotation=30, ha="right")
plt.title("Flye assembly summary metrics")
plt.ylabel("Value")
plt.tight_layout()
plt.savefig("figures/assembly_metrics.png", dpi=300)
plt.close()

# -----------------------------
# 2. Prokka counts
# -----------------------------
prokka_txt = find_one("annotations/prokka_flye/*.txt")

prokka_counts = {"CDS": 0, "rRNA": 0, "tRNA": 0}
with open(prokka_txt) as f:
    for line in f:
        for key in prokka_counts:
            if line.startswith(f"{key}:"):
                m = re.search(r"(\d+)", line)
                if m:
                    prokka_counts[key] = int(m.group(1))

# -----------------------------
# 3. Bakta counts from GFF3
# -----------------------------
bakta_gff = find_one("annotations/bakta_flye/*.gff3")

bakta_counts = {
    "CDS": count_gff_feature(bakta_gff, "CDS"),
    "rRNA": count_gff_feature(bakta_gff, "rRNA"),
    "tRNA": count_gff_feature(bakta_gff, "tRNA"),
}

ann_df = pd.DataFrame({
    "Feature": ["CDS", "rRNA", "tRNA"],
    "Prokka": [prokka_counts["CDS"], prokka_counts["rRNA"], prokka_counts["tRNA"]],
    "Bakta": [bakta_counts["CDS"], bakta_counts["rRNA"], bakta_counts["tRNA"]],
})

ann_long = ann_df.melt(id_vars="Feature", var_name="Tool", value_name="Count")

plt.figure(figsize=(8, 5))
for tool in ann_long["Tool"].unique():
    sub = ann_long[ann_long["Tool"] == tool]
    plt.bar(
        [f"{feat}\n{tool}" for feat in sub["Feature"]],
        sub["Count"]
    )
plt.title("Annotation comparison: Prokka vs Bakta")
plt.ylabel("Feature count")
plt.tight_layout()
plt.savefig("figures/annotation_comparison.png", dpi=300)
plt.close()

# -----------------------------
# 4. Prokka functional coverage
# -----------------------------
prokka_tsv = find_one("annotations/prokka_flye/*.tsv")
tsv = pd.read_csv(prokka_tsv, sep="\t")

# be defensive about column names
if "ftype" in tsv.columns:
    cds = tsv[tsv["ftype"] == "CDS"].copy()
else:
    cds = tsv.copy()

if "product" in cds.columns:
    annotated = cds["product"].fillna("").str.strip().ne("") & cds["product"].fillna("").str.lower().ne("hypothetical protein")
    annotated_count = int(annotated.sum())
    hypothetical_count = int((~annotated).sum())
else:
    annotated_count = 0
    hypothetical_count = len(cds)

plt.figure(figsize=(5, 5))
plt.pie(
    [annotated_count, hypothetical_count],
    labels=["Annotated", "Hypothetical"],
    autopct="%1.1f%%"
)
plt.title("Prokka functional annotation coverage")
plt.tight_layout()
plt.savefig("figures/functional_coverage.png", dpi=300)
plt.close()

print("Done. Figures written to:")
print("  figures/assembly_metrics.png")
print("  figures/annotation_comparison.png")
print("  figures/functional_coverage.png")
