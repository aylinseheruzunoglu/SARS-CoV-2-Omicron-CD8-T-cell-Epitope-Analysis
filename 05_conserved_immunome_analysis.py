# ============================================================
# VARIANT-SPECIFIC EPITOPE–ALLELE LOSS ANALYSIS + BAR CHART
# ============================================================

import pandas as pd
import os
import glob
import matplotlib.pyplot as plt

# ============================================================
# 1. EXTRACT LOST EPITOPE–ALLELE PAIRS PER VARIANT
# ============================================================

ancestral_file = "ancestral.tsv"
variant_files = [f for f in glob.glob("*.tsv") if not f.startswith("ancestral") and not f.endswith("_lost_pairs.tsv")]

if not os.path.exists(ancestral_file):
    raise FileNotFoundError("Missing 'ancestral.tsv' file in the directory.")

# Load ancestral data
ancestral_df = pd.read_csv(ancestral_file, sep="\t", dtype=str)
if not {"peptide", "allele"}.issubset(ancestral_df.columns):
    raise ValueError("Ancestral file must contain 'peptide' and 'allele' columns.")

ancestral_pairs = set(zip(ancestral_df["peptide"], ancestral_df["allele"]))

# Loop through each variant file and identify lost pairs
for file in variant_files:
    variant = os.path.basename(file).replace(".tsv", "")
    print(f"Processing variant: {variant}")

    try:
        variant_df = pd.read_csv(file, sep="\t", dtype=str)
        if not {"peptide", "allele"}.issubset(variant_df.columns):
            print(f"Skipping {variant} (missing 'peptide' or 'allele' columns).")
            continue

        variant_pairs = set(zip(variant_df["peptide"], variant_df["allele"]))
        lost_pairs = ancestral_pairs - variant_pairs

        lost_df = pd.DataFrame(list(lost_pairs), columns=["peptide", "allele"])
        lost_df.to_csv(f"{variant}_lost_pairs.tsv", sep="\t", index=False)
        print(f"→ Saved {variant}_lost_pairs.tsv ({len(lost_df)} lost pairs)")

    except Exception as e:
        print(f"Error processing {variant}: {e}")

print("\nAll lost-pairs files generated successfully.\n")

# ============================================================
# 2. VARIANT-SPECIFIC LOSS COUNTS (BAR CHART)
# ============================================================

lost_files = glob.glob("*_lost_pairs.tsv")
loss_data = []

for file in lost_files:
    variant = os.path.basename(file).replace("_lost_pairs.tsv", "")
    try:
        df = pd.read_csv(file, sep="\t", dtype=str)
        lost_count = len(df)
    except Exception:
        lost_count = 0
    loss_data.append({"Variant": variant, "Lost Pairs": lost_count})

loss_df = pd.DataFrame(loss_data)

if not loss_df.empty:
    loss_df = loss_df.sort_values("Variant")

    plt.figure(figsize=(8, 5))
    bars = plt.bar(loss_df["Variant"], loss_df["Lost Pairs"], color="#007994", zorder=3)
    plt.ylabel("Number of Ancestral Epitope–Allele Pairs Lost")
    plt.title("Variant-Specific Losses from Ancestral Spike Epitope Set")
    plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)

    # Label bars with counts
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height - 0.05*height,
                 f"{int(height)}", ha='center', va='top', color='white', fontsize=10)

    plt.tight_layout()
    plt.show()
else:
    print("No '*_lost_pairs.tsv' files found — skipping variant-specific loss counts plot.")

Kod 6: CD8+ T hücresi epitop-HLA çiftlerinin varyantlar arası korunum (conservation) analizi, farklı bağlanma eşik değerlerine göre filtrelenmesi ve allel/lokus bazlı dağılımlarının görselleştirilmesi 
# ============================================================
# CD8 T CELL CONSERVED EPITOPE–HLA PAIR ANALYSIS PIPELINE
# ============================================================
# ANALYSES
# 1. Strict conserved pairs (without threshold filtering)
#    - Load raw variant data
#    - Identify conserved epitope–HLA pairs across all variants
#    - Save results (TSV + Excel)
#    - Plot per-allele and per-locus conserved counts
#
# 2. Threshold-based conserved analysis
#    - Filter by NetMHCIIpan percentile thresholds
#    - Identify conserved epitope–HLA pairs across all variants
#    - Optional heatmap visualization
#    - Plot per-allele and per-locus conserved counts
# ============================================================

# ============================================================
# 1. IMPORTS
# ============================================================
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import os # Import os to handle file paths if needed

# ============================================================
# 2. INPUT FILES
# ============================================================
# Ensure these files exist in the same directory as the script
# or provide full paths to them.
variant_files = {
    "Ancestral": "ancestral.tsv",
    "BA.2.86": "ba286.tsv",
    "JN.1": "jn1.tsv",
    "LP.8.1": "lp81.tsv",
    "NB.1.8.1": "nb181.tsv",
    "XEC": "xec.tsv"
}

# ============================================================
# 3. HELPER FUNCTION TO CHECK FILES
# ============================================================
def check_input_files():
    """Checks if all input files exist before running."""
    print("Checking for input files...")
    all_files_found = True
    missing_files = []
    for variant, filepath in variant_files.items():
        if not os.path.exists(filepath):
            print(f"  [ERROR] File not found for {variant}: {filepath}")
            missing_files.append(filepath)
            all_files_found = False
        else:
            print(f"  [FOUND] File for {variant}: {filepath}")
            
    if not all_files_found:
        print("\n*** ABORTING SCRIPT ***")
        print("The following input files are missing:")
        for f in missing_files:
            print(f" - {f}")
        print("Please make sure all .tsv files are in the correct location.")
    return all_files_found

# ============================================================
# BLOCK 1: Strict conserved pairs (no threshold filtering)
# ============================================================

def run_strict_analysis():
    print("\n==============================================")
    print("Running BLOCK 1: Strict Conserved Pairs")
    print("==============================================")
    
    # --- Load epitope–HLA pairs for each variant ---
    variant_pairs = {}
    for variant, filepath in variant_files.items():
        try:
            df = pd.read_csv(filepath, sep="\t")
            df["peptide"] = df["peptide"].str.strip()
            df["allele"] = df["allele"].str.strip()
            variant_pairs[variant] = set(zip(df["peptide"], df["allele"]))
            print(f"  Loaded {len(variant_pairs[variant])} pairs for {variant}.")
        except Exception as e:
            print(f"  [ERROR] Could not process file {filepath}: {e}")
            return # Stop if a file fails to load

    # --- Find conserved pairs across all variants ---
    conserved_pairs = set.intersection(*variant_pairs.values())
    print(f"\nNumber of strictly conserved epitope–HLA pairs: {len(conserved_pairs)}")
    
    if not conserved_pairs:
        print("No conserved pairs found. Skipping plots for Block 1.")
        return

    # --- Save conserved pairs (TSV + Excel) ---
    conserved_df = pd.DataFrame(list(conserved_pairs), columns=["peptide", "allele"])
    conserved_df.to_csv("conserved_epitope_HLA_pairs_strict.tsv", sep="\t", index=False)
    conserved_df.to_excel("conserved_epitope_HLA_pairs_strict.xlsx", index=False)
    print("  Saved strict conserved pairs to .tsv and .xlsx files.")

    # --- FIGURE 1: Per-allele conserved counts ---
    allele_counts = {}
    for pep, allele in conserved_pairs:
        allele_counts[allele] = allele_counts.get(allele, 0) + 1

    allele_df = pd.DataFrame(list(allele_counts.items()), columns=["allele", "count"])
    allele_df = allele_df.sort_values(by="count", ascending=False)
    allele_df.to_csv("conserved_epitope_counts_per_allele_strict.tsv", sep="\t", index=False)

    plt.figure(figsize=(10, 6))
    bars = plt.barh(
        allele_df["allele"],
        allele_df["count"],
        color="#007994",  
        edgecolor="dimgray"
    )
    plt.gca().invert_yaxis()
    plt.xlabel("Number of Conserved Epitopes", fontsize=12)
    plt.ylabel("HLA Allele", fontsize=12)
    plt.title("Conserved Epitope–HLA Pairs per Allele (no threshold)", fontsize=12, color="black")
    plt.yticks(fontsize=10)
    plt.grid(axis="x", linestyle="--", alpha=0.7, color="dimgray", zorder=0)
    for bar in bars:
        bar.set_zorder(3)
    plt.tight_layout()
    plt.savefig("figure_1_strict_per_allele.png", dpi=300)
    plt.show()
    print("  Generated and saved strict per-allele plot.")

    # --- FIGURE 2: Per-locus conserved counts ---
    locus_counts = {}
    for pep, allele in conserved_pairs:
        try:
            # Assumes format HLA-A*... or HLA-B*...
            locus = allele.split("-")[1][0]  # <-- CHANGED LOGIC for HLA-A, HLA-B
            locus_counts[locus] = locus_counts.get(locus, 0) + 1
        except IndexError:
            print(f"  Warning: Could not parse locus from allele '{allele}'. Skipping.")

    locus_df = pd.DataFrame(list(locus_counts.items()), columns=["locus", "count"])
    locus_df = locus_df.sort_values(by="count", ascending=False)
    locus_df.to_csv("conserved_epitope_counts_per_locus_strict.tsv", sep="\t", index=False)

    plt.figure(figsize=(6, 5))
    bars = plt.bar(
        locus_df["locus"],
        locus_df["count"],
        color="#007994",  
        edgecolor="dimgray"
    )

    # Add numbers inside bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            plt.text(
                bar.get_x() + bar.get_width()/2., height - (0.25 * height),
                f"{int(height)}", ha="center", va="bottom",
                color="white", fontsize=11
            )

    plt.ylabel("Number of Conserved Epitopes", fontsize=12)
    plt.xlabel("HLA Locus", fontsize=12)
    plt.title("Per-locus Conserved Epitope Counts (no threshold)", fontsize=12, color="black")
    plt.grid(axis="y", linestyle="--", alpha=0.7, color="dimgray", zorder=0)
    for bar in bars:
        bar.set_zorder(3)
    plt.tight_layout()
    plt.savefig("figure_2_strict_per_locus.png", dpi=300)
    plt.show()
    print("  Generated and saved strict per-locus plot.")

# ============================================================
# BLOCK 2: Threshold-based conserved analysis
# ============================================================

def analyze_conserved(threshold, label, chunk_size=50, show_heatmap=False):
    """
    Analyze conserved epitope–HLA pairs across variants
    with a given NetMHCIIpan percentile threshold.
    """

    print(f"\n==============================================")
    print(f"Running BLOCK 2: Threshold Analysis {label}")
    print(f"==============================================")

    # --- Load filtered epitope–HLA sets ---
    variant_pairs = {}
    for variant, filepath in variant_files.items():
        try:
            df = pd.read_csv(filepath, sep="\t")
            
            # Check if the required column exists
            if "netmhcpan_el percentile" not in df.columns:
                print(f"  [ERROR] Column 'netmhcpan_el percentile' not found in {filepath}.")
                print("  Make sure the column name matches exactly.")
                return # Stop this analysis
                
            df["peptide"] = df["peptide"].str.strip()
            df["allele"] = df["allele"].str.strip()
            
            # Apply threshold filter
            df = df[df["netmhcpan_el percentile"] <= threshold]
            
            variant_pairs[variant] = set(zip(df["peptide"], df["allele"]))
            print(f"  Loaded {len(variant_pairs[variant])} pairs for {variant} (Threshold {label}).")
        except Exception as e:
            print(f"  [ERROR] Could not process file {filepath}: {e}")
            return # Stop if a file fails to load
            
    # --- Conserved across all variants ---
    conserved_pairs = set.intersection(*variant_pairs.values())
    print(f"\nNumber of conserved epitope–HLA pairs ({label}): {len(conserved_pairs)}")

    if not conserved_pairs:
        print(f"No conserved pairs found at threshold {label}. Skipping plots.")
        return
        
    # --- 1. Heatmap of conserved pairs ---
    if show_heatmap:
        print("  Generating heatmap...")
        heatmap_df = pd.DataFrame(index=[f"{pep}_{hla}" for pep, hla in conserved_pairs])
        for variant, pairs in variant_pairs.items():
            heatmap_df[variant] = heatmap_df.index.map(
                lambda x: tuple(x.split("_")) in pairs
            )
        heatmap_data = heatmap_df.astype(int)

        n_chunks = math.ceil(len(heatmap_data) / chunk_size)
        print(f"  Heatmap will be split into {n_chunks} chunk(s).")
        
        for i in range(n_chunks):
            chunk = heatmap_data.iloc[i*chunk_size:(i+1)*chunk_size]
            if chunk.empty:
                continue
            plt.figure(figsize=(10, max(5, 0.25*len(chunk)))) # Ensure min height
            sns.heatmap(
                chunk,
                cmap=["#007994", "lightgray"],
                cbar=False,
                linewidths=0.3,
                linecolor="white"
            )
            plt.title(f"Conserved Epitope–HLA Pairs ({label} percentile) [Part {i+1}/{n_chunks}]")
            plt.xlabel("Variants")
            plt.ylabel("Peptide–HLA Pairs")
            plt.yticks(rotation=0, fontsize=6)
            plt.tight_layout()
            plt.savefig(f"figure_3_heatmap_{label.replace('≤','')}_{i+1}.png", dpi=300)
            plt.show()
        print("  Saved heatmap plot(s).")

    # --- 2. Per-allele conserved counts ---
    allele_counts = {}
    for pep, allele in conserved_pairs:
        allele_counts[allele] = allele_counts.get(allele, 0) + 1

    allele_df = pd.DataFrame(list(allele_counts.items()), columns=["allele", "count"])
    allele_df = allele_df.sort_values(by="count", ascending=False)
    # Save allele counts for this threshold
    allele_df.to_csv(f"conserved_epitope_counts_per_allele_{label.replace('≤','')}.tsv", sep="\t", index=False)

    plt.figure(figsize=(10, 6))
    plt.barh(allele_df["allele"], allele_df["count"], color="#007994", edgecolor="dimgray", zorder=2) 
    plt.gca().invert_yaxis()
    plt.xlabel("Number of Conserved Epitopes", fontsize=12)
    plt.ylabel("HLA Allele", fontsize=12)
    plt.title(f"Conserved Epitope–HLA Pairs per Allele ({label} percentile)", fontsize=12)
    plt.yticks(fontsize=10)
    plt.grid(axis="x", linestyle="--", alpha=0.6, color="dimgray", zorder=0)
    plt.tight_layout()
    plt.savefig(f"figure_4_per_allele_{label.replace('≤','')}.png", dpi=300)
    plt.show()
    print(f"  Generated and saved per-allele plot ({label}).")

    # --- 3. Per-locus conserved counts ---
    locus_counts = {}
    for pep, allele in conserved_pairs:
        try:
            locus = allele.split("-")[1][0]  # <-- CHANGED LOGIC for # HLA-A, HLA-B
            locus_counts[locus] = locus_counts.get(locus, 0) + 1
        except IndexError:
            print(f"  Warning: Could not parse locus from allele '{allele}'. Skipping.")

    locus_df = pd.DataFrame(list(locus_counts.items()), columns=["locus", "count"])
    locus_df = locus_df.sort_values(by="count", ascending=False)
    # Save locus counts for this threshold
    locus_df.to_csv(f"conserved_epitope_counts_per_locus_{label.replace('≤','')}.tsv", sep="\t", index=False)

    plt.figure(figsize=(6, 5))
    bars = plt.bar(locus_df["locus"], locus_df["count"], color="#007994", edgecolor="dimgray", zorder=2) 

    # Add numbers inside bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            plt.text(
                bar.get_x() + bar.get_width()/2., height - (0.3 * height),
                f"{int(height)}", ha="center", va="bottom",
                color="white", fontsize=11
            )

    plt.ylabel("Number of Conserved Epitopes", fontsize=12)
    plt.xlabel("HLA Locus", fontsize=12)
    plt.title(f"Per-locus Conserved Epitope Counts ({label} percentile)", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.6, color="dimgray", zorder=0)
    plt.tight_layout()
    plt.savefig(f"figure_5_per_locus_{label.replace('≤','')}.png", dpi=300)
    plt.show()
    print(f"  Generated and saved per-locus plot ({label}).")

# ============================================================
# 4. RUN ANALYSES
# ============================================================
def main():
    if check_input_files():
        # --- Run Block 1 ---
        run_strict_analysis()
        
        # --- Run Block 2 with new thresholds ---
        analyze_conserved(0.5, "≤0.5", show_heatmap=False)
        analyze_conserved(0.2, "≤0.2", show_heatmap=False)
        analyze_conserved(0.1, "≤0.1", show_heatmap=True)
        analyze_conserved(0.05, "≤0.05", show_heatmap=True)
        
        print("\n==============================================")
        print("All analyses complete.")
        print("==============================================")

if __name__ == "__main__":
    main()
