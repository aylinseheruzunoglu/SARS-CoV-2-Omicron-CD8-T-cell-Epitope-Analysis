# ============================================================
# CD8 T CELL EPITOPE–HLA ALLELE PAIR ANALYSIS PIPELINE
# ============================================================
# ANALYSES
# 1. Variant-specific loss counts (bar chart)
# 2. Pairs preserved in JN.1/BA.2.86 but lost in all newer variants (heatmap)
# 3. Alternative check: BA.2.86/JN.1 preserved but lost in newer variants
# 4. Positional mapping of lost epitopes (JN.1 reference)
# 5. Domain mapping & summary plots for lost epitopes (mapping to a Spike reference)
# ============================================================

# ============================================================
# 1. IMPORTS
# ============================================================
import os
import glob
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Instaling Required Packages
try:
    # Install Biopython if not already installed
    %pip install biopython
    from Bio import SeqIO
except Exception as e:
    raise ImportError("Biopython is required for sequence parsing. Install it with: pip install biopython") from e

# ============================================================
# 2. FILE PATHS FOR VARIANT DATA
# ============================================================
variant_files = {
    "Ancestral": "ancestral.tsv",
    "BA.2.86": "ba286.tsv",
    "JN.1": "jn1.tsv",
    "LP.8.1": "lp81.tsv",
    "NB.1.8.1": "nb181.tsv",
    "XEC": "xec.tsv"
}

# Validate variant files exist
for v, p in variant_files.items():
    if not Path(p).exists():
        print(f"Warning: variant file for {v} not found at '{p}'. Make sure the file exists before running that section.")

# ============================================================
# 3. VARIANT-SPECIFIC LOSS COUNTS (BAR CHART)
# ============================================================
# Collect all *_lost_pairs.tsv files (these are expected outputs from other parts of the workflow)
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

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height - 0.05*height,
                 f"{int(height)}", ha='center', va='top', color='white', fontsize=10)
    plt.tight_layout()
    plt.show()
else:
    print("No '*_lost_pairs.tsv' files found — skipping variant-specific loss counts plot.")

# ============================================================
# 4. PAIRS PRESERVED IN JN.1 / BA.2.86 BUT LOST IN ALL NEWER VARIANTS (HEATMAP)
# ============================================================
# Load per-variant peptide/HLA pairs (as tuples)
variant_pairs = {}
for variant, filepath in variant_files.items():
    if Path(filepath).exists():
        df = pd.read_csv(filepath, sep="\t", dtype=str).rename(columns=str.strip)
        # Expect columns named 'peptide' and 'allele'
        if "peptide" not in df.columns or "allele" not in df.columns:
            print(f"Warning: '{filepath}' does not contain both 'peptide' and 'allele' columns.")
            variant_pairs[variant] = set()
        else:
            df["peptide"] = df["peptide"].astype(str).str.strip()
            df["allele"] = df["allele"].astype(str).str.strip()
            variant_pairs[variant] = set(zip(df["peptide"], df["allele"]))
    else:
        variant_pairs[variant] = set()

# Define reference set as union of JN.1 and BA.2.86
reference_pairs = variant_pairs.get("JN.1", set()).union(variant_pairs.get("BA.2.86", set()))

# Find pairs present in reference but absent in LP.8.1, NB.1.8.1, XEC
lost_in_all_newer = [
    pair for pair in reference_pairs
    if pair not in variant_pairs.get("LP.8.1", set())
    and pair not in variant_pairs.get("NB.1.8.1", set())
    and pair not in variant_pairs.get("XEC", set())
]

print(f"Number of epitope–HLA pairs preserved in JN.1/BA.2.86 but lost in all newer variants: {len(lost_in_all_newer)}")

if lost_in_all_newer:
    heatmap_df = pd.DataFrame(index=[f"{pep}_{hla}" for pep, hla in lost_in_all_newer])
    variant_order = ["Ancestral", "BA.2.86", "JN.1", "LP.8.1", "NB.1.8.1", "XEC"]
    for variant in variant_order:
        pairs = variant_pairs.get(variant, set())
        heatmap_df[variant] = heatmap_df.index.map(lambda x: tuple(x.split("_")) in pairs)
    heatmap_data = heatmap_df.astype(int)

    plt.figure(figsize=(10, max(1, 0.25*len(heatmap_data))))
    sns.heatmap(heatmap_data, cmap=["#ffffd9", "#007994"], cbar=False, linewidths=0.4, linecolor='#1a517d')
    plt.title("Epitope–HLA Pairs: Preserved in JN.1/BA.2.86 but Lost in All Newer Variants")
    plt.xlabel("Variants")
    plt.ylabel("Peptide–HLA Pairs")
    plt.yticks(rotation=0, fontsize=6)
    plt.tight_layout()
    plt.show()

    # Save list for downstream use
    lost_df = pd.DataFrame([tuple(x.split("_")) for x in heatmap_df.index], columns=["peptide", "allele"])
    lost_df.to_csv("lost_in_all_newer.tsv", sep="\t", index=False)
    print("Saved lost_in_all_newer.tsv")
else:
    print("No pairs found for the 'preserved in JN.1/BA.2.86 but lost in all newer variants' analysis.")

# ============================================================
# 5. ALTERNATIVE CHECK: BA.2.86/JN.1 PRESERVED BUT LOST IN NEWER VARIANTS
# ============================================================
variant_pairs_str = {}
for variant, filepath in variant_files.items():
    if Path(filepath).exists():
        df = pd.read_csv(filepath, sep="\t", dtype=str).rename(columns=str.strip)
        if "peptide" in df.columns and "allele" in df.columns:
            variant_pairs_str[variant] = set(df["peptide"].astype(str).str.strip() + "_" + df["allele"].astype(str).str.strip())
        else:
            variant_pairs_str[variant] = set()
    else:
        variant_pairs_str[variant] = set()

preserved_old = variant_pairs_str.get("BA.2.86", set()).intersection(variant_pairs_str.get("JN.1", set()))
lost_newer = preserved_old - (variant_pairs_str.get("LP.8.1", set()) | variant_pairs_str.get("NB.1.8.1", set()) | variant_pairs_str.get("XEC", set()))

if len(lost_newer) == 0:
    print("No epitope–HLA pairs were found that are preserved in BA.2.86/JN.1 but lost in newer variants.")
else:
    heatmap_df2 = pd.DataFrame(index=sorted(lost_newer), columns=variant_files.keys())
    for variant in variant_files.keys():
        heatmap_df2[variant] = heatmap_df2.index.isin(variant_pairs_str.get(variant, set()))
    heatmap_data2 = heatmap_df2.astype(int)

    plt.figure(figsize=(10, max(6, 0.3 * len(heatmap_data2))))
    sns.heatmap(heatmap_data2, cmap=["#ffffd9", "#007994"], cbar=False, linewidths=0.5, linecolor='#1a517d')
    plt.title("Epitope–HLA Pairs Preserved in BA.2.86/JN.1 but Lost in Newer Variants")
    plt.xlabel("Variants")
    plt.ylabel("Epitope–HLA Pairs")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()

# ============================================================
# 6. POSITIONAL MAPPING OF LOST EPITOPE–HLA PAIRS (JN.1 REFERENCE)
# ============================================================
variant_pairs_pos = {}
for variant, filepath in variant_files.items():
    if Path(filepath).exists():
        df = pd.read_csv(filepath, sep="\t", dtype=str).rename(columns=str.strip)
        if "peptide" in df.columns and "allele" in df.columns:
            df["Peptide_HLA"] = df["peptide"].astype(str).str.strip() + "_" + df["allele"].astype(str).str.strip()
            variant_pairs_pos[variant] = set(df["Peptide_HLA"])
        else:
            variant_pairs_pos[variant] = set()
    else:
        variant_pairs_pos[variant] = set()

# This logic is now consistent with Section 4/7:
# (JN.1 OR BA.2.86) - (LP.8.1 OR NB.1.8.1 OR XEC)
lost_pairs = (variant_pairs_pos.get("JN.1", set()) | variant_pairs_pos.get("BA.2.86", set())) - \
             (variant_pairs_pos.get("LP.8.1", set()) | variant_pairs_pos.get("NB.1.8.1", set()) | \
              variant_pairs_pos.get("XEC", set()))

print(f"Number of lost pairs (JN.1 or BA.2.86, but not in newer): {len(lost_pairs)}")

# Attempt to read JN.1 reference fasta
jn1_fasta = "JN.1_reference.fasta"
if Path(jn1_fasta).exists():
    jn1_record = SeqIO.read(jn1_fasta, "fasta")
    jn1_seq = str(jn1_record.seq)
else:
    print(f"Warning: '{jn1_fasta}' not found. Positional mapping section will be skipped unless you provide the file.")
    jn1_seq = None

lost_positions = []
if jn1_seq:
    # === MODIFIED SECTION ===
    # Load BOTH ancestral and BA.2.86 data to find positions, as
    # the lost_pairs set is a union of epitopes from both sources.

    ancestral_pos_df = pd.DataFrame()
    ba286_pos_df = pd.DataFrame()

    required_cols = ["peptide", "allele", "start", "end"]

    if Path("ancestral.tsv").exists():
        try:
            ancestral_pos_df = pd.read_csv("ancestral.tsv", sep="\t", dtype=str).rename(columns=str.strip)
            if not all(col in ancestral_pos_df.columns for col in required_cols):
                print("Warning: 'ancestral.tsv' is missing required coordinate columns. Will not be used for mapping.")
                ancestral_pos_df = pd.DataFrame() # Reset if bad
        except Exception as e:
            print(f"Error reading 'ancestral.tsv': {e}")
            ancestral_pos_df = pd.DataFrame()
    else:
        print("Warning: 'ancestral.tsv' not found. Cannot map ancestral-derived pairs.")

    if Path("ba286.tsv").exists():
        try:
            ba286_pos_df = pd.read_csv("ba286.tsv", sep="\t", dtype=str).rename(columns=str.strip)
            if not all(col in ba286_pos_df.columns for col in required_cols):
                print("Warning: 'ba286.tsv' is missing required coordinate columns. Will not be used for mapping.")
                ba286_pos_df = pd.DataFrame() # Reset if bad
        except Exception as e:
            print(f"Error reading 'ba286.tsv': {e}")
            ba286_pos_df = pd.DataFrame()
    else:
        print("Warning: 'ba286.tsv' not found. Cannot map BA.2.86-derived pairs.")

    # Combine the position data
    combined_pos_df = pd.concat([ancestral_pos_df, ba286_pos_df])

    if combined_pos_df.empty or not all(col in combined_pos_df.columns for col in required_cols):
         print(f"Warning: No valid position data found in 'ancestral.tsv' or 'ba286.tsv'. Cannot perform positional mapping.")
         lost_pos_data = pd.DataFrame() # Create empty DataFrame to skip mapping
    else:
        # Ensure only columns we know exist are kept
        combined_pos_df = combined_pos_df[required_cols]
        combined_pos_df["Peptide_HLA"] = combined_pos_df["peptide"].astype(str).str.strip() + "_" + combined_pos_df["allele"].astype(str).str.strip()
        # Keep the first match for any pair (in case it's in both)
        combined_pos_df.drop_duplicates(subset="Peptide_HLA", keep="first", inplace=True)

        # Filter for the identified lost pairs and get their positions
        lost_pos_data = combined_pos_df[combined_pos_df["Peptide_HLA"].isin(lost_pairs)].copy()
    # === END OF MODIFIED SECTION ===

    if not lost_pos_data.empty:
        # Use the start and end columns from the combined data
        lost_positions_df = lost_pos_data[["peptide", "allele", "start", "end"]].copy()
        lost_positions_df.rename(columns={"peptide": "Peptide", "allele": "HLA", "start": "Start", "end": "End"}, inplace=True)
        lost_positions_df["Start"] = pd.to_numeric(lost_positions_df["Start"], errors='coerce')
        lost_positions_df["End"] = pd.to_numeric(lost_positions_df["End"], errors='coerce')
        lost_positions_df.dropna(subset=["Start", "End"], inplace=True) # Drop rows with invalid positions

        if not lost_positions_df.empty:
            lost_positions_df.to_excel("Lost_Epitopes_Positions_JN1.xlsx", index=False)
            print("Positional mapping saved to Lost_Epitopes_Positions_JN1.xlsx")

            # Simple positional plot
            plt.figure(figsize=(12, 2))
            plt.hlines(1, 1, len(jn1_seq), colors="dimgray", linestyles="dashed")
            for _, row in lost_positions_df.iterrows():
                plt.hlines(1, row["Start"], row["End"], colors="#007994", linewidth=10)
            plt.title("Positional Mapping of Lost Epitopes (JN.1/BA.2.86 Preserved, Newer Lost)", fontsize=14)
            plt.xlabel("Spike Protein Position (JN.1 reference)")
            plt.yticks([])
            plt.tight_layout()
            plt.show()
        else:
            print("No lost pairs with valid positions found for positional mapping.")
    else:
        print("No lost pairs found in combined ancestral/BA.2.86 data for positional mapping.")
else:
    print("Skipping JN.1 positional plotting (no JN.1 reference sequence available).")

# ============================================================
# 7. DOMAIN MAPPING & SUMMARY PLOTS FOR LOST EPITOPE–HLA PAIRS
# ============================================================
# Use a lost-pairs input for mapping:
# Prefer the pipeline-created 'lost_in_all_newer.tsv' if present; otherwise expect 'lost_pairs.tsv'.
if Path("lost_in_all_newer.tsv").exists():
    lost_pairs_path = "lost_in_all_newer.tsv"
elif Path("lost_pairs.tsv").exists():
    lost_pairs_path = "lost_pairs.tsv"
else:
    lost_pairs_path = None
    print("No lost pairs input file found ('lost_in_all_newer.tsv' or 'lost_pairs.tsv'). Skipping domain mapping section.")

# Domain definitions for Spike (Wuhan-Hu-1 canonical coords; 1-based inclusive)
domains = {
    "Signal peptide": (1, 13),
    "NTD": (14, 305),
    "RBD": (319, 541),
    "RBM": (437, 508),
    "Furin cleavage / S1/S2": (681, 686),
    "S2 (ecto)": (686, 1208),
    "TM": (1213, 1237),
    "CT": (1238, 1273)
}

def read_fasta_sequence(path):
    seq_lines = []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)

if lost_pairs_path:
    # Attempt to find a spike reference fasta to map against. Prefer spike_reference.fasta, fallback to JN.1_reference.fasta
    if Path("spike_reference.fasta").exists():
        spike_fasta_path = "spike_reference.fasta"
    elif Path("JN.1_reference.fasta").exists():
        spike_fasta_path = "JN.1_reference.fasta"
    else:
        spike_fasta_path = None
        print("No spike FASTA found (checked 'spike_reference.fasta' and 'JN.1_reference.fasta'). Provide a FASTA or set 'spike_seq' variable inline.")

    if spike_fasta_path:
        spike_seq = read_fasta_sequence(spike_fasta_path).upper()
        print(f"Loaded spike sequence from {spike_fasta_path} (length {len(spike_seq)} aa).")
    else:
        spike_seq = None

    # Load lost pairs table
    if spike_seq is not None:
        mp_df = pd.read_csv(lost_pairs_path, sep="\t", dtype=str).rename(columns=str.strip)
        if "peptide" not in mp_df.columns:
            raise ValueError("Input lost pairs file must contain a 'peptide' column.")
        if "allele" not in mp_df.columns:
            mp_df["allele"] = ""
        mp_df["peptide"] = mp_df["peptide"].astype(str).str.strip().str.upper()

        # mapping helper
        def find_all_matches(seq, peptide):
            matches = []
            start_idx = 0
            while True:
                idx = seq.find(peptide, start_idx)
                if idx == -1:
                    break
                matches.append((idx+1, idx+len(peptide)))  # 1-based
                start_idx = idx + 1
            return matches

        rows = []
        has_coords = ("start" in mp_df.columns) and ("end" in mp_df.columns)
        for _, row in mp_df.iterrows():
            pep = row["peptide"].upper()
            allele = row.get("allele", "")
            matches = []

            # If coords present and valid, use them
            if has_coords and pd.notna(row.get("start")) and pd.notna(row.get("end")):
                try:
                    s = int(row["start"])
                    e = int(row["end"])
                    if 1 <= s <= e <= len(spike_seq):
                        matches = [(s, e)]
                except Exception:
                    matches = []

            # Otherwise search sequence
            if not matches:
                matches = find_all_matches(spike_seq, pep)

            if not matches:
                rows.append({
                    "peptide": pep,
                    "allele": allele,
                    "matches_count": 0,
                    "matched_positions": "",
                    "domain_list": "",
                    "notes": "not_found_in_reference"
                })
            else:
                pos_strs = []
                domain_set = set()
                for s, e in matches:
                    pos_strs.append(f"{s}-{e}")
                    for dname, (ds, de) in domains.items():
                        if not (e < ds or s > de):
                            domain_set.add(dname)
                rows.append({
                    "peptide": pep,
                    "allele": allele,
                    "matches_count": len(matches),
                    "matched_positions": ";".join(pos_strs),
                    "domain_list": ";".join(sorted(domain_set)),
                    "notes": ""
                })

        mapped_df = pd.DataFrame(rows)
        # Merge original columns (if any extras) for convenience
        mapped_df = pd.concat([mp_df.reset_index(drop=True), mapped_df.drop(columns=["peptide", "allele"])], axis=1)

        # Save mapping outputs
        mapped_df.to_csv("lost_pairs_mapped.tsv", sep="\t", index=False)
        mapped_df.to_excel("lost_pairs_mapped.xlsx", index=False)
        print(f"Saved mapped results to lost_pairs_mapped.tsv and lost_pairs_mapped.xlsx (rows: {len(mapped_df)})")

        # Create summary plots (if any mapped)
        found = mapped_df[mapped_df["matches_count"].astype(int) > 0].copy()
        if not found.empty:
            found["first_start"] = found["matched_positions"].str.split(";").str[0].str.split("-").str[0].astype(int)

            # Histogram of start positions
            plt.figure(figsize=(10, 5))
            counts, bins, patches = plt.hist(found["first_start"], bins=30, edgecolor="black", color="#007994")
            plt.xlabel("Spike amino acid position (reference, 1-based)")
            plt.ylabel("Number of epitopes")
            plt.title("Distribution of mapped epitope start positions")
            plt.grid(axis="y", linestyle="--", color="gray", alpha=0.7, zorder=0)
            for patch in patches:
                patch.set_zorder(3)
            # Add numeric labels for non-zero bars (center of bins)
            bin_width = bins[1] - bins[0] if len(bins) > 1 else 1
            for idx, cnt in enumerate(counts):
                if cnt > 0:
                    center = bins[idx] + bin_width / 2
                    plt.text(center, cnt, str(int(cnt)), ha="center", va="bottom", color="white", fontsize=9, fontweight="bold")
            plt.tight_layout()
            plt.show()

            # Domain summary bar chart
            exploded = found.assign(domain=found["domain_list"].str.split(";")).explode("domain")
            domain_counts = exploded["domain"].value_counts().sort_values(ascending=False)
            if not domain_counts.empty:
                plt.figure(figsize=(8, 4))
                domain_counts.plot(kind="bar", color="#007994", edgecolor="black")
                plt.ylabel("Number of mapped peptides")
                plt.title("Mapped peptides per Spike domain (overlap allowed)")
                plt.tight_layout()
                plt.show()
        else:
            print("No peptides were mapped to the reference sequence; no plots generated.")
    else:
        print("Spike reference not available; domain mapping skipped.")
# End of pipeline
