"""
Simple Linear Plasmid Mutation Map
Draws a horizontal backbone with GenBank features above and mutations below.
No circular plotting, no coordinate remapping, no complex rotations.
"""

import sys
import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO


# ============================================================
# SETTINGS
# ============================================================
FEATURE_COLORS = {
    "CDS": "#4daf4a",
    "gene": "#4daf4a",
    "promoter": "#377eb8",
    "rep_origin": "#ff7f00",
    "ori": "#ff7f00",
    "misc_feature": "#984ea3",
    "primer_bind": "#a65628",
    "protein_bind": "#f781bf",
    "regulatory": "#999999",
    "terminator": "#e41a1c",
    "RBS": "#ffff33",
}

FEATURE_TYPES = [
    "CDS", "gene", "promoter", "rep_origin", "ori",
    "misc_feature", "primer_bind", "protein_bind",
    "regulatory", "terminator", "RBS",
]

MUTATION_COLOR = "#e41a1c"
BACKBONE_COLOR = "#333333"


def read_genbank(genbank_path):
    """Read GenBank file, return record."""
    if not os.path.exists(genbank_path):
        raise FileNotFoundError(f"GenBank file not found: {genbank_path}")
    record = SeqIO.read(genbank_path, "genbank")
    print(f"[OK] GenBank loaded: {record.id} ({len(record.seq):,} bp)")
    return record


def read_mutations(csv_path):
    """Read mutation_report.csv, return DataFrame."""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    
    df = pd.read_csv(csv_path)
    required = ["mutation_id", "event_type", "reference_start", "reference_end"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")
    
    print(f"[OK] Mutations loaded: {len(df)}")
    return df


def extract_features(record, max_features=20):
    """Extract GenBank features as simple dicts."""
    features = []
    for f in record.features:
        if f.type not in FEATURE_TYPES:
            continue
        
        start = int(f.location.start) + 1
        end = int(f.location.end)
        
        label = None
        if "label" in f.qualifiers:
            label = f.qualifiers["label"][0]
        elif "gene" in f.qualifiers:
            label = f.qualifiers["gene"][0]
        elif "product" in f.qualifiers:
            label = f.qualifiers["product"][0][:30]
        
        features.append({
            "type": f.type,
            "start": start,
            "end": end,
            "label": label,
            "color": FEATURE_COLORS.get(f.type, "#cccccc"),
        })
    
    # Limit and sort by start
    features = sorted(features, key=lambda x: x["start"])[:max_features]
    return features


def draw_linear_map(record, features, mutations, output_prefix, title_suffix=""):
    """Draw simple linear plasmid map."""
    plasmid_length = len(record.seq)
    plasmid_name = record.id
    
    # ============================================================
    # FIGURE SETUP
    # ============================================================
    if len(mutations) > 0:
        fig, (ax_map, ax_table) = plt.subplots(
            2, 1,
            figsize=(16, 10),
            gridspec_kw={"height_ratios": [3, 2]},
        )
    else:
        fig, ax_map = plt.subplots(1, 1, figsize=(16, 6))
        ax_table = None
    
    # ============================================================
    # MAP AXES
    # ============================================================
    ax_map.set_xlim(0, plasmid_length * 1.05)
    ax_map.set_ylim(-3, 4)
    
    # Backbone line
    ax_map.axhline(y=0, color=BACKBONE_COLOR, linewidth=3, zorder=1)
    
    # Start and end markers
    ax_map.plot(1, 0, 'o', color=BACKBONE_COLOR, markersize=8, zorder=2)
    ax_map.plot(plasmid_length, 0, 'o', color=BACKBONE_COLOR, markersize=8, zorder=2)
    
    # Start label
    ax_map.text(1, -0.5, "1", fontsize=8, ha='center', va='top', color='#555555', fontweight='bold')
    
    # End label
    ax_map.text(plasmid_length, -0.5, str(plasmid_length), fontsize=8, ha='center', va='top', color='#555555', fontweight='bold')
    
    # Tick marks every 500 bp
    for pos in range(500, plasmid_length, 500):
        ax_map.plot([pos, pos], [-0.2, 0.2], color=BACKBONE_COLOR, linewidth=1.5)
        ax_map.text(pos, -0.5, str(pos), fontsize=7, ha='center', va='top', color='#555555')
    
    # ============================================================
    # GENBANK FEATURES (above backbone)
    # ============================================================
    y_levels = [0.5, 1.2, 1.9, 2.6]
    
    for feat in features:
        y = y_levels[0]
        height = 0.5
        x = feat["start"]
        width = feat["end"] - feat["start"]
        
        if width < plasmid_length * 0.005:
            width = plasmid_length * 0.005
        
        rect = mpatches.FancyBboxPatch(
            (x, y - height/2), width, height,
            boxstyle="round,pad=0.05",
            facecolor=feat["color"],
            edgecolor='white',
            linewidth=0.5,
            alpha=0.85,
            zorder=3,
        )
        ax_map.add_patch(rect)
        
        # Label
        if feat["label"] and width > plasmid_length * 0.02:
            ax_map.text(
                x + width/2, y, feat["label"],
                fontsize=6, ha='center', va='center',
                color='white', fontweight='bold', zorder=4,
            )
    
    # ============================================================
    # MUTATION MARKERS (below backbone) - only if mutations exist
    # ============================================================
    if len(mutations) > 0:
        for _, mut in mutations.iterrows():
            mid = (int(mut["reference_start"]) + int(mut["reference_end"])) / 2
            
            # Vertical line
            ax_map.plot([mid, mid], [-0.4, -1.2], color=MUTATION_COLOR, linewidth=2, zorder=5)
            
            # Diamond marker
            ax_map.plot(mid, -0.3, 'D', color=MUTATION_COLOR, markersize=8,
                       markeredgecolor='darkred', markeredgewidth=1.5, zorder=6)
            
            # Label
            event = mut["event_type"]
            if event == "mismatch":
                label = f"M{mut['mutation_id']} mismatch\npos {mut['reference_start']}"
            elif event == "deletion":
                label = f"M{mut['mutation_id']} deletion\n{mut['reference_start']}-{mut['reference_end']}"
            elif event == "insertion":
                label = f"M{mut['mutation_id']} insertion\nafter {mut['reference_start']}"
            else:
                label = f"M{mut['mutation_id']} {event}"
            
            ax_map.text(mid, -1.5, label, fontsize=7, color=MUTATION_COLOR,
                       ha='center', va='top', fontweight='bold', zorder=6,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                edgecolor=MUTATION_COLOR, alpha=0.9))
    
    # ============================================================
    # LEGEND
    # ============================================================
    legend_patches = []
    seen_types = set()
    for feat in features:
        if feat["type"] not in seen_types:
            seen_types.add(feat["type"])
            legend_patches.append(mpatches.Patch(color=feat["color"], label=feat["type"]))
    
    if len(mutations) > 0:
        legend_patches.append(mpatches.Patch(color=MUTATION_COLOR, label="Mutation"))
    
    ax_map.legend(
        handles=legend_patches,
        loc='upper right',
        fontsize=7,
        framealpha=0.9,
        ncol=2,
    )
    
    # Title
    ax_map.set_title(
        f"{plasmid_name} {title_suffix} ({plasmid_length:,} bp)",
        fontsize=13, fontweight='bold', pad=15,
    )
    
    ax_map.set_yticks([])
    ax_map.spines['top'].set_visible(False)
    ax_map.spines['right'].set_visible(False)
    ax_map.spines['left'].set_visible(False)
    
    # ============================================================
    # TABLE AXES (mutation summary) - only if mutations exist
    # ============================================================
    if ax_table is not None and len(mutations) > 0:
        ax_table.axis('off')
        
        table_data = []
        for _, mut in mutations.iterrows():
            table_data.append([
                str(mut["mutation_id"]),
                mut["event_type"],
                str(mut["reference_position"]),
                str(mut["length"]),
                mut["note"][:60],
            ])
        
        col_labels = ["ID", "Type", "Position", "Length (bp)", "Note"]
        
        table = ax_table.table(
            cellText=table_data,
            colLabels=col_labels,
            cellLoc='left',
            loc='center',
            colWidths=[0.05, 0.08, 0.12, 0.08, 0.42],
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 1.5)
        
        # Style header
        for i in range(len(col_labels)):
            table[0, i].set_facecolor('#333333')
            table[0, i].set_text_props(color='white', fontweight='bold')
        
        # Alternate row colors
        for row in range(1, len(table_data) + 1):
            for col in range(len(col_labels)):
                if row % 2 == 0:
                    table[row, col].set_facecolor('#f5f5f5')
                else:
                    table[row, col].set_facecolor('#ffffff')
        
        ax_table.set_title("Mutation Summary Table", fontsize=11, fontweight='bold', pad=10)
    
    # ============================================================
    # SAVE
    # ============================================================
    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    png_path = f"{output_prefix}.png"
    pdf_path = f"{output_prefix}.pdf"
    
    plt.tight_layout()
    fig.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[OK] PNG saved: {png_path}")
    fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"[OK] PDF saved: {pdf_path}")
    plt.close(fig)

def main():
    """Main entry point."""
    print("=" * 60)
    print("  Simple Linear Plasmid Mutation Map")
    print("=" * 60)
    print()
    
    genbank_path = input("Path to GenBank file (e.g., data/reference/reference_puc19.gb): ").strip()
    mutations_path = input("Path to mutation report CSV (e.g., outputs/tables/mutation_report.csv): ").strip()
    
    print()
    print("-" * 60)
    
    try:
        record = read_genbank(genbank_path)
        mutations_df = read_mutations(mutations_path)
        features = extract_features(record)
        print(f"[OK] Features extracted: {len(features)}")
        
        # Extract base name from GenBank file (without path and extension)
        genbank_base = os.path.splitext(os.path.basename(genbank_path))[0]
        
        # Create output directory
        Path("outputs/figures").mkdir(parents=True, exist_ok=True)
        
        # ============================================================
        # MAP 1: Original plasmid (no mutations)
        # ============================================================
        original_prefix = f"outputs/figures/{genbank_base}_original_map"
        print(f"\n[1/2] Drawing original plasmid map...")
        draw_linear_map(
            record=record,
            features=features,
            mutations=pd.DataFrame(),
            output_prefix=original_prefix,
            title_suffix="— Original Plasmid (Reference)",
        )
        
        # ============================================================
        # MAP 2: Mutation-marked map
        # ============================================================
        mutation_prefix = f"outputs/figures/{genbank_base}_mutation_map"
        print(f"\n[2/2] Drawing mutation-marked map...")
        draw_linear_map(
            record=record,
            features=features,
            mutations=mutations_df,
            output_prefix=mutation_prefix,
            title_suffix=f"— Mutation Map ({len(mutations_df)} mutations)",
        )
        
        print()
        print("=" * 60)
        print("  Done! Both maps generated.")
        print(f"  Original: {original_prefix}.png/.pdf")
        print(f"  Mutations: {mutation_prefix}.png/.pdf")
        print("=" * 60)
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    except ValueError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()