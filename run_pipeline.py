"""
run_pipeline.py

Master pipeline script for plasmid-mutation-qc.

Runs the complete workflow from a single config.yml file:
1. Simulate mutations → mutated FASTA
2. Compare sequences → mutation_report.csv + alignment.txt
3. Draw linear plasmid maps → original + mutation-marked PNG/PDF
4. Draw circular plasmid map → PNG/PDF
5. Generate summary report → Markdown + PDF

Usage:
    python run_pipeline.py
    python run_pipeline.py --config config.yml
"""

import sys
import os
import subprocess
from pathlib import Path
from datetime import datetime

import yaml


def print_header(title):
    """Print a formatted section header."""
    print()
    print("=" * 60)
    print(f"  {title}")
    print("=" * 60)
    print()


def print_step(step_num, description):
    """Print a step indicator."""
    print(f"\n[{step_num}/6] {description}...")
    print("-" * 40)


def run_command(cmd, description=""):
    """Run a shell command and handle errors."""
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  [ERROR] {description} failed!")
        print(f"  stdout: {result.stdout}")
        print(f"  stderr: {result.stderr}")
        return False
    
    if result.stdout:
        for line in result.stdout.strip().split('\n'):
            print(f"  {line}")
    
    return True


def load_config(config_path):
    """Load and validate config.yml."""
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Validate required fields
    required = ['reference_genbank', 'reference_fasta', 'mutated_fasta', 'mutations']
    missing = [k for k in required if k not in config]
    if missing:
        print(f"[ERROR] Config missing required fields: {', '.join(missing)}")
        sys.exit(1)
    
    # Set defaults
    config.setdefault('project_name', 'plasmid-mutation-qc')
    config.setdefault('topology', 'circular')
    config.setdefault('output_folder', 'outputs')
    
    return config


def get_genbank_base(config):
    """Get base name from GenBank file for output naming."""
    gb_path = config['reference_genbank']
    return os.path.splitext(os.path.basename(gb_path))[0]


def main():
    """Main pipeline runner."""
    start_time = datetime.now()
    
    print("=" * 60)
    print("  PLASMID MUTATION QC — FULL PIPELINE")
    print("=" * 60)
    print(f"  Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    
    # ============================================================
    # Load config
    # ============================================================
    config_path = "config.yml"
    if len(sys.argv) > 2 and sys.argv[1] == "--config":
        config_path = sys.argv[2]
    
    print(f"\n  Config file: {config_path}")
    config = load_config(config_path)
    gb_base = get_genbank_base(config)
    
    print(f"  Project: {config['project_name']}")
    print(f"  Plasmid: {gb_base}")
    print(f"  Topology: {config['topology']}")
    print(f"  Mutations defined: {len(config['mutations'])}")
    
    # ============================================================
    # Create output directories
    # ============================================================
    Path("outputs/tables").mkdir(parents=True, exist_ok=True)
    Path("outputs/figures").mkdir(parents=True, exist_ok=True)
    Path("outputs/reports").mkdir(parents=True, exist_ok=True)
    Path("data/mutated").mkdir(parents=True, exist_ok=True)
    
    # ============================================================
    # Step 1: Simulate mutations
    # ============================================================
    print_step(1, "Simulating mutations")
    
    if not run_command(
        ["python", "src/simulate_mutation.py", "--config", config_path],
        "Mutation simulation"
    ):
        sys.exit(1)
    
    # ============================================================
    # Step 2: Compare sequences → mutation_report.csv
    # ============================================================
    print_step(2, "Comparing sequences & generating mutation report")
    
    ref_fasta = config['reference_fasta']
    mut_fasta = config['mutated_fasta']
    report_csv = f"outputs/tables/{gb_base}_mutation_report.csv"
    alignment_txt = f"outputs/reports/{gb_base}_alignment.txt"
    
    if not run_command(
        [
            "python", "src/compare_sequences.py",
            "--reference", ref_fasta,
            "--sample", mut_fasta,
            "--output", report_csv,
            "--alignment-output", alignment_txt,
        ],
        "Sequence comparison"
    ):
        sys.exit(1)
    
    # ============================================================
    # Step 3: Draw linear plasmid maps
    # ============================================================
    print_step(3, "Drawing linear plasmid maps")
    
    genbank_path = config['reference_genbank']
    original_prefix = f"outputs/figures/{gb_base}_original_map"
    mutation_prefix = f"outputs/figures/{gb_base}_mutation_map"
    
    # Check if linear map script exists
    linear_script = "src/draw_linear_map.py"
    if os.path.exists(linear_script):
        # Run linear map script with original (no mutations)
        # We'll use a temporary approach: call the function directly
        print("  Generating original linear map...")
        try:
            import pandas as pd
            from Bio import SeqIO
            
            # Import the linear map function
            sys.path.insert(0, 'src')
            from draw_linear_map import read_genbank, extract_features, draw_linear_map
            
            record = read_genbank(genbank_path)
            features = extract_features(record)
            
            # Original map (no mutations)
            draw_linear_map(
                record=record,
                features=features,
                mutations=pd.DataFrame(),
                output_prefix=original_prefix,
                title_suffix="— Original Plasmid (Reference)",
            )
            
            # Mutation-marked map
            mutations_df = pd.read_csv(report_csv)
            draw_linear_map(
                record=record,
                features=features,
                mutations=mutations_df,
                output_prefix=mutation_prefix,
                title_suffix=f"— Mutation Map ({len(mutations_df)} mutations)",
            )
            print("  [OK] Linear maps generated")
        except Exception as e:
            print(f"  [WARNING] Linear map generation skipped: {e}")
    else:
        print(f"  [SKIP] Linear map script not found: {linear_script}")
    
    
    # ============================================================
    # Step 5: Generate summary report
    # ============================================================
    print_step(5, "Generating summary report")
    
    map_png = f"outputs/figures/{gb_base}_circular_map.png"
    if not os.path.exists(map_png):
        map_png = f"outputs/figures/{gb_base}_mutation_map.png"
    
    md_report = f"outputs/reports/{gb_base}_summary_report.md"
    pdf_report = f"outputs/reports/{gb_base}_summary_report.pdf"
    
    if not run_command(
        [
            "python", "src/generate_summary.py",
            "--genbank", genbank_path,
            "--mutations", report_csv,
            "--alignment", alignment_txt,
            "--map-image", map_png,
            "--output-md", md_report,
            "--output-pdf", pdf_report,
            "--topology", config['topology'],
            "--project-name", config['project_name'],
        ],
        "Summary report generation"
    ):
        print("  [WARNING] Summary report generation had issues")
    
    # ============================================================
    # Step 6: Print summary
    # ============================================================
    print_header("PIPELINE COMPLETE")
    
    end_time = datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    print(f"  Duration: {duration:.1f} seconds")
    print()
    print("  Generated files:")
    print(f"    data/mutated/{os.path.basename(config['mutated_fasta'])}")
    print(f"    {report_csv}")
    print(f"    {alignment_txt}")
    print(f"    {original_prefix}.png / .pdf")
    print(f"    {mutation_prefix}.png / .pdf")
    print(f"    outputs/figures/{gb_base}_circular_map.png / .pdf")
    print(f"    {md_report}")
    print(f"    {pdf_report}")
    print()
    print("=" * 60)


if __name__ == "__main__":
    main()