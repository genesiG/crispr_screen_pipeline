#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_5_color_pdb.py - Color PDB structures based on CRISPR tiling scores
Steps:
1. Load and process CRISPR tiling scores data
2. Map scores to PDB residues
3. Recolor PDB structure in PyMOL
4. Save colored PyMOL session
"""

import os
import argparse
import pandas as pd
import numpy as np
import pymol
from pymol import cmd
import config


def load_and_process_data(data_path, peptide_id):
    """Load CRISPR data and filter for specific peptide"""
    df = pd.read_csv(data_path)

    # Filter for target peptide and valid scores
    df = df[df['ensembl_peptide_id'] == peptide_id]
    df = df[['position', 'LFC_LOESS']].dropna()

    # Average scores for positions with multiple guides
    score_df = df.groupby('position')['LFC_LOESS'].mean().reset_index()
    return score_df

def map_scores_to_pdb(score_df, pdb_file, chain_id):
    """Map scores to PDB residues using residue numbering"""
    # Create score dictionary
    score_dict = dict(zip(score_df['position'], score_df['LFC_LOESS']))

    # Load PDB structure and extract residue info
    cmd.load(pdb_file, "target")
    cmd.remove(f"not chain {chain_id} and not resn HOH")  # Keep only target chain and waters

    # Get residue positions in PDB
    residues = []
    cmd.iterate("target and name CA",
                "residues.append((chain, resv, resn))", 
                space={'residues': residues})

    # Map scores to PDB residues
    mapping = []
    for chain, resi, resn in residues:
        try:
            pos = int(resi)
            score = score_dict.get(pos, None)
            mapping.append((chain, resi, resn, score))
        except ValueError:
            # Skip residues with insertion codes (e.g., "100A")
            continue

    return mapping

def create_color_gradient():
    """Create blue-white-red color gradient"""
    return {
        "blue_neg4": [0.0, 51/255, 153/255],
        "blue_neg2": [0.5, 0.5, 1.0],
        "white_zero": [1.0, 1.0, 1.0],
        "red_pos2": [1.0, 0.5, 0.5],
        "red_pos4": [102/255, 0.0, 153/255]
    }

def color_pdb_structure(mapping):
    """Color PDB structure based on scores"""
    # Calculate score range for normalization
    scores = [s for _,_,_,s in mapping if s is not None]
    if not scores:
        print("No scores found for coloring!")
        return

    max_score = max(scores)
    min_score = min(scores)
    max_abs = max(abs(max_score), abs(min_score))

    # Apply colors directly to residues
    for chain, resi, resn, score in mapping:
        sel = f"chain {chain} and resi {resi}"

        if score is None:
            cmd.color("gray", sel)
        else:
            # Normalize score to [-4, 4] range
            # norm_score = 4 * (score / max_abs) if max_abs > 0 else 0
            norm_score = score

            # Determine color based on normalized score
            if norm_score < -1.0:
                color_name = "blue_neg4"
            # elif norm_score < -0.8:
            #     color_name = "blue_neg2"
            # elif norm_score < 0.8:
            #     color_name = "white_zero"
            # elif norm_score < 2.4:
            #     color_name = "red_pos2"
            else:
                # color_name = "red_pos4"
                color_name = "white_zero"

            cmd.color(color_name, sel)


def main(input_folder, pdb, chain_id):
    """Main workflow to color PDB structure"""

    # Directory to save colored PDB files
    output_dir = os.path.join(config.CRISPRODIR,
                              input_folder,
                              "colored_pdb")  

    structures_dir = os.path.join(config.CRISPRODIR,
                                  input_folder,
                                  "structures",
                                  pdb)

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(structures_dir, exist_ok=True)

    # Validate input PDB file
    pdb_lowercase = pdb.lower()
    pdb_file = os.path.join(structures_dir, f"{pdb_lowercase}.pdb")
    if not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Get peptide ID
    summary_dir = os.path.join(config.CRISPRODIR, input_folder, "score_summary")
    summary_file = os.path.join(summary_dir, "crispro_stats_LFC.csv")
    df = pd.read_csv(summary_file)
    peptide_id = df['ensembl_peptide_id'].iloc[0]

    # Initialize PyMOL
    pymol.finish_launching(['pymol', '-cQ'])
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_flat_sheets", 0)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.show("cartoon")
    cmd.remove("resn HOH")  # Remove water molecules
    # cmd.show("surface", "all")

    # Process data
    data_path = os.path.join(config.CRISPRODIR,
                             input_folder,
                            "crispro_scores_merged_pdb.csv")
    score_df = load_and_process_data(data_path, peptide_id)
    residue_mapping = map_scores_to_pdb(score_df, pdb_file, chain_id)

    # Create color gradient
    color_dict = create_color_gradient()

    # Register colors in PyMOL
    for name, rgb in color_dict.items():
        cmd.set_color(name, rgb)

    # Color structure
    color_pdb_structure(residue_mapping)

    # cmd.set("surface_ramp_above_mode", 0)
    # Critical: Set surface color mode to "atomic" to use residue colors
    # cmd.set("surface_color", "atomic", "all")
    # Create surface representation without electrostatic coloring
    cmd.show("surface", "all")
    # Make surface semi-transparent
    # cmd.set("transparency", 0.5, "all")


    # Add visual enhancements
    # cmd.util.cnc()
    cmd.orient()
    cmd.zoom("visible", 5)

    # Save session
    output_file = os.path.join(output_dir, f"{pdb_file}_colored_neg1.0.pse")
    cmd.save(output_file)
    print(f"Saved colored session to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Color PDB structure based on CRISPR scores")
    parser.add_argument("--input_folder",
                        required=True,
                        help="Path to results folder from CRISPRO")
    parser.add_argument("--pdb",
                        required=True,
                        help="PDB ID")
    parser.add_argument("--chain",
                        default="A",
                        help="Chain ID to color. Default is 'A'")

    args = parser.parse_args()

    main(
        input_folder=args.input_folder,
        pdb=args.pdb,
        chain_id=args.chain
    )
