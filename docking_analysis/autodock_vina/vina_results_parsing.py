#!/usr/bin/env python3
import pandas as pd
import re
import os

def parsing_vina(vinaFF) -> pd.DataFrame:
    with open(vinaFF) as f:
        vina_output = [line.strip() for line in f.readlines()]
        
        start_index = None
        ligand = receptor = ""
        for i, line in enumerate(vina_output):
            if line.startswith("mode"):
                start_index = i
                break
            elif line.startswith("Ligand") and re.search(r"\bLigand\b", line):
                ligand = line.split()[1].split(".")[0]
            elif line.startswith("Ligand") and re.search(r"\bLigands\b", line):
                ligands = vina_output[i+1:i+3]
                ligand = [ligand.split()[-1].split(".")[0] for ligand in ligands]
                ligand = "_".join(ligand)
            elif line.startswith("Rigid receptor"):
                receptor = line.split()[2].split("-")[2]
           
        docking_data = []
        if start_index is not None:
            aff_binding = vina_output[start_index + 3:]
            for record in aff_binding:
                parts = record.split()
                mode = parts[0]
                affinity = parts[1]
                rmsd_lb = parts[2]
                rmsd_ub = parts[3]
                if ligand and receptor is not None:
                    docking_data.append({'mode': mode, "receptor": receptor, "ligand": ligand, 'affinity': affinity, 'rmsd_lb': rmsd_lb, 'rmsd_ub': rmsd_ub})
        docking_df = pd.DataFrame(docking_data)

        return docking_df

# Receptor: WT, PRPP, PRPP with 5FU
receptor_ligand_docking = []
for vina_output in os.listdir("./input/"):
    vina_output = os.path.join("input", vina_output)
    docking_output = parsing_vina(vina_output)
    receptor_ligand_docking.append(docking_output.loc[:, ["receptor", "ligand", "affinity"]])

# Concatenate all docking output in the list row-wise
merged_df = pd.concat(receptor_ligand_docking, ignore_index = True)
merged_df = merged_df.sort_values(by = ["receptor", "ligand"])
merged_df.reset_index(drop = True, inplace = True)

# Matching Internal_ID to Strains_ID
merged_df.replace(
    {
        "192R4": "R4",
        "R214T": "R1"
    }, inplace = True)

# Save in csv format
merged_df.to_csv("binding_affinity.csv", index = False)