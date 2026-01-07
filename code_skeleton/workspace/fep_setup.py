import os
import subprocess
import argparse
import shutil
import glob

def run_command(cmd, check=True):
    """Run shell command with basic error handling."""
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=check)

def setup_fep(protein_pdb, ligand_pdb, charge_method="gas", net_charge=0):
    """
    Automates the complex system setup:
    1. Parameterize Ligand (Acpype)
    2. Parameterize Protein (pdb2gmx)
    3. Merge Structures and Topologies
    """
    
    # --- 1. Ligand Parameterization ---
    print("\n--- Step 1: Ligand Parameterization (Acpype) ---")
    # Outputs: ligand.acpype/ligand_GMX.gro, ligand_GMX.itp
    run_command(f"acpype -i {ligand_pdb} -c {charge_method} -n {net_charge}")
    
    # Find output directory (usually *.acpype)
    dirs = glob.glob("*.acpype")
    if not dirs:
        raise ValueError("Acpype failed to generate output directory.")
    lig_dir = dirs[0]
    
    # Copy essential files
    shutil.copy(os.path.join(lig_dir, f"{os.path.splitext(os.path.basename(ligand_pdb))[0]}_GMX.gro"), "ligand.gro")
    shutil.copy(os.path.join(lig_dir, f"{os.path.splitext(os.path.basename(ligand_pdb))[0]}_GMX.itp"), "ligand.itp")
    # Copy posre if exists, or ignore
    
    # --- 2. Protein Parameterization ---
    print("\n--- Step 2: Protein Parameterization (pdb2gmx) ---")
    # Outputs: protein.gro, topol.top
    # Force CHARMM27 and TIP3P for consistency with ligand (GAFF/Amber) - wait, inconsistency risk?
    # Acpype uses GAFF (General Amber Force Field).
    # Best practice: Use Amber99sb-ildn for protein if using GAFF for ligand.
    # Selecting option by number/name is tricky in headless mode without -ff.
    # Let's use amber99sb-ildn (often included). If not, charmm27 is acceptable for MVP testing.
    # We use -ff amber99sb-ildn -water tip3p
    
    # Pre-process Protein PDB to remove HETATMs (Ligands/Waters) that cause pdb2gmx to fail
    clean_protein = "clean_protein.pdb"
    with open(protein_pdb, "r") as fin, open(clean_protein, "w") as fout:
        for line in fin:
            if not line.startswith("HETATM") and "AN2" not in line: # Explicitly strip AN2 if it appears in ATOM records errantly
                fout.write(line)
    
    # Use clean protein
    run_command("gmx pdb2gmx -f " + clean_protein + " -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh")
    
    # --- 3. Merge Structures (.gro) ---
    print("\n--- Step 3: Merge Structures ---")
    # We must append the ligand atoms to the protein.gro file.
    # protein.gro has:
    # Title
    # Count
    # Line 1
    # Line N
    # Box params
    
    with open("protein.gro", "r") as f:
        p_lines = f.readlines()
        
    with open("ligand.gro", "r") as f:
        l_lines = f.readlines()
    
    p_count = int(p_lines[1].strip())
    l_count = int(l_lines[1].strip())
    total_count = p_count + l_count
    
    # helper to parse last resid and atomid from protein lines
    # .gro fixed format: 
    # resid(0-5) resname(5-10) atomname(10-15) atomid(15-20) x(20-28) y(28-36) z(36-44)
    last_p_line = p_lines[-2] # last atom line
    try:
        current_resid = int(last_p_line[:5].strip())
        current_atomid = int(last_p_line[15:20].strip())
    except ValueError:
        current_resid = 1
        current_atomid = p_count

    # Construct complex.gro
    with open("complex.gro", "w") as f:
        f.write("Protein-Ligand Complex\n")
        f.write(f" {total_count}\n")
        f.write("".join(p_lines[2:-1])) # Protein atoms
        
        # Renumber Ligand lines
        new_lig_lines = []
        ligand_offset_resid = current_resid 
        
        # Ligand usually starts at resid 1. We want to increment it.
        # But Acpype output is usually 1 residue.
        # We just add to the offset.
        
        for line in l_lines[2:-1]:
            # Simple parsing (fixed width)
            # We assume standard 8.3 format from acpype
            resid_str = line[:5]
            resname = line[5:10]
            atomname = line[10:15]
            atomid_str = line[15:20]
            coords = line[20:]
            
            # Increment
            current_atomid += 1
            # For residue, we assume all ligand atoms belong to the next residue
            # Or we check if input resid changed? typically ligand is 1 residue.
            # Let's just force it to current_resid + 1
            this_resid = ligand_offset_resid + 1
            
            # Format %5d%5s%5s%5d
            # Note: GROMACS supports hex for >99999 but python format needs care.
            # For MVP < 99999 atoms, %5d is fine.
            new_line = f"{this_resid:5d}{resname[:5]}{atomname[:5]}{current_atomid%100000:5d}{coords}"
            f.write(new_line)
            
        f.write(p_lines[-1]) # Box vector
        
    print(f"Created complex.gro with {total_count} atoms.")
    
    # --- 4. Merge Topologies (.top) ---
    print("\n--- Step 4: Merge Topologies ---")
    # We need to insert `#include "ligand.itp"` and update `[ molecules ]`
    
    with open("topol.top", "r") as f:
        top_lines = f.readlines()
        
    with open("topol.top", "w") as f:
        written_include = False
        in_molecules = False
        
        for line in top_lines:
            # Insert ligand.itp after forcefield.itp (or near top)
            if "; Include forcefield parameters" in line and not written_include:
                f.write(line)
                # Find the actual #include line following this comment usually
                continue
                
            if ".itp" in line and not written_include:
                f.write(line)
                f.write('\n; Include Ligand Parameters\n#include "ligand.itp"\n')
                written_include = True
                continue
            
            # Update Molecules section
            if "[ molecules ]" in line:
                in_molecules = True
                f.write(line)
                continue
            
            if in_molecules and line.strip() and not line.startswith(";"):
                # Usually writes: Protein_chain_A     1
                f.write(line)
                # After writing protein, we want to add ligand? 
                # Actually, wait until the end?
                # Usually best to append "Ligand   1" at the very end of file.
                pass
            else:
                f.write(line)
        
        # Determine Ligand Moleculetype name from itp
        # Robust parsing for moleculetype name
        molname = "MOL" # Default fallback
        found_directive = False
        with open("ligand.itp", "r") as f_itp:
            for line in f_itp:
                sline = line.strip()
                if "[ moleculetype ]" in sline:
                    found_directive = True
                    continue
                if found_directive:
                    # Skip empty lines and comments
                    if not sline or sline.startswith(";"):
                        continue
                    # The first real word is the molecule name
                    molname = sline.split()[0]
                    break
            
        f.write(f"{molname}             1\n")
        
    print(f"Updated topol.top with molecule {molname}")

    # --- 5. Solvate and Ionize (Standard GROMACS) ---
    print("\n--- Step 5: Box, Solvate, Ionize ---")
    run_command("gmx editconf -f complex.gro -o box.gro -bt cubic -d 1.0 -c")
    run_command("gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top")
    
    # Ions
    run_command("gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1")
    # Replace SOL with ions. 
    # Notes: -neutral automatically adds Na/Cl to neutralize.
    # We need a basic ions.mdp. We can create it temporary.
    
    run_command("echo 'SOL' | gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral")
    
    print("\n--- Setup Complete: system.gro and topol.top are ready ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", required=True)
    parser.add_argument("--ligand", required=True)
    args = parser.parse_args()
    
    # Create temp ions.mdp
    with open("ions.mdp", "w") as f:
        f.write("integrator = steep\nnsteps = 1000\n")
        
    setup_fep(args.protein, args.ligand)
