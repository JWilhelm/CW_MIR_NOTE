import os
import subprocess

# Define the directory and file path
directory = "06_dist_3A"
file_name = "cp2k.out"
struc_name = "struc.xyz"
file_path = os.path.join("../../"+directory, file_name)
struc_path = os.path.join("../../"+directory,struc_name)

# Function to read lines from 3rd to 1011 following "Mulliken" and perform calculations
def process_mulliken_charges(file_path, directory, xyz_file):

    num_atoms, atomic_symbols, atomic_coordinates = xyz_file

    z_avg_tip_A = 0
    z_avg_surface_A = 0

    n_atom_tip = 352

    for i_atom in range(num_atoms):
        
       x, y, z = atomic_coordinates[i_atom]
       if i_atom <= n_atom_tip:
          z_avg_tip_A += z/n_atom_tip
       else:
          z_avg_surface_A += z/(num_atoms-n_atom_tip)

    print("z_avg_tip_A = "+str(z_avg_tip_A))
    print("z_avg_surface_A = "+str(z_avg_surface_A))

    subprocess.run("echo '' > "+directory+"_charge_tip", shell=True, check=True)
    subprocess.run("echo '' > "+directory+"_charge_surface", shell=True, check=True)
    subprocess.run("echo '' > "+directory+"_dipole_tip", shell=True, check=True)
    subprocess.run("echo '' > "+directory+"_dipole_surface", shell=True, check=True)
    subprocess.run("echo '' > "+directory+"_dipole_between_tip_and_surface", shell=True, check=True)
    subprocess.run("echo '' > "+directory+"_dipole", shell=True, check=True)

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Initialize variables to store the sums
        charge_tip = 0.0
        charge_surface = 0.0
        dipole_tip = 0.0
        dipole_surface = 0.0
        dipole_between_tip_and_surface = 0.0
        dipole = 0.0

        read_lines = False

        charge_tip_init = 0.0
        charge_surface_init = 0.0
        dipole_tip_init = 0.0
        dipole_surface_init = 0.0
        dipole_between_tip_and_surface_init = 0.0
        dipole_init = 0.0

        i_t = 0
        for line in lines:
            if "Mulliken" in line:

                if i_t == 1:
                  charge_tip_init = charge_tip
                  charge_surface_init = charge_surface
                  dipole_tip_init = dipole_tip
                  dipole_surface_init = dipole_surface
                  dipole_between_tip_and_surface_init = dipole_between_tip_and_surface
                  dipole_init = dipole

                print(f"q_tip = {charge_tip - charge_tip_init:.3f} q_surf = {charge_surface - charge_surface_init:.3f} " \
                      f"d_tip = {dipole_tip - dipole_tip_init:.3f} d_surf = {dipole_surface - dipole_surface_init:.3f} " \
                      f"d_tot = {dipole - dipole_init:.3f} " \
                      f"dipole_between_tip_and_surface = {dipole_between_tip_and_surface - dipole_between_tip_and_surface_init:.3f}")

                subprocess.run(f"echo '{charge_tip - charge_tip_init:.5f}' >> "        +directory+"_charge_tip", shell=True, check=True)
                subprocess.run(f"echo '{charge_surface - charge_surface_init:.5f}' >> "+directory+"_charge_surface", shell=True, check=True)
                subprocess.run(f"echo '{dipole_tip - dipole_tip_init:.5f}' >> "        +directory+"_dipole_tip", shell=True, check=True)
                subprocess.run(f"echo '{dipole_surface - dipole_surface_init:.5f}' >> "+directory+"_dipole_surface", shell=True, check=True)
                subprocess.run(f"echo '{dipole_between_tip_and_surface - dipole_between_tip_and_surface_init:.5f}' >> "+directory+"_dipole_between_tip_and_surface", shell=True, check=True)
                subprocess.run(f"echo '{dipole - dipole_init:.5f}' >> "                +directory+"_dipole", shell=True, check=True)

                # Reset sums when a new "Mulliken" occurrence is found
                charge_tip = 0.0
                charge_surface = 0.0
                dipole_tip = 0.0
                dipole_surface = 0.0
                dipole_between_tip_and_surface = 0.0
                dipole = 0.0

                read_lines = True
                i = 0
                i_t += 1
            elif read_lines:
                i += 1
                # Extract the 5th entry from the line
                entries = line.split()
                if len(entries) >= 5 and i > 2 and i < 1013:
                    entry_5 = float(entries[4])

                    x, y, z = atomic_coordinates[i-3]

                    # Update sums based on entry position
                    if i < 353:
                        charge_tip += entry_5
                        # 1 Angström = 1/0.529 atomic units
                        dipole_tip += (z-z_avg_tip_A)/0.529*entry_5
                        dipole_between_tip_and_surface += z_avg_tip_A/0.529*entry_5
                    else:
                        charge_surface += entry_5
                        # 1 Angström = 1/0.529 atomic units
                        dipole_surface += (z-z_avg_surface_A)/0.529*entry_5
                        dipole_between_tip_and_surface += z_avg_surface_A/0.529*entry_5

                    dipole += z/0.529*entry_5

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return

def read_xyz_file(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        num_atoms = int(lines[0].strip())
        atomic_symbols = []
        atomic_coordinates = []

        for line in lines[2:2 + num_atoms]:
            parts = line.split()
            if len(parts) >= 4:
                atomic_symbols.append(parts[0])
                atomic_coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])

        return num_atoms, atomic_symbols, atomic_coordinates

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None


xyz_file = read_xyz_file(struc_path)

# Call the function to process Mulliken data
process_mulliken_charges(file_path, directory, xyz_file)
