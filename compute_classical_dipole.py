import os

# Define the directory and file path
directory_path = "../../06_dist_3A"
file_name = "cp2k.out"
struc_name = "struc.xyz"
file_path = os.path.join(directory_path, file_name)
struc_path = os.path.join(directory_path,struc_name)

# Function to read lines from 3rd to 1011 following "Mulliken" and perform calculations
def process_mulliken_charges(file_path, xyz_file):

    num_atoms, atomic_symbols, atomic_coordinates = xyz_file

    z_avg_tip = 0
    z_avg_surface = 0

    n_atom_tip = 352

    for i_atom in range(num_atoms):
        
       x, y, z = atomic_coordinates[i_atom]
       if i_atom <= n_atom_tip:
          z_avg_tip += z/n_atom_tip/0.529
       else:
          z_avg_surface += z/(num_atoms-n_atom_tip)/0.529

    print("z_avg_tip = "+str(z_avg_tip))
    print("z_avg_surface = "+str(z_avg_surface))

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Initialize variables to store the sums
        charge_tip = 0.0
        charge_surface = 0.0

        read_lines = False

        charge_tip_init = 0.0

        i_t = 0
        for line in lines:
            if "Mulliken" in line:

                if i_t == 1:
                  charge_tip_init = charge_tip

                print(f"Charge of tip: {charge_tip-charge_tip_init}")
#                print(f"Charge of surface: {charge_surface}")

                # Reset sums when a new "Mulliken" occurrence is found
                charge_tip = 0.0
                charge_surface = 0.0
                dipole_tip = 0.0
                dipole_surface = 0.0

                read_lines = True
                i = 0
                i_t += 1
            elif read_lines:
                i += 1
                # Extract the 5th entry from the line
                entries = line.split()
                if len(entries) >= 5 and i > 2 and i < 1013:
                    entry_5 = float(entries[4])

                    # Update sums based on entry position
                    if i >= 353:
                        charge_tip += entry_5
                        x, y, z = atomic_coordinates[i-3]
                        # 1 Angström = 1/0.529 atomic units
                        dipole_tip += z/0.529*entry_5
                    else:
                        charge_surface += entry_5


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
process_mulliken_charges(file_path, xyz_file)
