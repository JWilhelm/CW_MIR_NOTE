import os

# Define the directory and file path
directory_path = "../../01_dist_10A"
file_name = "cp2k.out"
file_path = os.path.join(directory_path, file_name)

# Function to read lines from 3rd to 1011 following "Mulliken" and perform calculations
def process_mulliken_data(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Initialize variables to store the sums
        sum1_to_351 = 0.0
        sum352_to_1010 = 0.0

        # Flag to indicate when to start reading lines
        read_lines = False

        for line in lines:
            if "Mulliken" in line:

                print(f"Sum of entries 1 to 351: {sum1_to_351}")
                print(f"Sum of entries 352 to 1010: {sum352_to_1010}")

                # Reset sums when a new "Mulliken" occurrence is found
                sum1_to_351 = 0.0
                sum352_to_1010 = 0.0
                read_lines = True
                i = 0
            elif read_lines:
                i += 1
                # Extract the 5th entry from the line
                entries = line.split()
                if len(entries) >= 5 and i > 2 and i < 1013:
                    entry_5 = float(entries[4])

                    # Update sums based on entry position
                    if i >= 353:
                        sum1_to_351 += entry_5
                    else:
                        sum352_to_1010 += entry_5

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return

# Call the function to process Mulliken data
process_mulliken_data(file_path)
