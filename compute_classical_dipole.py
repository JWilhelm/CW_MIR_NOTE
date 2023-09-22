import os

# Define the directory and file path
directory_path = "../01_dist_10A"
file_name = "cp2k.out"
file_path = os.path.join(directory_path, file_name)

# Function to read lines from 3rd to 1013th following "Mulliken"
def read_lines_after_mulliken(file_path, start_line=3, num_lines=1011):
    found_mulliken = False
    result = []

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if "Mulliken" in line:
                    found_mulliken = True
                    result.extend(lines[i + start_line - 1:i + start_line + num_lines - 1])
                    break

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None

    if found_mulliken:
        return result
    else:
        print("Mulliken not found in the file.")
        return None

# Call the function to read lines
mulliken_lines = read_lines_after_mulliken(file_path)

# Print the lines if found
if mulliken_lines:
    for line in mulliken_lines:
        print(line.strip())
