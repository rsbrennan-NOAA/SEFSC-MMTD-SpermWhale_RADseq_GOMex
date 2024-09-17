
import sys

def process_file(input_file, output_file):
    chr_mapping = {} # create dictionary. empty
    chr_counter = 1 # start with number 1

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:

            # first separate each line into its parts
            parts = line.strip().split()
            if len(parts) != 4:
                continue

            # name these parts
            chr_name, _, _, position = parts

            # add chr_name and assigned value (number) to dictionary if it isn't in there
            if chr_name not in chr_mapping:
                chr_mapping[chr_name] = chr_counter
                chr_counter += 1

            chr_num = chr_mapping[chr_name]
            # use fstring to allow expressions
            new_second_column = f"{chr_num}:{position}"

            outfile.write(f"{chr_num} {new_second_column} 0 {position}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: GONE_map_format.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_file(input_file, output_file)
    print(f"File converson done. Output written to {output_file}")
