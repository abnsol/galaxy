import sys

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    with open(input_file) as fin, open(output_file, "w") as fout:
        for line in fin:
            fout.write(line.upper())  # Example: convert to uppercase

if __name__ == "__main__":
    main()