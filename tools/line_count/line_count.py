import argparse

def count_lines(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    return len(lines)

def main():
    parser = argparse.ArgumentParser(description="Count lines in a file")
    parser.add_argument('--input', required=True, help="Input text file")
    parser.add_argument('--output', required=True, help="Output file for line count")
    args = parser.parse_args()

    line_count = count_lines(args.input)
    with open(args.output, 'w') as f:
        f.write(f"Line count: {line_count}\n")

if __name__ == "__main__":
    main()