import argparse

def count_words(filename):
    with open(filename) as f:
        text = f.read()
    words = text.split()
    return len(words)

def main():
    parser = argparse.ArgumentParser(description="Count words in a file")
    parser.add_argument('--input', required=True, help="Input text file")
    parser.add_argument('--output', required=True, help="Output file for word count")
    args = parser.parse_args()

    word_count = count_words(args.input)
    with open(args.output, 'w') as f:
        f.write(f"Word count: {word_count}\n")

if __name__ == "__main__":
    main()
