#!/usr/bin/env python3

import os
import re
import sys


def extract_accession(filename):
    """
    Extract first genome accession from filename.
    Supports GCA_ and GCF_.
    """
    match = re.search(r'(GC[AF]_\d+\.\d+)', filename)
    if match:
        return match.group(1)
    return None


def extract_metrics(filepath):
    identity = similarity = gaps = None

    with open(filepath) as f:
        for line in f:
            if line.startswith("# Identity:"):
                m = re.search(r'\(\s*([\d\.]+)%\)', line)
                if m:
                    identity = m.group(1)

            elif line.startswith("# Similarity:"):
                m = re.search(r'\(\s*([\d\.]+)%\)', line)
                if m:
                    similarity = m.group(1)

            elif line.startswith("# Gaps:"):
                m = re.search(r'\(\s*([\d\.]+)%\)', line)
                if m:
                    gaps = m.group(1)

    return identity, similarity, gaps



def main(input_folder, output_file):

    with open(output_file, "w") as out:
        out.write("gen\tIdentity\tSimilarity\tGaps\n")

        for file in os.listdir(input_folder):
            if not file.endswith(".txt"):
                continue

            filepath = os.path.join(input_folder, file)

            gen = extract_accession(file)
            identity, similarity, gaps = extract_metrics(filepath)

            if gen:
                out.write(
                    f"{gen}\t"
                    f"{identity or 'NA'}\t"
                    f"{similarity or 'NA'}\t"
                    f"{gaps or 'NA'}\n"
                )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_needle_folder.py <input_folder> <output.tsv>")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_file = sys.argv[2]

    main(input_folder, output_file)
