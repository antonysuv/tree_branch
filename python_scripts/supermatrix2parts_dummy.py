
##########################UNFINISHED
import argparse
import os
import sys

from Bio import AlignIO


def get_coord(part_file):
    coords = []
    with open(part_file, "r") as f:
        for line in f:
            line = line.strip()
            if "charset" in line:
                line = line.split(",")[1].split(";")[0].strip().split("-")
                coords.append(line)
    return coords


def unlink_loci(supermatrix_file, coord_pairs):
    a = AlignIO.read(supermatrix_file, "fasta")
    l_msa = []
    for i in coord_pairs:
        l_msa += [a[:, int(i[0]) - 1 : int(i[1])]]
    AlignIO.write(l_msa, supermatrix_file + ".splt", "fasta")


def main():
    parser = argparse.ArgumentParser(description="Split supermatrix by partitions")
    parser.add_argument(
        "--fasta", "-f", help="Supermatrix in FASTA", dest="FASTA", required=True
    )
    parser.add_argument(
        "--part", "-p", help="Partition file in NEXUS", dest="PART", required=True
    )
    args = parser.parse_args()

    my_coords = get_coord(args.PART)
    unlink_loci(args.FASTA, my_coords)


if __name__ == "__main__":
    main()