import csv
import io
import math
import os
import pickle
import sys
import time
from itertools import groupby
from operator import itemgetter

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from PIL import Image


def reverse_compliment(seq: str) -> str:
    """Return the reverse compliment of a nucleotide sequence."""
    return ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[base] for base in reversed(seq)])


# https://www.biostars.org/p/710/
def fasta_iter(fh: io.TextIOWrapper) -> dict:
    """Yields a dict of header, sequence for each fasta record in a fasta file"""
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield {"header": header, "seq": seq}


def load_data() -> list:
    """Open bed and fa file, extract desired splice locations, and extract them from the fa file."""
    # trans_dict is used for changing the given names into standardized names.
    trans_dict = {"chr1": "1", "chr2": "2", "chr3": "3", "chr4": "4", "chr5": "5", "chr6": "6", "chr7": "7",
                  "chr8": "8", "chr9": "9", "chr10": "10", "chr11": "11", "chr12": "12", "chr13": "13", "chr14": "14",
                  "chr15": "15", "chr16": "16", "chr17": "17", "chr18": "18", "chr19": "19", "chrx": "x", "chry": "y"}
    # This try statement catches user error.
    try:
        with open(sys.argv[1]) as bed_file, open(sys.argv[2]) as fasta_file:
            fasta_records = []
            # Opens the bed file and splits into lists
            bed_file = list(csv.reader(bed_file, delimiter='\t'))
            # Changes the names of the chromosomes in bed file, does some light rearranging and formatting.
            bed_file = [[trans_dict[record[0].lower()], record[1], record[3][0:record[3].index(
                '\'')]] for record in bed_file]
            # Sorts the desired indices by chromosome, then by index in the chromosome.
            bed_file = sorted(bed_file, key=itemgetter(1))
            bed_file = sorted(bed_file, key=itemgetter(0))
            # This stores the desired indexes for each chromosome.
            indexable_bed_records = {'1': [], '2': [], '3': [], '4': [], '5': [], '6': [], '7': [], '8': [], '9': [],
                                     '10': [], '11': [], '12': [], '13': [], '14': [], '15': [], '16': [], '17': [],
                                     '18': [], '19': [], 'x': [], 'y': []}
            # Put each desired index into it's appropriate chromosome list.
            for record in bed_file:
                indexable_bed_records[record[0]].append([record[2], record[1]])
            # Loops over fasta records in the supplied fasta file
            for fasta_record in fasta_iter(fasta_file):
                # grabs the chromosome id
                chrom_id = fasta_record["header"][:fasta_record["header"].index(' ')].lower()
                # Some chromosomes are not desired, skip them.
                if chrom_id not in indexable_bed_records.keys():
                    continue
                # Grabs the indexes we want to extract from the chromosome.
                indexes = indexable_bed_records[chrom_id]
                # Grabs each index+/-10 from the sequence
                for index in indexes:
                    fasta_records.append([index[0], fasta_record["seq"][int(index[1]) - 10:int(index[1]) + 10]])
        # Returns a list of lists of format [5'/3',splice site sequence]
        return fasta_records
    # Catches user error.
    except (FileNotFoundError, IndexError) as e:
        if type(e) is IndexError:
            sys.stderr.write("Usage: {} bed_file fasta_file\n\tbed_file: The appropriate bed file. \n\t"
                             "fasta_file: The appropriate fasta file.\n".format(os.path.basename(__file__)))
        elif type(e) is FileNotFoundError:
            sys.stderr.write("One of the specified files was not found.\n")
        sys.exit(1)


def process_data(data: list) -> (list, list):
    """Splits a list of 3' and 5' splice sites into separate lists. Reverse-compliments all 3' splice sequences."""
    data = [[int(datum[0]), datum[1]] for datum in data]
    data_3 = [[datum[0], reverse_compliment(datum[1])] for datum in data if datum[0] is 3]
    data_5 = [datum for datum in data if datum[0] is 5]
    return data_3, data_5


def count_frequency(data: list) -> list:
    """Counts how many of each base appears in each column of the provided list of sequences."""
    counts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for _ in range(len(data[0][1]))]
    for i in range(len(data[0][1])):
        for datum in [datum[1] for datum in data]:
            counts[i][datum[i]] += 1
    return counts


def information(counts: list) -> list:
    """Calculates the bits of information and height of each character in the logo."""
    heights = []
    # magic
    e = (1 / math.log(2)) * ((4 - 1) / (2 * sum([counts[1][base] for base in "ACGT"])))
    for column_count in counts:
        relative_frqs = {base: column_count[base] / sum(column_count.values()) for base in "ACGT"}
        H = -1 * sum([relative_frqs[base] * math.log2(relative_frqs[base]) for base in "ACGT"])
        R = math.log2(4) - (H + e)
        heights.append({base: relative_frqs[base] * R for base in "ACGT"})
    # end magic
    return heights


def place_images(panel: plt.Axes, height_data: list) -> None:
    """Places images of A C T G on a panel in accordance with their information and height provided by fasta file."""
    # left right bottom top
    if "A_small.png" in os.listdir():
        pass
    else:
        for file in "ACGT":
            img = Image.open("{}.png".format(file))
            img = img.resize((img.size[0] // 4, img.size[1] // 4))
            img.save("{}_small.png".format(file))
    pic_dic = {'A': mpimg.imread('A_small.png'), 'C': mpimg.imread('C_small.png'),
               'G': mpimg.imread('G_small.png'), 'T': mpimg.imread('T_small.png')}
    for i, height_dict in enumerate(height_data):
        heights_sorted = sorted([[key, height_dict[key]] for key in height_dict.keys()], key=itemgetter(1))
        for j, base_value in enumerate(heights_sorted):
            if j is 0:
                bottom = 0
            else:
                bottom = sum(height_sort[1] for height_sort in heights_sorted[:j])
            panel.imshow(pic_dic[base_value[0].upper()], extent=[i - 10, i - 9, bottom, bottom + base_value[1]],
                         aspect="auto")

    return


def check_pickle() -> list:
    """Check to see if pickled data exists, if it does not, create it."""
    try:
        with open("data.pkl", mode='r+b') as open_pickle:
            data = pickle.load(open_pickle)
    except FileNotFoundError as _:
        data = load_data()
        with open("data.pkl", mode='w+b') as open_pickle:
            pickle.dump(data, open_pickle)
    return data


def main():
    """Load data, split data, calculate heights, create/modify plot and panels, save plot."""
    start = time.time()
    data = check_pickle()
    data_3, data_5 = process_data(data)
    three_counts = count_frequency(data_3)
    five_counts = count_frequency(data_5)
    three_letter_heights = information(three_counts)
    five_letter_heights = information(five_counts)
    plt.style.use('BME163.mplstyle')
    plt.figure(figsize=(6, 3))
    panel_left = plt.axes([.6 / 6, .9 / 3, 2.4 / 6, 1 / 3])
    panel_left.set(title="5'SS", xlim=(-10, 10), ylim=(0, 2), xlabel="Distance to\nSplice Site", ylabel="Bits")
    panel_left.vlines(0, panel_left.get_ylim()[0], panel_left.get_ylim()[1], linewidth=0.5, color="black")
    panel_right = plt.axes([3.3 / 6, 0.9 / 3, 2.4 / 6, 1 / 3])
    panel_right.set(title="3'SS", xlim=(-10, 10), ylim=(0, 2), xlabel="Distance to\nSplice Site", yticks=[])
    panel_right.vlines(0, panel_right.get_ylim()[0], panel_right.get_ylim()[1], linewidth=0.5, color="black")
    place_images(panel_left, five_letter_heights)
    place_images(panel_right, three_letter_heights)
    plt.savefig('McCreath_Benjamin_BME263_Week5.png', dpi=600)
    print("Time to complete: {}s".format(time.time() - start))


if __name__ == "__main__":
    main()
