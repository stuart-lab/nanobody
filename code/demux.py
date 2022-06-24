import pysam
from Bio import SeqIO
import os
import gzip
from pathlib import Path
from argparse import ArgumentParser
import sys

# demux using i5 index only

# parse args
parser = ArgumentParser(description='Read demultiplexer')
parser.add_argument('--read1', help='Path to read1')
parser.add_argument('--read_i5', help='Path to read2')
parser.add_argument('--read2', help='Path to read3')
parser.add_argument('--tn5', help='Path to Tn5 i5 index FASTA file')
parser.add_argument('--output', help='Path to output directory')
args = parser.parse_args()


def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def open_fastq(f):
    if os.path.isfile(f):
        if f.endswith(".gz"):
            f_open = gzip.open(f, "rt") # rt read text mode (decodes binary gzip)
        else:
            f_open = open(f, "r")
        return(f_open)
    else:
        raise Exception("File not found")


def read_sample_barcodes(inpath, maxlen=8):
    all_bc = dict()
    for i in SeqIO.parse(open(inpath),'fasta'):
        all_bc[str(i.seq[:maxlen])] = i.id
    return(all_bc)


def get_entry(f):
    return([f.readline(), f.readline(), f.readline(), f.readline()])


def extract_barcodes_simple(sequence, bc1_len=8, spacer_len=14):
    # version with fixed-length barcodes
    sequence = sequence.strip("\n")
    tn5_bc = sequence[:bc1_len]
    cell_bc = sequence[bc1_len+spacer_len:]
    if len(cell_bc) != 16:
        print(cell_bc)
    return((tn5_bc, cell_bc))

tn5_i5_barcodes = read_sample_barcodes(args.tn5, maxlen=8) # using first 8 bases of barcode

r1_read = open_fastq(f=args.read1)
i5_read = open_fastq(f=args.read_i5)
r2_read = open_fastq(f=args.read2)

# create dictionary with file handles
outf = dict()
outpath = Path(args.output)
if not outpath.exists():
    os.mkdir(outpath)

unique_i5_marks = [tn5_i5_barcodes[key] for key in tn5_i5_barcodes.keys()]
unique_i5_marks.append('unknown')
for i5 in unique_i5_marks:
    # outputting uncompressed fastq is ~10x faster
    fname_1 = open(outpath / (i5 + ".R1.fastq"), "w+")
    fname_2 = open(outpath / (i5 + ".R2.fastq"), "w+")
    outf[i5] = (fname_1, fname_2)

x = 0
while True:
    r1_entry = get_entry(r1_read)
    r2_entry = get_entry(r2_read)
    i5_entry = get_entry(i5_read)

    if r1_entry[0] == '':
        break

    # get barcodes
    tn5_barcode, cell_barcode = extract_barcodes_simple(sequence=i5_entry[1], bc1_len=8)

    if tn5_barcode in tn5_i5_barcodes.keys():
        mark = tn5_i5_barcodes[tn5_barcode]
    else:
        # compute hamming distance
        hams = [hamming_distance(tn5_barcode, x) < 3 for x in tn5_i5_barcodes.keys()]
        if sum(hams) == 1:
            mark = tn5_i5_barcodes[list(tn5_i5_barcodes.keys())[hams.index(True)]]
        else:
            mark = "unknown"

    # add barcodes to r1 and r2 genomic
    bc_combination = "@" + cell_barcode + ":" + tn5_barcode + "+"
    r1_entry[0] = bc_combination + r1_entry[0][1:]
    r2_entry[0] = bc_combination + r2_entry[0][1:]

    if mark == "unknown":
        r1_outf = outf['unknown'][0]
        r2_outf = outf['unknown'][1]
    else:
        # write to correct output files based on barcode combination
        r1_outf = outf[mark][0]
        r2_outf = outf[mark][1]

    r1_outf.write("".join(r1_entry))
    r2_outf.write("".join(r2_entry))
    x += 1
    if x % 1e6 == 0:
        print("Processed " + str(int(x/1e6)) + " million reads", file=sys.stderr, end="\r")

# close all files
for i in outf.keys():
    outf[i][0].close()
    outf[i][1].close()
