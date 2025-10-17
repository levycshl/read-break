
from read_break import FastqWriter, FastqReader

import os
from pathlib import Path



def convert_fasta(input_path, output_path,  head_name, from_base="C", to_base="T"):
    reader = FastqReader(input_path / f"{head_name}.R1.fastq.gz",
                         input_path / f"{head_name}.R2.fastq.gz")
    writer = FastqWriter(output_path, f"{head_name}_{from_base}{to_base}")

    for read_id, seq1, qual1, seq2, qual2 in reader:
        seq1 = seq1.replace(from_base, to_base)
        seq2 = seq2.replace(from_base, to_base)
        writer.write((read_id, seq1, qual1, seq2, qual2))
    writer.close()



input_dir = Path( "Z:/mnt/wigtop1/data/safe/levy/read_break/2025_09_19_alu_capture_deaminate" )
output_dir = Path( "Z:/mnt/wigtop1/data/safe/levy/read_break/2025_09_19_alu_capture_deaminate/converted" )
print( os.listdir(input_dir) )
head_names = [filename.stem.split(".", 1)[0] for filename in Path(input_dir).glob("*.fastq.gz")]
print(head_names)

for head_name in head_names:
    print(f"converting {head_name}, C->T")
    convert_fasta(input_dir, output_dir, head_name, from_base="C", to_base="T")
    print(f"converting {head_name}, G->A")
    convert_fasta(input_dir, output_dir, head_name, from_base="G", to_base="A")

print("done")
