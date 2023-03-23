import argparse
from multiprocessing.sharedctypes import Value
import subprocess
import glob
import os


parser = argparse.ArgumentParser(
    description="Processing multi-sample QC and quanlification for raw reads of RNA-seq.")
parser.add_argument("-i", "--input", type=str,
                    help="Specify your input directory including raw fastq files")
parser.add_argument("-f", "--suffix1", type=str, default=".r1",
                    help="Specify suffix of reads1, such as '.r1' for sample1_r1.fq. <PE model only>")
parser.add_argument("-b", "--suffix2", type=str, default=".r2",
                    help="Specfiy suffix of reads2, such as '.r2' for sample1_r2.fq. <PE model only>")
parser.add_argument("-o", "--output", type=str, default=".",
                    help="Sepecfiy output directory. <default = ./>")
parser.add_argument("-t", "--threads", type=int, default=1,
                    help="Specify threads you want to use. <default = 1>")
parser.add_argument("-z", action="store_true",
                    help="Specify your fq files are compressed")
parser.add_argument("-p", default=False, action="store_true",
                    help="Specify reads is paired")
parser.add_argument("-s", default=False, action="store_true",
                    help="Specfiy reads is single")
args = parser.parse_args()
input_d = args.input
out_d = args.output
threads = args.threads
suffix1 = args.suffix1
suffix2 = args.suffix2


_out = ["are", "is"]
print(
    f"Your input directory is {input_d}, output directory is {out_d}, threads {_out[0] if threads >1 else _out[1]} {threads}.")
if not input_d:
    raise TypeError("Please specify input directory!")
if args.s and (suffix1 != ".r1" or suffix2 != ".r2"):
    raise TypeError("Only PE model needs to specify -f and -b")
if args.p and not args.s:
    print("Model: paired reads")
elif args.s and not args.p:
    print("Model: single reads")
elif args.p and args.s:
    raise TypeError("Please specify your reads type!")
if not os.path.isdir(out_d):
    os.mkdir(out_d)


def runQC():
    qcdir = os.path.join(out_d, "QC")
    if not os.path.isdir(qcdir):
        os.mkdir(qcdir)
    for sample in glob.glob(f"{input_d}/*.fq*"):
        gzsuffix = ".fq.gz" if sample.endswith(".fq.gz") else ".fq"
        if args.s:
            sample_name = os.path.basename(sample).split(f".fq.gz")[0] if sample.endswith(
                ".fq.gz") else os.path.basename(sample).split(".fq")[0]
            outfile = os.path.join(qcdir, sample_name+".clean"+gzsuffix)
            print(f"Sample file is {sample}, outfile is {outfile}")
            # subprocess.Popen(f"fastp -i {sample} -o {outfile}", shell=True)
            os.system(f"fastp -i {sample} -o {outfile}")

        elif args.p:
            if suffix1 in sample:
                sample_name = os.path.basename(sample).split(f"{suffix1}.fq.gz")[0] if sample.endswith(
                    ".fq.gz") else os.path.basename(sample).split(f"{suffix1}.fq")[0]
                reads1_outfile = os.path.join(
                    qcdir, sample_name+suffix1+".clean"+gzsuffix)
                reads2_infile = os.path.join(os.path.dirname(
                    sample), sample_name+suffix2+gzsuffix)
                reads2_outfile = os.path.join(
                    qcdir, sample_name+suffix2+".clean"+gzsuffix)
                print(
                    f"Reads1 file is {sample}, outfile is {reads1_outfile}, reads2 file is {reads2_infile}, outfile is {reads2_outfile}")
                # subprocess.Popen(f"fastp -i {sample} -o {reads1_outfile} -I {reads2_infile} -O {reads2_outfile}", shell=True)
                os.system(
                    f"fastp -i {sample} -o {reads1_outfile} -I {reads2_infile} -O {reads2_outfile}")


if __name__ == "__main__":
    runQC()
