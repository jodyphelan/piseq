#! /usr/bin/python
import os
import argparse

scriptDir = os.path.dirname(os.path.realpath(__file__))

def trim(args):
    sample = args.sample
    threads = args.threads
    fastq1File = "fastq/" + sample + "_1.fastq.gz"
    fastq2File = "fastq/" + sample + "_2.fastq.gz"
    os.system("java -jar %s/trimmomatic.jar PE -threads %s -phred33 %s %s %s_1_trimmed_paired.txt %s_1_trimmed_unpaired.txt %s_2_trimmed_paired.txt %s_2_trimmed_unpaired.txt LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % (scriptDir,threads,fastq1File,fastq2File,sample,sample,sample,sample))

def bwa(args):
    sample = args.sample
    fastq1File = ""
    fastq2File = ""
    if args.f1==None and args.f2==None:
        fastq1File = "%s_1_trimmed_paired.txt" % sample
        fastq2File = "%s_2_trimmed_paired.txt" % sample
    else:
        fastq1File = args.f1
        fastq2File = args.f2
    refFile = args.ref
    threads = args.threads
    bamFile = sample + ".bam"
    mappingCMD = "bwa mem -t %s %s %s %s | samtools view -Sb - | samtools sort - -o %s" % (threads,refFile,fastq1File,fastq2File,bamFile)
    os.system("bwa aln -t %s %s %s > %s_1.sai" % (threads,refFile,fastq1File,sample))
    os.system("bwa aln -t %s %s %s > %s_2.sai" % (threads,refFile,fastq2File,sample))
    os.system("bwa sampe %s %s_1.sai %s_2.sai %s %s > %s.unsorted.sam" % (refFile,sample,sample,fastq1File,fastq2File,sample))
    os.system("samtools view -Sb %s.unsorted.sam | samtools sort -m 600M - -o %s.bam" % (sample,sample))
    os.system("rm %s.unsorted.sam %s*sai" % (sample,sample))

def samtoolsCall(args):
    refFile = args.ref
    vcfFile = args.sample + ".vcf.gz"
    bamFile = args.sample + ".bam"
    ploidy = args.ploidy
    os.system("samtools mpileup -ugf %s %s | bcftools call --ploidy %s -vmO z -o %s" % (refFile,bamFile,ploidy,vcfFile))



parser = argparse.ArgumentParser(description='Next-generation sequencing processing pipeline for the rasberry pi',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_bwa = subparsers.add_parser('bwa', help='Use BWA to map reads to a reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_bwa.add_argument('ref',help='Reference')
parser_bwa.add_argument('sample',help='Sample name')
parser_bwa.add_argument('-f1',help='Forward fastq')
parser_bwa.add_argument('-f2',help='Reverse fastq')
parser_bwa.add_argument('-nt','--threads', type=int, default='1', help='Number of threads')
parser_bwa.set_defaults(func=bwa)

parser_vcf = subparsers.add_parser('vcf', help='Use samtools to call variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_vcf.add_argument('ref',help='Reference')
parser_vcf.add_argument('sample',help='Sample name')
parser_vcf.add_argument('--ploidy',default=1,help='Ploidy')
parser_vcf.set_defaults(func=samtoolsCall)

parser_trim = subparsers.add_parser('trim', help='Use trimmomatic to trim reads',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_trim.add_argument('-nt','--threads', type=int, default='1', help='Number of threads')
parser_trim.add_argument('sample',help='Sample name')
parser_trim.set_defaults(func=trim)


args = parser.parse_args()
print(args)
args.func(args)
