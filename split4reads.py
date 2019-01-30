#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @File    : biopython.py
# @Date    : 2019-01-24 10:16:58
# @Author  : luyang(luyang@novogene.com)
from Bio import SeqIO
import argparse
import gzip
import sys
from multiprocessing import Pool
from random import randint

# Global variables
gzip_files = []
files = []
read_format = ''
cutoff = 0
expected_bases = 0



def index_read(file):
    idx = file + '.idx'
    reads = SeqIO.index_db(idx, file, read_format)
    length = []
    for read in reads:
        length.append(len(read.seq))
    return length


def choose_read(file):
    idx = file + '.idx'
    for read in SeqIO.index_db(idx, file, read_format):
        if len(read.seq) > cutoff:
            gzip_files[randint(1, len(gzip_files))].write(read.format(read_format))


def index_read_gz(file):
    length = []
    with gzip.open(file, 'rt') as handle:
        for read in SeqIO.parse(handle, read_format):
            length.append(len(read.seq))
    return length


def choose_read_gz(file):
    with gzip.open(file, 'rt') as handle:
        for read in SeqIO.parse(handle, read_format):
            if len(read.seq) > cutoff:
                gzip_files[randint(1, len(gzip_files))].write(read.format(read_format))


def choose_type(file):
    file_type = file.split('.')[-2:]
    gzip_set = ('gz', 'gzip')
    bgzf_set = ('bgz')
    fasta_set = ('fasta', 'fa')
    fastq_set = ('fastq', 'fq')
    if file_type[1] in gzip_set:
        if file_type[0] in fasta_set:
            return index_read_gz, choose_read_gz, 'fasta'
        elif file_type[0] in fastq_set:
            return index_read_gz, choose_read_gz, 'fastq'
        else:
            sys.exit('Unknown filetype')
    elif file_type[1] in bgzf_set:
        if file_type[0] in fasta_set:
            return index_read, choose_read, 'fasta'
        elif file_type[0] in fastq_set:
            return index_read, choose_read, 'fastq'
        else:
            sys.exit('Unknown filetype')
    elif file_type[1] in fasta_set:
        return index_read, choose_read, 'fasta'
    elif file_type[1] in fastq_set:
        return index_read, choose_read, 'fastq'
    else:
        sys.exit('Unknown filetype')


def determine_cutoff(length_total):
    length_total.sort(reverse=True)
    length_sum = sum(length_total)
    tmp = 0
    for cutoff in length_total:
        tmp += cutoff
        if tmp > length_sum * 0.6:
            break
    return cutoff


def run(samtools, kbm2, dataDir, alignDir, inputList, dataType, sizes, numsOfCopies, threads):

    # Expected output bases
    # sizes of each copy * total number of copies
    expected_bases = sizes * numsOfCopies

    # read the input list
    with open(sys.argv[1]) as f:
        for line in f.readlines():
            files.append(line.strip())

    for i in range(1, numsOfCopies+1):
        gzip_files.append(gzip.open(str(i) + '.gz', 'at'))

    # genome_size = int(sys.argv[2])
    genome_size = 0

    # 判断文件类型,fastq/fasta,plain/gz/bgz,返回方法
    index_fun, choose_fun, read_format = choose_type(files[0])
    length_total = []

    with Pool(threads) as pool:
        jobs = []
        for file in files:
            jobs.append(pool.apply_async(index_fun, args=(file,)))
        for job in jobs:
            length_total.extend(job.get())
        pool.close()
        pool.join()
    cutoff = determine_cutoff(length_total)

    print(cutoff)
    print('all sort')

    with Pool(threads) as pool:
        for file in files:
            pool.apply_async(choose_fun, args=(file,))
        pool.close()
        pool.join()

    for _ in gzip_files:
        _.close()

def parse_args(argv):
    
    description = 'Split reads into equal fractions to run the auto_wtdbg using multiple threads'
    parser = argparse.ArgumentParser(
            description=description,
    )
    parser.add_argument(
            '--samtools', required=True, type=str,
            help='The path of samtools')
    parser.add_argument(
            '--kbm2', required=True, type=str,
            help='The path of kbm2')
    parser.add_argument(
            '--dataDir', required=True, type=str,
            help='The dir of data output')
    parser.add_argument(
            '--alignDir', required=True, type=str,
            help='The dir of align output')
    parser.add_argument(
            '--inputList', required=True, type=str,
            help='The lists of inputs(fa/fasta/fasta.gz, fq/fastq/fastq.gz, bam)')
    parser.add_argument(
            '--dataType', required=True, type=str,
            help='The datatype of inputs, PACB or ONT')
    parser.add_argument(
            '--sizes', required=True, type=int,
            help='The sizes for each output')
    parser.add_argument(
            '--numsOfCopies', required=True, type=int,
            help='The total number of copies for total output')
    parser.add_argument(
            '--threads', required=True, type=int,
            help='The number of cpus for this program')

    return parser.parse_args(argv[1:])


def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))

if __name__ == "__main__":
    main(sys.argv)

