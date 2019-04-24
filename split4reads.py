#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @File    : biopython.py
# @Date    : 2019-01-24 10:16:58
# @Author  : luyang(luyang@novogene.com)


from __future__ import print_function, division

import argparse
import os
import sys
import time
from multiprocessing import Manager, Pool
from random import randint

import pysam


def index_ont_single(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    name_length = {read.name: len(read.sequence) for read in pysam.FastxFile(file) if len(read.sequence) > cutoff_min}
    origin_bases = sum(name_length.values())
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_ont_single(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.FastxFile(file):
        name = read.name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.sequence, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def index_ont_gz_single(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    name_length = {read.name: len(read.sequence) for read in pysam.FastxFile(file) if len(read.sequence) > cutoff_min}
    origin_bases = sum(name_length.values())
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_ont_gz_single(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.FastxFile(file):
        name = read.name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.sequence, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def index_ont(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    name_length = {read.name: len(read.sequence) for read in pysam.FastxFile(file) if len(read.sequence) > cutoff_min}
    origin_bases = sum(name_length.values())
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_ont(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.FastxFile(file):
        name = read.name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.sequence, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def index_ont_gz(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    name_length = {read.name: len(read.sequence) for read in pysam.FastxFile(file) if len(read.sequence) > cutoff_min}
    origin_bases = sum(name_length.values())
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_ont_gz(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.FastxFile(file):
        name = read.name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.sequence, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def index_pacb(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    zmw_rec = {}
    name_length = {}
    origin_bases = 0
    for read in pysam.FastxFile(file):
        length = len(read.sequence)
        origin_bases += length
        if length < cutoff_min:
            continue
        name = read.name
        zmw = name.split('/')[-2]
        if zmw in zmw_rec:
            if name_length[zmw_max] > length:
                pass
            else:
                name_length.pop(zmw_max)
                zmw_max = name
                name_length[name] = length
        else:
            zmw_max = name
            zmw_rec[zmw] = ''
            name_length[name] = length
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_pacb(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.FastxFile(file):
        name = read.name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.sequence, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def index_pacb_gz(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    zmw_rec = {}
    name_length = {}
    origin_bases = 0
    for read in pysam.FastxFile(file):
        length = len(read.sequence)
        origin_bases += length
        if length < cutoff_min:
            continue
        name = read.name
        zmw = name.split('/')[-2]
        if zmw in zmw_rec:
            if name_length[zmw_max] > length:
                pass
            else:
                name_length.pop(zmw_max)
                zmw_max = name
                name_length[name] = length
        else:
            zmw_max = name
            zmw_rec[zmw] = ''
            name_length[name] = length
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_pacb_gz(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.FastxFile(file):
        name = read.name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.sequence, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def index_pacb_bam(file):
    print('#####\t{0}\tStart indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    zmw_rec = {}
    name_length = {}
    origin_bases = 0
    for read in pysam.Samfile(file, 'rb', check_sq=False):
        length = read.query_length
        name = read.query_name
        origin_bases += length
        if length < cutoff_min:
            continue
        zmw = name.split('/')[-2]
        if zmw in zmw_rec:
            if name_length[zmw_max] > length:
                pass
            else:
                name_length.pop(zmw_max)
                zmw_max = name
                name_length[name] = length
        else:
            zmw_max = name
            zmw_rec[zmw] = ''
            name_length[name] = length
    print('#####\t{0}\tFinish indexing {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    return name_length, origin_bases


def choose_pacb_bam(file, name_length_filter, numsOfCopies):
    print('#####\t{0}\tStart choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()
    for read in pysam.Samfile(file, 'rb', check_sq=False):
        name = read.query_name
        if name in name_length_filter:
            queues[randint(0, numsOfCopies - 1)].put(''.join(['>', name, '\n', read.seq, '\n']))
    print('#####\t{0}\tFinish choosing in {1}'.format(time.asctime(time.localtime(time.time())), file), end='\n')
    sys.stdout.flush()


def write_fasta(index):
    print('#####\t{0}\tStart writing reads.{1}.fa'.format(time.asctime(time.localtime(time.time())), index + 1), end='\n')
    sys.stdout.flush()
    while True:
        read = queues[index].get()
        if not read:
            break
        out_files[index].write(read)
    out_files[index].close()
    print('#####\t{0}\tFinish writing reads.{1}.fa'.format(time.asctime(time.localtime(time.time())), index + 1), end='\n')
    sys.stdout.flush()


def choose_type(file, dataType):
    file_type = file.split('.')[-2:]
    gz_set = {'gz'}
    bgz_set = {'bgz'}
    fasta_set = {'fasta', 'fa'}
    fastq_set = {'fastq', 'fq'}
    bam_set = {'bam'}
    if dataType == 'PACB':
        if file_type[1] in fasta_set:
            return 'pacb_fun', 'fasta'
        elif file_type[1] in fastq_set:
            return 'pacb_fun', 'fastq'
        elif file_type[1] in bam_set:
            return 'pacb_bam_fun', 'bam'
        elif file_type[1] in bgz_set:
            if file_type[0] in fasta_set:
                return 'pacb_fun', 'fasta'
            elif file_type[0] in fastq_set:
                return 'pacb_fun', 'fastq'
            else:
                sys.exit('Error at filetype:\t{0}'.format('.'.join(file_type)))
        elif file_type[1] in gz_set:
            if file_type[0] in fasta_set:
                return 'pacb_gz_fun', 'fasta'
            elif file_type[0] in fastq_set:
                return 'pacb_gz_fun', 'fastq'
            else:
                sys.exit('Error at filetype:\t{0}'.format('.'.join(file_type)))
        else:
            sys.exit('Error at filetype:\t{0}'.format('.'.join(file_type)))
    elif dataType == 'ONT':
        if file_type[1] in fasta_set:
            if len(file) > 1:
                return 'ont_fun', 'fasta'
            else:
                return 'ont_single_fun', 'fasta'
        elif file_type[1] in fastq_set:
            if len(file) > 1:
                return 'ont_fun', 'fastq'
            else:
                return 'ont_single_fun', 'fastq'
        elif file_type[1] in bgz_set:
            if file_type[0] in fasta_set:
                if len(file) > 1:
                    return 'ont_fun', 'fasta'
                else:
                    return 'ont_single_fun', 'fasta'
            elif file_type[0] in fastq_set:
                if len(file) > 1:
                    return 'ont_fun', 'fastq'
                else:
                    return 'ont_single_fun', 'fastq'
            else:
                sys.exit('Error at filetype:\t{0}'.format('.'.join(file_type)))
        elif file_type[1] in gz_set:
            if file_type[0] in fasta_set:
                if len(file) > 1:
                    return 'ont_gz_fun', 'fasta'
                else:
                    return 'ont_gz_single_fun', 'fasta'
            elif file_type[0] in fastq_set:
                if len(file) > 1:
                    return 'ont_gz_fun', 'fastq'
                else:
                    return 'ont_gz_single_fun', 'fastq'
            else:
                sys.exit('Error at filetype:\t{0}'.format('.'.join(file_type)))
        else:
            sys.exit('Error at filetype:\t{0}'.format('.'.join(file_type)))
    else:
        sys.exit('Error at dataType:\t{0}'.format(dataType))


def determine_cutoff(name_length_total, expected_bases):
    length_sorted = sorted([v for k, v in name_length_total.items()], reverse=True)
    tmp = 0
    sub_sample_reads_flag = 0
    cutoff = 0
    for cutoff in length_sorted:
        tmp += cutoff
        if tmp >= expected_bases:
            sub_sample_reads_flag = 1
            break
    return cutoff, sub_sample_reads_flag


def write_shell(kbm2, dataDir, alignDir, dataType, numsOfCopies):
    wtdbg_K = int(1000 / numsOfCopies)
    wtdbg_n = int(1000 / numsOfCopies)

    shell_f = open(alignDir + '/00.shells/run1_kbm.sh.old', 'w')
    list_f = open(alignDir + '/kmaps.list', 'w')
    falist_f = open(dataDir + '/reads.list', 'w')

    for i in range(numsOfCopies):
        rname = dataDir + '/reads.' + str(i + 1) + '.fa'
        if dataType == 'PACB':
            print('{0} -m 300 -K {1} -n {2} -p 0 -k 15 -S 2 -t 10 -c -i {3} -fo {4}/reads.kbmap.{5}.{5} 2> alignDir/reads.kbmap.{5}.{5}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname, alignDir, i + 1), end='\n', file=shell_f)
            rsub_f = open(alignDir + '/00.shells/00.rsubs/rsub' + str(i + 1) + '.' + str(i + 1) + '.sh', 'w')
            print('{0} -m 300 -K {1} -n {2} -p 0 -k 15 -S 2 -t 10 -c -i {3} -fo {4}/reads.kbmap.{5}.{5} 2> alignDir/reads.kbmap.{5}.{5}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname, alignDir, i + 1), end='\n', file=rsub_f)
            rsub_f.close()
        elif dataType == 'ONT':
            print('{0} -m 300 -K {1} -n {2} -p 17 -k 0 -S 2 -t 10 -c -i {3} -fo {4}/reads.kbmap.{5}.{5} 2> alignDir/reads.kbmap.{5}.{5}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname, alignDir, i + 1), end='\n', file=shell_f)
            rsub_f = open(alignDir + '/00.shells/00.rsubs/rsub' + str(i + 1) + '.' + str(i + 1) + '.sh', 'w')
            print('{0} -m 300 -K {1} -n {2} -p 17 -k 0 -S 2 -t 10 -c -i {3} -fo {4}/reads.kbmap.{5}.{5} 2> alignDir/reads.kbmap.{5}.{5}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname, alignDir, i + 1), end='\n', file=rsub_f)
            rsub_f.close()
        else:
            sys.exit('Error at dataType:\t{0}'.format(dataType))

    for i in range(numsOfCopies - 1):
        for j in range(i + 1, numsOfCopies):
            rname1 = dataDir + '/reads.' + str(i + 1) + '.fa'
            rname2 = dataDir + '/reads.' + str(j + 1) + '.fa'
            if dataType == 'PACB':
                print('{0} -m 300 -K {1} -n {2} -p 0 -k 15 -S 2 -t 10 -c -i {3} -d {4} -fo {5}/reads.kbmap.{6}.{7} 2> alignDir/reads.kbmap.{6}.{7}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname1, rname2, alignDir, i + 1, j + 1), end='\n', file=shell_f)
                rsub_f = open(alignDir + '/00.shells/00.rsubs/rsub' + str(i + 1) + '.' + str(j + 1) + '.sh', 'w')
                print('{0} -m 300 -K {1} -n {2} -p 0 -k 15 -S 2 -t 10 -c -i {3} -d {4} -fo {5}/reads.kbmap.{6}.{7} 2> alignDir/reads.kbmap.{6}.{7}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname1, rname2, alignDir, i + 1, j + 1), end='\n', file=rsub_f)
                rsub_f.close()
            elif dataType == 'ONT':
                print('{0} -m 300 -K {1} -n {2} -p 17 -k 0 -S 2 -t 10 -c -i {3} -d {4} -fo {5}/reads.kbmap.{6}.{7} 2> alignDir/reads.kbmap.{6}.{7}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname1, rname2, alignDir, i + 1, j + 1), end='\n', file=shell_f)
                rsub_f = open(alignDir + '/00.shells/00.rsubs/rsub' + str(i + 1) + '.' + str(j + 1) + '.sh', 'w')
                print('{0} -m 300 -K {1} -n {2} -p 17 -k 0 -S 2 -t 10 -c -i {3} -d {4} -fo {5}/reads.kbmap.{6}.{7} 2> alignDir/reads.kbmap.{6}.{7}.log'.format(kbm2, wtdbg_K, wtdbg_n, rname1, rname2, alignDir, i + 1, j + 1), end='\n', file=rsub_f)
                rsub_f.close()
            else:
                sys.exit('Error at dataType:\t{0}'.format(dataType))

    for i in range(numsOfCopies):
        print('{0}/reads.kbmap.{1}.{1}'.format(alignDir, i + 1), end='\n', file=list_f)
        print('{0}/reads.{1}.fa'.format(dataDir, i + 1), end='\n', file=falist_f)
        for j in range(i + 1, numsOfCopies):
            print('{0}/reads.kbmap.{1}.{2}'.format(alignDir, i + 1, j + 1), end='\n', file=list_f)
    shell_f.close()
    list_f.close()
    falist_f.close()


def run(samtools, kbm2, dataDir, alignDir, inputList, dataType, sizes, numsOfCopies, threads):
    # prepare folder, inputList index tmp
    if not os.path.exists(alignDir + '/00.shells/00.rsubs'):
        os.makedirs(alignDir + '/00.shells/00.rsubs')
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)

    # Expected output bases
    # sizes of each copy * total number of copies
    expected_bases = sizes * numsOfCopies
    name_length_total = {}

    # read the input list
    print('#####\t{0}\tReading {1}'.format(time.asctime(time.localtime(time.time())), inputList), end='\n')
    with open(inputList) as f:
        for line in f.readlines():
            if line.strip():
                files.append(line.strip())

    # prepare output files handle
    print('#####\t{0}\tSpliting into {1} files'.format(time.asctime(time.localtime(time.time())), numsOfCopies), end='\n')
    for i in range(numsOfCopies):
        out_files.append(open(dataDir + '/reads.' + str(i + 1) + '.fa', 'w'))
        queues.append(Manager().Queue())

    # 判断reads是PACB(多文件)/ONT(单文件),fastq/fasta,plain/gz/bgz,返回方法.
    print('#####\t{0}\tChoosing process function'.format(time.asctime(time.localtime(time.time()))), end='\n')
    fun, read_format = choose_type(files[0], dataType)
    print('#####\t{0}\tReading the reads for the first time'.format(time.asctime(time.localtime(time.time()))), end='\n')
    origin_bases_total = 0

    if fun == 'pacb_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_pacb, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_pacb, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    elif fun == 'pacb_gz_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_pacb_gz, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_pacb_gz, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    elif fun == 'pacb_bam_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_pacb_bam, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_pacb_bam, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    elif fun == 'ont_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_ont, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_ont, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    elif fun == 'ont_gz_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_ont_gz, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_ont_gz, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    elif fun == 'ont_single_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_ont_single, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_ont_single, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    elif fun == 'ont_gz_single_fun':
        jobs = []
        pool = Pool(threads)
        for file in files:
            jobs.append(pool.apply_async(index_ont_gz_single, args=(file,)))
        for job in jobs:
            name_length, origin_bases = job.get()
            name_length_total.update(name_length)
            origin_bases_total += origin_bases
        pool.close()
        pool.join()
        # 确定cutoff,删减候选字典
        cutoff, sub_sample_reads_flag = determine_cutoff(name_length_total, expected_bases)
        name_length_filter = {k: v for k, v in name_length_total.items() if v >= cutoff}
        actual_bases = sum(name_length_filter.values())
        print('####\t{0}\tProcess for outputs'.format(time.asctime(time.localtime(time.time()))), end='\n')
        print('Total output numbers:\t{0}'.format(len(name_length_filter)), end='\n')
        print('#####\t{0}\tReading the reads for the second time'.format(time.asctime(time.localtime(time.time()))), end='\n')
        consumePool = Pool(numsOfCopies)
        for i in range(numsOfCopies):
            consumePool.apply_async(write_fasta, args=(i,))
        producePool = Pool(threads)
        for file in files:
            producePool.apply_async(choose_ont_gz_single, args=(file, name_length_filter, numsOfCopies,))
        producePool.close()
        producePool.join()
        for queue in queues:
            queue.put('')
        consumePool.close()
        consumePool.join()
    else:
        sys.exit('Error at dataType:\t{0}'.format(dataType))

    print('#####\t{0}\tWriting shells for kbm alignment'.format(time.asctime(time.localtime(time.time()))), end='\n')
    write_shell(kbm2, dataDir, alignDir, dataType, numsOfCopies)

    # print log记录
    print('#####\t{0}\tOutput the log file'.format(time.asctime(time.localtime(time.time()))), end='\n')
    with open(dataDir + '/split.log', 'w') as f:
        print('[data_log]', end='\n\n', file=f)
        print('input={0}'.format(inputList), end='\n\n', file=f)
        print('expected_index={0}, and actual_index={1}'.format(numsOfCopies, numsOfCopies), end='\n\n', file=f)
        print('original_total_bases={0}'.format(origin_bases_total), end='\n\n', file=f)
        print('reads_length_cut={0}'.format(cutoff), end='\n\n', file=f)
        print('expected_alignment_base={0}, about {1:.2f} of original total bases was expected used!'.format(expected_bases, expected_bases / origin_bases_total), end='\n\n', file=f)
        print('actual_alignment_base={0}, about {1:.2f} of original total bases was expected used!'.format(actual_bases, actual_bases / origin_bases_total), end='\n\n', file=f)
        print('alignment_used_reads_count={0}'.format(len(name_length_filter)), end='\n\n', file=f)
        print('alignment_used_reads_aveLen={0}'.format(int(actual_bases / len(name_length_filter))), end='\n\n', file=f)
        print('sub_sample_flag={0}'.format(sub_sample_reads_flag), end='\n\n', file=f)


def parse_args(argv):
    description = 'Split reads into equal fractions to run the auto_wtdbg using multiple threads'
    parser = argparse.ArgumentParser(description=description, )
    parser.add_argument('--samtools', required=True, type=str, help='The path of samtools')
    parser.add_argument('--kbm2', required=True, type=str, help='The path of kbm2')
    parser.add_argument('--dataDir', required=True, type=str, help='The dir of data output')
    parser.add_argument('--alignDir', required=True, type=str, help='The dir of align output')
    parser.add_argument('--inputList', required=True, type=str, help='The lists of inputs(fa/fasta/fasta.gz, fq/fastq/fastq.gz, bam)')
    parser.add_argument('--dataType', required=True, type=str, help='The datatype of inputs, PACB or ONT')
    parser.add_argument('--sizes', required=True, type=int, help='The sizes for each output')
    parser.add_argument('--numsOfCopies', required=True, type=int, help='The total number of copies for total output')
    parser.add_argument('--threads', required=True, type=int, help='The number of cpus for this program')
    return parser.parse_args(argv[1:])


def main(argv=sys.argv):
    args = parse_args(argv)
    run(**vars(args))


if __name__ == "__main__":
    # Global variables
    files = []
    out_files = []
    queues = []
    cutoff_min = 5000
    main(sys.argv)
