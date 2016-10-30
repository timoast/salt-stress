#! /usr/bin/env python

from __future__ import division
import pysam
import os


def check_bam(bam, p, make_new_index=False):
    """
    Sort and index bam file
    returns dictionary of chromosome names and lengths
    """
    # check if sorted
    test_head = pysam.AlignmentFile(bam, 'rb')
    chrom_sizes = {}
    p = str(p)
    for i in test_head.header['SQ']:
        chrom_sizes[i['SN']] = int(i['LN'])
    try:
        test_head.header['HD']['SO']
    except KeyError:
        print '  sorting bam file'
        pysam.sort('-@', p, bam, 'sorted.temp')
        os.remove(bam)
        os.rename('sorted.temp.bam', bam)
    else:
        if test_head.header['HD']['SO'] == 'coordinate':
            pass
        else:
            print '  sorting bam file'
            pysam.sort('-@', p, bam, 'sorted.temp')
            os.remove(bam)
            os.rename('sorted.temp.bam', bam)
    test_head.close()
    # check if indexed
    if '{}.bai'.format(bam) in os.listdir('.') and make_new_index is False:
        pass
    else:
        print '  indexing bam file'
        pysam.index(bam)
    return chrom_sizes


def get_coverages(chrom, start, stop, bam, chrom_sizes):
    """
    find average coverage in given region
    """
    coverage = 0
    if chrom not in chrom_sizes.keys():
        raise Exception('Chromosome names do not match bam file')
    else:
        pass
    stop = chrom_sizes[chrom] if stop > chrom_sizes[chrom] else stop
    interval_length = stop - start
    for read in bam.pileup(chrom, start, stop):
        coverage += read.n
    data_string = "\t".join([chrom, str(start), str(stop)]) + "\t"  + str(coverage / interval_length)
    return data_string


def main(options):
    # first check bam file is sorted and has an index
    chrom_sizes = check_bam(options.file)
    bam = pysam.AlignmentFile(option.file, "rb")
    # now find the mean coverage in each bin and write to file
    with open(options.output, "w+") as outfile:
        for chrom in range(len(chrom_sizes)):
            start_position = 1
            while start_position < chrom_sizes[chrom]:
                interval_data = get_coverages(str(chrom+1),
                                              start_position,
                                              start_position + options.size,
                                              bam, chrom_sizes)
                outfile.write(interval_data + "\n")
                start_position += options.size
    bam.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="get_cov.py -- find coverage in X kb bins for each chromosome")
    parser.add_argument("-f", "--file", help="input bam file, must be sorted")
    parser.add_argument("-o", "--output", help="name of output file")
    parser.add_argument("-s", "--size", help="size of bins")
    options = parser.parse_args()

    main(options)
