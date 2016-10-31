#! /usr/bin/env python

from __future__ import division
import pysam
import os


def check_bam(bam, p=1, make_new_index=False):
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


def get_coverages(chrom, start, stop, bam, chrom_sizes, scaling):
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
    if scaling[chrom] > 0:
        norm_coverage = ((coverage / interval_length) / scaling[chrom]) * 10**6
    else:
        norm_coverage = 0
    counts = bam.count(chrom, start, stop)
    if counts > 0:
        rp100kpm = ((counts / interval_length) / scaling[chrom]) * 10**6
    else:
        rp100kpm = 0
    data_string = "\t".join([chrom, str(start), str(stop), str(coverage / interval_length), str(norm_coverage), str(rp100kpm)])
    return data_string


def normalization(options):
    """
    find total number of mapped reads for each chromosome
    to be used as a scaling factor
    when comparing between samples
    """
    stats = pysam.idxstats(options.file)
    norms = {}
    for i in stats:
        norms[i.rsplit()[0]] = int(i.rsplit()[2])        
    return norms


def main(options):
    # first check bam file is sorted and has an index
    chrom_sizes = check_bam(options.file)
    bam = pysam.AlignmentFile(options.file, "rb")
    scaling_factors = normalization(options)
    # now find the mean coverage in each bin and write to file
    with open(options.output, "w+") as outfile:
        outfile.write("chromosome\tstart\tstop\tcoverage\tnormalized_coverage\trpkm\n")
        for chrom in chrom_sizes.keys():
            start_position = 0
            while start_position < chrom_sizes[chrom]:
                interval_data = get_coverages(chrom,
                                              start_position,
                                              start_position + int(options.size),
                                              bam, chrom_sizes, scaling_factors)
                outfile.write(interval_data + "\n")
                start_position += int(options.size)
    bam.close()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="get_cov.py -- find coverage in X kb bins for each chromosome")
    parser.add_argument("-f", "--file", help="input bam file, must be sorted")
    parser.add_argument("-o", "--output", help="name of output file")
    parser.add_argument("-s", "--size", help="size of bins")
    options = parser.parse_args()

    main(options)
