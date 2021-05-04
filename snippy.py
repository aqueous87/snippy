from __future__ import print_function

import os

#from kerneltree import IntervalTree
from intervaltree import Interval, IntervalTree
import argparse
import pysam
import sys
import uuid
import logging
import mappy as mp

__version__ = "1.0.0"


name2id = {}


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-V", "--version",
                        action="version",
                        version="%(prog)s " + __version__)
    parser.add_argument("-l", "--log",
                        dest="log",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL',
                                 'debug', 'info', 'warning', 'error', 'critical'],
                        default="INFO",
                        help="Print debug information")
    parser.add_argument("-p", "--primary-only",
                        action="store_true",
                        dest="PRIMARY_ONLY",
                        help="Report only primary alignments")
    parser.add_argument("-o", "--output",
                        dest="OUT",
                        type=str,
                        required=True,
                        help="Output BAM file")
    parser.add_argument("--tsv",
                        dest="TSV_PATH",
                        type=str,
                        default=None,
                        help="Output TSV file")
    parser.add_argument("TARGET_IDX",
                        type=str,
                        nargs=1,
                        help="Target/reference minimap2 index file")
    parser.add_argument("QUERY_FA",
                        type=str,
                        nargs="*",
                        help="Query files")

    args = parser.parse_args(argv)

    return args


def _revcomp(seq):
    """
    Compute reverse complement. Supports full FASTA alphabet and assumes
    that seq does not contain non FASTA characters.

    :param seq: str
    :return: str
    """
    if seq is not None:
        return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd',
                                           'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]
    else:
        return None


def _parse_reads_info(comment):
    """
    Parse ONT FASTQ header

    :param comment: str
    :return: dict
    """
    single_read_info = {'runid': None}
    if comment:
        for pair in comment.split(" "):
            cols = pair.split('=')
            if len(cols) == 2:
                single_read_info[cols[0]] = cols[1]
    return single_read_info


def print_read(samfile, name, comment, quals, seq, hit):
    samfile.write(
        _hit2sam(uuid.uuid4().hex, seq, len(seq), quals, comment, name,
                 hit=hit))


def filter_intervals(intervals):
    it = IntervalTree()

    intervals_filtered = []
    for start, end in intervals:
        #if it.search(start, end):
        if it.overlap(start, end):
            pass
        else:
            it.addi(start, end, 1)
            #it.add(start, end, 1)
            intervals_filtered.append((start, end))
    return sorted(intervals_filtered, key=lambda tup: tup[0])


def map_intervals(aligner, samfile, name, comment, qual, seq, intervals, only_primary):
    reads_written = 0
    for start, end in intervals:
        sub_seq = seq[start:end]
        sub_qual = qual[start:end]
        hits = list(aligner.map(sub_seq))
        for hit in hits:
            if not only_primary or hit.is_primary:
                print_read(samfile, name, comment, sub_qual, sub_seq, hit)
                reads_written += 1
    return reads_written


def get_splits(aligner, seq, intervals, offset=0, depth=0):
    if len(seq) < 100:
        return
    read_count = 0
    hits = list(aligner.map(seq))
    for hit in hits:  # traverse alignments
        q_st = hit.q_st
        q_en = hit.q_en
        intervals.append((q_st + offset, q_en + offset))


def _hit2sam(name, seq, seq_len, qual, comment, orig_name, hit=None):
    """
    Convert minimap2 hit object to pysam SAM record

    :param name: Read name
    :type name: str
    :param seq:  Read sequence
    :type seq: str
    :param seq_len: Read length
    :type seq_len: int
    :param qual: Read qualities
    :type qaul: list
    :param comment: Read comment
    type comment: str
    :param hit: Minimap2 hit
    :return: AlignedSegment
    """
    segment = pysam.AlignedSegment()
    segment.query_name = name

    segment.query_sequence = seq
    #     print(qual)
    if qual:
        segment.query_qualities = qual

    segment.set_tag("CO", comment)
    segment.set_tag("ID", orig_name)

    if not hit:
        segment.is_unmapped = True
        return segment

    segment.mapping_quality = hit.mapq
    segment.reference_id = name2id[hit.ctg]
    segment.reference_start = hit.r_st
    if not hit.is_primary:
        segment.is_supplementary = True

    if hit.strand < 0:
        segment.is_reverse = True
        segment.query_sequence = _revcomp(seq)
        if qual:
            segment.query_qualities = qual[::-1]

    # Add clipping to CIGAR string
    c1 = ""
    if hit.q_st > 0:
        c1 = "{}S".format(hit.q_st)
    c2 = ""
    if hit.q_en < seq_len:
        c2 = "{}S".format(hit.q_en - seq_len)

    if segment.is_reverse:
        segment.cigarstring = "{}{}{}".format(c2, hit.cigar_str, c1)
    else:
        segment.cigarstring = "{}{}{}".format(c1, hit.cigar_str, c2)

    segment.set_tag("NM", hit.NM)
    #     segment.set_tag("MD", str(hit.md))

    return segment


def print_tsv_header(tsv):
    print('id',
          'info',
          'barcode',
          'length',
          'subread_number',
          'sum_subread_length',
          'avg_subread_length',
          file=tsv,
          sep='\t')


def print_tsv(read_stats, tsv):
    print(read_stats['id'],
          read_stats['info'],
          read_stats['barcode'],
          read_stats['length'],
          read_stats['subread_number'],
          read_stats['sum_subread_length'],
          read_stats['avg_subread_length'],
          file=tsv,
          sep='\t')


def map(args):

    # Check input files
    if not os.path.exists(args.TARGET_IDX[0]):
        logging.error("Could not open: {}".format(args.TARGET_IDX[0]))

    if not args.QUERY_FA:
        args.QUERY_FA = ['/dev/stdin']

    aligner = mp.Aligner(args.TARGET_IDX[0],
                         preset="map-ont")  # load or build index
    if not aligner: raise Exception("ERROR: failed to load/build index")

    references = []
    lengths = []
    for i, name in enumerate(aligner.seq_names):
        name2id[name] = i
        references.append(name)
        lengths.append(len(aligner.seq(name)))

    # Open output BAM
    samfile = pysam.AlignmentFile("{}_tmp.bam".format(args.OUT), "wb",
                                  reference_names=references,
                                  reference_lengths=lengths)

    tsv = None
    if args.TSV_PATH:
        tsv = open(args.TSV_PATH, "w")
        print_tsv_header(tsv)

    read_n = 0
    written_n = 0
    sum_subreads_n = 0
    failed_n = 0

    for query_fa in args.QUERY_FA:
        logging.info("Processing {}".format(query_fa))
        with pysam.FastxFile(query_fa) as fh:

            for entry in fh:
                try:
                    name = entry.name
                    seq = entry.sequence
                    comment = entry.comment
                    qual = entry.get_quality_array()

                    intervals = []
                    # Split concatenated read into individual reads
                    get_splits(aligner, seq, intervals)
                    intervals_filtered = filter_intervals(intervals)

                    total_length = 0
                    for start, end in intervals_filtered:
                        total_length += (end - start)

                    sum_subreads_n += len(intervals_filtered)

                    written_n += map_intervals(aligner, samfile, name, comment, qual, seq, intervals_filtered, args.PRIMARY_ONLY)
                    read_n += 1

                    # Print TSV stats per read
                    if tsv:
                        info = _parse_reads_info(comment)

                        read_stats = {"length": len(seq),
                                      "subread_number": len(intervals_filtered),
                                      "sum_subread_length": total_length,
                                      "avg_subread_length": total_length / max(1, len(intervals_filtered)),
                                      "id": name,
                                      "info": comment,
                                      "barcode": info.get('barcode', "NA")
                                      }

                        print_tsv(read_stats, tsv)
                except RecursionError as e:
                    logging.warning("Ups something went wrong with read {} "
                                    "(length {}). Ignore if it is not happening "
                                    "too often".format(name, len(seq)))
                    logging.warning(e)
                    failed_n += 1

    samfile.close()
    if tsv:
        tsv.close()

    logging.info("{} concatenated reads found in FASTQ file".format(read_n))
    logging.info("{} subreads reads written to BAM file".format(written_n))
    logging.info("{} failed".format(failed_n))

    logging.info("Found on average {} subreads per concatenated read".format(round(sum_subreads_n / read_n, 2)))

    logging.info("Sorting BAM file")
    pysam.sort("-o", args.OUT,
               "{}_tmp.bam".format(args.OUT))
    logging.info("Deleting temp file")
    os.unlink("{}_tmp.bam".format(args.OUT))
    logging.info("Indexing BAM file")
    pysam.index(args.OUT)


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log.upper())
    logging.basicConfig(level=numeric_level, format='%(message)s')

    map(args)


if __name__ == '__main__':
    main()


