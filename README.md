# snippy
Install
=======

```
pip install pysam mappy intervaltree
```

Run
===
```
./snippy.py --tsv out.tsv -o out.bam reference.minimap2.idx barcode01.fastq 

```

Output
======

- Indexed BAM file
- TSV file (optional)

Example:
```
id                                    info                                                                                                                                                                        barcode  length  subread_number  sum_subread_length  avg_subread_length
65061d97-5003-4a5a-a64a-3a00a9ea7c9e  runid=626e32a86d1a696e7c31d8a37dd46c7f25a4f63f read=19 ch=481 start_time=2019-03-07T10:22:29Z flow_cell_id=FAK65965 protocol_group_id=snip_enrich sample_id=snip_enrich     NA       1364    7               1154                164.85714285714286
3ec6ffb2-ddb1-4968-9a95-118ad6774884  runid=626e32a86d1a696e7c31d8a37dd46c7f25a4f63f read=45 ch=19 start_time=2019-03-07T10:22:38Z flow_cell_id=FAK65965 protocol_group_id=snip_enrich sample_id=snip_enrich      NA       253     1               152                 152.0

```

Parameters
==========
```
usage: snippy.py [-h] [-V]
                 [-l {DEBUG,INFO,WARNING,ERROR,CRITICAL,debug,info,warning,error,critical}]
                 [-p] -o OUT [--tsv TSV_PATH]
                 TARGET_IDX [QUERY_FA [QUERY_FA ...]]

Command line interface to telemap

positional arguments:
  TARGET_IDX            Target/reference minimap2 index file
  QUERY_FA              Query files

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  -l {DEBUG,INFO,WARNING,ERROR,CRITICAL,debug,info,warning,error,critical}, --log {DEBUG,INFO,WARNING,ERROR,CRITICAL,debug,info,warning,error,critical}
                        Print debug information
  -p, --primary-only    Report only primary alignments
  -o OUT, --output OUT  Output BAM file
  --tsv TSV_PATH        Output TSV file
```
