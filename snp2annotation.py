'''
snps2annotation.py - assign SNPs to overlapping enriched annotations
===================================================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Uses the output from GAT to assign SNPs to overlapping
annotations using bedtools

Input
-----

SNP file: BED4 format, 4th column is the SNP ID

GAT results: table containing merged results from GAT.  This should contain column headers::
  <celltype_<statistic>
  Row names are the annotations tested

Annotation directory: location of BED files containing the features tested by GAT


Usage
-----

.. Example use case

Example::

   python snps2anotation.py

Type::

   python snps2annotation.py --help

for command line help.

Command line options
--------------------

Requirements:
* bedtools >= 2.2

TODO:
* documentation
* cleaner output
* split by SNP?
* summarise over annotations for each SNP
* include enrichment statistics?
* output/allow for SNPs that don't intersect any features

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pandas as pd
import re
import itertools
import subprocess
import os


def findSigResults(input_file, qval_thresh=0.01,
                   l2fold_thresh=2.0):
    '''
    Parse a merged results table from GAT,
    contains qval, pval, fold change and
    log2 fold change for each cell type (columns),
    and annotations (rows).

    Annotations not measured in a given cell type
    contain NA values.

    Arguments
    ---------
    input_file: string
      PATH to merged GAT results

    qval_thresh: float
      threshold for declaring statistically
      significant enriched annotations

    l2fold_thresh: float
      threshold for declaring biologically
      significantly enriched annotations

    Returns
    -------
    enriched_annots: dict
      dict object containing enriched annotations
      for each cell type (where relevant)    
    '''

    qval_re = re.compile("qval")
    l2fold_re = re.compile("l2")

    _df = pd.read_table(input_file, sep="\t",
                        header=0, index_col=0)
    qval_cols = [qv for qv in _df.columns if re.search(qval_re, qv)]
    l2f_cols = [fv for fv in _df.columns if re.search(l2fold_re, fv)]

    qval_df = _df.loc[:, qval_cols]
    l2_df = _df.loc[:, l2f_cols]
    
    # Pandas often reads scientific notation in as str
    # convert to float. Can't use pd.to_numeric over
    # data frame, need to apply to each columns series

    qval_df = qval_df.apply(lambda x: pd.to_numeric(x,
                                                    errors='coerce'))
    l2_df = l2_df.apply(lambda x: pd.to_numeric(x,
                                                errors='coerce'))

    # selct annotations breaching thresholds
    # separately, then merge them - makes it easier
    # to read rather than chaining operations

    sig_qvals = lambda q: q <= qval_thresh
    sig_l2fold = lambda l: l >= l2fold_thresh

    qval_df["annot"] = qval_df.index
    l2_df["annot"] = l2_df.index
    
    q_melt = pd.melt(qval_df, id_vars="annot")
    l_melt = pd.melt(l2_df, id_vars="annot")

    q_sigs = q_melt["value"].apply(sig_qvals)
    l_sigs = l_melt["value"].apply(sig_l2fold)

    q_sig_melt = q_melt.loc[q_sigs]
    l_sig_melt = l_melt.loc[l_sigs]

    # split stat from cell type for merging results
    q_sig_melt.loc[:, "cell_type"] = [qc.split("_")[0] for qc in q_sig_melt["variable"]]
    l_sig_melt.loc[:, "cell_type"] = [lc.split("_")[0] for lc in l_sig_melt["variable"]]

    merge_sig = pd.merge(q_sig_melt, l_sig_melt,
                         on=["annot", "cell_type"],
                         how='inner')

    merge_sig.columns = ["annot", "varX", "qvalue", "cell_type",
                         "varY", "l2fold"]
    merge_sig.drop(["varX", "varY"], inplace=True, axis=1)

    cells = set(merge_sig["cell_type"])

    enriched_annots = {}

    for cell in cells:
        cell_df = merge_sig[merge_sig["cell_type"] == cell]
        enriched_annots[cell] = cell_df["annot"].values
        

    return enriched_annots


def parseBedtools(bedout):
    '''
    Parse the output from bedtools intersect
    Expected format:

    contigA StartA EndA NameA contigB StartB EndB NameB
    .       .      .    .     .       .      .    . 
    .       .      .    .     .       .      .    .
    .       .      .    .     .       .      .    . 

    
    Arguments
    ---------
    bedout: string
      a single output line from bedtools intersect

    Returns
    -------
    snp_dict: dict
      dict object mapping annotations onto SNPs 

    '''

    snp_dict = {"snp": [],
                "chr": [],
                "start": [],
                "end": [],
                "annot_chr": [],
                "annot_start": [],
                "annot_end": [],
                "annot_name": []}

    rows = bedout.split("\n")
    for row in rows:
        if len(row):
            components = row.split("\t")
            snp_chr, snp_str, snp_end, snp_id = components[:4]
            annot_chr, annot_str, annot_end, annot_name = components[4:]

            snp_dict["snp"].append(snp_id)
            snp_dict["chr"].append(snp_chr)
            snp_dict["start"].append(snp_str)
            snp_dict["end"].append(snp_end)
            snp_dict["annot_chr"].append(annot_chr)
            snp_dict["annot_start"].append(annot_str)
            snp_dict["annot_end"].append(annot_end)
            snp_dict["annot_name"].append(annot_name)

    return snp_dict


def intersectSnpWithAnnotation(snp_bed, annotation_bed,
                               snp_dict):
    '''
    Using bedtools intersect to get the snps
    that overlap that enriched annotation

    Arguments
    ---------
    snp_bed: string
      PATH to BED4 format file containing SNP positions

    annotation_bed: string
      PATH to BED format file containing enriched annotations

    snp_dict: dict
      dict object containing SNPs as keys with container of
      overlapping annotations as values. Iteratively
      updated by successive calls to this function

    Returns
    -------
    snp_annotation: dict
      dict object with overlapping annotations for
      that SNP
    '''

    command = "bedtools intersect -a {} -b {} -wb "
    
    process = subprocess.Popen(command.format(snp_bed,
                                              annotation_bed),
                               shell=True, executable="/bin/bash",
                               stdout=subprocess.PIPE)

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise OSError(
            "-------------------------------------------\n"
            "Child was terminated by signal %i: \n"
            "The stderr was \n%s\n%s\n"
            "-------------------------------------------" %
            (-process.returncode, stderr, statement))

    return parseBedtools(stdout)
        

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--snp-bed", dest="snp_bed", type="string",
                      help="path to file containing SNP BED file")

    parser.add_option("--annotation-dir", dest="annot_dir", type="string",
                      help="path to annotations directory")

    parser.add_option("--q-value-threshold", dest="q_thresh",
                      type="float", help="q-value threshold to use "
                      "for selecting enriched annotations")

    parser.add_option("--l2fold-threshold", dest="l2_thresh",
                      type="float", help="log2 fold change "
                      "threshold to select enriched annotations")

    parser.add_option("--l2fold-direction", dest="l2_direct",
                      type="choice", choices=["up", "down",
                                              "absolute"],
                      help="direction of log2 fold change. If "
                      "set to absolute will ignore direction, "
                      "else will use enrich/depeleted only")

    parser.set_defaults(q_thresh=0.01,
                        l2_thresh=2.0,
                        l2_direct="up")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # infile is the merged GAT results
    infile = argv[-1]

    annots = findSigResults(infile)

    # also contains log files
    # make sure these are excluded
    dirlist = os.listdir(options.annot_dir)

    abs_path = os.path.abspath(options.annot_dir)

    snp_dict = {}
    counter = 0
    for cell in annots.keys():
        cell_re = re.compile(r"^{}.bed".format(cell))
        cell_files = [cx for cx in dirlist if re.match(cell_re, cx)]

        # ignore log files
        cell_annot = [os.path.join(abs_path,
                                   an) for an in cell_files if not re.search(r"log", an)]
        snp_annot = intersectSnpWithAnnotation(snp_bed=options.snp_bed,
                                               annotation_bed=cell_annot[0],
                                               snp_dict=snp_dict)

        cell_df = pd.DataFrame(snp_annot)
        cell_df["cell_type"] = cell
        
        if counter == 0:
            _df = cell_df
            counter += 1
        else:
            _df = _df.append(cell_df)

    # output in SNP position order and order headers
    # SNP first, then annotation
    out_df = _df.sort(["chr","start"])
    headers = ["chr", "start", "end", "snp", "annot_chr",
               "annot_start", "annot_end", "annot_name", "cell_type"]
    out_df = out_df[headers]

    out_df.to_csv(options.stdout, sep="\t",
                  index=None, header=None)
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
