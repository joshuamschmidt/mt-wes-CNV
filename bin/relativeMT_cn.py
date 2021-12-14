#!/usr/bin/env python
import argparse
import os
import sys
import textwrap
import mmap
import numpy as np
import pandas as pd
from difflib import SequenceMatcher

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='relativeMTcn',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        Per Sample MosDepth regions to reative MTDNA CN
                        ------------------------------------------------
                        Generates a table of sample + MTDNA CN estimates

                        Column headers are sample names derived from
                        the list of input files and the common file
                        {suffix} e.g. SAMPLE_1{.regions.bed.gz}.
                        Works just as well when data is mean coverage.

                        !!!!!Python needs to be >=3.9!!!!!!
                        '''),
                        add_help=False,
                        epilog="Questions, bugs etc?\njoshmschmidt1@gmail.com\ngithub.com/joshuamschmidt")
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

# Add back help
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=argparse.SUPPRESS,
    help='show this help message and exit'
)

parser.add_argument('files', metavar='FS', type=str, nargs='+',
                    help='files to convert')

required.add_argument('--suffix', type=str, dest='file_suffix',
                    help='common file suffix')
parser.set_defaults(file_suffix=".mosdepth.summary.txt")

optional.add_argument('--threads', type=int, dest='threads',
                       help='multithreaded mode with int threads:NOT IMPLEMENTED')
optional.add_argument('--out', type=argparse.FileType('w'),
                       dest='out_file',
                       help='output file (or stdout if not set)',
                       default=sys.stdout)

# '''class for counts/coverage files'''
class Counts():
    def __init__(self, files, file_suffix: str):
        self.files = files
        self.file_suffix = file_suffix
        self.n_files = len(self.files)
        self.sample_info_list = []
        self.fill_sample_info_list()
    def fill_sample_info_list(self):
        for i, file in enumerate(self.files):
            sample_name = file.removesuffix(file_suffix)
            countDF = pd.read_csv(file, sep='\t',header=None,names= ["chr", "start","end", "annotation", "counts"], index_col=False)
            total_reads = countDF['counts'].sum()
            scaling_factor = total_reads / 1e6
            countDF['width'] = (countDF['end'] - countDF['start']) / 1000
            countDF['scaled_counts'] = countDF['counts'] / scaling_factor
            countDF['fkpm'] = countDF['scaled_counts'] / countDF['width']
            fpkm_by_chr_annotation = countDF.groupby(['chr', 'annotation'])['fkpm'].mean().reset_index()
            autosome_target_fkpm = fpkm_by_chr_annotation[((fpkm_by_chr_annotation['chr']!="chrY") & (fpkm_by_chr_annotation['chr']!="chrX") & (fpkm_by_chr_annotation['chr']!="chrM") & (fpkm_by_chr_annotation['annotation']=="T"))]['fkpm'].values
            autosome_off_target_fkpm = fpkm_by_chr_annotation[((fpkm_by_chr_annotation['chr']!="chrY") & (fpkm_by_chr_annotation['chr']!="chrX") & (fpkm_by_chr_annotation['chr']!="chrM") & (fpkm_by_chr_annotation['annotation']=="O"))]['fkpm'].values
            mt_fkpm = fpkm_by_chr_annotation[fpkm_by_chr_annotation['chr']=="chrM"]['fkpm'].values
            CN_ratio_target_mean = np.mean(mt_fkpm / autosome_target_fkpm)
            CN_ratio_off_target_mean = np.mean(mt_fkpm / autosome_off_target_fkpm)
            CN_ratio_target_sd = np.std(mt_fkpm / autosome_target_fkpm)
            CN_ratio_off_target_sd = np.std(mt_fkpm / autosome_off_target_fkpm)
            mean_autosomal_target_fkpm = np.mean(autosome_target_fkpm)
            sample_list = [sample_name, mt_fkpm[0], mean_autosomal_target_fkpm, CN_ratio_target_mean,CN_ratio_target_sd, CN_ratio_off_target_mean, CN_ratio_off_target_sd, total_reads]
            self.sample_info_list.append(sample_list)
        self.sampleDF=pd.DataFrame(self.sample_info_list,columns= ["sample_id", "mtFKPM","autosomalFKPM", "CN_ratio_target","CN_ratio_target_sd", "CN_ratio_off_target","CN_ratio_off_target_sd","total_reads"])

def main():
    args = parser.parse_args()
    assert all(args.file_suffix in file for file in args.files), "files have different suffixes"
    samples = Counts(args.files, file_suffix=args.file_suffix)
    samples.sampleDF.to_csv(path_or_buf=args.out_file, sep='\t', encoding='utf-8',index=False, float_format='%3f',header=True)

if __name__ == '__main__':
    main()

