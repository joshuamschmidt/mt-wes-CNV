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
# '''class for counts/coverage files'''
class Counts():
    def __init__(self, files, count_column: int, file_suffix: str):
        self.files = files
        self.file_suffix = file_suffix
        self.count_column = count_column
        self.n_files = len(self.files)
        self.sample_names = []
        self.n_rows = None
        self.first_key=None
        self.get_n_rows()
        self.count_array = np.empty([self.n_rows,self.n_files], dtype=np.int64)
        self.fill_count_array()
        #self.sample_names = [file.removesuffix(file_suffix) for file in self.files]

    def get_n_rows(self):
        file=self.files[0]
        first_fileDF=pd.read_csv(file, sep='\t',header=None)
        lines = len(first_fileDF.index)
        self.n_rows = lines

    def fill_count_array(self):
        count_index = self.count_column - 1
        for i, file in enumerate(self.files):
            countDF = pd.read_csv(file, sep='\t',header=None)
            count_array=countDF.iloc[:,count_index].values
            #check dtype of count file. change type of self.count_array
            if i==0:
                self.count_array = self.count_array.astype(count_array.dtype)
                self.first_key = str(countDF.iloc[0,0]) + "_" + str(countDF.iloc[0,1])
            assert self.n_rows == np.size(count_array), f"File {file} has an unexpected number of rows. Expected {self.n_rows}. Observed {np.size(count_array)}"
            self.count_array[:,i] = count_array
            self.sample_names.append(file.removesuffix(self.file_suffix))

    def convert_to_fkpm(self,bedObject):
        assert self.first_key == bedObject.first_key, f"count files are not sorted in the same order as the bed file. File key: {self.first_key}, Bed key: {bedObject.first_key}"
        assert self.n_rows == bedObject.n_rows, f"count files and bed file differ in nummber of features. File: {self.n_rows}, Bed: {bedObject.n_rows}"
        perM_scaling_factors = self.count_array.sum(axis=0) / 1e6
        fpm_scaled_count_array = self.count_array / perM_scaling_factors
        fpkm_array = np.array(fpm_scaled_count_array / bedObject.kb_lengths[:,None])
        self.count_array = fpkm_array

    def array_to_df(self):
        self.count_array = pd.DataFrame(self.count_array)
        self.count_array.columns = self.sample_names

    def append_bed(self, bed):
       assert isinstance(self.count_array, pd.DataFrame), "you need to convert to DF before this! (array_to_df)"
       self.count_array = pd.concat([bed.bed_data,self.count_array],axis=1)


def main():
    args = parser.parse_args()
    assert all(args.file_suffix in file for file in args.files), "files have different suffixes"
    samples = coverageSummaries(args.files, file_suffix=args.file_suffix)
    samples.summaries.to_csv(path_or_buf=args.out_file, sep='\t', encoding='utf-8',index=False, float_format='%3f',header=True)

if __name__ == '__main__':
    main()

