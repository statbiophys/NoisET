# -*- coding: utf-8 -*-
"""Command line script to learn the detection model and spot responding clones to a stimulus from two samples taken at two different time points 
for second function of NoisET (2)

   Copyright (C) 2021 Meriem Bensouda Koraichi

   This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from argparse import ArgumentParser
import os

import numpy as np
import pandas as pd

from noisets.noisettes import Data_Process
from noisets.noisettes import ExpansionModel

def main():

    """
    Detect responding clones
    noise_model = 0 : Negative Binomial + Poisson
    noise_model = 1 : Negative Binomial
    noise_model = 2 : Poisson


    pval : choose pval limit to detect clones
    smedthreshold: choose logfoldchange limit to detect clones - see Methods of README for more conceptual details
    """

    parser = ArgumentParser(conflict_handler='resolve')

    #Choice of Noise Model
    parser.add_argument("--NBPoisson", action = 'store_true', dest = 'Negative_bino_poisson', default= False, help = 'Choose the negative binomial + poisson noise model' )
    parser.add_argument("--NB", action = 'store_true', dest = 'Negative_bino', default= False, help = 'Choose negative binomial noise model' )
    parser.add_argument("--Poisson", action = 'store_true', dest = 'Poisson', default= False, help = 'Choose Poisson noise model' )

    #Input and output
    parser.add_argument('--path', type = str, metavar='PATH/TO/FILE', dest = 'path',  help='set path name to file ', required=True)
    parser.add_argument('--f1', type = str, metavar='filename1', dest = 'filename1',  help='name of first sample', required=True)
    parser.add_argument('--f2', type = str, metavar='filename2', dest = 'filename2', help='name of second sample', required=True)
    parser.add_argument('--output', type = str, metavar='detection_filename', dest = 'detect_filename',  help='set path name to file ', required=True)


    #Characteristics of data-sets
    parser.add_argument('--specify', action = 'store_true', dest = 'give_details', default= False, help = 'give names of data columns')
    parser.add_argument('--freq', type = str, metavar='frequencies', dest = 'freq_label', help=' clone fractions - data - label ')
    parser.add_argument('--counts', type = str, metavar='counts', dest = 'count_label', help=' clone counts - data - label ')
    parser.add_argument('--ntCDR3', type = str, metavar='ntCDR3label', dest = 'ntCDR3_label', help=' CDR3 nucleotides sequences - data - label ')
    parser.add_argument('--AACDR3', type = str, metavar='AACDR3label', dest = 'AACDR3_label', help=' CDR3 Amino acids - data - label ')

    # Null model parameters     
    parser.add_argument('--nullpara1', type = str, metavar='nullpara1filename', dest = 'nullparastime1',  help= 'Null parameters at first time point', required=True)
    parser.add_argument('--nullpara2', type = str, metavar='nullpara2filename', dest = 'nullparastime2', help= 'Null parameters at second time point', required=True)


    # Thresholds parameters 
    parser.add_argument('--pval', type = float, metavar='pvaluethreshol', dest = 'pval',  help='set pval threshold value for detection of responding clones')
    parser.add_argument('--smedthresh', type = float, metavar='smedthreshold', dest = 'smed',  help='set  mediane threshold value  for detection of responding clones ')

    arguments = parser.parse_args()

    bool_sum = arguments.Negative_bino_poisson + arguments.Negative_bino + arguments.Poisson
    if bool_sum == 0:
        raise ValueError('A noise model has not been selected.')
    if bool_sum > 1:
        raise ValueError('Only one noise model must be selected.')

    ## Code 
    ## Load data
    path = arguments.path
    filename1 = arguments.filename1 # first biological replicate
    filename2 = arguments.filename2 # second biological replicate
    if arguments.give_details:
        colnames1 = [arguments.freq_label, arguments.count_label, arguments.ntCDR3_label, arguments.AACDR3_label] #colnames that will change if you work with a different data-set
        colnames2 = [arguments.freq_label, arguments.count_label, arguments.ntCDR3_label, arguments.AACDR3_label]

    else:
        colnames1 = ['Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3'] #colnames that will change if you work with a different data-set
        colnames2 = ['Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3']

    # check 
    cl_S1 = Data_Process(path, filename1, filename2, colnames1,  colnames2)
    print("First Filename is : " , cl_S1.filename1)
    print("Second Filename is : ",  cl_S1.filename2)
    print("Name of columns of first file are : ", cl_S1.colnames1)
    print("Name of columns of second file are : ", cl_S1.colnames2)

    df_1 = pd.read_csv(arguments.nullparastime1, sep ='\t')
    paras_1 = np.array(df_1['value'])
    print('paras time 1: ' + str(paras_1))

    df_2 = pd.read_csv(arguments.nullparastime2, sep ='\t')
    paras_2 = np.array(df_2['value'])
    print('paras time 2: ' + str(paras_2))


    # Create dataframe

    n, df = cl_S1.import_data()

    if arguments.Negative_bino_poisson:
        noise_model = 0

    elif arguments.Negative_bino:
        noise_model = 1

    elif arguments.Poisson:
        noise_model = 2

    ## Expansion Model

    expansion = ExpansionModel()
    pval_threshold = arguments.pval
    smed_threshold = arguments.smed
    print(df.columns)
    outdf = expansion.expansion_table(df, noise_model, paras_1, paras_2,
                                      count_1_col='Clone_count_1',
                                      count_2_col='Clone_count_2')
    if pval_threshold is not None:
        pval_expand_subset = outdf['pval_expand'] <= pval_threshold
        pval_contract_subset = outdf['pval_contract'] <= pval_threshold
        outdf = outdf[pval_expand_subset | pval_contract_subset]
    if smed_threshold is not None:
        smed_expand_subset = outdf['s_med_expand'] >= smed_threshold
        smed_contract_subset = outdf['s_med_contract'] <= -smed_threshold
        outdf = outdf[smed_expand_subset | smed_contract_subset]

    outdf.to_csv(arguments.detect_filename)

if __name__ == '__main__': main()
