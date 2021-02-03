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


import os
import numpy as np
import pandas as pd

from optparse import OptionParser
from noisets.noisettes import Data_Process
from noisets.noisettes import Expansion_Model

def main():

    """
    Detect responding clones
    noise_model = 0 : Negative Binomial + Poisson
    noise_model = 1 : Negative Binomial
    noise_model = 2 : Poisson


    pval : choose pval limit to detect clones
    smedthreshold: choose logfoldchange limit to detect clones - see Methods of README for more conceptual details
    """

    parser = OptionParser(conflict_handler='resolve')

    #Choice of Noise Model
    parser.add_option("--NBPoisson", action = 'store_true', dest = 'Negative_bino_poisson', default= False, help = 'Choose the negative binomial + poisson noise model' )
    parser.add_option("--NB", action = 'store_true', dest = 'Negative_bino', default= False, help = 'Choose negative binomial noise model' )
    parser.add_option("--Poisson", action = 'store_true', dest = 'Poisson', default= False, help = 'Choose Poisson noise model' )

    #Input and output
    parser.add_option('--path', type = str, metavar='PATH/TO/FILE', dest = 'path',  help='set path name to file ')
    parser.add_option('--f1', type = str, metavar='filename1', dest = 'filename1',  help='name of first sample')
    parser.add_option('--f2', type = str, metavar='filename2', dest = 'filename2', help='name of second sample')


    #Characteristics of data-sets
    parser.add_option('--specify', action = 'store_true', dest = 'give_details', default= False, help = 'give names of data columns')
    parser.add_option('--freq', type = str, metavar='frequencies', dest = 'freq_label', help=' clone fractions - data - label ')
    parser.add_option('--counts', type = str, metavar='counts', dest = 'count_label', help=' clone counts - data - label ')
    parser.add_option('--ntCDR3', type = str, metavar='ntCDR3label', dest = 'ntCDR3_label', help=' CDR3 nucleotides sequences - data - label ')
    parser.add_option('--AACDR3', type = str, metavar='AACDR3label', dest = 'AACDR3_label', help=' CDR3 Amino acids - data - label ')

    # Null model parameters     
    parser.add_option('--nullpara1', type = str, metavar='nullpara1filename', dest = 'nullparastime1',  help= 'Null parameters at first time point')
    parser.add_option('--nullpara2', type = str, metavar='nullpara2filename', dest = 'nullparastime2', help= 'Null parameters at second time point')


    # Thresholds parameters 
    parser.add_option('--pval', type = float, metavar='pvaluethreshol', dest = 'pval',  help='set pval threshold value for detection of responding clones')
    parser.add_option('--smedthresh', type = float, metavar='smedthreshold', dest = 'smed',  help='set  mediane threshold value  for detection of responding clones ')

    parser.add_option('--output', type = str, metavar='detection_filename', dest = 'detect_filename',  help='set path name to file ')


    (options, args) = parser.parse_args()

    ## Code 
    ## Load data
    path = options.path
    filename1 = options.filename1 # first biological replicate
    filename2 = options.filename2 # second biological replicate
    if options.give_details:
        colnames1 = [options.freq_label, options.count_label, options.ntCDR3_label, options.AACDR3_label] #colnames that will change if you work with a different data-set
        colnames2 = [options.freq_label, options.count_label, options.ntCDR3_label, options.AACDR3_label]

    else:
        colnames1 = ['Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3'] #colnames that will change if you work with a different data-set
        colnames2 = ['Clone fraction','Clone count', 'N. Seq. CDR3', 'AA. Seq. CDR3'] 

    # check 
    cl_S1 = Data_Process(path, filename1, filename2, colnames1,  colnames2)
    print("First Filename is : " , cl_S1.filename1)
    print("Second Filename is : ",  cl_S1.filename2)
    print("Name of columns of first file are : ", cl_S1.colnames1)
    print("Name of columns of second file are : ", cl_S1.colnames2)

    df_1 = pd.read_csv(options.nullparastime1, sep ='\t')
    paras_1 = np.array(df_1['value'])
    print('paras time 1: ' + str(paras_1)) 

    df_2 = pd.read_csv(options.nullparastime2, sep ='\t')
    paras_2 = np.array(df_2['value'])
    print('paras time 2: ' + str(paras_2)) 


    # Create dataframe

    n, df = cl_S1.import_data()

    if options.Negative_bino_poisson:
        noise_model = 0

    elif options.Negative_bino:
        noise_model = 1

    elif options.Poisson:
        noise_model = 2

    ## Expansion Model

    expansion = Expansion_Model()
    pval_threshold = options.pval 
    smed_threshold = options.smed

    outpath = options.detect_filename + filename1 + filename2

    expansion.expansion_table(outpath, paras_1, paras_2, df, noise_model, pval_threshold, smed_threshold)


if __name__ == '__main__': main()

    