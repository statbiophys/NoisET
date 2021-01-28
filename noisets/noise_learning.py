#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line script to learn the null model parameters for first function of NoisET (1)

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
from optparse import OptionParser
from noisets.noisettes import Data_Process
from noisets.noisettes import Noise_Model

def main():

    """
    Learn noise parameters
    noise_model = 0 : Negative Binomial + Poisson
    noise_model = 1 : Negative Binomial
    noise_model = 2 : Poisson
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

    # Create dataframe

    n, df = cl_S1.import_data()

    # Learn Noise Model parameters

    init_paras_arr = [ np.asarray([ -2.046736,    1.539405,    1.234712,    6.652190,  -9.714225]), \
                    np.asarray([-2.02192528,   0.45220384,   1.06806274, -10.18866972]), \
                     np.asarray([-2.15206189,  -9.46699067])
                 ]

    if options.Negative_bino_poisson:
        noise_model = 0
        init_paras = init_paras_arr[noise_model]

    elif options.Negative_bino:
        noise_model = 1
        init_paras = init_paras_arr[noise_model]

    elif options.Poisson:
        noise_model = 2
        init_paras = init_paras_arr[noise_model]

    null_model = Noise_Model() 
    null_model.learn_null_model(df, noise_model, init_paras)



if __name__ == '__main__': main()

    

