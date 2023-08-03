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
from argparse import ArgumentParser
import os

import numpy as np

from noisets.noisettes import Data_Process
from noisets.noisettes import NoiseModel

def main():

    """
    Learn noise parameters
    noise_model = 0 : Negative Binomial + Poisson
    noise_model = 1 : Negative Binomial
    noise_model = 2 : Poisson
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

    #Characteristics of data-sets

    parser.add_argument('--specify', action = 'store_true', dest = 'give_details', default= False, help = 'give names of data columns')
    parser.add_argument('--freq', type = str, metavar='frequencies', dest = 'freq_label', help=' clone fractions - data - label ')
    parser.add_argument('--counts', type = str, metavar='counts', dest = 'count_label', help=' clone counts - data - label ')
    parser.add_argument('--ntCDR3', type = str, metavar='ntCDR3label', dest = 'ntCDR3_label', help=' CDR3 nucleotides sequences - data - label ')
    parser.add_argument('--AACDR3', type = str, metavar='AACDR3label', dest = 'AACDR3_label', help=' CDR3 Amino acids - data - label ')


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

    # Create dataframe

    n, df = cl_S1.import_data()

    # Learn Noise Model parameters

    init_paras_arr = [np.asarray([ -2.046736,    1.539405,    1.234712,    6.652190,  -9.714225]),
                      np.asarray([-2.02192528,   0.45220384,   1.06806274, -10.18866972]),
                      np.asarray([-2.15206189,  -9.46699067])
                     ]

    if arguments.Negative_bino_poisson:
        noise_model = 0
        init_paras = init_paras_arr[noise_model]

    elif arguments.Negative_bino:
        noise_model = 1
        init_paras = init_paras_arr[noise_model]

    elif arguments.Poisson:
        noise_model = 2
        init_paras = init_paras_arr[noise_model]

    null_model = NoiseModel()
    null_model.learn_null_model(df, noise_model, init_paras)

if __name__ == '__main__': main()
