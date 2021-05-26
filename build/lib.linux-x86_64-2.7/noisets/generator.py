# -*- coding: utf-8 -*-
"""Command line to generate biological replicates counts from same patients same day

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
from noisets.noisettes import Generator


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


    # Null model parameters     
    parser.add_option('--nullpara', type = str, metavar='nullparafilename', dest = 'nullparas',  help= 'path to null parameters')


    # Samples parameters 
    parser.add_option('--NreadsI', type = float, metavar='NrI', dest = 'NreadsI',  help='Number of reads in first replicate')
    parser.add_option('--NreadsII', type = float, metavar='NrII', dest = 'NreadsII',  help='Number of reads in second replicate  ')
    parser.add_option('--Nclones', type = float, metavar='Nclones', dest = 'Nclones',  help='Number  of clones present in replicate 1 OR in replicate 2 ')

    parser.add_option('--output', type = str, metavar='synthetic_data_filename', dest = 'synth_filename', help = 'filename of generated data')


    (options, args) = parser.parse_args()

    ## Code 
   

    # check 
    cl_rep_gen_gDNA = Generator()
    
    #paras = np.load(options.nullparas)
    df = pd.read_csv(options.nullparas, sep ='\t')
    paras = np.array(df['value'])
    print('paras: ' + str(paras)) 

    NreadsI = int(options.NreadsI)
    NreadsII = int(options.NreadsII)
    Nsamp = int(options.Nclones)

    if options.Negative_bino_poisson:
        noise_model = 0

    elif options.Negative_bino:
        noise_model = 1

    elif options.Poisson:
        noise_model = 2



    # Create dataframe

    f_samples_gDNA, pair_samples_gDNA = cl_rep_gen_gDNA.gen_synthetic_data_Null(paras, noise_model, NreadsI, NreadsII, Nsamp)

    pair_samples_gDNA.to_csv('synthetic_'+  options.synth_filename + '.csv', sep = '\t')


if __name__ == '__main__': main()