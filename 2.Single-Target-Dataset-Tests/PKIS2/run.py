#!/bin/python3.9.4
#Ian Outhwaite, 2021

import JS_singletarget_020222_no_dilutions as JS
import pandas as pd
import numpy as np
import sys
import openpyxl

if __name__ == "__main__":

    dataset=pd.read_excel("./PKIS2_dataset.xlsx", engine='openpyxl', nrows=645)

    print(dataset)

    prior_type = 1
    bdist_alpha = 3
    bdist_beta = 1
    mu = 400
    size = 100000
    prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size))
    noise_variance = 2.5
    max_combination_iter = 5
    influence = 0.1
    output_name = ''
    num_p = 60

    if len(sys.argv) > 1:
        print('Number of processes to use: ' + str(sys.argv[1]))
        num_p = int(sys.argv[1])


    ##################### Variable Scanning Tests: Repeat each test 5x 
    
    ##################### TEST 1: Vary Prior

    for x in range(0,5):
    
        excel_prefix = '020222_PKIS2_full_90thresh_nodilutions_max5_'+str(x)

        for new_mu in [700,200]:
            output_name = excel_prefix + '_poisson_' + str(new_mu)
            new_prior_type = 2
            new_prior_settings = (new_prior_type, (bdist_alpha, bdist_beta, new_mu, size))
            JS.main_function(dataset,new_prior_settings,noise_variance,max_combination_iter,influence,output_name,num_p, 90, 1)
    
