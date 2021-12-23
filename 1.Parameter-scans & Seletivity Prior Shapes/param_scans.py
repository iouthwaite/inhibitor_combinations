#!/bin/python3.9.4
#Ian Outhwaite, 2021
#082521
#param_scans.py
#runs tests for the JS distance protocol using the PKIS2 dataset

import JS_distance_drugpotency_111821 as JS
import pandas as pd
import numpy as np
import sys
import openpyxl

if __name__ == "__main__":

    dataset=pd.read_excel("./PKIS2_dataset_no_promks.xlsx", engine='openpyxl', nrows=645)

    print(dataset)

    prior_type = 1
    bdist_alpha = 3
    bdist_beta = 1
    mu = 400
    size = 10000
    prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size))
    noise_variance = 5
    max_combination_iter = 3
    influence = 0.0
    output_name = ''
    num_p = 5

    if len(sys.argv) > 1:
        print('Number of processes to use: ' + str(sys.argv[1]))
        num_p = int(sys.argv[1])


    ##################### Variable Scanning Tests: Repeat each test 5x 
    
    ##################### TEST 1: Vary Prior Settings 

    for x in range(0,5):
    
        excel_prefix = '111921_PKIS2_nopromks_90thresh_NoInfluence_x'+str(x)+'_vary_prior_'
        
        for new_bdist_alpha in np.arange(2,8,0.5):
            output_name = excel_prefix + 'beta_' + str(new_bdist_alpha)+'_1'
            prior_type = 1
            new_prior_settings = (prior_type, (new_bdist_alpha, bdist_beta, mu, size))
            JS.main_function(dataset,new_prior_settings,noise_variance,max_combination_iter,influence,output_name,num_p)

        for new_mu in np.arange(100,1300,100):
            output_name = excel_prefix + 'poisson_' + str(new_mu)
            new_prior_type = 2
            new_prior_settings = (new_prior_type, (bdist_alpha, bdist_beta, new_mu, size))
            JS.main_function(dataset,new_prior_settings,noise_variance,max_combination_iter,influence,output_name,num_p)

    

    ##################### TEST 2: Vary Gaussian Noise Added to Measurments

    for x in range(0,5):

        excel_prefix = '111921_PKIS2_nopromks_90thresh_NoInfluence_x'+str(x)+'_vary_noise_'

        for new_noise in np.arange(0,6,1):
            output_name = excel_prefix+str(new_noise)
            JS.main_function(dataset,prior_settings,new_noise,max_combination_iter,influence,output_name,num_p)
    
    ##################### TEST 3: Vary the weight of the on-target measurment versus the on-target prior 

    for x in range(0,5):

        excel_prefix = '111921_PKIS2_nopromks_90thresh_NoInfluence_x'+str(x)+'_vary_influence_'

        for new_influence in np.arange(0, 0.6, 0.1):
            output_name = excel_prefix+str(new_influence)
            JS.main_function(dataset,prior_settings,noise_variance,max_combination_iter,new_influence,output_name,num_p)
    
    ##################### END VARIABLE SCANNING TESTS


    
