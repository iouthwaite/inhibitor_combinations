#!/bin/python3.9.4
#Ian Outhwaite, 2021
#082521
#param_scans.py
#runs tests for the JS distance protocol using the PKIS2 dataset

import JS_multipletargets_120621 as JS
import pandas as pd
import numpy as np
import sys
import os
import openpyxl

if __name__ == "__main__":

    ##################### GET DATASETS

    datasets = []

    dataset = pd.read_excel('PKIS2_dataset.xlsx', engine='openpyxl',  nrows=645)
    dataset.name = 'PKIS2_'

    datasets.append(dataset)

    ###################### SETTINGS
    
    prior_type = 2  #1=use a beta distribution, 2=use a poisson distribution
    bdist_alpha = 3 #The alpha value for the beta distribution
    bdist_beta = 1 #The beta value  for the beta distribution
    mu = 700 #The mu value for the poisson distribution. We find 200 and 700 to be good values to use; use both in different tests to compare selectivity profiles
    size = 100000 #The number of points to sample to generate the selectivity prior. We find 100,000 to be a decent number, it facilitates reproucibility without being too large.
    on_target_inhibition_threshold = 90 #the minimum on-target activity required for inhibitors. Try 90, 80, 70 as desired
    number_top_UW_scores_to_maximize = 100 #the top number_top_UW_scores_to_maximize unweighed combinations to select for the subsequent inhibitor concentration optimization step. We suggest using a larger number here than in the single target program since ony a single target set (of multiple targets) is under consideration - depends upon user preferance, but a number in the tens (we use 20) is generally sufficient. If a dataset contains a large number (hundreds) of possible combinations one might consider bumping this up, with the only downside of added processing time
    noise_variance = 2.5 #the variance of the gaussian noise added to measurments
    max_combination_iter = 5 #the maximum number of inhibitor combinations to try (ex: 3 means using 3 inhibitors together in one mix)
    influence = 0.1 #the weight that is used to add an additional penalty to the highest off-target effects (in the [95, 100] bin). We advise using a value of 0.1.
    output_name = '' #Name of the output file
    num_p = 1 #number of proccesses to call during multiprocessing

    target_kinases = ['AKT1', 'MTOR'] #Names of the target kinases - must match exactly with the characters in the dataset (be aware of spaces after the names!)

    ####################### USER DEFINED # OF PROC

    prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size)) #package the prior settings for the program

    if len(sys.argv) > 1:
        print('Number of processes to use: ' + str(sys.argv[1]))
        num_p = int(sys.argv[1])

    ##################### Experiment: Repeat 5 technical replicates for each condition

    for x in range(0,5):

        for dataset in datasets:

            for mu in [200,700]:

                output_name = str(dataset.name) + 'AKT1-MTOR_mu_' + str(mu) +'_' + str(x)

                print('Working on...')
                print(output_name)
                
                new_prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size))
                JS.main_function(dataset, new_prior_settings, noise_variance, max_combination_iter, influence, output_name, num_p, target_kinases, on_target_inhibition_threshold, number_top_UW_scores_to_maximize)
    
