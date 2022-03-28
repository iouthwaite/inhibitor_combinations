#!/bin/python3.9.4
#Ian Outhwaite, 2021

import JS_singletarget_120621 as JS
import pandas as pd
import numpy as np
import sys
import openpyxl

if __name__ == "__main__":

    dataset=pd.read_excel("./PKIS2_dataset.xlsx", engine='openpyxl', nrows=645)

    print(dataset)

    prior_type = 2  #1=use a beta distribution, 2=use a poisson distribution
    bdist_alpha = 3 #The alpha value for the beta distribution
    bdist_beta = 1 #The beta value  for the beta distribution
    mu = 700 #The mu value for the poisson distribution. We find 200 and 700 to be good values to use; use both in different tests to compare selectivity profiles
    size = 100000 #The number of points to sample to generate the selectivity prior. We find 100,000 to be a decent number, it facilitates reproucibility without being too large.
    on_target_inhibition_threshold = 90 #the minimum on-target activity required for inhibitors. Try 90, 80, 70 as desired
    number_top_UW_scores_to_maximize = 3 #the top number_top_UW_scores_to_maximize unweighed combinations to select for the subsequent inhibitor concentration optimization step. We suggest using 3 or 5 - only in very rare instances would a combination that doesn't perform well in unweighed scoring perform exceptionally better in weighed scoring, and reducing the number of combinations to optimize speeds up processing time at this step linearly.
    noise_variance = 2.5 #the variance of the gaussian noise added to measurments
    max_combination_iter = 3 #the maximum number of inhibitor combinations to try (ex: 3 means using 3 inhibitors together in one mix)
    influence = 0.1 #the weight that is used to add an additional penalty to the highest off-target effects (in the [95, 100] bin). We advise using a value of 0.1.
    output_name = '' #Name of the output file
    num_p = 60 #number of proccesses to call during multiprocessing

    ####################### RUN

    prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size)) #package the prior settings for the program

    if len(sys.argv) > 1:
        print('Number of processes to use: ' + str(sys.argv[1]))
        num_p = int(sys.argv[1])

    for x in range(0,5):
    
        excel_prefix = '120721_PKIS2_full_90thresh_Inf0.1_optsteps_R1s5_R2s1p1_noise2p5_x'+str(x)

        for new_mu in [200,700]:
            output_name = excel_prefix + '_poisson_' + str(new_mu)
            new_prior_type = 2
            new_prior_settings = (new_prior_type, (bdist_alpha, bdist_beta, new_mu, size))

            JS.main_function(dataset, new_prior_settings, noise_variance, max_combination_iter, influence, output_name, num_p, on_target_inhibition_threshold, number_top_UW_scores_to_maximize)
    
