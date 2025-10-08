# (C) Ian Outhwaite & Sukrit Singh, 2023
# Death by a Thousand Cuts – Combining Kinase Inhibitors for Selective Target Inhibition and Rational Polypharmacology
# Ian R. Outhwaite, Sukrit Singh, Benedict-Tilman Berger, Stefan Knapp, John D. Chodera, Markus A. Seeliger
# Run JSD scoring of inhibitor combinations
# python version 3.9.13

import MMS_methods_100825 as JS
import pandas as pd
import sys
import openpyxl
import warnings
import itertools

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

if __name__ == "__main__":

    ###########################################################################################################
    ######################################### Primary User Settings ###########################################

    ################ single or multiple target analysis?

    single=False #set to False for multiple target analysis

    ################ user defined dataset and number of inhibitors to use
    #                can pre-sort the datasets and set the number of inhibitors to use only a subset of the data
    #                otherwise, leave num_inhibitors to be the full number in the dataset
    
    dataset_path = './Klaeger_2017_MET_converted_subset.csv'
                   
    num_inhibitors = 3
    
    #645 Inhibitors considered for inclusion in PKIS2
    #72 Davis 2011
    #20 Fabian 2005
    #38 Karaman 2008
    #243 Klaeger 2017

    ################ user defined targets for multiple target analysis. For single target analysis this can be left blank.
    #                Ensure that case matches input dataset (and trailing whitespace in names in some cases)
    #                If conducting multiple-target analysis, each set of kinases must be formatted as a list within a list
    #                Ex: [ [kinase1, kinase2] , [kinase3, kinase4] ]
    #                Set to 'ALL' to query all kinases in a dataset as individual targets
    
    data = pd.read_csv(dataset_path)

    target_kinases = ['MET(WT)','MET(ex14del)','MET(Y1234E/Y1235E)','MET(V1092I)','MET(H1094R)','MET(H1094Y)','MET(L1195V)','MET(F1200I)','MET(D1228H)','MET(D1228N)','MET(D1228V)','MET(Y1230A)','MET(Y1230C)','MET(Y1230D)','MET(Y1230H)','MET(Y1235D)','MET(M1250T)']

    ################ user defined kinases not to include in off-target distributions (ex: disease relevant mutants in a dataset that wouldn't be relevent when considering off-target effects)

    non_off_target_kinases = ['MET(WT)','MET(ex14del)','MET(Y1234E/Y1235E)','MET(V1092I)','MET(H1094R)','MET(H1094Y)','MET(L1195V)','MET(F1200I)','MET(D1228H)','MET(D1228N)','MET(D1228V)','MET(Y1230A)','MET(Y1230C)','MET(Y1230D)','MET(Y1230H)','MET(Y1235D)','MET(M1250T)']
        
    ################ critical off-targets that the user does NOT want to inhibit, as well as a weight to assign to each. The weight is roughly equal to one equivalent off-target (ex: w=10 would be simmilar to inhibiting 10 off-targets by the assigned amount)
    
    critical_off_targets = []
    critical_off_target_weights = []
    
    ################ activity threshold for considering inhibitors for on-target inhibition. Input from 1->99. 90, 80, or 70 is suggested, depends on data and reference concentration frame.
    
    min_potency = 0
    
    ################ activity threshold for combinations to reach (ex: must reach 90% on-target activity) at indicated concentrations. Input from 1->99. 90, 80, or 70 is suggested, depends on data and reference concentration frame.

    on_target_activity_threshold = 90

    ################ maximum number of inhibitors to try in a combination (ex: up to 4 inhibitors in a single combination)

    max_number_inhibs_per_combo = 2

    ################ the number of combinations to consider for concentration optimization and then scoring
    #                -for single-target scoring, this can be left relatively low to speed up calculations (ex: 10 is probably fine)
    #                -for multiple-target scoring, this number should be higher depending upon the expected number of combinations
    #                 you are testing to cover possible combinations.
    #                -Set to 'ALL' to cover all possibilities

    num_UW_to_test = 'ALL'

    ################ number of processes to call. Can also be set by the first command-line argument which will overwrite the variable (ex: python JS_run.py 25)

    num_p = 1

    ################ number of technical replicates to conduct. This is helpful for judging the significance of observed improvements in JSD score.

    n_replicates = 1

    ################ prefix for output files

    outfile_prefix = '041025_multi_target_combination_screen_inhibs_50uM_upper_limit_R1_1p001_R2_1p00001_nonoise_targ90_inf0p5'

    ################ lower limit for compound dilutions given input reference frame. Ex: given a reference frame of 1uM, 0.001 would mean using compounds at a lower limit of 1 nM

    lower_limit = 0.0000000000001

    ################ upper limit for concentration of compounds given input reference frame. Ex: given a reference frame of 1uM, 100 would mean using compounds at an upper limit of 100 uM

    upper_limit = 50

    ################ maximum number of combinations to report per target kinase or set of target kinases
    #                set this to an int, or 'ALL' to report all possible combinations/inhibitors

    max_num_to_report = 10

    ###########################################################################################################
    ########################################## Additional Settings ############################################

    ################ use a penalty prior with a base shape of a poisson or beta distribution? Poisson is suggested
    #                set to 2 for poisson, 1 for beta
    
    prior_type = 2

    ################ define the parameters for the penalty prior distribution
    #                for a beta distribution, set alpha / beta.
    #                for poisson, set the mu
    #                200, 700, 1200 corresponds to using all three of the high, medium, and broad penalty priors used in the associated work

    bdist_alpha = 3
    bdist_beta = 1

    mu_range = [200,700]

    ################ set the high off-target penalty (0.1->0.3 is suggested)

    influence = 0.5

    ################ set the noise added to off-target measurments when building off-target distributions
    #                2.5 is suggested, set to 0 if you want to improve reproducibility of results
    #                this value corresponds to the variance of the normal distributions that are used for sampling each off-target calculation

    noise_variance = 0

    ################ set the step sizes for calculating optimal concentrations of inhibitors
    #                The R1 step refers to how much inhibitor pairs are increased or decreased by (a factor of R1) during initial branching
    #                The R2 step refers to the factor that concentrations are decreased by at each step in order to obtain minimum on-target activity
    #                Increasing the step sizes will increase the speed of the optimization steps, but will reduce the accuracy of obtained concentrations
    #                We reccomend using larger step sizes for initial screens (ex: R1=5, R2=1.2), and then decreasing the step sizes for smaller analysis using drug-target pairs

    R1_step = 1.001
    R2_step = 1.00001

    ################ set the time limit per single target-combination in the optimization concentrations protocol (in seconds)

    time_lim = 1000

    ################ set the number of points samples from the prior distribution used to build the penalty prior.

    size = 100000

    ###########################################################################################################
    ################################################### Run ###################################################

    print('\nInitiating JSD Analysis\nLoading dataset..')
    
    dataset=pd.read_csv(dataset_path)#, engine='openpyxl', nrows=num_inhibitors)

    print('dataset loaded\n')
    
    if len(sys.argv) > 1:
        print('Number of procs initiated: ' + str(sys.argv[1])+'\n')
        num_p = int(sys.argv[1]) #overwrite the number preset by the user

    if len(non_off_target_kinases) != 0:
        print('Excluding the following from off-target distributions:')
        for k in non_off_target_kinases:
            print(k)
        print()
        
    if len(critical_off_targets) != 0:
        print('Critical Off-Targets MMS will attempt to minimize inhibiton for:')
        for k in critical_off_targets:
            print(k)
        print()

    if single:
        print('Conducting single-target analysis')
        if target_kinases != []:
            if target_kinases == 'ALL':
                target_kinases = list(dataset.columns.values)[1:]
            print('Target kinases: ')
            for tk in target_kinases:
                print(tk)
            for x in range(n_replicates):
                print('------Conducting technical replicate round '+str(x+1)+' of ' + str(n_replicates)+'------')
                for mu in mu_range:
                    outname = outfile_prefix + '_mu_'+str(mu) + '_rep_' + str(x)
                    prior_settings = (2, (bdist_alpha, bdist_beta, mu, size))
                    JS.main_function(
                        dataset,
                        prior_settings,
                        noise_variance,
                        max_number_inhibs_per_combo,
                        influence,
                        outname,
                        num_p,
                        on_target_activity_threshold,
                        num_UW_to_test,
                        target_kinases,
                        non_off_target_kinases,
                        R1_step,
                        R2_step,
                        time_lim,
                        lower_limit,
                        upper_limit,
                        single,
                        max_num_to_report,
                        critical_off_targets,
                        critical_off_target_weights,
                        min_potency
                    )

    if not single:
        print('Conducting multiple-target analysis')
        print('Target kinases:')
        for k in target_kinases:
            print(k)
        
        for x in range(n_replicates):
            print('------Conducting technical replicate round '+str(x+1)+' of ' + str(n_replicates)+'------')
            for mu in mu_range:
                outname = outfile_prefix + '_mu_'+str(mu) + '_rep_' + str(x)
                prior_settings = (2, (bdist_alpha, bdist_beta, mu, size))
                JS.main_function(
                    dataset,
                    prior_settings,
                    noise_variance,
                    max_number_inhibs_per_combo,
                    influence,
                    outname,
                    num_p,
                    on_target_activity_threshold,
                    num_UW_to_test,
                    target_kinases,
                    non_off_target_kinases,
                    R1_step,
                    R2_step,
                    time_lim,
                    lower_limit,
                    upper_limit,
                    single,
                    max_num_to_report,
                    critical_off_targets,
                    critical_off_target_weights,
                    min_potency
                )
