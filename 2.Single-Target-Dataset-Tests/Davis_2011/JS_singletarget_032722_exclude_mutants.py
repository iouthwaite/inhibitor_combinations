#!/bin/python3.9.4
#Ian Outhwaite & Sukrit Singh, 2021

#################### import

import pandas as pd
import numpy as np
import xlsxwriter
import openpyxl
import itertools
from itertools import repeat
import random
import statistics
from math import log2
from scipy.spatial import distance
from scipy.stats import beta
from scipy.stats import poisson
import multiprocessing
from multiprocessing import Process
from progressbar import ProgressBar, SimpleProgress
import time

#################### functions

# Get kinase names, inhibitor names, and a 2D array of data from the input dataset
def import_from_dataset(dataset):
    kinase_names = list(dataset.columns.values)[7:]
    compound_names = dataset['Compound'].tolist()
    array_values = dataset[dataset.columns[7:]].to_numpy()
    return kinase_names, compound_names, array_values

# Generates a left-skewed distribution using a Poisson prior (as a basis for the shape) with the mean (max) at 100, only selected from [0, 100] 
def get_Poisson_Prior(mu, _size):
    r = poisson.rvs(mu, size = _size*2)
    mean = np.mean(r)
    if mean > 100:
        r = [x - (mean-100.0) for x in r]
    else:
        r = [x + (100.0 - mean) for x in r]
    prior = []
    for x in r:
        if x <= 100:
            prior.append(x)
    return prior

# Generates a beta distribution scaled to [0, 100]
def get_Beta_Prior(a, b, _size):
    prior = beta.rvs(a, b, scale=1, size=_size)
    prior = [x * 100.0 for x in prior]
    return prior

# Depending upon the prior settings, generate a set of values and fill probability bins to generate the selectivity prior
# Return the bins and the number of entries used to fill them (to help normalize values later on) 
def generate_prior(settings):
    prior = []
    size = settings[1][3]
    if settings[0] == 1:
        print('Generating a Beta Distribution Prior')
        print('Alpha='+str(settings[1][0])+', Beta='+str(settings[1][1])+', size='+str(size)) 
        prior = get_Beta_Prior(settings[1][0], settings[1][1], size)
    elif settings[0] == 2:
        print('Generating a Poisson Prior')
        print('Mu='+str(settings[1][2])+', size='+str(size))
        prior = get_Poisson_Prior(settings[1][2], size)
    bins = [0 for i in range(0,20)]
    for x in range(0, len(prior)):
        temp = int(prior[x]/5)
        if temp < 20 and temp >= 0:
            bins[temp] += 1.0
        elif temp >= 20:
            bins[len(bins)-1] += 1.0
        else:
            bins[0] += 1.0
    return (bins,size)

# Returns normally distributed random values
def generate_noise(noise_variance):
    noise_v = np.random.normal(0, noise_variance, 100)
    return noise_v

# Returns the Jensen-Shannon distance metric
def get_score(on, off):
    JS_score = distance.jensenshannon(on, off, 2.0)
    return JS_score

# Given a % inhibition (activity) value in a dataset, returns the Ki of the inhibitor.
# Ki will be in uM assuming the inhibitor is used at a concentration of 1uM (as in PKIS2)
def get_Ki(percent_inhibition):
    if percent_inhibition == 100:
        percent_inhibition = 99.5
    return ((1.0/(percent_inhibition/100.0))-1)

# Returns the concentration of a drug needed to reach a certain threshold (in % activity). Returns the concentration in uM, requires the Ki of the drug.
# ( [drug] + Ki ) * activity = [drug]
# ( [drug] * activity + Ki * activity ) = [drug]
# if we use concentrations relative to 1uM, then
# Ki * activity = [drug] * (1-activity)
# (Ki * activity) / (1-activity) = [drug]
# activity must be in 0->1 range, not 0->100
def get_concentration(Ki, threshold):
    t = threshold / 100.0
    c = (Ki * t) / (1 - t)
    return c

#returns the indicies of targets that should be excluded from the off-target distribution (ex: disease-mutants that wouldn't normally be present)
def get_indicies_to_exclude(names_to_exclude, dataset):
    names = list(dataset.columns.values)[7:]
    inds = []
    for name in names_to_exclude:
        inds.append(names.index(name))
    return inds
    
# Given a set of inhibition values and inhibitor concentrations against a single target, returns the total inhibition at a single target
# for each drug, let u = [drug]/Ki_drug
# cumulative % inhibition (activity) = (sum(u for all u)) * 100% / (sum(u for all u) + 1)
def get_t_inhib_weights(inhibitor_values, concentrations):
    if len(inhibitor_values) != len(concentrations):
        print('ERROR')
    sum_score = 0
    if len(inhibitor_values) > 1:
        for x in range(0, len(inhibitor_values)):
            if inhibitor_values[x] != 0:
                sum_score += (float(concentrations[x]) / get_Ki(inhibitor_values[x]))
        t_i = sum_score/(1+sum_score)
        return (t_i*100)
    else:
        if inhibitor_values[0] != 0:
             sum_score += (float(concentrations[0]) / get_Ki(inhibitor_values[0]))
             t_i = sum_score/(1+sum_score)
             return (t_i*100)
        else:
            return 0
             
# Returns the inhibitor indicies that target each kinase
def get_i_target_k(data,otit):
    res = []
    for k in range(0,len(data[0])):
        i_data = data[:,k]
        target_i = [i for i,val in enumerate(i_data) if val >= otit]
        res.append(target_i)
    return res

# Given a slice of the data matrix corresponding to the correct inhibitor data, recalculate the activity values and pass them back
# conc is a value in uM calculated by the get_concentration() method
# This method is not used in the current iteration of the program, kept here for legacy reasons & to facilitate user hacks if so desired 
def calculate_activity(conc, data_slice):
    if len(data_slice) >1:
        res = []
        for col in data_slice:
            if col > 0:
                Ki = get_Ki(col)
                activity = ((conc/Ki)/((conc/Ki)+1)*100)
                res.append(activity)
            else:
                res.append(0)
        return res
    else:
        if data_slice[0] > 0:
            Ki = get_Ki(data_slice[0])
            activity = ((conc/Ki)/((conc/Ki)+1)*100)
            return activity
        else:
            return 0

# Gets all combinations, and for each combination also pulls a copy of the on-target values and the off-target values
# Also returns the concentration of each inhibitor needed to reach (activity_thresh) target inhibition. Ex: 90%.
def get_combinations(i_target_k, cmbn, data, otit, ind_to_exclude):
    res = []
    for k in range(0, len(i_target_k)):
        inhibs = i_target_k[k]
        for c in itertools.combinations(inhibs,cmbn):
            c = list(c)
            on = []
            off = []
            off_ind = []
            concentrations = []
            for i in c:
                target = data[i,k]
                Ki = get_Ki(target)
                conc = get_concentration(Ki, otit)
                concentrations.append(conc)
                on.append(target)
                off_target_selector = [x for x in range(data.shape[1]) if x != k and x not in ind_to_exclude]
                off.append(data[i, off_target_selector])
                off_ind.append(off_target_selector)
            on = np.array(on)
            off = np.array(off)
            res.append((k,c,on,off,concentrations,off_ind))
    return res

# Generates an on-target distribution
# penalty represents the penalty (from 0->1) added against the highest off-target effects. This weight corresponds to percent of the initial (not final) prior size.
# Ex: 0.2 would add an aditional penalty weight of 20% of the initial distribution (or 16.6% of the final distribution: 20/120) added to the highest off-target effects in the (95,100] bin.
# User defined, in practice a value around ~0.1 is helpful if you want to penalize the highest off-target effects.
# This useage is also helpful to facilitate proper scoring of the results; we find that a value of ~0.1 tends to be helpful in this regard.
# Values that are too high put too much emphasis on the highest off-target effects
def make_on_t_dist(values, concentrations, prior,  noise_variance, influence, extended_output):
    bins = prior[0].copy()
    prior_size = prior[1]
    num_to_add = int(prior_size * influence)
    bins[len(bins)-1] += num_to_add
    total_num_points = float(sum(bins))
    bins[:] =  [(x/total_num_points) for x in bins]
    if extended_output == False:
        return bins
    else:
        ot_inhib = get_t_inhib_weights(values, concentrations)
        return (bins, ot_inhib)

# Generates the off-target distribution, bin size=5
# Using 100 randomly sampled points from a normal distribution, replicate these values at each off-target measurment.
# The off-target distribution is the sum of these noise distributions, normalized to 1
def make_off_t_dist(values, weights, noise_variance, extended_output):
    off_t_dist = []
    off_t_inhib_vals = []
    noise_dist = np.random.normal(0, noise_variance, 100) #sample 100 points of gaussian noise to add at each measurment
    numzeros = 0
    for x in range(0, len(values[0])):
        p_vals = values[:, x]
        if np.sum(p_vals) != 0: # for every off-target measurment that isn't 0...
            offt_inhib = get_t_inhib_weights(p_vals, weights)
            off_t_inhib_vals.append(offt_inhib)
            noise_vals = [offt_inhib + noise_dist[x] for x in range(0,100)] #move the noise to be centered at the measurment
            off_t_dist += noise_vals #add the noise to the off-target distributin
        else:
            off_t_inhib_vals.append(0) #otherwise, take note of how many were 0
            numzeros += 1
    nvals = len(off_t_dist)
    bins = [0 for i in range(0,20)] #bin size=5 for the final off-target distribution
    for x in range(0, nvals): #for every point in the set of sample noise (100 per measrument), add it to the distribution
        temp = int(off_t_dist[x]/5)
        if temp < 20 and temp >= 0:
            bins[temp] += 1.0
        elif temp >= 20:
            bins[len(bins)-1] += 1.0
        else:
            bins[0] += 1.0
    for y in range(0,numzeros): #this will ensure that we normalize by the same number - add 100 "points of noise" for each off-target value of 0
        bins[0] += 100.0 #if we don't add 100 points of noise above per measurment, change this value as well
    normfactor = np.sum(bins)
    bins[:] = [(x/normfactor) for x in bins] #normalize the off-target distribution to 1
    if extended_output == False:
        return bins
    else:
        avg_offt = np.mean(off_t_inhib_vals)
        num_offt = np.count_nonzero(off_t_inhib_vals)
        return (bins, avg_offt, num_offt, off_t_inhib_vals)


# JS unweighed inhibitor scoring
# For each combination, pass in equal divisions of each inhibitor (ex: 50%-50% at combination #=2), then score the on-target versus off-target distributions
def JS_UW(inpt):

    combination = inpt[0]
    prior = inpt[1]
    noise_variance = inpt[2]
    influence = inpt[3]
    n_iter = inpt[4]
    t_iter = inpt[5]
    prior_t = inpt[6]

    k_ind = combination[0]
    inhibitors = combination[1]
    on_target_vals = combination[2]
    off_target_vals = combination[3]
    concentrations = combination[4]
    off_target_ind = combination[5]

    if n_iter % 1000 == 0 and n_iter != 0:
        #if multiprocessing is set to return threads in order, the next line can be used to estimate the time remaining in the program
        time_r = (((time.time() - prior_t) * t_iter) / n_iter) - (time.time() - prior_t) 
        print('Working on iter:', n_iter, '/',  t_iter, flush=True)

    # Reduce the concentration of each inhibitor by the number of inhibitors used.
    # Ex: for 2 inhibitors, use each at 50% to maintain a total of (ex: 90%) on-target inhibition
    new_conc = [(concentrations[x]/len(inhibitors)) for x in range(0, len(concentrations))]

    on_target = make_on_t_dist(on_target_vals, new_conc, prior, noise_variance, influence, False)
    off_target = make_off_t_dist(off_target_vals, new_conc, noise_variance, False)
    score = get_score(on_target, off_target)

    return (k_ind, inhibitors, score, on_target_vals, off_target_vals, new_conc, off_target_ind)

# Generate all possible variations of +/- inhibitor weights given a certain step size
# Will make branches no more than 100X greater or 1000X less than the initial inhibitor concentrations
def branch(nodes,step,init_concentrations):
    branches = []
    for n in nodes:
        for x in range(0, len(n)):
            for y in range(0, len(n)): #paired changes to "get off the ground" as it were - increase the [] of one, decrease the [] of another
                if x != y:
                    n1 = n.copy()
                    n1[x] = n1[x] * step
                    n1[y] = n1[y] * (1.0/step)
                    #make sure the concentrations of the drugs are realistic - stay within 100X (greater) 1000X (less) of initial concentrations
                    if (n1[x] < 100 * init_concentrations[x]) and (n1[x] * 1000 > init_concentrations[x]):
                        if (n1[y] < 100 * init_concentrations[y]) and (n1[y] * 1000 > init_concentrations[y]):
                            branches.append(n1)
    return branches

# Decrease (-) all inhibitor weights given a certain step size
def reduce_conc(node,step,init_concentrations):
    #reduce ALL inhibitors by the step size
    n1 = node.copy()
    for x in range(0, len(n1)):
        n1[x] = n1[x] * (1.0/step)
    return n1

# opt_w() Returns optimal weights

# Score optimization using alternating optimization sub-rounds
# Round 1: all varyations of +/- inhibitor pairs are generated
# Round 2: all concentrations are diluted by equal percentages, such that minimum threshold on-target activity (ex: 90%) is maintained

# During R1 we examine the concentrations of the inhibitors - those that are more than 100-fold different from their starting concentration represent anomalous results (ex: a case where 1 inhibitor is better than 2, so the concentration of the 2nd is being radically reduced). We exclude these possibilities by comparing the concentrations to the initial concentrations. Due to the nature of the R1/R2 rounds (with minimization of all concentrations in R2), we find that using a lower limit that is more forgiving than the upper limit (ex: no more than 100x greater, 1,000x less) improves program performance
# This test is performed in the R1 branching method - otherwise, we might get more than minimum threshold on-target inhibition if we perform the test during subsequent steps as well
# This means that the final concentrations could be less than 1000X following R2 minimization, but none will be more than 100X.

# Combinations are then scored, all conditions (concentration combinations) that exceed the top-scoring condition from the prior round are used as the basis for the next round
# If a round yields no scores that are better than the prior round, the method is complete
# The step for each round represents the fold-change in inhibitor concentration. Ex: step=2 would double/halve the concentration of inhibitors
# If maximization is taking a long time, a cutoff (by time, in seconds) can be set, or the step sizes can be increased

# We find an R1 step of 5 (fivefold increase / decrease) to be reasonable maximum for the first step, matches what might be performed experimentally & improves program speed.
# Reduce to R1 step to 2 if you want to slightly improve chance of retrieving the best estimations (no guarantee that this will improve results, however)
# Ex (with R1 step=5): Given 100nM inhibitor, will try 500nM, 2.5uM (stepping up) and 20nM, 4nM (stepping down). Actual concentrations will vary due to the dilution step in R2.
# We find a much smaller step size (1.1 for example) to be reasonable for R2
# Too large an R2 size may result in on-target activity undesirably greater than threshold activity
# All step sizes should be set to greater than 1
def opt_w(on_target_values, off_target_values, prior, noise_variance, influence, num_w, concentrations, otit, R1, R2):

    R1_step = R1
    R2_step = R2

    init_time = time.time()
    init_concentrations = concentrations.copy()

    on_target = make_on_t_dist(on_target_values, concentrations, prior, noise_variance, influence, False)
    off_target = make_off_t_dist(off_target_values, concentrations, noise_variance, False)

    #set our initial top score and weights to those we already have
    top_score = get_score(on_target, off_target)
    top_concentrations = concentrations

    nodes = []
    nodes.append(concentrations)

    while len(nodes) > 0:

        #if we have been running the method for too long, then break
        curr_t = time.time()
        if (curr_t - init_time) > 300: # if we have been running for > 5 minutes
            break

        #R1: first generate +/- combinations of the inhibitors: don't score them yet
        R1_variations = branch(nodes,R1_step,init_concentrations)

        #R2: lower the concentration of the inhibitors in each combination so that we get as close to threshold on-target inhibition (ex: >= 90% activity) as possible, reducing off-target effects
        #Note, we still do not score the combinations
        
        inhibition_threshold = otit
        R2_variations = []
        
        for combination in R1_variations:
            
            oti_test = False
            concs = combination.copy()
            on_target_inhibition = (make_on_t_dist(on_target_values, concs, prior, noise_variance, influence, True))[1]
            
            if on_target_inhibition >= inhibition_threshold: #only return a value later in the method if we get at least threshold inhibition
                oti_test = True

            if oti_test:
                
                while True:
                
                    new_concs = reduce_conc(concs, R2_step, init_concentrations)

                    #if the reduction was reasonable (that is, we don't get an empty list because one of the concentrations was absurdly small)
                    if new_concs == []:
                        break
                
                    else:
                        new_on_target_inhibition = (make_on_t_dist(on_target_values, new_concs, prior, noise_variance, influence, True))[1]
                        #as long as we stay above threshold, keep decreasing the [] of the inhibitors
                        if new_on_target_inhibition >= inhibition_threshold:
                            concs = new_concs
                            #on_target_inhibition = new_on_target_inhibition
                        else:
                            break
                
            if oti_test: #as long as our initial set was good enough, so we can be sure we have at least 1 reasonable result
                R2_variations.append(concs)

        #Now, for every combination (with close to threshold on-target inhibition) we can score and evaluate the on-versus off targets using JSD

        #reset nodes to nothing for now
        nodes = []
        scores = []

        for combo in R2_variations:
            on_target_bins = make_on_t_dist(on_target_values, combo, prior, noise_variance, influence, False)
            off_target_bins = make_off_t_dist(off_target_values, combo, noise_variance, False)
            score = get_score(on_target_bins, off_target_bins)
            if score > top_score:
                nodes.append(combo)
                scores.append(score)

        #Now, set our threshold score for the next round to be the best score from our current round as long as we have at least 1 good result
        if len(nodes) > 0:
            top_score = max(scores)
            top_concentrations = nodes[scores.index(top_score)]

    #the variable top_score represents our maximum-scoring results, and the variable top_concentrations represents the concentrations that yielded this score
    return top_concentrations
     
# Jensen-Shannon weighed inhibitor scoring method
# First calculates optimal weights, then scores the combination
def JS_WE(inpt):

    k_ind = inpt[0][0]
    inhibitors = inpt[0][1]
    on_target_values = inpt[0][3]
    off_target_values = inpt[0][4]
    off_target_ind = inpt[0][6]

    #unweighted concentrations in uM
    prior_concentrations = inpt[0][5]

    ot_prior = inpt[1]
    noise_variance = inpt[2]
    influence = inpt[3]

    n_iter = inpt[4]
    t_iter = inpt[5]
    prior_t = inpt[6]

    R1 = inpt[8]
    R2 = inpt[9]

    on_target_inhibition_threshold = inpt[7]

    if n_iter % 20 == 0 and n_iter != 0:
        #if multiprocessing is set to return threads in order, the next line can be used to estimate the time remaining in the program
        time_r = (((time.time() - prior_t) * t_iter) / n_iter) - (time.time() - prior_t)
        print('Working on iter:', n_iter, '/',  t_iter, flush=True)

    final_concentrations = opt_w(on_target_values, off_target_values, ot_prior, noise_variance, influence, len(inhibitors), prior_concentrations, on_target_inhibition_threshold, R1, R2)
    
    on_target_res = make_on_t_dist(on_target_values, final_concentrations, ot_prior, noise_variance, influence, True)
    off_target_res = make_off_t_dist(off_target_values, final_concentrations, noise_variance, True)
    
    on_target_distribution = on_target_res[0]
    on_target_percent_inhib = on_target_res[1]
    off_target_distribution = off_target_res[0]
    off_target_average_percent_inhib = off_target_res[1]
    num_off_target_kinases = off_target_res[2]
    off_target_percent_inhib_values = off_target_res[3]
    score = get_score(on_target_distribution, off_target_distribution)

    return (k_ind, inhibitors, score, on_target_percent_inhib, on_target_distribution, off_target_average_percent_inhib, num_off_target_kinases, off_target_distribution, off_target_percent_inhib_values, final_concentrations, off_target_ind)

#################### main function

# This function runs the program
# Given a dataset, will return the highest-scoring inhibitor combination for each of ALL targets (kinases) in the dataset, up to the maximum user-defined combination number of inhibitors
# Operates on a single target kinase at a time, although it will eventually proccess all targets in the dataset

def main_function(dataset,prior_settings,noise_variance,max_combination_iter,influence,output_name,num_p,otit,number_top_UW_scores_to_maximize,names_to_exclude,R1,R2):

    #import our data from the dataset
    kinase_names, compound_names, data = import_from_dataset(dataset)

    #Generate the selectivity prior that will be used for all analyses
    on_target_prior = generate_prior(prior_settings)

    #Get the indicies (columns) of the targets to exclude from the off-target distribution (ex: disease mutants in a dataset)
    ind_to_exclude = get_indicies_to_exclude(names_to_exclude, dataset)

    #For each combination # of inhibitors... (ex: i=1, i=2...)
    for x in range(1,max_combination_iter+1):
        
        curr_t = time.time()
        print('Initializing inhibitor combinations at combination number: ' + str(x))
        
        #build a list of the indicies to score
        i_target_k = get_i_target_k(data,otit)

        combos = get_combinations(i_target_k, x, data, otit, ind_to_exclude)
                
        print('Time: "',time.time() - curr_t)
              
        print('Calculating unweighed initial scores at combination number: ' + str(x))

        print('Total number of combinations: ' + str(len(combos)))
              
        curr_t = time.time()
        p_pool = multiprocessing.Pool(processes=num_p)
        num_r = range(0,len(combos))
        t_r = len(combos)
        p_input = list(zip(combos, repeat(on_target_prior), repeat(noise_variance), repeat(influence), num_r, repeat(t_r), repeat(curr_t)))
        UW_JS_scores = p_pool.map(JS_UW, p_input)
        p_pool.close()
        p_pool.join()
        
        temp_UW = [ [] for x in range(0,len(kinase_names))]
        for y in range(0, len(UW_JS_scores)):
            res = UW_JS_scores[y]
            temp_UW[res[0]].append(res)

        top_scores = []
        for k in temp_UW:
            sorted_UW = sorted(k, key=lambda res: res[2], reverse=True)
            for res in sorted_UW[:number_top_UW_scores_to_maximize]:  #HOW MANY TOP SCORES TO MAXIMIZE
                top_scores.append(res)

        print('\nTime: "',time.time() - curr_t)
                
        print('Optimizing Inhibitor Weights at combination number: ' + str(x)+' (only maximizes top '+str(number_top_UW_scores_to_maximize)+' combinations per kinase)')

        print('Total number of combinations: ' + str(len(top_scores)))
        
        curr_t = time.time()
        num_r = range(0, len(top_scores))
        t_r = len(top_scores)        

        p_pool = multiprocessing.Pool(processes=num_p)
        p_input = list(zip(top_scores, repeat(on_target_prior), repeat(noise_variance), repeat(influence), num_r, repeat(t_r), repeat(curr_t), repeat(otit), repeat(R1), repeat(R2)))
        WE_JS_scores = p_pool.map(JS_WE, p_input)
        p_pool.close()
        p_pool.join()

        temp_WE = [ [] for x in range(0, len(kinase_names))]
        for y in range(0, len(WE_JS_scores)):
            res = WE_JS_scores[y]
            temp_WE[res[0]].append(res)
            
        top_scores = []
        for k in temp_WE:
            sorted_WE = sorted(k, key=lambda res: res[2], reverse=True)
            if len(sorted_WE) >= 1:
                top_scores.append(sorted_WE[0])

        print('Time: "',time.time() - curr_t)

        print('Organizing results for output...')

        #for each best score, get the kinase name, inhibitor names, and output the results to an excel file
        output_to_excel = []
        for k_target in top_scores:
            kinase_name = kinase_names[k_target[0]]
            inhibitor_names = [compound_names[i] for i in k_target[1]]
            score = k_target[2]
            on_target_percent_inhib = k_target[3]
            on_target_distribution = k_target[4]
            off_target_average_percent_inhib = k_target[5]
            num_off_target_kinases = k_target[6]
            off_target_distribution = k_target[7]
            off_target_percent_inhib_values = k_target[8]
            inhibitor_weights = k_target[9]
            off_target_ind = k_target[10][0]
            off_target_names = [kinase_names[x] for x in off_target_ind]
            output_to_excel.append([kinase_name,
                                    str(x),
                                    inhibitor_names,
                                    score,
                                    on_target_percent_inhib,
                                    off_target_average_percent_inhib,
                                    num_off_target_kinases,
                                    inhibitor_weights,
                                    on_target_distribution,
                                    off_target_distribution,
                                    off_target_percent_inhib_values,
                                    off_target_names])
                                    
        output_columns = ['kinase target',
                   'number of inhibitors',
                   'inhibitors',
                   'JS Distance Score',
                   'Target kinase %inhib',
                   'Average off-target kinase %inhib',
                   'Number of off-target kinases',
                   'inhibitor concentrations (micromolar uM)',
                   'on-target probability distribution',
                   'off-target probability distribution',
                   'off-target %inhibition values',
                   'off-target names'] 

        #print(output_to_excel)
        o_name = output_name+'.xlsx'
        if x==1:
            output = pd.ExcelWriter(o_name)
            pdframe = pd.DataFrame(output_to_excel,columns=output_columns)
            pdframe.to_excel(output, sheet_name=str(x))
            output.save()
        else:
            output = pd.ExcelWriter(o_name, engine = 'openpyxl', mode='a')
            pdframe = pd.DataFrame(output_to_excel,columns=output_columns)
            pdframe.to_excel(output, sheet_name=str(x))
            output.save()

    
################### main

if __name__ == '__main__':

    ######################## USER DEFINED DATASET, should note correct number of rows to use (cols too if desired)

    # Note, dataset must have the same formatting as the excel file for the PKIS2 dataset.
    # User should add null values for some of the initial daata descriptions (ex: chemotype family) if not present in the dataset to ensure overal formatting is consistent. 

    dataset=pd.read_excel("./110421_Karaman_2008_dataset.xlsx", engine='openpyxl', nrows=39, usecols=range(0,328))

    print(dataset)

    ######################## USER DEFINED SETTINGS
    
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
    output_name = 'my_new_results' #Name of the output file
    num_p = 1 #number of proccesses to call during multiprocessing

    ####################### RUN

    prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size)) #package the prior settings for the program
    main_function(dataset, prior_settings, noise_variance, max_combination_iter, influence, output_name, num_p, on_target_inhibition_threshold, number_top_UW_scores_to_maximize)
 
