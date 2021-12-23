#!/bin/python3.9.4
#Ian Outhwaite, 2021

#################### import

import pandas as pd
import numpy as np
import xlsxwriter
import openpyxl
import itertools
from itertools import repeat
import operator
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
from time import sleep

#################### functions

#get kinase names, inhibitor names, and a 2D array of data from PKIS2
def import_from_PKIS2(dataset):
    kinase_names = list(dataset.columns.values)[7:]
    compound_names = dataset['Compound'].tolist()
    array_values = dataset[dataset.columns[7:]].to_numpy()
    return kinase_names, compound_names, array_values

#generates a left-skewed distribution using a Poisson distribution as a basis for the shape scaled to 100 with the mean at 100
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

#generates a beta distribution scaled to 100
def get_Beta_Prior(a, b, _size):
    prior = beta.rvs(a, b, scale=1, size=_size)
    prior = [x * 100.0 for x in prior]
    return prior

#depending upon the prior settings, generate a set of values and fill probability bins. Return the bins and the number of entries used to fill them (to help normalize values later on) 
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

#returns normally distributed random values
def generate_noise(noise_variance):
    noise_v = np.random.normal(0, noise_variance, 100)
    return noise_v

#returns Jensen-Shannon distance metric
def get_score(on, off):
    JS_score = distance.jensenshannon(on, off, 2.0)
    return JS_score

#Given a % inhibition (activity) value in a dataset, returns the Ki of the inhibitor assuming the inhibitor is used at a concentration of 1uM
def get_Ki(percent_inhibition):
    if percent_inhibition == 100:
        percent_inhibition = 99.5
    return ((1.0/(percent_inhibition/100.0))-1)

#returns the concentration of a drug needed to reach a certain threshold (in % activity). Returns the concentration in uM, requires the Ki of the drug.
# ( [drug] + Ki ) * activity = [drug]
# ( [drug] * activity + Ki * activity ) = [drug]
# if we use concentrations relative to 1uM, then
# Ki * activity = [drug] * (1-activity)
# (Ki * activity) / (1-activity) = [drug]
#activity must be in 0->1 range, not 0->100
def get_concentration(Ki, threshold):
    t = threshold / 100.0
    c = (Ki * t) / (1 - t)
    return c
    
#given a set of inhibition values and inhibitor concentrations against a single target, returns the total inhibition at a single target
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
             
#returns the inhibitor indicies that target each kinase
def get_i_target_k(data):
    res = []
    for k in range(0,len(data[0])):
        i_data = data[:,k]
        target_i = [i for i,val in enumerate(i_data) if val > 90]
        res.append(target_i)
    return res

#given a slice of the data matrix corresponding to the correct inhibitor data, recalculate the activity values and pass them back
#conc is a value in uM calculated by the get_concentration() method
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

#gets all combinations, and for each combination also pulls a copy of the on-target values and the off-target values
#The values are adjusted given the amount of inhibitor needed to reach a potency of 90% against the target kinase
def get_combinations(i_target_k, cmbn, data):
    res = []
    for k in range(0, len(i_target_k)):
        inhibs = i_target_k[k]
        for c in itertools.combinations(inhibs,cmbn):
            c = list(c)
            on = []
            off = []
            concentrations = []
            activity_thresh = 90
            for i in c:
                target = data[i,k]
                Ki = get_Ki(target)
                conc = get_concentration(Ki, activity_thresh)
                concentrations.append(conc)
                on.append(target)
                off_target_selector = [x for x in range(data.shape[1]) if x != k]
                off.append(data[i, off_target_selector])
            on = np.array(on)
            off = np.array(off)
            res.append((k,c,on,off,concentrations))
    return res

#generates an on-target distribution
#bin size = 5
def make_on_t_dist(values, concentrations, prior, noise_variance, influence, extended_output):
    prior_size = prior[1]
    size_ot = prior_size * influence
    relative_weight = size_ot/100.0
    ot_inhib = get_t_inhib_weights(values, concentrations)
    noise_dist = np.random.normal(0, noise_variance, 100)
    noise = [ot_inhib + noise_dist[x] for x in range(0, 100)]
    bins = prior[0].copy()
    for x in range(0, len(noise)):
        temp = int(noise[x]/5)
        if temp < 20 and temp >= 0:
            bins[temp] += relative_weight
        elif temp >= 20:
            bins[len(bins)-1] += relative_weight
        else:
            bins[0] += relative_weight
    ot_dist_size = prior_size + size_ot
    bins[:] = [(x/ot_dist_size) for x in bins]
    if extended_output == False:
        return bins
    else:
        return (bins, ot_inhib)
        
#generates the off-target distribution, bin size=5
def make_off_t_dist(values, weights, noise_variance, extended_output):
    off_t_dist = []
    off_t_inhib_vals = []
    noise_dist = np.random.normal(0, noise_variance, 100)
    numzeros = 0
    for x in range(0, len(values[0])):
        p_vals = values[:, x]
        if np.sum(p_vals) != 0:
            offt_inhib = get_t_inhib_weights(p_vals, weights)
            off_t_inhib_vals.append(offt_inhib)
            noise_vals = [offt_inhib + noise_dist[x] for x in range(0,100)]
            off_t_dist += noise_vals
        else:
            off_t_inhib_vals.append(0)
            numzeros += 1
    nvals = len(off_t_dist)
    bins = [0 for i in range(0,20)]
    for x in range(0, nvals):
        temp = int(off_t_dist[x]/5)
        if temp < 20 and temp >= 0:
            bins[temp] += 1.0
        elif temp >= 20:
            bins[len(bins)-1] += 1.0
        else:
            bins[0] += 1.0
    for y in range(0,numzeros): #this will ensure that we normalize by the same number
        bins[0] += 100.0 #if we don't add 100 points of noise above per measurment, change this value as well
    normfactor = np.sum(bins)
    bins[:] = [(x/normfactor) for x in bins]
    if extended_output == False:
        return bins
    else:
        avg_offt = np.mean(off_t_inhib_vals)
        num_offt = np.count_nonzero(off_t_inhib_vals)
        return (bins, avg_offt, num_offt, off_t_inhib_vals)


#JS unweighed inhibitor scoring
#for each combination, pass in equal divisions of each inhibitor (ex: 50-50), then score the on-target versus off-target distributions
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

    if n_iter % 1000 == 0 and n_iter != 0:
        time_r = (((time.time() - prior_t) * t_iter) / n_iter) - (time.time() - prior_t)
        print('Working on iter:', n_iter, '/',  t_iter, flush=True)

    new_conc = [(concentrations[x]/len(inhibitors)) for x in range(0, len(concentrations))]

    on_target = make_on_t_dist(on_target_vals, new_conc, prior, noise_variance, influence, False)
    off_target = make_off_t_dist(off_target_vals, new_conc, noise_variance, False)
    score = get_score(on_target, off_target)

    return (k_ind, inhibitors, score, on_target_vals, off_target_vals, new_conc)

#generate all possible variations of +/- inhibitor weights given a certain step size
def branch(nodes,step):
    branches = []
    for n in nodes:
        for x in range(0, len(n)):
            for y in range(0, len(n)):
                if x != y:
                    n1 = n.copy()
                    n1[x] = n1[x] * step
                    n1[y] = n1[y] * (1.0/step)
                    if (n1[x] <= 1000 and n1[x] >= 0.0001):
                        if (n1[y] <= 1000 and n1[y] >= 0.0001):
                            branches.append(n1)
    return branches

#returns optimal weights - simple maximization, if the landscape is too complicated this may get caught at local maxes. Prioritizes speed.
def opt_w(on_target_values, off_target_values, prior, noise_variance, influence, num_w, concentrations):
    step = 2 #% interval to try. Ex: 2 +/- a factor of 2 (double, half) of an inhibitor at each step
    nodes = []
    nodes.append(concentrations)
    top_s = 0
    top_w = []
    while len(nodes) > 0:
        good_n = []
        good_s = []
        for n in range(0, len(nodes)):
            cur_n = nodes[n]
            on_target_res = make_on_t_dist(on_target_values, cur_n, prior, noise_variance, influence, True)
            if on_target_res[1] >= 90:
                on_target = on_target_res[0]
                off_target = make_off_t_dist(off_target_values, cur_n, noise_variance, False)
                score = get_score(on_target, off_target)
                if score > top_s:
                    good_s.append(score)
                    good_n.append(cur_n)
        if len(good_s) > 0:
            top_s = max(good_s)
            top_w = good_n[good_s.index(top_s)]
        nodes = branch(good_n,step)
    return top_w
                
#Jensen-Shannon weighed inhibitor scoring
def JS_WE(inpt):

    k_ind = inpt[0][0]
    inhibitors = inpt[0][1]
    on_target_values = inpt[0][3]
    off_target_values = inpt[0][4]

    #unweighted concentrations in uM
    prior_concentrations = inpt[0][5]

    ot_prior = inpt[1]
    noise_variance = inpt[2]
    influence = inpt[3]

    n_iter = inpt[4]
    t_iter = inpt[5]
    prior_t = inpt[6]

    if n_iter % 20 == 0 and n_iter != 0:
        time_r = (((time.time() - prior_t) * t_iter) / n_iter) - (time.time() - prior_t)
        print('Working on iter:', n_iter, '/',  t_iter, flush=True)

    final_concentrations = opt_w(on_target_values, off_target_values, ot_prior, noise_variance, influence, len(inhibitors), prior_concentrations)
    
    on_target_res = make_on_t_dist(on_target_values, final_concentrations, ot_prior, noise_variance, influence, True)
    off_target_res = make_off_t_dist(off_target_values, final_concentrations, noise_variance, True)
    
    on_target_distribution = on_target_res[0]
    on_target_percent_inhib = on_target_res[1]
    off_target_distribution = off_target_res[0]
    off_target_average_percent_inhib = off_target_res[1]
    num_off_target_kinases = off_target_res[2]
    off_target_percent_inhib_values = off_target_res[3]
    score = get_score(on_target_distribution, off_target_distribution)

    return (k_ind, inhibitors, score, on_target_percent_inhib, on_target_distribution, off_target_average_percent_inhib, num_off_target_kinases, off_target_distribution, off_target_percent_inhib_values, final_concentrations)

def main_function(dataset,prior_settings,noise_variance,max_combination_iter,influence,output_name,num_p):
    kinase_names, compound_names, data = import_from_PKIS2(dataset)
    on_target_prior = generate_prior(prior_settings)
    for x in range(1,max_combination_iter+1):
        curr_t = time.time()
        print('Initializing inhibitor combinations at combination number: ' + str(x))
        
        #build a list of the indicies to score
        i_target_k = get_i_target_k(data)

        combos = get_combinations(i_target_k, x, data)
                
        print('Time:"',time.time() - curr_t)
              
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
            for res in sorted_UW[:3]:  #HOW MANY TOP SCORES TO MAXIMIZE
                top_scores.append(res)

        print('\nTime:"',time.time() - curr_t)
                
        print('Optimizing Inhibitor Weights at combination number: ' + str(x)+' (only maximizes top 3 combinations per kinase)')

        print('Total number of combinations: ' + str(len(top_scores)))
        
        curr_t = time.time()
        num_r = range(0, len(top_scores))
        t_r = len(top_scores)        

        p_pool = multiprocessing.Pool(processes=num_p)
        p_input = list(zip(top_scores, repeat(on_target_prior), repeat(noise_variance), repeat(influence), num_r, repeat(t_r), repeat(curr_t)))
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

        print('Time:"',time.time() - curr_t)

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
                                    off_target_percent_inhib_values])
                                    
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
                   'off-target %inhibition values'] 

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

    dataset=pd.read_excel("./110421_Karaman_2008_dataset.xlsx", engine='openpyxl', nrows=39, usecols=range(0,328))

    print(dataset)
    
    prior_type = 2
    bdist_alpha = 3
    bdist_beta = 1
    mu = 700
    size = 10000
    prior_settings = (prior_type, (bdist_alpha, bdist_beta, mu, size))
    noise_variance = 2.5
    max_combination_iter = 3
    influence = 0.0
    output_name = 'TEST'
    num_p = 1

    main_function(dataset,prior_settings,noise_variance,max_combination_iter,influence,output_name,num_p)
