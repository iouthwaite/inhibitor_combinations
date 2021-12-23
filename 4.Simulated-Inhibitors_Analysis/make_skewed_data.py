#!/bin/python3.9
#Ian Outhwaite, 2021

import pandas as pd
from scipy.stats import skewnorm
import string
import random 

def main_function_make_skew(r):

    for ninhib in r:

        output = pd.ExcelWriter('C_simulated_'+str(ninhib)+'_inhibitors_100_kinases_binary2.xlsx')

        fake_data = []

        for i in range(0,ninhib): #number of inhibitors goes here

            letters = string.ascii_letters
            fake_regno = (''.join(random.choice(letters) for i in range(10)))
            fake_inhibitor_name = (''.join(random.choice(letters) for i in range(10)))
            fake_chemotype = (''.join(random.choice(letters) for i in range(10)))
            fake_smiles = (''.join(random.choice(letters) for i in range(10)))

            fake_inhib = []

            binmax = [x for x in range(0,21)]

            weights = [90,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5] #binary-inhibitors

            #weights = [0,15000,10000,9000,8000,7000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,6000,1500,3500] simulate even less specific inhibitors
            #weights = [30000,25000,25000,20000,15000,13000,10000,7000,5000,2700,2500,2300,2000,2000,2000,2000,2000,2000,2000,1500,3500] #to simulate less specific inhibitors
            #weights = [74966,11571,14045,12227,9963,7622,5722,4247,3071,2246,1800,1443,1169,1012,939,929,882,865,1005,1293,3588] #these are the frequencies from the PKIS2 dataset

            my_list = []
            for h in range(0,21):
                to_add = [binmax[h] for x in range(0,weights[h])]
                my_list += to_add

            for j in range(0,100): #number of kinase targets goes here (first of two locations)

                inhibition_val_upper_bound = random.choice(my_list) * 5
                fake_value = None

                if inhibition_val_upper_bound == 0:
                    fake_inhib.append(0)
                elif inhibition_val_upper_bound == 5:
                    fake_inhib.append(random.uniform(0.01,4.99)) 
                elif inhibition_val_upper_bound == 100:
                    fake_inhib.append(random.uniform(95,100))
                else:
                    fake_inhib.append(random.uniform(inhibition_val_upper_bound-5,inhibition_val_upper_bound-0.01))
            temp = [fake_regno,fake_inhibitor_name,fake_chemotype,fake_smiles,'temp','temp','temp'] + fake_inhib
            fake_data.append(temp)

        temp_names = []

        for i in range(0,100): #number of kinases again (second of two locations)
            temp_names.append(str(i))

        col = ['Regno','Compound','Chemotype','Smiles','>90','>80','>70'] + temp_names

        pdframe = pd.DataFrame(fake_data,columns=col)
        pdframe.to_excel(output)
        output.save()
        print ('wrote fake data')

if __name__ == '__main__':

    range_init = [200,300,400,500,600]
    
    main_function_make_skew(range_init)      

#EOF
