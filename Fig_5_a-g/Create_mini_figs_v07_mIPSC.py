#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 01:18:24 2024

@author: wirrbel
"""

import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, fft, signal

def load_raw_mIPSC_from_txt(txtfile_path):

    #load txt file. Assumes that file is organized into sweeps with a single 's' timing column at the beginning 
    #raw_rec = pd.read_table('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/Amygdala_Data/2020-10-12 WT P13 Amygdala mIPSCs/Minis_MN_121020_003.txt')
    raw_rec = pd.read_table(txtfile_path,skiprows=1,dtype=float,header=None,)


    #drop first column and first row (ms timing and headers)
    #raw_rec = raw_rec.drop(labels='s',axis=1)
    raw_rec = raw_rec.drop([0],axis=1)
    #raw_rec = raw_rec.drop([0],axis=0)

    #stack, will multi-index the numeric index and sweep name 
    raw_rec = raw_rec.stack(level=0)

    #sort according to numeric index 
    #otherwise the index for each sweep would be next to each other (index 1, sweep 1, sweep 2... Not good.)
    raw_rec = raw_rec.reset_index()

    #sort by frame so that everything from frame 1 is subsequent
    raw_rec = raw_rec.sort_values(by=['level_1','level_0'])    

    #clean up indices that now are column names 
    raw_rec = raw_rec.drop(labels=['level_0','level_1'],axis=1)

    #reset index once more to account for 300-step timing 
    raw_rec = raw_rec.reset_index()
    raw_rec = raw_rec.drop(labels=['index'],axis=1)

    return raw_rec

#load representative WT mIPSC sample
minis_wt = load_raw_mIPSC_from_txt("/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/mIPSCs_Amygdala/Amygdala_Data/2020-10-28 WT P13 Amygdala mIPSCs/Minis_MN_281020_000.txt")

#load representative Het mIPSC sample
minis_het = load_raw_mIPSC_from_txt("/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/mIPSCs_Amygdala/Amygdala_Data/2020-10-29_m2_P14_Amygdala_Het/Minis_MN_291020_004.txt")


#define example mini
example_sec_wt = minis_wt[2307000:2407000]

#smooth (1ms moving average) and prepare for graph
ex_wt_mean = example_sec_wt.rolling(20).mean()
ex_wt_mean = ex_wt_mean - ex_wt_mean.mean()
#ex_wt_mean = ex_wt_mean.drop(columns=["index"])
ex_wt_mean = ex_wt_mean.reset_index(drop=True)

#add to plot 
plt.plot(ex_wt_mean,c="black")

#save plot
plt.savefig("/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/mIPSCs_Amygdala/Amygdala_Data/example_second_WT_201020_000_1msec_moving_average.svg")

#define example mini
example_sec_het = minis_het[2307000:2407000]

#smooth (1ms moving average) and prepare for graph
ex_het_mean = example_sec_het.rolling(20).mean()
ex_het_mean = ex_het_mean - ex_het_mean.mean()
#ex_het_mean = ex_het_mean.drop(columns=["index"])
ex_het_mean = ex_het_mean.reset_index(drop=True)

#add to plot 
plt.plot(ex_het_mean,c="red")
plt.savefig("/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/mIPSCs_Amygdala/Amygdala_Data/example_second_Het_291020_004_1msec_moving_average.svg")