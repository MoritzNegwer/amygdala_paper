#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 11:41:50 2020

@author: wirrbel
"""

import pandas as pd
import os 
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, fft, signal

#import reference table 
#Don't forget to change if using a different atlas
#referencetable = pd.read_csv("/media/kwaakbuntu/MN_2/atlas_p21/Regions IDs_backup.csv", encoding = "ISO-8859-1", usecols=[0,1,6,7])

#Create overview dataframes 
Time_overview = pd.DataFrame()
Amp_overview = pd.DataFrame()
IEI_overview = pd.DataFrame()
minis_overview = pd.DataFrame()
#merged = merged.append(referencetable)


#This is where all brain folders are stored. 
#Assumes Clearmap setup, i.e. one folder per brain, with an "Annotated_counts.csv" inside
rootdir = '/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/Amygdala_Data/'


def walk_through_folders(rootdir, Time_overview, Amp_overview, IEI_overview,minis_overview):
    #iterate 
    for subdir, dirs, files in os.walk(rootdir, followlinks=True):
        for file in files:
            if file.endswith("xlsx"):
                #extract name of the brain
                Brain_name = (os.path.split(subdir)[-1]+'_'+file)
                #print(os.path.split(subdir)[-1])
                
                #read Annotated_counts.csv into dataframe
                #data = pd.read_excel(os.path.join(subdir,file), header=None,usecols=[1,2],names=["ID",Brain_name])
                
                #Extract event time
                Time = pd.read_excel(os.path.join(subdir,file), header=0,usecols=[1],names=[Brain_name])
                #sort ascending
                Time = Time.sort_values(by=Brain_name)
                #add to overview table
                Time_overview = pd.concat([Time_overview,Time[Brain_name]],axis=1)  
                print('file found:', file, ' entries: ',len(Time))
                '''
                #extract event Amp
                Amplitude = pd.read_excel(os.path.join(subdir,file), header=0,usecols=[2],names=[Brain_name])
                #add to overview table
                #Amp_overview = pd.merge(Amp_overview,Amplitude[Brain_name],how='right',left_index=True,right_index=True)
                Amp_overview = pd.concat([Amp_overview,Amplitude[Brain_name]],axis=1)  
                            
                #extract ISI
                IEI = Time.diff()
                #add to overview table

                IEI_overview = pd.concat([IEI_overview,IEI[Brain_name]],axis=1)
                '''
                #===average minis:===
                #find txt file 
                txtfile = file.replace("xlsx","txt")
                txtfile_path = os.path.join(subdir,txtfile)
                #print("txtfile_path = ",txtfile_path)
                
                #get raw_rec 
                raw_rec = load_raw_mIPSC_from_txt(txtfile_path)
                
                #define minis_list 
                minis_list = pd.DataFrame()
                mean_minis = pd.DataFrame()
                minis_list_temp = np.empty((20000,1))
                
                #collect minis 
                for event in Time.iterrows():
                    #assuming a sampling rate of 20Khz 
                    event_sample = int(event[1][0]*20)
                    print ('event sampled',event[1][0])
                    #safety check if event_sample is not before 10k samples 
                    if event_sample < 10000 or event_sample > 11090000:
                        pass
                    else:                                      
                        minis_list_temp = add_mIPSCs(raw_rec,event_sample,minis_list_temp)
                    
                
                minis_list = pd.DataFrame(minis_list_temp)
                
                #troubleshooting: save all 300 minis in a single table
                #minis_list.to_excel(os.path.join('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/',"mIPSC_Amygdala_minis_overview_v06"+Brain_name+"_allvalues.xlsx"))

                '''
                #average minis 
                minis_list = minis_list.mean(axis=1)
                
                #remove 50Hz noise and prominent 100Hz/200Hz harmonics
                plot_noise_signals(minis_list)
                
                minis_list = filter_noise(minis_list,50,20000,30)
                minis_list = filter_noise(minis_list,100,20000,30)
                minis_list = filter_noise(minis_list,200,20000,30)
                plot_noise_signals(minis_list)
                
                #transfer back into panda df 
                mean_minis[Brain_name] = minis_list
                                
                #add to minis_overview (one averaged mini per cell)
                minis_overview = pd.concat([minis_overview,mean_minis[Brain_name]],axis=1)
                
                #troubleshooting: Export every mini average file 
                #mean_minis.to_excel(os.path.join('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/',"mIPSC_Amygdala_minis_overview_v06"+Brain_name+".xlsx"))
                print("average minis added! ",len(mean_minis))
                '''

    #output to xlsx
    #Time_overview.to_excel(os.path.join('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/',"mIPSC_Amygdala_Time_overview_v06.xlsx"))
    #Amp_overview.to_excel(os.path.join('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/',"mIPSC_Amygdala_Amp_overview_v06.xlsx"))
    #IEI_overview.to_excel(os.path.join('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/',"mIPSC_Amygdala_IEI_overview_v06.xlsx"))
    #minis_overview.to_excel(os.path.join('/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/',"mIPSC_Amygdala_minis_overview_v06.xlsx"))
    return raw_rec, minis_list_temp, minis_list
    

#====Histograms====
def sort_into_genotypes(input_dataframe,genotype1,genotype2):
    #Define which dataframe you want a histogram of 
    Histogram_dataframe = input_dataframe
    
    WT_Hist = pd.DataFrame()
    Het_Hist = pd.DataFrame()
    
    #sort into genotypes 
    for column_name in Histogram_dataframe:
        #check if name contains WT. If so, copy to WT dataframe
        if genotype1 in column_name:
            WT_Hist = pd.concat([WT_Hist,Histogram_dataframe[column_name]],axis=1)  
        #same, but with Het 
        elif genotype2 in column_name:
            Het_Hist = pd.concat([Het_Hist,Histogram_dataframe[column_name]],axis=1) 
        #If there are neither genotype names in there, 
        #something went wrong! Not throwing error, just notifying user 
        else:
            print ("Warning: No genotype info found in filename of file: ", column_name)
    
    #return both sorted dataframes
    return WT_Hist,Het_Hist


def generate_cumulative_plot(input_dataframe,range_min,range_max,save_dir,save_filename):
    #example: generate_cumulative_plot(IEI_overview,0,3000,'/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/','IEI_v06')
    
    #sort by genotypes 
    WT_Hist, Het_Hist = sort_into_genotypes(input_dataframe, 'WT', 'Het')
        
    #make histogram
    
    #max over all values (used in earlier version)
    #global_max = pd.concat([WT_Hist.stack(),Het_Hist.stack()],axis=0).max()
    
    #Generate histogram for WT and Het 
    #(Note: .stack() makes a single series out of the dataframe, super useful) 
    WT_values, WT_base = np.histogram(WT_Hist.stack(),bins=100,range=(range_min,range_max))
    Het_values, Het_base = np.histogram(Het_Hist.stack(),bins=100,range=(range_min,range_max))
    
    #generate a cumulative plot of the bins 
    WT_cumulative = np.cumsum(WT_values)/len(WT_Hist.stack())
    Het_cumulative = np.cumsum(Het_values)/len(Het_Hist.stack())
    
    # plot the cumulative function
    plt.plot(WT_base[:-1], WT_cumulative, c='black')
    plt.plot(Het_base[:-1], Het_cumulative, c='red')
    
    #show plot 
    #plt.show()
    
    #export as SVG 
    plt.savefig(os.path.join(save_dir,save_filename+'.svg'))
    
    #clear plot
    plt.clf()

def run_KS_test(input_dataframe):
    #example: run_KS_test(IEI_overview)
    
    #sort by genotypes 
    WT_Hist, Het_Hist = sort_into_genotypes(input_dataframe, 'WT', 'Het')
    
    #stack and sort dataframes 
    WT_Hist = WT_Hist.stack().sort_values()
    Het_Hist = Het_Hist.stack().sort_values()
    
    #run 2-sample Kolmogorov-Smirnov test
    KS_stat, p_value = stats.ks_2samp(WT_Hist,Het_Hist)
    
    print("2-Sample KS test, KS Stat: ", KS_stat, ", p value: ",p_value)
    
    return KS_stat, p_value


#==== Average mIPSC plot =====

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

def add_mIPSCs (raw_rec, event_ms,minis_list):
    #for test purposes 
    #event_ms=30306
    
    event_df = pd.DataFrame()
   
    #return values 0.5s before +0.5s after peak 
    event_df = raw_rec.iloc[int(event_ms-10000):int(event_ms+10000)]
    
    #reset index 
    event_df = event_df.reset_index()
    event_df = event_df.drop(labels=['index'],axis=1)
    
    #normalize for average 0.4s before peak 
    average= event_df.iloc[0:8000].mean()[0]
    
    #substract from all values to normalize 
    event_df = event_df-average
    
    #rolling average 
    #event_df = event_df.rolling(250,min_periods=1).mean()
    #event_df.ewm(alpha=0.01,adjust=False).mean()
    
    #convert to numpy for greater speed 
    event_df = event_df.to_numpy()
    #print(event_df.shape)
    #add event to list 
    #minis_list = pd.concat([minis_list, event_df],axis=1,ignore_index=True)
    
    #print('minis_list_temp=',minis_list.shape)
    minis_list = np.concatenate([minis_list,event_df],axis=1)

        
    #event_df = event_df.to_numpy()
    
    return minis_list

#========
#frequency filtering 
def plot_noise_signals(minis_df):
    #assumes that minis_df contains a single column 
    
    #check if minis_df is a pandas df 
    if not type(minis_df) is np.ndarray:
        minis_df = minis_df.to_numpy()
    
    #perform FFT
    yf = fft.fft(minis_df)
    x = fft.fftfreq(yf.size, 1/20000) #assuming 20Khz sampling 
       
    #plot frequency 0-300 Hz
    plt.plot(np.abs(x[:300]),np.abs(yf[:300]))
    plt.show()
    #plt.clf()

def filter_noise(minis_df,freq_to_filter,sampling_rate,quality_setting):
    #check if minis_df is a pandas df 
    if not type(minis_df) is np.ndarray:
        minis_df = minis_df.to_numpy()

    
    #creates a notch filter and filters minis_df 
    b,a = signal.iirnotch(freq_to_filter,quality_setting,sampling_rate)
    
    #filter minis_df
    minis_df = signal.lfilter(b,a,minis_df)
    
    return minis_df
    
#actual commands    
raw_rec, minis_list_temp, minis_list = walk_through_folders(rootdir, Time_overview, Amp_overview, IEI_overview, minis_overview)

'''
#generate plots 
generate_cumulative_plot(IEI_overview,0,3000,'/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/','IEI_v06')
generate_cumulative_plot(Amp_overview,6,30,'/media/wirrbel/Moritz_Negwer/mIPSCs_Amygdala/','Amp_v06')

#run KS tests 
#(Note that the results are not automatically saved 
run_KS_test(IEI_overview)
run_KS_test(Amp_overview)
'''

