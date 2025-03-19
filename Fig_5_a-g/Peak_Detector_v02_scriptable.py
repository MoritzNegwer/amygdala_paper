# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import csv
import numpy as np
import os
import pandas as pd

'''
def input_range (min_range,max_range,list_of_measurements):
    #creates a temporary dict to store the values 
    output_list_raw = dict.fromkeys(list_of_measurements[0].keys(),[])
    
    #This is necessary to make sure each element of output_list_raw is a separate list. 
    #By default, dict.fromkeys generates a single object with several keys, which we don't want 
    for k,_ in output_list_raw.items():
        output_list_raw[k] = []

    #collect all the values in the range in output_list_raw
    for i in range(min_range,max_range):
        current_line = list_of_measurements[i]
        
        for key,value in current_line.items():
            output_list_raw[key].append(float(value))
            
    return output_list_raw


def average_over_range (min_range,max_range,list_of_measurements):
    
    #output_list_raw = input_range(min_range,max_range,list_of_measurements)
    
    #this dict (same categories as in list_of_measurements) will finally be output
    output_list_averaged = dict.fromkeys(list_of_measurements[0].keys())
    
    #Per key (time, frame #, etc), create an average 
    for key,j in output_list_raw.items():
        output_list_averaged[key] = np.average(j)
    
    #remove 's' time column, not necessary here 
    del output_list_averaged['s']
        
    #output the averages 
    return output_list_averaged

def find_maximum_over_range (min_range,max_range,list_of_measurements,list_of_baselines):
    output_list_raw = input_range(min_range,max_range,list_of_measurements)
    
    #remove 's' time column, not necessary here 
    del output_list_raw['s']
    
    #define the list (=dictionary) of minimal values, with the same Frames (=keys) as the original measurements file 
    output_list_maximum = dict.fromkeys(list_of_measurements[0].keys())
    #remove 's' time column, not necessary here 
    del output_list_maximum['s']
    
    # same, but for amplitudes (baseline - peak) 
    output_list_amplitudes = dict.fromkeys(list_of_measurements[0].keys())
    #remove 's' time column, not necessary here 
    del output_list_amplitudes['s']
    
    #Per key (time, frame #, etc), find the maximum 
    for key,j in output_list_raw.items():
        #calculate maximum
        peak_maximum = np.amax(j)
        
        #put into list 
        output_list_maximum[key] = peak_maximum
        
        #calculate amplitude and put into separate list 
        output_list_amplitudes[key] = peak_maximum - list_of_baselines[key]
           
    #output the averages 
    return output_list_maximum, output_list_amplitudes

def find_time_to_peak (min_range,max_range,list_of_measurements,peaks_list,current_stim_position):
    output_list_raw = input_range(min_range,max_range,list_of_measurements)
    
    #this dict (same categories as in list_of_measurements) will finally be output
    output_list_times = dict.fromkeys(list_of_measurements[0].keys())
    #remove 's' time column, not necessary here 
    del output_list_times['s']
    
    #Per key (time, frame #, etc), find the maximum 
    for key,j in peaks_list.items():
        index = np.where(output_list_raw[key] == j)
       
        timing_values = output_list_raw['s']
        
        index = index [0][0]
        
        peak_timing = timing_values[index]
        
        stim_to_peak_time = peak_timing - current_stim_position
        
        output_list_times[key] = stim_to_peak_time
            
    #output the averages 
    return output_list_times

def write_to_csv(output_file_path,output_list,i):
    # add time to output list 
    output_list['time'] = i
    
    #define list of Frames (=Field names)
    fieldnames = []
    for k,_ in output_list.items():
        fieldnames.append(k)
    
    #Write current list to csv file 
    with open(output_file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if not(os.path.exists(output_file_path)):
            writer.writeheader()
        writer.writerow (output_list)
    
'''

def average_over_range (min_range,max_range,list_of_measurements,df_to_add):
    #find avg in specific time frame 
    output_list_averaged = list_of_measurements.loc[min_range:max_range].mean(axis=0)
    
    #add to dataframe 
    df_to_add = df_to_add.append(output_list_averaged,ignore_index = True)
    
    return df_to_add
    
def find_maximum_over_range (min_range,max_range,list_of_measurements,list_of_baselines,amp_df, peak_df,norm_df):
    #find max in specific time frame  
    peak_list = list_of_measurements.loc[min_range:max_range].max(axis=0)
    
    #calculate amplitude by substracting from the last 
    amp_list = peak_list - list_of_baselines.iloc[[-1]]
    
    #calculate normalized Amplitude (normalized to 1st peak)
    if amp_df.empty:
        #first peak is always 1
        norm_list = amp_list / amp_list
    else: 
        #if there is more than one peak already, calculate normalized 
        norm_list = amp_list / amp_df.iloc[0]     
    
    #add to dataframes 
    amp_df = amp_df.append(amp_list, ignore_index=True)
    peak_df = peak_df.append(peak_list, ignore_index = True)
    norm_df = norm_df.append(norm_list, ignore_index = True)
    
    return amp_df, peak_df, norm_df

def find_time_to_peak (min_range,max_range,list_of_measurements,peaks_list,current_stim_position,peaks_time_df,seconds):
    
    #create truncated list 
    output_list_raw = list_of_measurements.loc[min_range:max_range]
    
    #calculate peak timing 
    peak_indices = np.where(output_list_raw == peaks_list.iloc[-1])[0]
    
    #add min_range to get to a proper index
    peak_indices = peak_indices + min_range 
    
    #find seconds time 
    peak_timing = seconds[peak_indices]
    
    #calculate time elapsed between stim and peak 
    stim_to_peak_time = peak_timing - current_stim_position
    
    #add to df 
    peaks_time_df = peaks_time_df.append(stim_to_peak_time.reset_index()[0], ignore_index = True)
    
    #output the averages 
    return peaks_time_df

def normalize_avg_stdev_2min(Normalized_df, Amplitudes_df):
    #Compare between peak and subsequent peak 
    norm_2min = pd.DataFrame(columns=Amplitudes_df.columns,index=[0])
    
    for i,peak in Amplitudes_df.iloc[0].iteritems():
        if i < 1: 
            norm_2min.iat[0,i] = np.NaN
        elif i <10:
            norm_2min.iat[0,i] = Amplitudes_df.iloc[-1][i+1] / peak 

    
    #Add to normalized list 
    Normalized_df = Normalized_df.append(norm_2min,ignore_index = True)
    
    #Add means column 
    Normalized_df ['avg'] = Normalized_df.mean(axis=1)
    
    #Add SEM column 
    Normalized_df ['sem'] = Normalized_df.sem(axis=1)
    
    return Normalized_df
    

#===Main program body===


def main (input_path, input_file_name,output_path=None ):
    
    #define your input file path here. 
    #input_path = None
    #input_file_name = None

    #default, place results files in the same folder. Can be adapted though.
    if output_path is None: 
        output_path = input_path
        mouse_name = '' # if files are saved in the same folder, no need to add super long names
    else:
        #separate out the folder's name 
        #(assuming one folders with mouse ID containing 10Hz... files )
        mouse_name =  os.path.split(os.path.split(input_path)[0])[-1]

    #input_text_file = 'C:/Users/wirrbel/Desktop/10Hz_MN_151019_002_WT_2mMCa2+.txt'
    input_text_file = os.path.join(input_path,input_file_name)

    #Prefixes for output file paths. 
    output_text_file_amplitudes = os.path.join(output_path,'Amplitudes_'+mouse_name+'_'+input_file_name)
    output_text_file_normalized = os.path.join(output_path,'Normalized_'+mouse_name+'_'+input_file_name)
    output_text_file_peak_time = os.path.join(output_path,'Peak_time_'+mouse_name+'_'+input_file_name)
    
    '''
    measurement = []
    
    #read all lines into an array
    with open(input_text_file) as f:    
        #    lines = f.read().splitlines() #old-school reading entire line into one array, keeps the tabs as '\t'
        lines = csv.DictReader(f, delimiter='\t')
        
        #to fully load file, take this: 
        for line in lines:
            measurement.append(dict(line))
     '''
    #read file into pandas dataframe 
    measurement =  pd.read_csv(os.path.join(input_path,input_file_name),delim_whitespace=True,skiprows=1,header=None)
    
    #create a separate 'seconds' column 
    seconds = measurement[0]
    
    #remove 'seconds' column, assumed to be the first column 
    measurement = measurement.drop(0,axis=1)
    
    #List of pulses, assuming 100 pulses every 100 ms, starting at 500 ms 
    Pulses_ms = list(range(500,10500,100))

    #add another pulse @ 30s 
    Pulses_ms.append(29800)

    #define empty DFs 
    Baselines = pd.DataFrame()
    Peaks = pd.DataFrame()
    Amplitudes = pd.DataFrame()
    Normalized = pd.DataFrame()
    Peaks_timing = pd.DataFrame()
    
    #fill DFs
    for i in Pulses_ms:
    
        #define time / line number (assuming 20 kHz sampling ) 
        Baseline_line = (i-1)*20 # Measure baseline 1 ms before stimulus 
        Peak_Range_min = (i+4)*20     # Start searching peak 4 ms after the stimulus (x 20 to get the lines_)
        Peak_range_max = (i+20)*20    #search peak until 20 ms after stimulus (x20 to get the line number )
    
        # Now we start collecting baseline values 
        Baselines = average_over_range ((Baseline_line-20),Baseline_line,measurement, Baselines)
    
        #Determine amplitudes and write to csv
        Amplitudes, Peaks, Normalized = find_maximum_over_range(Peak_Range_min,Peak_range_max,measurement,Baselines, Amplitudes, Peaks, Normalized)
        
        #Determine timing of peaks (in ms) and write to a separate CSV
        Peaks_timing = find_time_to_peak(Peak_Range_min,Peak_range_max,measurement,Peaks,i,Peaks_timing, seconds)
    
    #calculate averages + stdev + 2 min peak-to-peak
    Normalized = normalize_avg_stdev_2min(Normalized, Amplitudes)
   
    #Save all 
    Amplitudes.to_excel(output_text_file_amplitudes+'.xlsx')
    Normalized.to_excel(output_text_file_normalized+'.xlsx')
    Peaks_timing.to_excel(output_text_file_peak_time+'.xlsx')
    
    #alternative: Return average and normalized plus name
    return Normalized['avg'], Normalized['sem'],str(mouse_name+'_'+input_file_name)
    
    
#debug    
input_path = 'D:/Amygdala AMPA-GABA/2021-02-27_m2_Amygdala_Het_AMPA-GABA+GABA_Rundown/'
input_file_name = '10Hz_MN_280221_001_exp.txt'
main(input_path, input_file_name,output_path='D:/Amygdala AMPA-GABA/2021-02-27_m2_Amygdala_Het_AMPA-GABA+GABA_Rundown/')
        

        
 