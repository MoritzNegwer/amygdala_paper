#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 20:56:56 2024

@author: wirrbel
"""
import numpy as np
import os
import glob
import pandas
import matplotlib.pyplot as plt
from matplotlib import colormaps 


#define input + output folder (make sure to omit trailing slash!)
input_folder_list = [
    "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-02-26_P12_Amygdala_Het_mIPSCs+AMPA-GABA+GABA_Rundown",
    "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-02-27_m1_Amygdala_WT_AMPA-GABA+GABA_Rundown",
    "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-02-27_m2_Amygdala_Het_AMPA-GABA+GABA_Rundown",
    "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-02-28_m1_Amygdala_WT_AMPA-GABA+GABA_Rundown",
    "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-02-28_m2_Amygdala_Het_AMPA-GABA+GABA_Rundown",
    "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-03-01_m1_Amygdala_WT_AMPA-GABA+GABA_Rundown",
    ]

#input_folder = "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2021-02-27_m2_Amygdala_Het_AMPA-GABA+GABA_Rundown"
output_folder = "/media/wirrbel/MN_4/Amygdala_paper_collection/Ephys_Amygdala/2024-02-18_10Hz_overview_plots/"

for input_folder in input_folder_list:
    #file parameters
    parent_dir = os.path.split(input_folder)[1]
    input_file_list = sorted(glob.glob(os.path.join(input_folder,"*10Hz*.txt")))

    #plotting parameters
    linewidth = 0.1
    #define colors 
    if "WT" in parent_dir:
        colormap = colormaps["gray"]
        colormap = colormap(np.linspace(0.0,1.0,256))
        colormap = np.array((*colormap[:,0],*colormap[:,1],*colormap[:,2]))
        colormap = np.reshape(colormap,(256,3),order="F")
    elif "Het" in parent_dir:
        colormap = colormaps["Reds"]
        colormap = colormap(np.linspace(0.0,1.0,256))
        colormap = np.array((*colormap[:,0],*np.linspace(0.0,0.0,256),*np.linspace(0.0,0.0,256)))
        colormap = np.reshape(colormap,(256,3),order="F")
    else:
        colormap = colormaps["inferno"].colors

    
    for input_file in input_file_list:
        #read file
        txt_file = pandas.read_csv(os.path.join(input_folder,input_file),sep="\t")
        
        #split input_filename
        input_filename = os.path.split(input_file)[1]
        #define output name
        output_filename = parent_dir+input_filename[:-3]+"_10Hz_graph_v01"
        
        '''
        #10Hz plot
        #plot + save
        plt.figure().set_figwidth(20)
        plt.plot(txt_file["s"],txt_file["Frame 1"],color=colormap[20],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 2"],color=colormap[40],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 3"],color=colormap[60],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 4"],color=colormap[80],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 5"],color=colormap[100],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 6"],color=colormap[120],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 7"],color=colormap[140],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 8"],color=colormap[160],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 9"],color=colormap[180],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 10"],color=colormap[200],linewidth=linewidth)
        plt.xlim(0,11000)
        plt.ylim(0,1000)
        #plt.show()
        plt.savefig(os.path.join(output_folder,"svg",output_filename+".svg"))
        plt.savefig(os.path.join(output_folder,"png",output_filename+".png"),dpi=600)
        '''
        
        #30s spike
        plt.figure().set_figwidth(1)
        plt.plot(txt_file["s"],txt_file["Frame 1"],color=colormap[20],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 2"],color=colormap[40],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 3"],color=colormap[60],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 4"],color=colormap[80],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 5"],color=colormap[100],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 6"],color=colormap[120],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 7"],color=colormap[140],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 8"],color=colormap[160],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 9"],color=colormap[180],linewidth=linewidth)
        plt.plot(txt_file["s"],txt_file["Frame 10"],color=colormap[200],linewidth=linewidth)
        plt.xlim(29700,30250)
        plt.ylim(0,1000)
        plt.savefig(os.path.join(output_folder,"svg",output_filename+"_recovery.svg"))
        plt.savefig(os.path.join(output_folder,"png",output_filename+"_recovery.png"),dpi=600)
        
        