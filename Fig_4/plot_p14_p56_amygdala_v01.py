#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 12:46:10 2025

@author: wirrbel
"""
import pandas as pd
import numpy as np
import os, glob
import tifffile
import cv2
from scipy.stats import ttest_ind
#define a random number generator to be used in permutation t-tests
rng = np.random.default_rng()
#first, load the files from the Allen SDK
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import matplotlib.pyplot as plt
import seaborn as sns

p14 = pd.read_excel("/media/wirrbel/MN_4/Amygdala_paper_collection/2020-08-18_P14_Amygdala/P14_sheet_v01.xlsx",sheet_name="Sheet1")

regions = p14["Region"].unique()
mice = p14["Mouse ID"].unique()
p14_convert = pd.DataFrame(data=[],index=pd.Index([*regions,"gt"]),columns=pd.Index(mice))
for mouse in mice:
    #extract mouse from general df 
    current_mouse = p14[p14["Mouse ID"] == mouse]
    p14_convert.at["gt",mouse] = current_mouse["Genotype"].unique()[0]
    #add new column to convert dataframe
    for region in regions:
        p14_convert.at[region,mouse] = current_mouse[current_mouse["Region"] == region]["Cells/mm²"].sum()

reduced_df["gt"]
reduced_df[["gt"]]
reduced_df.loc["gt"]
output_dir = "/media/wirrbel/MN_4/Amygdala_paper_collection/2025-03-14_replot_P14_P56_amygdala/"
output_name = "P14_Amygdala_PV"
reduced_df = p14_convert.T
#reattach to cleaned-up dataframe
reduced_df["gt"] = reduced_df["gt"].str.replace("WT","Ehmt1++")
reduced_df = reduced_df.sort_values("gt")

# transfer to "long form" for plotting https://stackoverflow.com/questions/65975715/how-to-add-a-hue-of-one-categorical-column-to-seaborn-barplot-of-each-column-me
longform_df = pd.melt(reduced_df,id_vars="gt")

#create plot 
fig,ax = plt.subplots()
ax = sns.violinplot(data=longform_df,x="variable",y="value",hue="gt",split=True,inner=None,density_norm="count",palette=["gray","red"])
ax.figure.set_size_inches(len(reduced_df.T)*0.5,5)
ax.set_ylim(0)
ax.tick_params(axis='x', labelrotation=45)
plt.savefig(os.path.join(output_dir,output_name+".svg"))