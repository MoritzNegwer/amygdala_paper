#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 14:04:04 2022

@author: wirrbel
"""

### === single areas, uncorrected t-tests === ### 
### List from excel, with fiber tracts removed for visibility ###

#import sys
#sys.path.append('/home/wirrbel/brainrender_update_2022-01/')
import os
import pandas as pd
import numpy as np
from shutil import copyfile
import in_place

#from vedo import embedWindow  # for more explanations about these two lines checks the notebooks workflow example
#embedWindow(None)

# Import variables
from brainrender import * # <- these can be changed to personalize the look of your renders

# Import brainrender classes and useful functions
from brainrender import Scene
from brainrender.actors import Points, PointsDensity, Volume
from brainrender import settings 

import myterial

from matplotlib import colormaps

brainrender.settings.SHOW_AXES = False
brainrender.settings.SCREENSHOT_SCALE = 2
brainrender.settings.SHADER_STYLE = "glossy" # [cartoon, metallic, plastic, shiny, glossy]

def convert_colorvalues(region_cmap,value):
    #extract a list of RGBA values from the colormap 
    colorvalues = tuple(region_cmap(float(value),bytes=True))
    #remove alpha value 
    colorvalues = colorvalues[0:-1]
    #convert to hex 
    hex_colorvalues = '#%02x%02x%02x' % colorvalues 
    
    return hex_colorvalues 
    
    


def render_scene(title,folder,region_list,single_hemisphere=False,camera=None,colormap=None,split_zscores=False):
    scene = Scene (title=None,screenshots_folder=folder,inset=None)
    
    if colormap is None:
        #if no colormap is supplied, assume that region_list has the structure {region:alpha}
        #regions will be colored by their ABA color scheme 
        for idx,row in region_list.iterrows():
            #separate out the region ID and intensity value
            region = int(row["id"])
            intensity = row["Het log2 FC for seismic colorscale"]
            #add to the scene, color-coded for percentage differences 
            scene.add_brain_region(region,alpha=intensity)
        
    else:
        #assume tiat the region_list has the structure {region:intensity} 
        #regions will be color-coded according to the colormap
        #get a colormap to color-code the areas with
        #see here for an overview: https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html 
        region_cmap = colormaps.get_cmap(colormap)
        
        #offset for seismic. 
        #Colormap scales from 0-1 with a center at 0.5, and the log2-FC centers at 0.     
        #This compresses the values (so that the scale reaches from -1 to +1) and shifts the center to 0.5
        if colormap == "seismic":
            region_list["Het log2 FC for seismic colorscale"] = region_list["Het log2 FC for seismic colorscale"]/2
            region_list["Het log2 FC for seismic colorscale"] = region_list["Het log2 FC for seismic colorscale"]+0.5

        
        #assuming that region_list is an item with the structure {region:intensity}
        for idx,row in region_list.iterrows():
            #separate out the region ID and intensity value
            region = int(row["id"])
            intensity = row["Het log2 FC for seismic colorscale"]
            #get intensity 
            intensity_color = convert_colorvalues(region_cmap,intensity)
            #add to the scene, color-coded for percentage differences 
            scene.add_brain_region(region,color=intensity_color,alpha=0.5)
        
    
    if single_hemisphere:
        plane = scene.atlas.get_plane(plane="sagittal", norm=[0, 0, -1])
        #set sagittal cutting plane ever so slightly off the center 
        #to prevent a render error with the "cartoon" render setting (see settings.py)
        #plane.mesh.SetPosition(7829.739797878185, 4296.026746612369, -5700.497969379508)
        scene.slice(plane,close_actors=True)
    
    
    
    #if region has a special camera angle attached to it, use this one. Othewise, use default (see brainrender settings)
    if camera is not None:
        scene.render(camera=camera, interactive=False,offscreen=True)
        scene.screenshot(name=title)
        scene.close()
    else:
        scene.render()
    
    #scene.render()
    
    ## manually save screenshot, then close ### 

###quasi-sagittal, I think 
techpaper_cam_03 = {
     "pos": (-14632, -20614, -35939),
     "viewup": (0, -1, 0),
     "clipping_range": (27940, 61510),
     "focalPoint": (6888, 3571, -5717),
     "distance": 44288
     }

#sagittal
dorsal_cam = {
    "pos": (4672, -40653, -6604),
    "viewup": (-1, 0, 0),
    "clipping_range": (35913, 55415),
    "focalPoint": (6888, 3571, -5717),
    "distance": 44288
    }

sagittal_cam = {
    "pos": (8512, 5167, -49946),
    "viewup": (0, -1, 0),
    "clipping_range": (35391, 49773),
    "focalPoint": (6888, 3571, -5717),
    "distance": 44288
    }

coronal_cam = {
    "pos": (-37325, 988, -5786),
    "viewup": (0, -1, 0),
    "clipping_range": (30001, 61730),
    "focalPoint": (6888, 3571, -5717),
    "distance": 44288
    }

#load plaque counts 
counts_table = pd.read_csv('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/2022-05-29_overview_Ehmt1_figs_v02/SST_significant_areas_v08_log2_FC.csv')
SST = counts_table[["id","Het log2 FC for seismic colorscale"]]

'''#extract WT into dict
SST = {}

for region in counts_table.iterrows():
    abbreviation = region[1]['id']
    percentage = region[1]['Het log2 FC for seismic colorscale']
    SST[abbreviation] = percentage
'''


#load plaque counts 
counts_table = pd.read_csv('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/2022-05-29_overview_Ehmt1_figs_v02/PV_P14_significant_v05_log2_FC.csv')
PV_P14 = counts_table[["id","Het log2 FC for seismic colorscale"]]

'''
#extract WT into dict
PV_P14 = {}

for region in counts_table.iterrows():
    abbreviation = region[1]['id']
    percentage = region[1]['Het log2 FC for seismic colorscale']
    PV_P14[abbreviation] = percentage
'''

#load plaque counts 
counts_table = pd.read_csv('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/2022-05-29_overview_Ehmt1_figs_v02/PV_P56_significant_v05_log2_FC.csv')
PV_P56 = counts_table[["id","Het log2 FC for seismic colorscale"]]

'''
#extract WT into dict
PV_P56 = {}

for region in counts_table.iterrows():
    abbreviation = region[1]['id']
    percentage = region[1]['Het log2 FC for seismic colorscale']
    PV_P56[abbreviation] = percentage
'''

#load plaque counts 
counts_table = pd.read_csv('/home/wirrbel/2021-03-29_brainrender_2_preprocessing/2022-05-29_overview_Ehmt1_figs_v02/VIP_significant_prelim_v02_log2_FC.csv')
VIP = counts_table[["id","Het log2 FC for seismic colorscale"]]

'''
#extract WT into dict
VIP = {}

for region in counts_table.iterrows():
    abbreviation = region[1]['id']
    percentage = region[1]['Het log2 FC for seismic colorscale']
    VIP[abbreviation] = percentage
'''

#result from general uncorrected t-tests (minus fiber tracts)
#list_of_regions_and_titles = [['uncorrected_t_test_c26_vs_nc26.png',['RSPv6a','PL5','SSp-n6a','VISp6a','VISal6a','ORBm5','VISl6a','ORBvl6a','TRN','VISpm2/3','VISam4','SF','BLAp','ICB','ORBvl2/3','VISam6a','AVP','VISpm5','VISpm6a','AUDp1','RSPv5','AUDv1','SSp-ll5','DN','LHA','ORBvl1','ORBvl5','VISam5','ECT6a','TEa6a','SLD','GRN','PL2/3','IPN','APr','VISam2/3','SSp-tr6a','ILA5','SPVC','ACAd5','IF','RSPd6a','SSs6a','SPA','SSp-n5','PMv','MOs6a','VISpl6a','VISpor6a','VISli6a','PRNr','VISpm4','och','SSp-bfd6a','MGd','SMT','ACAd6a','VISa6a','ACAv5']]]

target_folder = '/home/wirrbel/2021-03-29_brainrender_2_preprocessing/2022-05-29_overview_Ehmt1_figs_v02/'

os.makedirs(target_folder,exist_ok=True)

### ==== areas color-coded by intensity ====


#colormap_name = 'Reds'

#for level in [[SST,'WT'],[plaque_count_6_mon,'Het']]:
    
level_list = [[SST,'SST_Het FC','seismic'],
              [PV_P14,'PV_P14_Het FC','seismic'],
              [PV_P56,'PV_P56_Het FC','seismic'],
              [VIP,'VIP_Het FC','seismic'],
              #[SST,'SST_Het FC','Greens'],
              #[PV_P14,'PV_P14_Het FC','Oranges'],
              #[PV_P56,'PV_P56_Het FC','Reds'],
              #[VIP,'VIP_Het FC','Blues']
              ]

for level in level_list:
    prefix = 'significant_log2_FC_'
    colormap_name = level[2]
    
    #both higher + lower 
    title = str(prefix + colormap_name + '_hi_low_' + level[1] )
    region = level[0].copy()
    render_scene(title, target_folder, region, single_hemisphere = True, camera=techpaper_cam_03,colormap=colormap_name,split_zscores=True)
    
    #only high
    title = str(prefix + colormap_name + '_only_high_' + level[1] )
    region = level[0].copy()
    region = region[region["Het log2 FC for seismic colorscale"]>0]
    render_scene(title, target_folder, region, single_hemisphere = True, camera=techpaper_cam_03,colormap=colormap_name,split_zscores=True)
    
    #only low
    title = str(prefix + colormap_name + '_only_low_' + level[1] )
    region = level[0].copy()
    region = region[region["Het log2 FC for seismic colorscale"]<0]
    render_scene(title, target_folder, region, single_hemisphere = True, camera=techpaper_cam_03,colormap=colormap_name,split_zscores=True)
    

'''
### ==== area transparency adapted by percentage ====
prefix = '72_areas_fraction_single_hemisphere_transparency_mapped_angled_color_'

for level in [[SST,'SST'],[plaque_count_6_mon,'plaque_count_6_mon']]:
    title = str(prefix + level[1] )
    region = level[0]
    render_scene(title, target_folder, region, single_hemisphere = True, camera=techpaper_cam_03,colormap=None)
'''