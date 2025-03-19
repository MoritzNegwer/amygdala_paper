# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 23:14:50 2021

@author: Moritz Negwer
"""

import pandas as pd
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import multiprocessing as mp

#make sure that lmtanalysis can be imported 
os.chdir(r"C:\Users\Moritz Negwer\Documents\Work 2021\LMT Data\lmt-analysis-master\LMT")

from lmtanalysis.FileUtil import getFilesToProcess
from lmtanalysis.Animal import AnimalPool
from lmtanalysis.Measure import oneHour, oneMinute
import matplotlib.patches as mpatches
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from lmtanalysis.Measure import cornerCoordinates50x50Area, scaleFactor
from lmtanalysis.Point import Point


#from scripts/Compute_Measures_Identity_Profile_Dyadic
from lmtanalysis.Animal import *
from lmtanalysis.Event import *
from lmtanalysis.Measure import *
from lmtanalysis import BuildEventTrain3, BuildEventTrain4, BuildEventFollowZone, BuildEventRear5, BuildEventFloorSniffing,\
    BuildEventSocialApproach, BuildEventSocialEscape, BuildEventApproachContact,\
    BuildEventApproachRear, BuildEventGroup2, BuildEventGroup3, BuildEventGroup4,\
    BuildEventStop, BuildEventWaterPoint
    
from scripts.Rebuild_All_Event import *

#===== definitions ===== 

def modifiedGetDetectionTable(Animal):
    #modified from the AnimalPool class to also include front XYZ coordinates, to later be able to calculate turnarounds
    
    """
    Returns detections as pandas table for all animals in pool.

    * adds location also in cm
    * adds time column in pandas timedelta
    * adds column to indicate detections in center region

    Returns:
        pd.DataFrame: detections as pandas table
    """

    data = defaultdict(list)
    for animal in Animal.getAnimalList():
        for frame, detection in animal.detectionDictionnary.items():
            data["RFID"]         .append(f"{animal.name}_{animal.RFID}")
            data["name"]         .append(f"{animal.name}")
            data["genotype"]     .append(f"{animal.genotype}")
            data['frame']        .append(frame)
            data['sec']          .append(frame / oneSecond)
            data['x']            .append(detection.massX)
            data['y']            .append(detection.massY)
            data['z']            .append(detection.massZ)
            data['front_x']      .append(detection.frontX)
            data['front_y']      .append(detection.frontY)
            data['front_z']      .append(detection.frontZ)

    df = pd.DataFrame(data)


    df["x_cm"] = (df.x - cornerCoordinates50x50Area[0][0]) / (cornerCoordinates50x50Area[1][0] - cornerCoordinates50x50Area[0][0]) * ARENA_SIZE
    df["y_cm"] = (df.y - cornerCoordinates50x50Area[1][1]) / (cornerCoordinates50x50Area[2][1] - cornerCoordinates50x50Area[1][1]) * ARENA_SIZE

    df["front_x_cm"] = (df.front_x - cornerCoordinates50x50Area[0][0]) / (cornerCoordinates50x50Area[1][0] - cornerCoordinates50x50Area[0][0]) * ARENA_SIZE
    df["front_y_cm"] = (df.front_y - cornerCoordinates50x50Area[1][1]) / (cornerCoordinates50x50Area[2][1] - cornerCoordinates50x50Area[1][1]) * ARENA_SIZE
    
    return df.sort_values("frame").reset_index(drop=True)


def read_and_filter_database(database,folder,zone_offset,collection_df):
    tmin = folder['tmin'] 
    tmax = folder['tmax']
    
    #rebuild all events (only on first run)
    #process(file = os.path.join(folder['Folder'],file))
                
    #conversion to cm in cage coords 
    areax_min = ((folder['areax_min_px'] - cornerCoordinates50x50Area[0][0]) * scaleFactor) - zone_offset
    areax_max = ((folder['areax_max_px'] - cornerCoordinates50x50Area[0][0]) * scaleFactor) + zone_offset
    areay_min = ((folder['areay_min_px'] - cornerCoordinates50x50Area[0][1]) * scaleFactor) - zone_offset
    areay_max = ((folder['areay_max_px'] - cornerCoordinates50x50Area[0][1]) * scaleFactor) + zone_offset
    
    #define cage floor edges (as measured) into cm
    cagex_min = ((folder['cagex_min'] - cornerCoordinates50x50Area[0][0]) * scaleFactor) 
    cagex_max = ((folder['cagex_max'] - cornerCoordinates50x50Area[0][0]) * scaleFactor) 
    cagey_min = ((folder['cagey_min'] - cornerCoordinates50x50Area[0][1]) * scaleFactor) 
    cagey_max = ((folder['cagey_max'] - cornerCoordinates50x50Area[0][1]) * scaleFactor)
    
    #load into pool
    #load database
    connection = sqlite3.connect(os.path.join(folder['Folder'],database))
    
    # create an animalPool, which basically contains your animals
    animalPool = AnimalPool()
    
    # load infos about the animals
    animalPool.loadAnimals( connection )

    # load all detection (positions) of all animals for the first hour
    animalPool.loadDetection( start = folder['tmin'], end = folder['tmax'] )
    
    #extract detection table as pandas df
    allAnimals_pandas = modifiedGetDetectionTable(animalPool)
    
    #filter everything outside of the zone 
    #this will generate several animal readings, but only outside of the zone 
    for index,frame in allAnimals_pandas.iterrows():
        x = frame.x_cm
        y = frame.y_cm
        front_x = frame.front_x_cm
        front_y = frame.front_y_cm
        relative_frame = frame.frame-folder['tmin']
        
        #filter only detections inside the cage floor coords. Excludes e.g. wall jumps and faulty readings on the wall 
        if (cagex_min < x < cagex_max and cagey_min < y < cagey_max): # if the animal is within the arena 
            
            if not (areax_min < x < areax_max and areay_min < y < areay_max): 
                #if the animal is within the outside zone, add to df
                collection_df = collection_df.append({'frame': relative_frame,'x_cm':x, 'y_cm':y,'front_x':front_x,'front_y':front_y},ignore_index = True)
            
            
    # Remove all with duplicate frames (i.e. where more than one mouse was erroneously detected in the outside zone). 
    # This is a crude way of error-checking and will likely lead to slight underestimates in the outside_zone time, 
    # but it's the only one that is likely to work
    
    collection_df = collection_df.drop_duplicates(subset=['frame'],keep=False)
    
    return collection_df

def calculate_distance_vector_and_speed(collection_df,folder):
    #calculates the distance, vector (to center of minicage) and speed, each in a separate column 
    
    #calculate distance moved per frame 
    collection_df['distance_moved'] = np.hypot((collection_df.x_cm-collection_df.x_cm.shift(1)),(collection_df.y_cm-collection_df.y_cm.shift(1)))
    
    #center of mini-cage
    cage_center_x = (np.mean((folder['areax_min_px'],folder['areax_max_px'])) - cornerCoordinates50x50Area[0][0]) * scaleFactor
    cage_center_y = (np.mean((folder['areay_min_px'],folder['areay_max_px'])) - cornerCoordinates50x50Area[0][0]) * scaleFactor
    
    #calculate distance to center of mini-cage 
    collection_df['distance_cage'] = np.hypot((collection_df.x_cm-cage_center_x),(collection_df.y_cm-cage_center_y))
    
    #calculate orientation re: cage center in degrees for polar plot
    collection_df['orientation_cage'] = np.rad2deg(np.arctan2(*(((collection_df.y_cm-cage_center_y),(collection_df.x_cm-cage_center_x))[::-1])) % (2 * np.pi))
    #collection_df['orientation_cage'] = np.rad2deg(np.arctan2(*(((collection_df.x_cm-cage_center_x),(collection_df.y_cm-cage_center_y))[::-1])) % (2 * np.pi))
    
    #calculate mouse angle 
    collection_df['orientation_mouse'] = np.rad2deg(np.arctan2(*(((collection_df.front_y-collection_df.y_cm),(collection_df.front_x-collection_df.x_cm))[::-1])) % (2 * np.pi))
    
    '''
    #calculate turning (160-180 degrees within ca. 1s )
    outside['turning_deg'] = (outside['orientation_mouse']-outside['orientation_mouse'].shift(20)).abs()\
        .rolling(window=pd.api.indexers.FixedForwardWindowIndexer(window_size=10),min_periods=1).max() 
        # essentially, create a window of 10 elements and shift it back by 9 to re-align
        #then take max value 
    '''
    #filter for large changes in distance from cage (i.e. running towards or away) within ca 1s 
    collection_df['turning'] = collection_df[collection_df['distance_moved'] > 0.3]['distance_cage'].pct_change(30).\
        rolling(window=pd.api.indexers.FixedForwardWindowIndexer(window_size=60),min_periods=3).median().\
            apply(lambda x: True if  x < -0.2 or x > 0.2 else False)
    
    
    #return 
    return collection_df

def frame2index (outside_df):
    
    #make 18000 entry df 
    corrected_df = pd.DataFrame(np.empty((18000,0)))

    #cast framenumber as int, then set as index
    outside_df['frame'] = outside_df['frame'].astype(int)
    outside_df = outside_df.set_index('frame',drop=True)
            
    #merge both 
    corrected_df = pd.merge(corrected_df, outside_df, how='left',left_index=True,right_index=True)

    return corrected_df

def add_to_collection(collection_df,to_add_df,folder):
    #generate unique id, e.g. B1M3_d1_pt1
    measurement_id = (folder['mouse_id']+'_d'+str(folder['day'])+'_pt'+str(folder['pt']))
    
    #drop frame number from to_add_df (index will be equivalent to frame#)
    #to_add_df = to_add_df.drop('frame',1)
    
    #rename dataframe column heads to include measurement_id
    to_add_df = to_add_df.add_prefix(measurement_id+'_')
    
    #merge to existing collection 
    #collection_df = collection_df.merge(to_add_df,left_index=True,right_index=True)
    collection_df = pd.merge(collection_df, to_add_df, how='left',left_index=True,right_index = True)

    return collection_df


def calculate_binned_timing_distance(collection_df):
    #calcualates the distance travelled (or time in zone) per bin 
    
    #define binned dfs
    binned_count_df = pd.DataFrame(data=[],columns=collection_df.columns)
    binned_sum_df =  pd.DataFrame(data=[],columns=collection_df.columns)
    binned_average_df =  pd.DataFrame(data=[],columns=collection_df.columns)
    
    #divide 10 minutes into bins (default=10)
    bins = range(0,18001,1800)
    
    for n,m in enumerate(bins):
        if n>0:
            #print (bins[n-1],m)
            binned_count_df = binned_count_df.append(collection_df.iloc[bins[n-1]:m].count(axis=0),ignore_index=True)
            binned_sum_df = binned_sum_df.append(collection_df.iloc[bins[n-1]:m].sum(axis=0),ignore_index=True)
            binned_average_df = binned_average_df.append(collection_df.iloc[bins[n-1]:m].mean(axis=0),ignore_index=True)
    
    return binned_count_df, binned_sum_df, binned_average_df
        
def plot_1d_heatmap(data_to_plot,color_map='hot',vmin=0,vmax=700,x_label=None,add_colorbar=False,figure_filename='fig'):
    #plot all animals 
    #from Animal.AnimalPool.plotTrajectory 

    
    #determine rows
    nbRows = len(data_to_plot.columns)
    
    #specify plot distribution. Note the high-res (dpi=1200) and figsize to make the figure the right size (in inches)
    #current figsize settings (10, nbRows/2) work well for 18000 frame recordings, might need to adapt for binned recordings 
    fig, axes = plt.subplots( nrows = nbRows , ncols = 1 , sharex='all', sharey='all',  dpi=1200, figsize = (10, nbRows/2))
    axis = axes[0]
    
    #add horizontal spacing
    #adjust hspace to make the plots fit more snugly or less so 
    plt.subplots_adjust(bottom=0.1, top=0.9, hspace=2) 

    for mouse_no in data_to_plot.columns:
            #get binned trajectory
            binned_trajectory = data_to_plot[mouse_no]
                        
            #draw in the right position
            axis = axes[data_to_plot.columns.get_loc(mouse_no)]
            
            #generate heatmap. Change vmax manually to match values 
            img = axis.imshow(binned_trajectory[np.newaxis,:], cmap=color_map,aspect='auto',vmin=vmin,vmax=vmax)
            
            #reset y-ticks
            axis.set_yticks([])
            axis.set_aspect(1.0/axis.get_data_ratio()*0.05)
                                  
            #add the mouse ID on top of each graph
            axis.set_title(str(mouse_no),fontsize=8)
        
    #X label for the bins. Don't forget to adapt!
    axis.set_xticks(np.arange(0,len(data_to_plot)+1,len(data_to_plot)/10)) #settings for 10 min recordings
    axis.set_xticklabels(np.arange(0,11,1)) 
    
    #add colorbar 
    #Formatting optimized for 4 mice. 
    #Will probably need to be adjusted for <4 mice. 
    if add_colorbar == True:
        axins = inset_axes(axes[-1], # here using axis of the lowest plot
                   width="2%",  
                   height="145%",  
                   loc='lower left',
                   bbox_to_anchor=(1.025, -1, 1, 10),
                   bbox_transform=axes[-1].transAxes,
                   borderpad=0,
                   )
        fig.colorbar(img, cax=axins,ticks=[vmin,(vmax/2),vmax])
        
    #show (TODO: Implement saving option)
    #plt.show()
    plt.savefig(os.path.join(target_dir,figure_filename + '.svg'), format='svg')

###======= Main sequence starts here ===== 
if __name__ == "__main__":
    #read reference table for all batches
    #assumes there is a folder column (path to .sqlite file), x/y min/max from background in px, and mouse ID, day, part, genotype
    reftable_complete = pd.read_excel('C:/Users/Moritz Negwer/Documents/Work 2021/LMT Data/times_files_for_sniffing_v04_time_sorted_b6m3_pt2-2.xlsx')
    
    #debug
    #reftable = pd.read_excel('C:/Users/Moritz Negwer/Documents/Work 2021/LMT Data/times_tables_debug_v01.xlsx')
    
    #define genotype list 
    genotype_list = reftable_complete.genotype.unique().astype(str)
    genotype_list = genotype_list.tolist()
    #debug
    #genotype_list = ["WT"]
    
    #multiprocessing seems to be broken
    #pool = mp.Pool(processes = mp.cpu_count())        
    #genotype_list = [pool.apply_async(process_per_genotype,args=(genotype,))for genotype in genotype_list]

    
    
    #run once per genotype
    for genotype in genotype_list:
           #make one folder per genotype
        target_dir = os.path.join(r'C:\Users\Moritz Negwer\Documents\Work 2021\LMT Data\2023-11-26_turn_analysis' + '_' + genotype)
        try:
            os.mkdir(target_dir)
        except:
            pass #if the folder already exists, proceed (but overwrite contents)
        
        #subset reftable
        reftable = reftable_complete.loc[reftable_complete.genotype == genotype]
        
        #debug
        #target_dir = r'C:\Users\Moritz Negwer\Documents\Work 2021\LMT Data\2021-12-08_turn_troubleshooting'
        
        #define offset for zone size here (in cm) - the zone is mini-cage + this on every side
        zone_offset = 5
        
        #all collected works go here 
        t_outside_collection = pd.DataFrame(np.empty((18000,0)))
        
        
        #iterate through reftable, load files and calculate timings 
        for index,folder in reftable.iterrows():
            
            #check whether there is a database file 
            for file in os.listdir(folder['Folder']):
                if file.endswith('.sqlite'):
                    #if there is, process the database file 
                    #print("file found: ",file)
                    
                    #generate empty df 
                    outside = pd.DataFrame(columns=['frame','x_cm','y_cm','front_x','front_y'])
        
                    #read database and extract all detections that are outside of the mini-cage + offset zone 
                    #this assumes there are 2 animals (that may flip IDs at random intervals, as for Moritz's experiments)
                    outside = read_and_filter_database(file,folder,zone_offset,outside)
                    
                    #bring index to right frame
                    outside = frame2index(outside)
                    
                    #add distance and speed vectors
                    outside = calculate_distance_vector_and_speed(outside,folder)
                    
                    #add to collection 
                    t_outside_collection = add_to_collection(t_outside_collection,outside,folder)
                    
        
        #then, export to excel
        t_outside_collection.to_excel(os.path.join(target_dir,'t_outside_collection.xlsx'))
        
        #and to pickle for easy re-loading 
        t_outside_collection.to_pickle(os.path.join(target_dir,'t_outside_collection.pkl'))
        
        #calculate summed, average and means for each zone 
        outside_binned_count, outside_binned_sum, outside_binned_average = calculate_binned_timing_distance(t_outside_collection)
        
        #make selection for easier plotting
        select_outside_binned_sum_distance_moved = outside_binned_sum.filter(regex='distance_moved')
        #select_outside_binned_count = outside_binned_count.filter(regex='x_cm') #otherwise it will miss all the first rows of a movement block (those have x/y coords, but no movement distance)
        
        #choose the following regex when using the simple drop_duplicates method at the end of read_and_filter_database. 
        #That method produces approx 30-90 frames of loose single entries that are not useful tracking. 
        #Thus, using distance_moved regex eliminates all single trackings plus a bit more. 
        select_outside_binned_count = outside_binned_count.filter(regex='distance_moved') 
        
        #calculate number of turning events
        '''
        #extracct degrees 
        turning_deg = t_outside_collection.filter(regex='turning_deg')
        
        #calculate when events happen 
        #turning_bool = turning_deg.apply(lambda x: True if 160 < x < 200 else False,axis=0)
        turning_bool = pd.DataFrame(index = turning_deg.index)
        for column in turning_deg:
            #mark the turns with True / False
            iseventhappening = turning_deg[column].apply(lambda x: True if 140 < x < 220 else False)
            turning_bool[column] = iseventhappening
        '''
        #make turning heatmaps 
        turning_bool = t_outside_collection.filter(regex='turning')
        #turning_bool_plot = turning_bool.astype(int)
        #low-pass filter: only turns longer than 10 frames will be counted 
        turning_bool= turning_bool.rolling(window=pd.api.indexers.FixedForwardWindowIndexer(window_size=10),min_periods=3).min()
        #replace the 0s with NaN
        turning_bool.replace(0, np.nan, inplace=True)
        turning_bool.replace(True, 1, inplace=True)
        #make 1D heatmaps
        plot_1d_heatmap(turning_bool,color_map='Blues',vmax=1,figure_filename='turn_bool_int_plot')
        
        #count events
        turning_bool_count = turning_bool.replace(np.nan, 0)
        turning_bool_count = turning_bool_count.astype('bool')
        #mask by its own shifted self = count changes from True -> False 
        turning_bool_count = (turning_bool_count & (turning_bool_count != turning_bool_count.shift(1))).sum()
        #Cast as int, otherwise it's a float 
        turning_bool_count.astype(int) 
        #save as excel 
        turning_bool_count.to_excel(os.path.join(target_dir,'turn_counts' + '.xlsx'))
        
        #save to excel + pickle
        select_outside_binned_sum_distance_moved.to_excel(os.path.join(target_dir,'select_outside_binned_sum_distance_moved_' + '.xlsx'))
        select_outside_binned_count.to_excel(os.path.join(target_dir,'select_outside_binned_time' + '.xlsx'))
        
        #make separate cumulative sum table 
        cumulative_time = select_outside_binned_count.cumsum()
        cumulative_time.to_excel(os.path.join(target_dir,'select_outside_binned_time_cumulative' + '.xlsx'))
        
        #generate overview plots 
        plot_1d_heatmap(t_outside_collection.filter(regex='distance_moved'),vmax=2,figure_filename='distance_moved_all_mice_v01_vmax_2')
        
        
        ###=== per mouse list of binned results ===
        #create single_mouse subdirectory
        try:
            os.mkdir(os.path.join(target_dir,'single_mice'))
        except: 
            pass
        
        #make mouse list 
        mouse_list = list(reftable['mouse_id'])
        
        #make table for each mouse 
        for mouse in mouse_list:
            #make per-mouse selection for outside_binned_count
            mouse_selection = select_outside_binned_count.filter(regex=mouse)
            mouse_selection.to_excel(os.path.join(target_dir,'single_mice',mouse + '_outside_time' + '.xlsx'))
            
            #make per-mouse selection for sum_distance 
            mouse_selection = select_outside_binned_sum_distance_moved.filter(regex=mouse)
            mouse_selection.to_excel(os.path.join(target_dir,'single_mice',mouse + '_outside_distance' + '.xlsx'))
        
    
        
        