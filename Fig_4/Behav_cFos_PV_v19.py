#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 01:48:33 2022

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


def generate_collection(filelist):
    
    #define new collection 
    collection_df = pd.DataFrame()
    
    #extract list
    for csv in filelist:
        #define mouse name from file name 
        mouse_name = os.path.split(csv)[1][:10]
        
        try:
            #extract list only if it says median
            current_list = pd.read_csv(csv,usecols=['Median'])
            #reanme to mouse_name
            current_list = current_list.rename(columns={'Median': mouse_name})
            
        
            #find unique occurrences 
            current_areas = current_list.value_counts().sort_index()
            #add mouse name
            current_areas = current_areas.rename(mouse_name)
            
            #add to collection df 
            collection_df = pd.concat([collection_df,current_areas],axis=1)
            
        except:
            #if the read_csv does not work (i.e. the column is not named 'Median', then report and move on)
            print ("did not find a column called \'Median\' in the file: " + mouse_name + ". Please check whether it needs to be adapted")
            
    return collection_df

'''
#not needed anymore
def count_occurrences (collection_df):
    #generate new df 
    occurrences = pd.DataFrame()

    #count how often which area occurs 
    for column in collection_df:
        occurrences[column] = collection_df[column].value_counts()
    
    return occurrences
'''

def values_to_acronyms(collection_df):
    #loads values from AllenSDK 
    
    #See here for example: https://allensdk.readthedocs.io/en/latest/_static/examples/nb/reference_space.html#Constructing-a-structure-tree 
    
    reference_space_key = 'annotation/ccf_2017'
    resolution = 25
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    # ID 1 is the adult mouse structure graph
    tree = rspc.get_structure_tree(structure_graph_id=1)
    
    #from the example
    #tree.get_structures_by_id([962])
    
    try:
        #first, make an interpretable list out of the strange Pandas multi-index 
        collection_df_list = collection_df.index.tolist()
        collection_df_list_flat = [i for (i,) in collection_df_list]
    except:
        #if the df has already a normal list as index (i.e. pixel_count)
        collection_df_list_flat = collection_df.index.to_list()

    #get all the Allen SDK into 
    Atlas_info = tree.get_structures_by_id(collection_df_list_flat)  
    
    #extract only the acronyms
    #Note: Replace any value that does not occur in the reference atlas with "None" 
    #This is for example due to a weird quirk of FIJI's "Median" function when a ROI encompasses two areas, 
    #but has exactly the same pixel count in both. 
    acronyms_list = ['unknown' if d is None else d['acronym'] for d in Atlas_info]

    #set the acronym list as index 
    collection_df = collection_df.set_index(pd.Series(acronyms_list))
    
    #clean up all 'unknown' counts 
    try:
        collection_df = collection_df.drop('unknown')
    except:
        pass
    
    return collection_df

    



#crawl through folder, check whether tiffs are nearest-neighbor interpolated 
#(if they are not, they have weird intermediate values and this will give wrong results. 
#Solution: Re-annotate landmarks in BigWarp, then re-make the image, then re-extract point values. 
def check_atlas_files_for_nearest_neighbor_interpolation(source_dir):
    #make list of _atlas tiff files in source_dir 
    atlas_img_list = sorted(glob.glob(str(source_dir + '*_atlas.tif')))
    
    for image_location in atlas_img_list:
        #load image 
        img_float = tifffile.imread(image_location)
        
        #test whether the image has been rendered with nearest-neighbor interpolation 
        #(i.e. there are only whole numbers in the image) 
        img_int = img_float.astype(int)
        if np.all((img_float - img_int)==0):
            #print('all ok')
            interpolation_correct = True
        else:
            print('Image likely not interpolated correctly, and needs to be redone: ' \
                  + image_location )
            interpolation_correct = False
            

def mask_atlas_file_and_count_px(source_dir):
    #make new pandas df 
    area_counts = pd.DataFrame()
    
    #Go through folder, make a list of atlas_tif files 
    atlas_img_list = sorted(glob.glob(str(source_dir + '*_atlas.tif')))
    
    for image_location in atlas_img_list:
        #debug
        #print("processing: ",image_location)
        
        #load atlas image 
        atlas_img = tifffile.imread(image_location)
        
        #define image name
        img_name = os.path.split(image_location)[1].replace('_Probabilities_atlas.tif','')
        
        #extract the name of the original tif 
        corresponding_original_location = os.path.join(source_dir,img_name + '.tif')
        corresponding_original = tifffile.imread(corresponding_original_location)
        
        #create mask by thresholding the first channel of the original image above 600 (=outline of the slice)
        mask = cv2.threshold(corresponding_original[0,:,:],400,255,cv2.THRESH_BINARY)[1].astype('uint8')
        
        #debug
        #import matplotlib.pyplot as plt
        #plt.imshow(mask)
        
        #workaround for weird image issues 
        try:
            #mask the atlas image with the mask (=slice outline) 
            #to make sure that cut-off slice bits are not counted as empty atlas areas 
            thresholded = cv2.bitwise_and(atlas_img,atlas_img,mask=mask)
            
        except:
            #if it fails, report the image and don't count the contents
            print("failed masking, using non-masked instead: ", img_name)
            thresholded = atlas_img
        #optional: save as thresholded_atlas.tif
        #save_filename = os.path.join(source_dir,img_name + '_Probabilities_atlas_thresholded.tif')
        #tifffile.imsave(save_filename,thresholded)
        
        #count the areas (px value = area ID) and pixel number per area 
        unique_areas, area_px_counts = np.unique(thresholded,return_counts=True)
        #combine in pandas series 
        current_slice_count = pd.Series(data=area_px_counts,index=unique_areas.astype('int64'),name=img_name[:10])
        
        area_counts = pd.concat([area_counts,current_slice_count],axis=1)
        
    #at the end, return complete area_count
    return area_counts
    

def assemble_overview_table (source_dir, target_dir, celltype):
    #cFos data is in this folder 
    filelist = sorted(glob.glob(str(source_dir + "*" + celltype + "_ROIs.csv")))
    
    #extract collection of all cFos data 
    collection_df = generate_collection(filelist)
    
    #fix naming 
    collection_df = values_to_acronyms(collection_df)
    
    #sort index in alphabetical order 
    collection_df = collection_df.sort_index()
    
    #optional: Make pivot table
    #occurrences = occurrences.T
    #collection_df = collection_df.T
    
    #optional: save to excel file 
    #collection_df.to_excel(os.path.join(target_dir, celltype + "_occurrences_v04_pivot.xlsx"))
    
    return collection_df

def subset_and_save(cells_df,areas_list,target_dir,prefix='please define a name',collection_df=None,pivot=False):
    #creates a subset of the dataframe and saves as prefix + excel in target_dir
    
    #also preserve batch/mouse/slice numbers + genotypes 
    areas_list = areas_list+['batch_no','mouse_no','slice_no','genotype','time_outside_s','distance_moved_cm','turn_counts']
    
    #create subset 
    subset_df = cells_df[cells_df.index.isin(areas_list)]
    
    #save subset_df as excel in target_dir
    if pivot:
        subset_df.T.to_excel(os.path.join(target_dir,prefix+'.xlsx'))
    else: 
        subset_df.to_excel(os.path.join(target_dir,prefix+'.xlsx'))
        
    if collection_df is not None:
        collection_df = pd.merge(collection_df,cells_df,how='left',left_index=True,right_index=True)
    #return 
    return collection_df



def assign_genotypes(cells_df, genotype_list_file):
    
    #get list of mouse IDs and genotypes
    genotype_list = pd.read_excel(genotype_list_file,usecols=['mouse_id',
                                                              'genotype',
                                                              'time_outside_s',
                                                              'distance_moved_cm',
                                                              'turn_counts'])
    
    #extract index (test, this gets only first entry)
    #slice_id = cells_df.index[0]
    
    #transpose cells_df to make work with the code below
    cells_df = cells_df.T
    
    #extract batch, mouse and slice # from slice_id 
    ids = cells_df.index.str.extract(r'B(\d)_M(\d)_Sl(\d)').copy()
    ids = ids.rename(columns={0 : 'batch_no', 1 : 'mouse_no', 2 : 'slice_no'})
    ids = ids.set_index(cells_df.index)
    #merge back into cells_df
    cells_df = pd.merge(left=cells_df,right=ids,left_index=True,right_index=True)
    
    #extract batch + mouse ID from genotype_list 
    #(note this will only work on Ehmt1, as the pattern matches EB and nothing else)
    genotype_list = genotype_list[0:24]
    gt_ids = genotype_list['mouse_id'].str.extract('EB(\d)M(\d)')
    gt_ids = gt_ids.rename(columns={0 : 'batch_no', 1 : 'mouse_no'})
    genotype_list = pd.merge(left=genotype_list,right=gt_ids,left_index=True,right_index=True)
    
    #add extra genotype column to cells_df
    cells_df['genotype'] = ''
    cells_df['time_outside_s'] = ''
    cells_df['distance_moved_cm'] = ''
    cells_df['turn_counts'] = ''
    
    #walk through genotype_list and match to cells_df 
    for idx,mouse in genotype_list.iterrows():
        cells_df.loc[(cells_df['batch_no'] == mouse['batch_no']) & (cells_df['mouse_no'] == mouse['mouse_no']),'genotype'] = mouse['genotype']
        cells_df.loc[(cells_df['batch_no'] == mouse['batch_no']) & (cells_df['mouse_no'] == mouse['mouse_no']),'time_outside_s'] = mouse['time_outside_s']
        cells_df.loc[(cells_df['batch_no'] == mouse['batch_no']) & (cells_df['mouse_no'] == mouse['mouse_no']),'distance_moved_cm'] = mouse['distance_moved_cm']
        cells_df.loc[(cells_df['batch_no'] == mouse['batch_no']) & (cells_df['mouse_no'] == mouse['mouse_no']),'turn_counts'] = mouse['turn_counts']
    
    #add unique mouse ID (assuming this was done in blocks of 4)
    cells_df["mouse_id"] = (cells_df["batch_no"].astype(int)-1)*4+(cells_df["mouse_no"].astype(int))
    
    #transpose cells_df back 
    cells_df = cells_df.T
    
    return cells_df 

def reappend_info(cells_df,template_df):
    #re-assigns batch_no, mouse_no, etc
    
    #flip both dfs 
    cells_df = cells_df.T
    template_df = template_df.T
    #template_df = template_dft.reset_index() #to make sure mouse_id is present even if otherwise used as index 
    
    #re-assign genotypes, batch + mouse numbers 
    cells_df[['batch_no','mouse_no','slice_no','genotype','time_outside_s','distance_moved_cm','turn_counts']] = ''
    
    for mouse in cells_df.index:
        try:
            cells_df.at[mouse,'batch_no'] = template_df.at[mouse,'batch_no'].values[0]
        except:
            #print(template_df.at[mouse,'batch_no'])
            cells_df.at[mouse,'batch_no'] = template_df.at[mouse,'batch_no']
        try:
            cells_df.at[mouse,'mouse_no'] = template_df.at[mouse,'mouse_no'].values[0]
        except:
            cells_df.at[mouse,'mouse_no'] = template_df.at[mouse,'mouse_no']
        
        cells_df.at[mouse,'slice_no'] = 'all'
        
        try:
            cells_df.at[mouse,'genotype'] = template_df.at[mouse,'genotype'].values[0]
        except:
            cells_df.at[mouse,'genotype'] = template_df.at[mouse,'genotype']
        try:
            cells_df.at[mouse,'time_outside_s'] = template_df.at[mouse,'time_outside_s'].values[0]
        except:
            cells_df.at[mouse,'time_outside_s'] = template_df.at[mouse,'time_outside_s']
        try:
            cells_df.at[mouse,'distance_moved_cm'] = template_df.at[mouse,'distance_moved_cm'].values[0]
        except:
            cells_df.at[mouse,'distance_moved_cm'] = template_df.at[mouse,'distance_moved_cm']
        try:
            cells_df.at[mouse,'turn_counts'] = template_df.at[mouse,'turn_counts'].values[0]
        except:
            cells_df.at[mouse,'turn_counts'] = template_df.at[mouse,'turn_counts']


    #flip back
    cells_df = cells_df.T
    
    return cells_df
    

def summarize_per_mouse(cells_df):
    
    #flip dataframe 
    cells_dft = cells_df.T
    
    #create a copy for safekeeping
    cells_dftcopy = cells_dft.copy()
    cells_dftcopy = cells_dftcopy.set_index('mouse_id')
    
    #summarize per mouse 
    cells_dft_sum = cells_dft.set_index('mouse_id').loc[:, cells_dft.columns [0:-8]].groupby('mouse_id').sum()
    cells_dft_count = cells_dft.set_index('mouse_id').loc[:, cells_dft.columns [0:-8]].groupby('mouse_id').count().replace(0,np.nan)
    
    #normalize for slice counts 
    cells_dft = cells_dft_sum / cells_dft_count
        
    #reassign infos 
    cells_dft = reappend_info(cells_dft.T, cells_dftcopy.T).T
    
    #flip df back
    cells_df = cells_dft.T
    
    return cells_df

def calculate_averages_and_ttest(df,normalize_row=None):
    #define empty results list 
    resultslist = []
    
    #optional: normalize by row
    if normalize_row is not None:
        df[:-7].loc[:,df.columns[0:-5]] = df[:-7].loc[:,df.columns[0:-5]]/df.T[normalize_row][0:-5]
    
    #iterate through rows (i.e. calculate once per are )
    for idx,row in df.iterrows():
        #skip the last 4 rows
        if idx not in ["batch_no","mouse_no","slice_no","genotype",'time_outside_s','distance_moved_cm','turn_counts']:
            #extract WT + Het subsets
            Het = df.T.loc[df.T["genotype"] == "Het"][idx].dropna().values
            WT = df.T.loc[df.T["genotype"] == "WT"][idx].dropna().values
            #calculate ttest 
            ttest = ttest_ind(list(WT),list(Het),nan_policy="omit",permutations=10000,random_state=rng)
            #ttest = ttest_ind(list(WT),list(Het),nan_policy="omit")
            #add to resultslist 
            resultslist.append(ttest.pvalue)
    
    #debug
    #print("resultslist: ",len(resultslist))
    #print("df: ",df)
    
    #add group mean values per are to df            
    df["Het mean"] = df.T.loc[df.T["genotype"] == "Het"].T.iloc[:-7].T.mean()
    df["WT mean"] = df.T.loc[df.T["genotype"] == "WT"].T.iloc[:-7].T.mean()
    #add stdev 
    df["Het stdev"] = df.T.loc[df.T["genotype"] == "Het"].T.iloc[:-7].T.std()
    df["WT stdev"] = df.T.loc[df.T["genotype"] == "WT"].T.iloc[:-7].T.std()
    #add the ttest resultslist to the df
    df["ttest"] = ""
    df["ttest"].iloc[:-7] = resultslist
    
    return df 

def plot_multiplot(df,output_dir,output_name):
    #TODO: update to include seaborn swarmplot
    # http://seaborn.pydata.org/generated/seaborn.swarmplot.html#seaborn.swarmplot
    # https://stackoverflow.com/questions/8671808/avoiding-overlapping-datapoints-in-a-scatter-dot-beeswarm-plot

    '''    
    #extract genotypes
    Het = df.T.loc[df.T["genotype"] == "Het"].dropna().values
    WT = df.T.loc[df.T["genotype"] == "WT"].dropna().values

    #define line properties    
    whiskerprops = dict(linestyle = "-", linewidth=2.5, color='black')
    meanlineprops = dict(linestyle='-', linewidth=2.5, color='white')
    capprops = whiskerprops
    lineprops = dict(linestyle = "-", linewidth=2.5, color='black')

    #draw boxplots    
    fig, ax = plt.subplots()
    ax.set_ylabel('fruit weight (g)')
    bplot = ax.boxplot([WT, Het], patch_artist=True, medianprops=meanlineprops,whiskerprops=whiskerprops,capprops=capprops,boxprops=lineprops)
    # fill with colors
    for patch, color in zip(bplot['boxes'], ["gray","red"]):
        patch.set_facecolor(color)
    '''
    #cut off unnecessary bits (e.g. batch #, ttest, summaries etc )
    reduced_df = df.T[0:-5].T[0:-7].T

    #extract genotypes and rename so that WT is first 
    gt = df.T[0:-5].T.loc["genotype"]
    gt = gt.str.replace("WT","Ehmt1++")

    #reattach to cleaned-up dataframe
    reduced_df["gt"] = gt
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

def collapse_areas(cells_df, areas_list):
    
    #first, generate list of area names in cells_df
    colnames = cells_df.T.columns.to_list()
    
    #then, generate a subset of area names that contain the keyword (e.g. "ORB")
    for area in areas_list:
        areas_to_collapse = [x for x in colnames if area in x]
        
        #summarize the areas into the superset, then drop the original columns 
        cells_df = cells_df.T
        try:
            cells_df[area] = cells_df.T.loc[areas_to_collapse].sum()
            #at the end, remove the areas from the 
            cells_df.drop(areas_to_collapse,axis=1,inplace=True)
        except:
            print("areas not present, skipping: ",areas_to_collapse)
            pass
        cells_df = cells_df.T
   
    return cells_df

#####
#debug
#source_dirs = ['/media/wirrbel/MN_4/Amygdala_paper_collection/Soc_Exp_cFos_PV/2022-07-03_Behav_cFos_PV_atlas_alignment/Soc_Exp_Stitch_B1+2+3/test/']

'''
#6TB HDD
source_dirs = ['/media/wirrbel/MN_4/Amygdala_paper_collection/Soc_Exp_cFos_PV/2022-07-03_Behav_cFos_PV_atlas_alignment/Soc_Exp_Stitch_B1+2+3/B1-3_tiffs/',
               '/media/wirrbel/MN_4/Amygdala_paper_collection/Soc_Exp_cFos_PV/Soc_Exp_B4-5_overviews/',
               '/media/wirrbel/MN_4/Amygdala_paper_collection/Soc_Exp_cFos_PV/Soc_Exp_B6-12_overviews/B6/'
               ]

'''

# m2 SSD
source_dirs = ['/media/wirrbel/7024783D696D0AAE/amygdala_paper_data_temp/behav_cfos/B1-3_tiffs/',
               '/media/wirrbel/7024783D696D0AAE/amygdala_paper_data_temp/behav_cfos/Soc_Exp_B4-5_overviews/',
               '/media/wirrbel/7024783D696D0AAE/amygdala_paper_data_temp/behav_cfos/Soc_Exp_B6-12_overviews/B6/'
               ]


target_dir = '/home/wirrbel/2022-01-09_Behav_cFos_counts/v19_all/'

genotype_list_file = '/home/wirrbel/2022-01-09_Behav_cFos_counts/times_files_for_sniffing_v06.xlsx'
#v05 contains updated timing + turn count from 2023-09-30 

#create target_dir if not already present 
if not os.path.exists(target_dir):
    os.mkdir(target_dir)

#define areas of interest 

#uncollapsed
areas_amygdala = ['BLAa','BLAp','BLAv','BMA','CEA','COAp','COApl','EPd','EPv','IA','LA','LHA','MEA','PA']

areas_frontal = ['ILA1','ILA2/3','ILA5','ILA6a','ILA6b','MOp1','MOp2/3','MOp5','MOp6a','MOp6b',
                 'MOs1','MOs2/3','MOs5','MOs6a','MOs6b','ORBl1','ORBl2/3','ORBl5','ORBl6a','ORBl6b',
                 'ORBm1','ORBm2/3','ORBm5','ORBm6a','ORBvl1','ORBvl2/3','ORBvl5','ORBvl6a','ORBvl6b',
                 'PL1','PL2/3','PL5','PL6a','PL6b']
'''
#for area collapse
areas_amygdala = ['BLA','CEA','COA','EP','IA','LA','LHA','MEA','PA']

areas_frontal = ['ILA','MOp','MOs','ORB','PL']

all_collapsible_areas = [*areas_amygdala,*areas_frontal]
'''
#####

cfos_density_amygdala_collection = pd.DataFrame(index=[*areas_amygdala,"batch_no","mouse_no","slice_no","genotype",'time_outside_s','distance_moved_cm','turn_counts'])
cfos_density_frontal_collection = pd.DataFrame(index=[*areas_frontal,"batch_no","mouse_no","slice_no","genotype",'time_outside_s','distance_moved_cm','turn_counts'])

PV_density_amygdala_collection = pd.DataFrame(index=[*areas_amygdala,"batch_no","mouse_no","slice_no","genotype",'time_outside_s','distance_moved_cm','turn_counts'])
PV_density_frontal_collection = pd.DataFrame(index=[*areas_frontal,"batch_no","mouse_no","slice_no","genotype",'time_outside_s','distance_moved_cm','turn_counts'])

for source_dir in source_dirs:
    
    #extract source dir name (e.g. B1-3_tiffs)
    source_dir_name = os.path.split(os.path.dirname(source_dir))[-1]
    
    #count # of pixels 
    pixel_count = mask_atlas_file_and_count_px(source_dir)
    #rename to real area names 
    pixel_count = values_to_acronyms(pixel_count)
    
    #calculate mm²
    mm2_per_voxel = 0.000454**2
    mm2_count = pixel_count * mm2_per_voxel
    
    #collapse areas 
    #mm2_count = collapse_areas(mm2_count,all_collapsible_areas)
    
    #summarize per mouse 
    mm2_count = assign_genotypes(mm2_count, genotype_list_file)
    mm2_count = summarize_per_mouse(mm2_count)
    #replace zeros with NaN (for density calculations later)
    mm2_count[0:-8] = mm2_count[0:-8].replace({0:np.nan})
    
    #collect cFos cells 
    cFos_cells = assemble_overview_table (source_dir, target_dir, 'cFos')
    #collapse areas 
    #cFos_cells = collapse_areas(cFos_cells,all_collapsible_areas)
    #replace 0 by nan
    cFos_cells = cFos_cells.replace(0,np.nan)
    #assign genotypes
    cFos_cells = assign_genotypes(cFos_cells,genotype_list_file)
    #summarize per mouse 
    cFos_cells = summarize_per_mouse(cFos_cells)
    #calculate cFos density 
    cFos_density = cFos_cells[0:-8] / mm2_count[0:-8]
    #re-create genotypes etc 
    cFos_density = reappend_info(cFos_density,mm2_count)
    
    
    #export to excel 
    cFos_cells.T.to_excel(os.path.join(target_dir, "cFos_occurrences_v09_pivot"+source_dir_name+".xlsx"))
    cFos_density.T.to_excel(os.path.join(target_dir, "cFos_densities_v09_pivot"+source_dir_name+".xlsx"))
    cFos_cells.to_excel(os.path.join(target_dir, "cFos_occurrences_v09_"+source_dir_name+".xlsx"))
    cFos_density.to_excel(os.path.join(target_dir, "cFos_densities_v09_"+source_dir_name+".xlsx"))
    
    
    #collect PV cells 
    PV_cells = assemble_overview_table (source_dir, target_dir, 'PV')
    #collapse areas 
    #PV_cells = collapse_areas(PV_cells,all_collapsible_areas)
    #replace 0 by nan
    PV_cells = PV_cells.replace(0,np.nan)
    #assign genotypes
    PV_cells = assign_genotypes(PV_cells,genotype_list_file)
    #summarize per mouse 
    PV_cells = summarize_per_mouse(PV_cells)
    #calculate PV density 
    PV_density = PV_cells[0:-8] / mm2_count[0:-8]
    #re-create genotypes etc 
    PV_density = reappend_info(PV_density,mm2_count)
    
    #export to excel 
    PV_cells.T.to_excel(os.path.join(target_dir, "PV_occurrences_v09_pivot"+source_dir_name+".xlsx"))
    PV_density.T.to_excel(os.path.join(target_dir, "PV_densities_v09_pivot"+source_dir_name+".xlsx"))
    PV_cells.to_excel(os.path.join(target_dir, "PV_occurrences_v09_"+source_dir_name+".xlsx"))
    PV_density.to_excel(os.path.join(target_dir, "PV_densities_v09_"+source_dir_name+".xlsx"))
    
    #subset for OFC and amygdala
    subset_and_save(cFos_cells,areas_amygdala,target_dir,str('cFos_occurrences_v09_amygdala'+source_dir_name))
    subset_and_save(cFos_cells,areas_frontal,target_dir,str('cFos_occurrences_v09_OFC'+source_dir_name))
    cfos_density_amygdala_collection = subset_and_save(cFos_density,areas_amygdala,target_dir,str('cFos_density_v09_amygdala'+source_dir_name),cfos_density_amygdala_collection)
    cfos_density_frontal_collection = subset_and_save(cFos_density,areas_frontal,target_dir,str('cFos_density_v09_OFC'+source_dir_name),cfos_density_frontal_collection)
    
    #subset for OFC and amygdala
    subset_and_save(PV_cells,areas_amygdala,target_dir,str('PV_occurrences_v09_amygdala'+source_dir_name))
    subset_and_save(PV_cells,areas_frontal,target_dir,str('PV_occurrences_v09_OFC'+source_dir_name))
    PV_density_amygdala_collection = subset_and_save(PV_density,areas_amygdala,target_dir,str('PV_density_v09_amygdala'+source_dir_name),PV_density_amygdala_collection)
    PV_density_frontal_collection = subset_and_save(PV_density,areas_frontal,target_dir,str('PV_density_v09_OFC'+source_dir_name),PV_density_frontal_collection)
    
    #save mm2_count
    #mm2_count = assign_genotypes(mm2_count,genotype_list_file)
    subset_and_save(mm2_count,areas_amygdala,target_dir,str('mm2_count_v09_amygdala'+source_dir_name))
    subset_and_save(mm2_count,areas_frontal,target_dir,str('mm2_count_v09_OFC'+source_dir_name))

    #save pixel_count
    pixel_count = assign_genotypes(pixel_count,genotype_list_file)
    subset_and_save(pixel_count,areas_amygdala,target_dir,str('pixel_count_v09_amygdala'+source_dir_name))
    subset_and_save(pixel_count,areas_frontal,target_dir,str('pixel_count_v09_OFC'+source_dir_name))
    

#cleanup
try:
    #remove b4m1 from amygdala dataset, it only has 2 slices there 
    cfos_density_amygdala_collection.drop(13,axis=1,inplace=True)
    PV_density_amygdala_collection.drop(13,axis=1,inplace=True)
except:
     pass

try:
    #replace zeros with NaN
    cfos_density_amygdala_collection = cfos_density_amygdala_collection.replace(0,np.nan)
    cfos_density_frontal_collection = cfos_density_frontal_collection.replace(0,np.nan)
    PV_density_amygdala_collection = PV_density_amygdala_collection.replace(0,np.nan)
    PV_density_frontal_collection = PV_density_frontal_collection.replace(0,np.nan)
except:
    pass

try:
    #remove b4-m6 (since reassigned)
    cfos_density_amygdala_collection.drop("B4_M6_Sl6_",axis=1,inplace=True)
    cfos_density_frontal_collection.drop("B4_M6_Sl6_",axis=1,inplace=True)
    PV_density_amygdala_collection.drop("B4_M6_Sl6_",axis=1,inplace=True)
    PV_density_frontal_collection.drop("B4_M6_Sl6_",axis=1,inplace=True)
except:
    pass

#calculate ttest stats 
cfos_density_amygdala_collection = calculate_averages_and_ttest(cfos_density_amygdala_collection)
cfos_density_frontal_collection = calculate_averages_and_ttest(cfos_density_frontal_collection)
PV_density_amygdala_collection = calculate_averages_and_ttest(PV_density_amygdala_collection)
PV_density_frontal_collection = calculate_averages_and_ttest(PV_density_frontal_collection)

#save density collections 
cfos_density_amygdala_collection.to_excel(os.path.join(target_dir, "density_supercollection_cFos_amygdala_v09.xlsx"))
cfos_density_frontal_collection.to_excel(os.path.join(target_dir, "density_supercollection_cFos_frontal_v09.xlsx"))
PV_density_amygdala_collection.to_excel(os.path.join(target_dir, "density_supercollection_PV_amygdala_v09.xlsx"))
PV_density_frontal_collection.to_excel(os.path.join(target_dir, "density_supercollection_PV_frontal_v09.xlsx"))

#generate multi-violin plots
plot_multiplot(cfos_density_amygdala_collection,target_dir,"cfos_density_amygdala_collection")
plot_multiplot(cfos_density_frontal_collection,target_dir,"cfos_density_frontal_collection")
plot_multiplot(PV_density_amygdala_collection,target_dir,"PV_density_amygdala_collection")
plot_multiplot(PV_density_frontal_collection,target_dir,"PV_density_frontal_collection")

'''
#normalize cfos data to the behavioural stats
cfos_amygdala_normalized_time = calculate_averages_and_ttest(cfos_density_amygdala_collection.copy(),'time_outside_s')
cfos_frontal_normalized_time = calculate_averages_and_ttest(cfos_density_frontal_collection.copy(),'time_outside_s')

cfos_amygdala_normalized_distance = calculate_averages_and_ttest(cfos_density_amygdala_collection.copy(),'distance_moved_cm')
cfos_frontal_normalized_distance = calculate_averages_and_ttest(cfos_density_frontal_collection.copy(),'distance_moved_cm')

cfos_amygdala_normalized_turns = calculate_averages_and_ttest(cfos_density_amygdala_collection.copy(),'turn_counts')
cfos_frontal_normalized_turns = calculate_averages_and_ttest(cfos_density_frontal_collection.copy(),'turn_counts')

cfos_amygdala_normalized_time.to_excel(os.path.join(target_dir, "density_normalized_time_cFos_amygdala_v09.xlsx"))
cfos_frontal_normalized_time.to_excel(os.path.join(target_dir, "density_normalized_time_cFos_frontal_v09.xlsx"))
cfos_amygdala_normalized_distance.to_excel(os.path.join(target_dir, "density_normalized_distance_cFos_amygdala_v09.xlsx"))
cfos_frontal_normalized_distance.to_excel(os.path.join(target_dir, "density_normalized_distance_cFos_frontal_v09.xlsx"))
cfos_amygdala_normalized_turns.to_excel(os.path.join(target_dir, "density_normalized_turns_cFos_amygdala_v09.xlsx"))
cfos_frontal_normalized_turns.to_excel(os.path.join(target_dir, "density_normalized_turns_cFos_frontal_v09.xlsx"))
'''

'''
#define directory where files should go 
target_dir = "/home/wirrbel/2022-01-09_Behav_cFos_counts/"

#cFos data is in this folder 
cFos_filelist = sorted(glob.glob("/media/wirrbel/MN_2/Soc_Exp_cFos_counts/*cFos_ROIs.csv"))

#extract collection of all cFos data 
cFos_collection_df = generate_collection(cFos_filelist)

#fix naming 
cFos_collection_df = values_to_acronyms(cFos_collection_df)

#sort index in alphabetical order 
cFos_collection_df = cFos_collection_df.sort_index()

#count the number of cFos cells per area 
#cFos_occurrences = count_occurrences(cFos_collection_df)

#optional: Make pivot table
#cFos_occurrences = cFos_occurrences.T
cFos_collection_df = cFos_collection_df.T

#save to excel file 
cFos_collection_df.to_excel(os.path.join(target_dir,"cFos_occurrences_v04_pivot.xlsx"))


####
#PV data is in this folder 
PV_filelist = sorted(glob.glob("/media/wirrbel/MN_2/Soc_Exp_cFos_counts/*PV_ROIs.csv"))

#extract collection of all PV data 
PV_collection_df = generate_collection(PV_filelist)

#fix naming 
PV_collection_df = values_to_acronyms(PV_collection_df)

#sort index in alphabetical order 
PV_collection_df = PV_collection_df.sort_index()

#count the number of PV cells per area 
#PV_occurrences = count_occurrences(PV_collection_df)

#optional: Make pivot table
#PV_occurrences = PV_occurrences.T
PV_collection_df = PV_collection_df.T

#save to excel file 
PV_collection_df.to_excel(os.path.join(target_dir,"PV_occurrences_v04.xlsx"))
'''