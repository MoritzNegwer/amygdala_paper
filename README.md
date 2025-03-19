This is the code that belongs with the upcoming publication "Whole-Brain Clearing Reveals Region- and Cell-Specific Imbalances in Inhibitory Neurons in Kleefstra Syndrome". 

Please find the related dataset on Zenodo: https://doi.org/10.5281/zenodo.15048387 

Abstract:
Inhibition is crucial for balanced brain function and often disrupted in autism spectrum disorder (ASD). Inhibitory neurons are uniquely vulnerable to developmental insults due to their large diversity and long migration trajectories. However, the diversity of inhibitory neurons – cell types as well as their spatial distribution – makes it difficult to assess the effects of a mutation. Using an unbiased whole-brain clearing and light-sheet imaging approach, we mapped the spatial distribution of the three main inhibitory cell types – Parvalbumin (PV+), Somatostatin (SST+) and Vasoactive Intestinal Polypeptide (VIP+) across the mouse brain in a mouse model of Kleefstra Syndrome, a monogenetic Intellectual Disability syndrome with strong ASD component. We found changes in all three types, including a surprising early hypermaturation of PV+ neurons in the basolateral amygdala (BLA) in adolescence (P14). This leads to an enhanced inhibitory drive in the adolescent BLA, but leaves the excitatory/inhibitory ratio unchanged, indicating a circuit-level homeostatic mechanism. Consequently, we found no change in behaviour or amygdalar activity (using c-Fos+ density as a proxy) following a social novelty challenge. However, we discovered a reduction in diurnal activity in Ehmt1+/- males in a group-based social activity tracking setup, which fits with higher PV+ and VIP+ density in sleep-associated thalamic and midbrain regions we identified with our whole-brain screening approach. Here we demonstrate that monogenetic mutations in ASD model mice have diverse region- and cell-type specific effects, and that an unbiased whole-brain screening can point towards especially affected circuits that are promising targets for in-depth investigation.

The code is divided by figure: 
Figure 2: 
- iDISCO+ light-sheet image stacks for VIP+, PV+ and SST+ reporter mouse brain hemispheres, immunostained against the tdTomato reporter and imaged with a light-sheet microscope.
- The brains were processed with FriendlyClearMap ([https://github.com/moritzne](https://github.com/MoritzNegwer/FriendlyClearMap-scripts)
- The resulting regions were rendered with BrainRender, the script to do so is here.

Figure 3:
- Flat-map projection of the cortex in Fig. 2, using the https://github.com/MoritzNegwer/FriendlyClearMap-scripts/tree/main/cortical_flatmaps version of the https://github.com/AllenInstitute/ccf_streamlines code from the Allen Brain Institute.
- The script to generate the flatmaps used in Fig. 3 is placed here.

Figure 4 supplementary figure 2: 
- Slice immunostainings against PV and c-Fos.
- Behav_cFos_PV_v19.py: Python script that I used to calculate the cell densities from the Ilastik-detected PV+ and c-Fos+ cells and the slice-registered atlas image (a plane from the Allen Brain Atlas CCF3 warped with Fiji's BigWarp plugin (from https://github.com/saalfeldlab/bigwarp))
- plot_p14_p56_amygdala: generates the violin plots in Fig. 4a-b.

Figure 5 a-g:
- Single-cell electrophysiology measuring the inhibitory input to principal neurons in the BLA
- Create_mini_figs_v07_mIPSC.py: Script to process the miniIPSC traces from Fig. 5a
- Create_Overview_table_v07_mIPSC.py: Script to generate the figures in Fig. 5b-c
- Peak_Detector_v02_scriptable.py: Script to call the peaks in the 10Hz stimulus train and analyze their relative amplitudes
- 10Hz_overview_plot_v01.py: Generates the overview plot in Fig. 5d

Figure 5 h-k and supplementary fig. 3 a-b:
- Single mouse reaction to a social novelty challenge
- calculate_distance_social_uncluttering_v11.py: Python script used to calculate movement in the "away zone" using the functions from lmt-analysis software (https://github.com/fdechaumont/lmt-analysis). Used to generate the 1D-heatmap plots in Fig. 5j and Supp. Fig. 3 a-b

Supplementary Fig. 5c: 
- Activity measurements in a social setting
- R scripts used to generate the plots in Supplementary Fig. 3c, based on the results from lmt-analysis (https://github.com/fdechaumont/lmt-analysis)
