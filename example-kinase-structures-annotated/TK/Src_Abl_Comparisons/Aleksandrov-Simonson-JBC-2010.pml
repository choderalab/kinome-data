# This file loads a pymol input macro that compares the structures used in
# Aleksandrov & Simonson, JBC (2010)  
# The main aim is to look at DFG-out in Abl vs. Src with the following PDBs:
# 2HYY - human abl1 kinase domain; imatinib
# 2OIQ - Chicken Src kinase domain; imatinib
          
# load PDBs

load ../Abl1/2HYY.pdb
load ../src/2OIQ.pdb

# show only chain A (of A-D of 2HYY and of A and B of 2OIQ)

# define and color relevant motifs
# Kinase Domain
select 2HYY_kin=(2HYY & POL & chain A)
select 2OIQ_kin=(2OIQ & POL & chain A)
color gray70, 2HYY_kin
color gray30, 2OIQ_kin
# ligands
select 2HYY_IMA=(2HYY & resname STI & chain A)
select 2OIQ_IMA=(2OIQ & resname STI & chain A)
color palecyan, 2HYY_IMA
color tv_blue, 2OIQ_IMA

# DFG
select 2HYY_DFG=(2HYY_kin & i;381-383)
select 2OIQ_DFG=(2OIQ_kin & i;404-406)
color limon, 2HYY_DFG
color forest, 2OIQ_DFG

# P-Loop (LGGGQYGEVY or LGQGCFGEVW) 
# P for Purple
select 2HYY_Ploop=(2HYY_kin& POL & i;248-257)
select 2OIQ_Ploop=(2OIQ_kin& POL & i;273-282)
color violet, 2HYY_Ploop
color deeppurple, 2OIQ_Ploop
# C-Helix (EVEEFLKEAAVMKEI or MSPEAFLQEAEQVMKKL)
select 2HYY_CHelix=(2HYY_kin & POL & i;279-293)
select 2OIQ_CHelix=(2OIQ_kin & POL & i;302-317)
color paleyellow, 2HYY_CHelix
color brightorange, 2OIQ_CHelix
# A-Loop (LSRLMTGDTYTAHAGAKFP or LARLIEDNEYTARQGAKFP)
# Note A loop not resolved in this structure
select 2HYY_Aloop=(2HYY_kin & POL & i;384-402)
select 2OIQ_Aloop=(2OIQ_kin & POL & i;407-425)
color salmon, 2HYY_Aloop
color firebrick, 2OIQ_Aloop

# align two structures
align 2HYY_kin, 2OIQ_kin

# rotate to relevant view point
rotate x, -90
rotate z, -40
rotate y, -20

center 2HYY_kin

# set rendering conditions
hide all
show cartoon, 2HYY_kin or 2OIQ_kin
show sticks, 2HYY_IMA or 2OIQ_IMA
show sticks, 2HYY_DFG or 2OIQ_DFG
util.cnc("2HYY_DFG or 2OIQ_DFG")
util.cnc("2HYY_IMA or 2OIQ_IMA")

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue