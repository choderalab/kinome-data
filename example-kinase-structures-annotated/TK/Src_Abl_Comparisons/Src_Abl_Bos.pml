# This file loads a pymol input macro that compares the structures used in
# Src vs. Abl bound to BOS.
          
# load PDBs

load ../Abl1/3UE4.pdb
load ../src/4MXO.pdb

# define and color relevant motifs
# Kinase Domain
select 3UE4_kin=(3UE4 & POL & chain A)
select 4MXO_kin=(4MXO & POL & chain A)
color gray70, 3UE4_kin
color gray30, 4MXO_kin
# ligands
select 3UE4_BOS=(3UE4 & resname DB8 & chain A)
select 4MXO_BOS=(4MXO & resname DB8 & chain A)
color palecyan, 3UE4_BOS
color tv_blue, 4MXO_BOS

# DFG
select 3UE4_DFG=(3UE4_kin & i;381-383)
select 4MXO_DFG=(4MXO_kin & i;404-406)
color limon, 3UE4_DFG
color forest, 4MXO_DFG

# P-Loop (LGGGQYGEVY or LGQGCFGEVW) 
# P for Purple
select 3UE4_Ploop=(3UE4_kin& POL & i;248-257)
select 4MXO_Ploop=(4MXO_kin& POL & i;273-282)
color violet, 3UE4_Ploop
color deeppurple, 4MXO_Ploop
# C-Helix (EVEEFLKEAAVMKEI or MSPEAFLQEAEQVMKKL)
select 3UE4_CHelix=(3UE4_kin & POL & i;279-293)
select 4MXO_CHelix=(4MXO_kin & POL & i;302-317)
color paleyellow, 3UE4_CHelix
color brightorange, 4MXO_CHelix
# A-Loop (LSRLMTGDTYTAHAGAKFP or LARLIEDNEYTARQGAKFP)
# Note A loop not resolved in this structure
select 3UE4_Aloop=(3UE4_kin & POL & i;384-402)
select 4MXO_Aloop=(4MXO_kin & POL & i;407-425)
color salmon, 3UE4_Aloop
color firebrick, 4MXO_Aloop

# align two structures
align 3UE4_kin, 4MXO_kin

# rotate to relevant view point
rotate x, -90
rotate z, -40
rotate y, -20

center 3UE4_kin

# set rendering conditions
hide all
show cartoon, 3UE4_kin or 4MXO_kin
show sticks, 3UE4_BOS or 4MXO_BOS
show sticks, 3UE4_DFG or 4MXO_DFG
util.cnc("3UE4_DFG or 4MXO_DFG")
util.cnc("3UE4_BOS or 4MXO_BOS")

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue