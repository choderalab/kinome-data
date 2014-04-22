# This file loads a pymol input macro that compares the structures used in
# Shan...Shaw, PNAS (2009)  
# The main aim is to look at DFG-in vs DFG-out structures of Abl with the following PDBs:
# 1OPK - mouse SH3-SH2-kinase domain; 'Structural basis for auto-inhibition'; PD166326
# 2F4J - kinase domain with C-term up to 513 (straight helix)

# load PDBs

load ../Abl1/2F4J.pdb
load ../Abl1/1OPK.pdb

# show only a single chain

#     These structures are each only one chain, so need need to hide chains.
#     However, SH3-SH2 domain is present on one, but this will just be annotated
#     as a separate domain for the user to decide how they want to visualize.

# define and color relevant motifs
# Kinase Domain
select 2F4J_kin=(2F4J & POL & i;227-513)
select 1OPK_kin=(1OPK & POL & i;246-531)
color gray70, 2F4J_kin
color gray30, 1OPK_kin
# ligands
select VX6=(resname VX6)
select P16=(resname P16)
color palecyan, VX6
color tv_blue, P16

# DFG
select 2F4J_DFG=(2F4J & POL & i;381-383)
select 1OPK_DFG=(1OPK & POL & i;400-402)
color limon, 2F4J_DFG
color forest, 1OPK_DFG

# P-Loop (LGGGQYGEVY) 
# P for Purple
select 2F4J_Ploop=(2F4J & POL & i;248-257)
select 1OPK_Ploop=(1OPK & POL & i;267-276)
color violet, 2F4J_Ploop
color deeppurple, 1OPK_Ploop
# C-Helix (EVEEFLKEAAVMKEI)
select 2F4J_CHelix=(2F4J & POL & i;279-293)
select 1OPK_CHelix=(1OPK & POL & i;298-312)
color paleyellow, 2F4J_CHelix
color brightorange, 1OPK_CHelix
# A-Loop (LSRLMTGDTYTAHAGAKFP)
select 2F4J_Aloop=(2F4J & POL & i;384-402)
select 1OPK_Aloop=(1OPK & POL & i;403-421)
color salmon, 2F4J_Aloop
color firebrick, 1OPK_Aloop
# SH3_SH2 Domain
select 1OPK_SH2SH3=(1OPK & POL & i;83-245)

# align two structures
align 2F4J_kin, 1OPK_kin

# rotate to relevant view point
rotate x, 45
rotate z, 180
rotate y, 90

# set rendering conditions
hide all
show cartoon
show sticks, VX6 or P16 
show sticks, 2F4J_DFG or 1OPK_DFG
util.cnc("2F4J_DFG or 1OPK_DFG")
util.cnc("VX6 or P16")

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue