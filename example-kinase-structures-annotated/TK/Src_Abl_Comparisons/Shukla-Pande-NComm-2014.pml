# This file loads a pymol input macro that compares the structures used in
# Shukla...Pande, Nat. Comm. (2014)  
# The main aim is to look at inactive vs. actove of Src with the following PDBs:
# 1Y57 - SH3, SH2 and kinase domains; apo; (semi-)active extended conformation; Y530 not phosphorylated; ATP-competitive inhibitor
# 2SRC - Human Src residues 86-836, containing SH2, SH3, Kinase domain, and C-terminal tail; AMP-PNP

# load PDBs

load ../src/1Y57.pdb
load ../src/2SRC.pdb

# Each of these already only consists of just Chain A, but there are SH3-SH2 domains in both structures

# align two structures

# define and color relevant motifs
# Kinase Domain
select 1Y57_kin=(1Y57 & POL & i;258-522)
select 2SRC_kin=(2SRC & POL & i;258-522)
color gray70, 1Y57_kin
color gray30, 2SRC_kin
# ligands
select 1Y57_MPZ=(1Y57 & resname MPZ)
select 2SRC_ANP=(2SRC & resname ANP)
color palecyan, 1Y57_MPZ
color tv_blue, 2SRC_ANP

# DFG
select 1Y57_DFG=(1Y57_kin & i;404-406)
select 2SRC_DFG=(2SRC_kin & i;404-406)
color limon, 1Y57_DFG
color forest, 2SRC_DFG

# P-Loop (LGGGQYGEVY or LGQGCFGEVW) 
# P for Purple
select 1Y57_Ploop=(1Y57_kin& POL & i;273-282)
select 2SRC_Ploop=(2SRC_kin& POL & i;273-282)
color violet, 1Y57_Ploop
color deeppurple, 2SRC_Ploop
# C-Helix (EVEEFLKEAAVMKEI or MSPEAFLQEAEQVMKKL)
select 1Y57_CHelix=(1Y57_kin & POL & i;302-317)
select 2SRC_CHelix=(2SRC_kin & POL & i;302-317)
color paleyellow, 1Y57_CHelix
color brightorange, 2SRC_CHelix
# A-Loop (LSRLMTGDTYTAHAGAKFP or LARLIEDNEYTARQGAKFP)
# Note A loop not resolved in this structure
select 1Y57_Aloop=(1Y57_kin & POL & i;407-425)
select 2SRC_Aloop=(2SRC_kin & POL & i;407-425)
color salmon, 1Y57_Aloop
color firebrick, 2SRC_Aloop

# SH3_SH2 Domain
select 1Y57_SH2SH3=(1Y57 & POL & i;84-258)
select 2SRC_SH2SH3=(2SRC & POL & i;84-258)
# C-term loop
select 2SRC_Cterm=(2SRC & POL & i;522-533)

# align two structures
align 1Y57_kin, 2SRC_kin

# rotate to relevant view point
rotate y, 80
rotate z, -80
rotate y, 50
rotate x, 40



center 1Y57_kin

# set rendering conditions
hide all
show cartoon, 1Y57_kin or 2SRC_kin
show sticks, 1Y57_MPZ or 2SRC_ANP
show sticks, 1Y57_DFG or 2SRC_DFG
util.cnc("1Y57_DFG or 2SRC_DFG")
util.cnc("1Y57_MPZ or 2SRC_ANP")

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue