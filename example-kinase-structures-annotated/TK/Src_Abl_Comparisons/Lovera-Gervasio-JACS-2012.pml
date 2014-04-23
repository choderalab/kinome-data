# This file loads a pymol input macro that compares the structures used in
# Lovera...Gervasio, JACS (2012)
# The main aim is to look at DFG-in in Abl vs. Src with the following PDBs:
# 2G1T - human abl1 kinase domain; 'Src-like Inactive Conformation'
# 2SRC - Human Src residues 86-836, containing SH2, SH3, Kinase domain, and C-terminal tail; AMP-PNP

# load PDBs

load ../Abl1/2G1T.pdb
load ../src/2SRC.pdb

# show only chain A (of A-H for 2G1T and just A for 2SRC)
# 2SRC also has the SH2, SH3 domain and C-terminal tail

# define and color relevant motifs
# Kinase Domain
select 2G1T_kin=(2G1T & POL & chain A)
select 2SRC_kin=(2SRC & POL & i;258-522)
color gray70, 2G1T_kin
color gray30, 2SRC_kin
# ligands
select 2G1T_112=(2G1T & resname 112 & chain E)
select 2SRC_ANP=(2SRC & resname ANP)
color palecyan, 2G1T_112
color tv_blue, 2SRC_ANP

# DFG
select 2G1T_DFG=(2G1T_kin & i;381-383)
select 2SRC_DFG=(2SRC_kin & i;404-406)
color limon, 2G1T_DFG
color forest, 2SRC_DFG

# P-Loop (LGGGQYGEVY or LGQGCFGEVW) 
# P for Purple
select 2G1T_Ploop=(2G1T_kin& POL & i;248-257)
select 2SRC_Ploop=(2SRC_kin& POL & i;273-282)
color violet, 2G1T_Ploop
color deeppurple, 2SRC_Ploop
# C-Helix (EVEEFLKEAAVMKEI or MSPEAFLQEAEQVMKKL)
select 2G1T_CHelix=(2G1T_kin & POL & i;279-293)
select 2SRC_CHelix=(2SRC_kin & POL & i;302-317)
color paleyellow, 2G1T_CHelix
color brightorange, 2SRC_CHelix
# A-Loop (LSRLMTGDTYTAHAGAKFP or LARLIEDNEYTARQGAKFP)
# Note A loop not resolved in this structure
select 2G1T_Aloop=(2G1T_kin & POL & i;384-402)
select 2SRC_Aloop=(2SRC_kin & POL & i;407-425)
color salmon, 2G1T_Aloop
color firebrick, 2SRC_Aloop

# SH3_SH2 Domain
select 2SRC_SH2SH3=(2SRC & POL & i;84-258)
# C-term loop
select 2SRC_Cterm=(2SRC & POL & i;522-533)

# align two structures
align 2G1T_kin, 2SRC_kin

# rotate to relevant view point
rotate x, -90
rotate z, 180
rotate y, -40

center 2G1T_kin

# set rendering conditions

hide all
show cartoon, 2G1T_kin or 2SRC_kin
show sticks, 2G1T_112 or 2SRC_ANP
show sticks, 2G1T_DFG or 2SRC_DFG
util.cnc("2G1T_DFG or 2SRC_DFG")
util.cnc("2G1T_112 or 2SRC_ANP")

bg_color white
set ray_shadows, 0
set orthoscopic, on
set ray_opaque_background, off

unset specular
set ray_trace_gain, 0
set ray_trace_mode, 3
set ray_trace_color, black
unset depth_cue