load 2FO0.pdb

select sh3=(2fo0 & i;80-140)
select sh2=(2fo0 & i;146-236)
select kin=(2fo0 & i;261-512)
select myr=(2fo0 & resname;MYR)
select p16=(2fo0 & resname;P16)

hide all
show cartoon
show sticks, myr
show sticks, p16
color gray90, all
color red, kin
color cyan, sh2
color yellow, sh3
color gray50, myr

#background white
set seq_view, 1

