load 2F4J.pdb
load 2FO0.pdb

select 1kin=(2f4j & i;242-493)
select 2kin=(2fo0 & i;261-512)
select 1helI=(2f4j & i;500-513)
select 2helI=(2fo0 & i;519-532)
select 2sh2=(2fo0 & i;127-212)
select 2sh3=(2fo0 & i;80-140)
select 1lig=(2f4j & resname;VX6)
select 2lig=(2fo0 & resname;P16)
select myr=(2fo0 & resname;MYR)

align 2kin, 1kin

hide all
show cartoon
color gray90, all
color red, 1kin
color orange, 2kin
color forest, 1helI
color green, 2helI
color cyan, 2sh2
color yellow, 2sh3
show sticks, 1lig
show sticks, 2lig
util.cbag 1lig
util.cbag 2lig
show sticks, myr

#background white

center 2f4j

