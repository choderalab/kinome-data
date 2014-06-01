load 2FO0.pdb

select sh3=(2fo0 & i;80-140)
select sh2=(2fo0 & i;146-236)
select kin=(2fo0 & i;261-512)
select myr=(2fo0 & resname;MYR)
select p16=(2fo0 & resname;P16)

select starts=(i;247+248+249+256+257+259+260+261+262)
select ends=(i;509+512+513+515+519+530+531+532+534)
select starts_ca=(starts & name;CA)
select ends_ca=(ends & name;CA)

hide all
show cartoon
#show sticks, myr
#show sticks, p16
color gray90, all
color red, kin
color cyan, sh2
color yellow, sh3
color gray50, myr

show spheres, starts_ca
show spheres, ends_ca
#color cyan, starts_ca
#color cyan, ends_ca

#background white
set seq_view, 1

