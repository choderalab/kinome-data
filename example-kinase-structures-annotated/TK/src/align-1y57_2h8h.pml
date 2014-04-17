load 1Y57.pdb
load 2H8H.pdb

select ackin=(1y57 & i;267-520)
select acsh2=(1y57 & i;148-245)
select acsh3=(1y57 & i;81-142)
select inackin=(2h8h & i;267-520)
#select phosphates=(1y57 & i;416+527)  # 409 and 424 are not real phospho sites - these are just to mark rough location of 416 (in substrate binding cleft)
#select gatekeeper=(1y57 & i;338)

align inackin, ackin

hide all
show cartoon
color gray90, all
color red, ackin
color cyan, inackin
#color cyan, acsh2
#color yellow, acsh3
#show sticks, phosphates
#color magenta, phosphates
#color blue, gatekeeper

#background white

