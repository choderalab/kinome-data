load 1Y57.pdb

select 1kin=(1y57 & i;267-520)
select 1sh2=(1y57 & i;148-245)
select 1sh3=(1y57 & i;81-142)
select phosphates=(1y57 & i;416+527)  # 409 and 424 are not real phospho sites - these are just to mark rough location of 416 (in substrate binding cleft)
select gatekeeper=(1y57 & i;338)

hide all
show cartoon
color gray90, all
color red, 1kin
color cyan, 1sh2
color yellow, 1sh3
show sticks, phosphates
color magenta, phosphates
color blue, gatekeeper

#background white

