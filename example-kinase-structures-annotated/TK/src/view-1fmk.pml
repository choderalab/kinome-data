load 1FMK.pdb

select 1kin=(1fmk & i;267-520)
select 1sh2=(1fmk & i;148-245)
select 1sh3=(1fmk & i;81-142)
select phosphates=(1fmk & i;409+424+416+527)  # 409 and 424 are not real phospho sites - these are just to mark rough location of 416 (in substrate binding cleft)
select gatekeeper=(1fmk & i;338)

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

