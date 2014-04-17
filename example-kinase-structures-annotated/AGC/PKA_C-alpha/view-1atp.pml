load 1ATP.pdb

# annotations based on: Endicott Johnson Annu Rev Biochem 2012

select inhibitor_peptide=(c;I)

select atp=(resname;ATP)
select Mg=(resname;MN)

select nterm_extension=(i;15-39)
select nterm_lobe=(i;40-125)
select cterm_lobe=(i;126-299)
select cterm_extension=(i;300-350)

select gly_rich_loop=(i;50-56)

select chelix=(i;85-97)

select catalytic_aspartate=(i;166)

# Activation segment aka T-loop. often partially disordered in inactive kinases. Adoption of the catalytically competent conformation is triggered in many kinases by phosphorylation.
select activation_segment=(i;184-208)

# DFG loop. Highly conserved. In inactive conformation, this structure is often inverted, so that the F points toward the ATP instead of the D
select DFG=(i;184-186)

# APE loop. conserved, although less so than DFG
select APE=(i;206-208)

select other_conserved_atp_binding_residues=(i;72+91+127+170+171+184)

hide all
show cartoon
hide cartoon, inhibitor_peptide
show sticks, atp
show nb_spheres, Mg
show sticks, DFG
show sticks, APE
show sticks, other_conserved_atp_binding_residues
show sticks, catalytic_aspartate

color gray20, nterm_extension
color gray20, cterm_extension
color white, nterm_lobe
color yellow, cterm_lobe
color red, chelix
color pink, gly_rich_loop
color cyan, activation_segment

