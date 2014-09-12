# This script aligns all PDBs of BTK in pymol.
# Sonya Hanson
# September 11, 2014

#Right now this is made stupidly using ipython.

# BTK_structures = ["1aww", "1awx", "1b55", "1btk", "1bwn", "1k2p", "1qly", "2ge9", "2z0p", "3gen", "3k54", "3ocs", "3oct", "3p08", "3pix", "3piy", "3piz", "3pj1", "3pj2", "3pj3", "4nwm", "4ot5", "4ot6", "4otq", "4otr"]
# ref = BTK_structures[0]
# for struct in BTK_structures[1:]:
#    print 'align %s, %s \n' % (struct,ref)

# But actually we're only interested in the kinase domain.
# Checked that "1aww", "1awx", "1b55", "1btk", "1bwn", "1qly", "2ge9", "2z0p" are not kinases.

#REDO
# BTK_structures = ["1k2p", "3gen", "3k54", "3ocs", "3oct", "3p08", "3pix", "3piy", "3piz", "3pj1", "3pj2", "3pj3", "4nwm", "4ot5", "4ot6", "4otq", "4otr"]
# ref = BTK_structures[0]
# for struct in BTK_structures[1:]:
#    print 'align %s, %s \n' % (struct,ref)

# Then  'open -a MacPyMOL *.pdb' and @BTK_align.pml

# align PDBs

align 3gen, 1k2p 
align 3k54, 1k2p 
align 3ocs, 1k2p 
align 3oct, 1k2p 
align 3p08, 1k2p 
align 3pix, 1k2p 
align 3piy, 1k2p 
align 3piz, 1k2p 
align 3pj1, 1k2p 
align 3pj2, 1k2p 
align 3pj3, 1k2p 
align 4nwm, 1k2p 
align 4ot5, 1k2p 
align 4ot6, 1k2p 
align 4otq, 1k2p 
align 4otr, 1k2p 