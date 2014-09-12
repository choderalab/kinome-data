# This script downloads all BTK PDB structures.
# Sonya Hanson
# September 11, 2014

# Right now I'm just using uniprot 'cross-refs' PDBs from file Q06187_PDBs.txt and
# extracting the list of PDBs with:
# awk '{print $1","}' Q06187_PDBs.txt | xargs | tr '[:upper:]' '[:lower:]'

# This will be done in the future directly from the UNIPROT server or the UNIPROT XML File.


from Bio import PDB

BTK_structures = ["1aww", "1awx", "1b55", "1btk", "1bwn", "1k2p", "1qly", "2ge9", "2z0p", "3gen", "3k54", "3ocs", "3oct", "3p08", "3pix", "3piy", "3piz", "3pj1", "3pj2", "3pj3", "4nwm", "4ot5", "4ot6", "4otq", "4otr"]

for struct in BTK_structures:
    
    dir = struct[1:3]

    pdb1=PDB.PDBList()
    pdb1.retrieve_pdb_file(struct)

    parser = PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(struct,"%s/pdb%s.ent" % (dir,struct))
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save('%s.pdb' % struct)

ref = BTK_structures[0]

with open('align.pml') as f:
    for struct in BTK_structures[1:]:
        f.write('align %s, %s \n' % (struct,ref))