import os
from pylab import *

plasmid_alns = [
'../DFHCC-PlasmID/HIP-human_kinase_collection-pJP1520/aln.txt',
'../DFHCC-PlasmID/HIP-human_kinase_collection-pDNR-Dual-Ctermtagonly/aln.txt',
'../addgene/human-kinase-ORF-collection/aln.txt',
]

labels = [
'HIP-pJP1520',
'HIP-pDNR-Dual',
'addgene-ORF',
]

pctidentities_data = []
nconflicts_data = []
pctconflicts_data = []

for p in range(len(plasmid_alns)):

    plasmid_aln_path = plasmid_alns[p]
    pctidentities = []
    nconflicts = []
    with open(plasmid_aln_path, 'r') as infile:
        for l, line in enumerate(infile.readlines()):
            if l % 3 == 1:
                words = line.split()
                pctidentities.append( float(words[1]) )
                nconflicts.append( int(words[2].split('/')[0]) )

    pctconflicts = [100. - x for x in pctidentities]

    pctidentities = array(pctidentities)
    nconflicts = array(nconflicts)
    pctconflicts = array(pctconflicts)

    pctidentities_data.append(pctidentities)
    nconflicts_data.append(nconflicts)
    pctconflicts_data.append(pctconflicts)


    # summary statistics

    print '=== %s ===' % labels[p]

    print 'Number of plasmids with 0 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts == 0), len(nconflicts), (sum(nconflicts == 0)*100./len(nconflicts)))
    print 'Number of plasmids with 1 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts == 1), len(nconflicts), (sum(nconflicts == 1)*100./len(nconflicts)))
    print 'Number of plasmids with 2 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts == 2), len(nconflicts), (sum(nconflicts == 2)*100./len(nconflicts)))
    print 'Number of plasmids with > 0 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 0), len(nconflicts), (sum(nconflicts > 0)*100./len(nconflicts)))
    print 'Number of plasmids with > 1 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 1), len(nconflicts), (sum(nconflicts > 1)*100./len(nconflicts)))
    print 'Number of plasmids with > 2 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 2), len(nconflicts), (sum(nconflicts > 2)*100./len(nconflicts)))
    print 'Number of plasmids with > 3 conflicts: %d / %d (%.2f %%)' % (sum(nconflicts > 3), len(nconflicts), (sum(nconflicts > 3)*100./len(nconflicts)))
    print 'Number of plasmids with > 10%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 10.), len(nconflicts), (sum(pctconflicts > 10.)*100./len(nconflicts)))
    print 'Number of plasmids with > 25%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 25.), len(nconflicts), (sum(pctconflicts > 25.)*100./len(nconflicts)))
    print 'Number of plasmids with > 50%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 50.), len(nconflicts), (sum(pctconflicts > 50.)*100./len(nconflicts)))
    print 'Number of plasmids with > 75%% conflicts: %d / %d (%.2f %%)' % (sum(pctconflicts > 75.), len(nconflicts), (sum(pctconflicts > 75.)*100./len(nconflicts)))
    print ''


nbins = 10
#bins = [0,0.1,0.2,0.4,0.6,0.8,1,2,3,4,6,8,10,15,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]
bins = [0,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100]
#bins = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20]

hist(pctconflicts_data, nbins, normed=False, log=True, label=labels)
#hist(nconflicts_data, bins=bins, normed=True, log=True, label=labels)

xlabel('% conflicts')
ylabel('population')
ax = gca()
#ax.set_xscale('log')
#ax.set_yscale('log')
#xlim(0,100)
#ylim(0,100)
legend()

a = axes([.57,.35,.3,.3], axisbg='y')
plot = hist(nconflicts_data, bins=bins, normed=False, log=True, label=labels)
xlim(0,10)
xlabel('nconflicts', fontsize='large')
a.tick_params(axis='both', labelsize='medium')

savefig('pctconflicts-hist.png')

