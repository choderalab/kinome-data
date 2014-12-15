% title: Kinases
% subtitle: or...
% subtitle: Stuff I learned while writing Fellowship Applications.
% author: Sonya Hanson
% favicon: figures/ATP_1.png

---

<div>
<img height="350" class="center" src=figures/Fellowships.png />
</div>

---
title: Kinases!

Categories of things that were learned:

- What are disease-relevant kinases?
- What are the drugs that target them?
- What are the resistance mutants that develop?
- What are some things we need to worry about with specific kinases?
- What other significant computational work has been done on kinases?

---
subtitle:"Recent years have observed the emergence of the kinase family (Ser/Thr and Tyr kinase) as one of the most intensively pursued target classes because of its intimate involvement in oncogenic signal transduction pathways that present multiple physiological responses, tumour cell proliferation and cell survival." 
class: dark 

-Blanc, Geney, & Menet: Anti-Cancer Agents in Medicinal Chemistry (2013)  <a href="http://dx.doi.org/10.2174/1871520611313050008">DOI</a>. 


---
title: Scientifically why are we interested in kinases?

- Understanding selectivity is tremendously medically relevant.

- The relationship of conformational heterogeneity to ligand binding affinity
	* Our ability to probe this experimentally with a fortuitously fluorescent set of inhibitors
	* The tractibility of being able to address this computationally due to the number and similarity of structures

<div>
<img height="250" class="center" src=figures/lapatinib.png />
</div>

---
title: Kinases: Structure

<div>
<img height="500" class="center" src=figures/Labeled-kinase.png />
</div>

---
title: Kinases: Problems with structures

<img height="500" src=figures/Src_Abl_Ima.png />

---
title: Kinases: Families

<div>
<img height="500" class="center" src=figures/kinome_structures.png />
</div>

---
title: Disease-relevant kinases

<div style="float:left">

<p>  &emsp;   &emsp;   &emsp; <b> Breast Cancer </b> </p>
<p>  &emsp;   &emsp;   &emsp; - PI3K pathway (RAS/RAF/MEK/ERK)</p>
<p>  &emsp;   &emsp;   &emsp; - HER2 </p>
<p>  &emsp;   &emsp;   &emsp; <b> Leukemia and Lymphoma </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Abl </p>
<p>  &emsp;   &emsp;   &emsp; - Src </p>
<p>  &emsp;   &emsp;   &emsp; - BTK  </p>
<p>  &emsp;   &emsp;   &emsp; - SYK </p>
<p>  &emsp;   &emsp;   &emsp; <b> Gastrointestinal Stromal Tumors </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Kit </p>

</div>
<div style="float:right"><img height="500" class="right" src=figures/Families_circled.png /></div>
<div style="clear:both"/>

---
title: This PI3K pathway.

<div>
<img height="325" class="center" src=figures/her2.png />
</div>


---
title: What are some FDA approved kinase inhibitors?


<div style="float:left">

<p></p>
<p>  &emsp;   &emsp;   &emsp; - Imatinib  : <b> Abl/cKit </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Dasatinib : <b> Abl </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Ponatinib : <b> Abl </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Gefitinib : <b> EGFR </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Erlotinib : <b> EGFR </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Neratinib : <b> EGFR/HER2 </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Lapatinib : <b> HER2 </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Sorafenib : <b> VEGFR/BRAF </b> </p>
<p>  &emsp;   &emsp;   &emsp; - Ibrutinib : <b> BTK </b> </p>

</div>
<div style="float:right"><img height="500" class="right" src=figures/kinase_inhibitor_dorm_room_poster.png /></div>
<div style="clear:both"/>

---
title: What are the categories of kinase inhibitors?

<div>
<img height="150" class="center" src=figures/types.png />
</div>

<div>
<img height="325" class="center" src=figures/2HW0-covalent.png />
covalent Src inhibitor (2HW0)
</div>

<footer class="source"> Blanc, Geney, & Menet: Anti-Cancer Agents in Medicinal Chemistry (2013).  <a href="http://dx.doi.org/10.2174/1871520611313050008">DOI</a>. </footer>

<!---
title: What are some other categories of kinase inhibitors?

- Fluorescent inhibitors
- Allosteric inhibitors
- non-competitive inhibitors
- covalent inhibitors: ibrutinib (BTK), neratinib (HER2)
-->

---
title: Type I vs. Type II inhibitors

<div>
<img height="450" class="center" src=figures/type1vtype2.png />
</div>

<footer class="source"> Blanc, Geney, & Menet: Anti-Cancer Agents in Medicinal Chemistry (2013). <a href="http://dx.doi.org/10.2174/1871520611313050008">DOI</a>. </footer>

---
title: Fluorescent inhibitors

<div>
<img height="500" class="center" src=figures/DQA_compounds.png />
</div>


---
title: Type III or Allosteric inhibitors
subtitle: 4, 6-disubstituted pyrimidines (GNF-2 and GNF-5)

Abl crystal structure (3K5V)
<div>
<img height="425" class="center" src=figures/pyridines-2.png /> 
</div>

<footer class="source"> Zhang et al: Nature (2010). <a href="http://dx.doi.org/10.1038/nature08675">DOI</a>. </footer>


---
title: Where do we see resistance to TKI's develop?

<div>
<img height="500" class="center" src=figures/res_muts-2.png />
</div>

---
title: Ponatinib overcomes the T315 I gatekeeper mutant.

<div>
<img height="300" class="center" src=figures/ponatinib.png />
</div>

<footer class="source"> Blanc, Geney, & Menet: Anti-Cancer Agents in Medicinal Chemistry (2013). <a href="http://dx.doi.org/10.2174/1871520611313050008">DOI</a>. </footer>

---
title: But even ponatinib is not immune to resistance.

<div>
<img height="450" class="center" src=figures/ponres-0.png />
</div>

<footer class="source"> Zabriskie et al: Cancer Cell (2014). <a href="http://dx.doi.org/10.1016/j.ccr.2014.07.006">DOI</a>.</footer>

---
title: Computational work.
subtitle: We are not alone.
class: segue dark nobackground

---
title: DE Shaw - Dasatinib and Src
class: img-top-center

<div>
<video id="sampleMovie" class="center" src="figures/shaw-dasatanib-2.mov" loop=\"true\ autoPlay=\"true\  width="512" height="384"></video>
</div>

<footer class="source"> Shan et al: J. Am. Chem. Soc. (2011). <a href="http://dx.doi.org/10.1021/ja202726y">DOI</a>. </footer>
---
title: The  story of Src and Abl
class: img-top-center

<div style="float:left"><img height="500" style="PADDING-LEFT: 5px" class="right" src=figures/Src-IMA-xtal-ligplot-vert.png />Src + Imatinib (2OIQ)</div>
<div style="float:right"><img height="500" class="right" src=figures/Abl-IMA-xtal-ligplot-vert.png />Abl + Imatinib (1IEP)</div>
<div style="clear:both"/>

---
title: Comparing free energy calcs using crystal structures.
class: img-top-center

<div>
<img height="180" class="center" src=figures/delG-Aleksandrov.png />
</div>

<footer class="source"> Aleksandrov & Simonson: J. Biol. Chem (2010). <a href="http://dx.doi.org/10.1074/jbc.M110.109660">DOI</a>. </footer>

---
title: Free Energy Surfaces for the DFG flip from metadynamics sims.
class: img-top-center

<div>
<img height="250" class="center" src=figures/src-abl.png />
</div>

<footer class="source"> Lovera et al: J. Am. Chem. Soc. (2012). <a href="http://dx.doi.org/10.1021/ja210751t">DOI</a>. </footer>

---
title: Explaining why Gleevec is a specific and potent inhibitor of Abl kinase - Benoit Roux
class: img-top-center

<div>
<img height="400" class="center" src=figures/Roux_PNAS.png />
</div>

<footer class="source"> Lin et al: PNAS (2013). <a href="http://dx.doi.org/10.1073/pnas.1214330110">DOI</a>.</footer>

---
title: Where can we really contribute?

- Investigating further fluorescent compounds. (DQA scaffolds are common - are they all characterized?)
- Do kinases with similar compound affinity have similar conformational distributions?
- How do resistance mutants shift conformational distributions?
- How do allosteric drugs affect the binding site?
- Other complications: unstructured loops, ordered waters, pH dependence.

---
title: Potential Issues
class: img-top-center

<div style="float:left"><img height="350" style="PADDING-LEFT: 5px" class="right" src=figures/MEK-ring.png />A MEK model</div>
<div style="float:right"><img height="350" class="right" src=figures/ATP-MG.png />Magnesium in the ATP binding site</div>
<div style="clear:both"/>

---
title: Questions? / Discussion.
class: img-top-center

<div>
<video id="sampleMovie" class="center" src="figures/shaw-dasatanib-2.mov" loop=\"true\ autoPlay=\"true\  width="512" height="384"></video>
</div>
