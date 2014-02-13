# pH titration experiments

## Experiments

### 26-point pH titration 2.6-7.6 with 30 mM (?) bosutinib

Used pH buffer table from Sigma-Aldrich:
http://www.sigmaaldrich.com/life-science/core-bioreagents/biological-buffers/learning-center/buffer-reference-center.html#citric

Note that we had 0.1M citric acid and 0.1M dibasic sodium phosphate stock solutions.

Assay volume 100 uL.  Wells 1-26 (Tecanese) go from pH 2.6 to 7.6 and contain DMSO control dispensed to equal DMSO content in other wells.  Wells 96 down to 71 contain inversion symmetry version of pH titration (96 => pH 2.6) and 30 uM bosutinib.

Code in `robots` github repo, `robots/pH-array/create-pH-array-fast.py`.

#### Data
 
* `2014-02-12 09-43-41_plate_1.xml` - absorbance, fluorescence (top; bottom); significant evaporation observed during course of 4h experiment; ambient temperature

