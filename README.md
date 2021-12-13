# PROXYANC 
PROXYANC has two novel algorithms including the correlation between observed linkage disequilibrium in an admixed population and population genetic differentiation in ancestral populations, and an optimal quadratic programming based on the linear combination of population genetic distances (FST). PROXYANC was evaluated against other methods, such as the f3 statistic using a simulated 5-way admixed population as well as real data for the local South African Coloured (SAC) population, which is also 5-way admixed. The simulation results showed that PROXYANC was a significant improvement on existing methods for multi-way admixed populations. PROXYANC implements an approach to select the best proxy ancestral populations for admixed populations. It searches for the best combination of reference populations that can minimize the genetic distance between the admixed population and all possible synthetic populations, consisting of a linear combination from reference populations. PROXYANC also computes a proxy-ancestry score by regressing a statistic for LD (at short distance &lt; 0.25 Morgan) between a pair of SNPs in the admixed population against a weighted ancestral allele frequency differentiation. Download PROXYANC. PROXYANANC can select AIMs based on the relationship between the observed local multi-locus linkage disequilibrium in a recently admixed population and ancestral population difference in allele frequency and based on the Kernel principal component analysis (Kernel-PCA), which is the extension of the linear PCA. PROXYANC can identify possible unusual difference in allele frequency between pair-wise populations, as signal of natural selection. PROXYANC compute the expected maximum admixture LD from proxy ancestral populations of the admixed population. PROXYANC compute population pair-wise Fst (Genetic distance).
# PROXYANC program is free software and was developed by Dr. Emile R. Chimusa: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  GNU General Public License for more details. 
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

# Implemetation in PROXYANC:
--------------------------
1. PROXYANC implements an approach to select the best proxy ancestral populations for admixed populations. It searches for the best combination of reference populations that can minimize the genetic distance between the admixed population and all possible synthetic populations, consisting of a linear combination from reference populations. PROXYANC also computes a proxy-ancestry score by regressing a statistic for LD (at short distance < 0.25 Morgan) between a pair of SNPs in the admixed population against a weighted ancestral allele frequency differentiation.

2. PROXYANANC can select AIMs based on the relationship between the observed local multi-locus linkage disequilibrium in a recently admixed population and ancestral population difference in allele frequency and based on the Kernel principal component analysis (Kernel-PCA), which is the extension of the linear PCA.
3. PROXYANC can identifies possible unusual difference in allele frequency between pair-wise popualtions, as signal of natural selection.
4. PROXYANC compute the expected maximun admixture LD from proxy ancestral populations of the admixed population.
*********************Contact Emile R. Chimusa (echimusa@gmail.com).******************************************
To run the PROXYANC, do:
using the parameter file:
    ./PROXYANC.02.2013.py paranc.txt 

An example for a fake South African Coloured simulation over 3964 SNP on chromosome 1 is available by calling:
    ./PROXYANC.02.2013.py paranc.txt
Note: this call will produce output files based on the method used and one method can be specified parameter file at time.
Information about each of the parameter in paranc.txt follows.
--------------------------------------------------------------------------------
popDD --> Path to where are different folders of candidate populations data, the number of these folder should be equal to subpop parameter (see example parameter file, param.txt).
admixD --> Path to GENOTYPE data of the admixed population or population under study (see example parameter file, param.txt).
SNPd --> Path to SNP data of the admixed population or population under study (see example parameter file, param.txt).
windows --> An integer specifies the Windows side to compute the LD (default = 100).
subSNP -->  An integer speicifies the number of random subset of SNPs to be chosen for the analysis (Default = 500).
subpop -->  An integer specifies the number of admixture way or number of ancestry component in the admixed populations should equal to the number of groups (gathering in diffrent folder) of populations (default = 5).
maf -->  Float specifies the cutoff of Minor allele frequency to be considered (Default = 0.005)
LDcutoff --> Float specifies the cutoff of LD (Default = 0.5)
relatedInd --> Boolean (YES/NO) specifies whether the population samples under studies are related or unrelated (Default = related (YES)).
Fstscore --> Boolean (YES/NO) specifies whether the Fst Optimal Cone programming method to select best proxy ancestry. The Output of this method is "OptFST.out".
Frqscore --> Boolean (YES/NO) specifies whether the proxy-ancestry score method, to select best proxy ancestry. The Output of this method is "proxyanc.out".
pwFst --> Boolean (YES/NO) specifies whether to write popualtion pair-wise Fst, and the output is "PairwiseFST.out".
Ufrq --> Boolean (YES/NO) specifies whether to investigate popualtion pair-wise unusual difference in allele frequency as signal of natual selection. The output is "labelpop1-labelpop2.Unsual.Frq.out"
writeLD --> Boolean (YES/NO) specifies whether to write the LD in a file.
WeightedLD --> Boolean (YES/NO) specifies whether to consider a weighted model of LD
aims --> Boolean (YES/NO) specifies whether to select AIMS panel of the admixed population based on each single best proxy ancestral from each group. The file is "AIMS.score"
PCAaims --> Boolean (YES/NO) specifies whether to select AIMS panel using a Kernel PCA method (module not yet included)
Sampling --> An integer specifies bootstrap sampling (Default = 1000)
cMspace --> A float specifies physical space between adjant markers to limit the background LD (Default = 0.2 cM).
no-log -->  Boolean (YES/NO) specifies whether -log10(p) to be applied on the value in manhattan from unusual difference in allele frequencies (active if option Ufrq if TRUE) method.
colors --> Boolean (YES/NO) specifies the cycle through colors in manhattan from unusual difference in allele frequencies (active if option Ufrq if TRUE) method.
lines --> Boolean (YES/NO) whether manhattan to be lines mode (not recommanded)
Lexp --> Boolean (YES/NO) specifies whether to examine Admixed LD based on the expected maxiumun admixture LD from pair-wise proxy ancestral populations. The output is "labelpop1-labelpop2.exp"
outp --> Path to where direct all outputs.


#Information about inputs data.
--------------------------------------------------------------------------------
data-file:
  An example data file is in ~/PROXYANC/REF/ for candidate groups of proxy ancestral pops and data of fake 5-way admixed popuation in the same directory as the program PROXYANC.
All data inputs file must be in EIGENSTRAT format.

snp-file:
  A Reich lab format SNP file specifying the genetic map to be used.  The
  map distance is **required to be in Morgans.**
