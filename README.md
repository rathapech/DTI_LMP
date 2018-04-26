# Drug-Target-Interaction-based-on-Low-Rank-Matrix-Projection (LMP)
The codes for predicting interactions between drugs and targets based on low-rank matrix projection on heterogeneous biological data. 

The code is written by Ratha Pech. For comments, suggestions and questions, please contact me at ratha.pech(at)gmail.com 

Please cite the paper if you use this code: 
"A generalized method toward drug-target interaction prediction via low-rank matrix projection "

There are 5 benchmark datasets in this implementation:
1. Matador: contains only interaction information
2. Enzyme 
3. Ion Channel
4. GPCR
5. Nuclear receptors

The enzyme, ion channel, GPCR and nuclear receptors contain three types of information 
e.g.,interaction information, drug similarity and target similarity  

Using the above the datasets, please cite the original papers. 

For Matador, please cite:
Günther S, Kuhn M, Dunkel M, Campillos M, Senger C, Petsalaki E, Ahmed J, Urdiales EG, Gewiess A, Jensen LJ, 
Schneider R, Skoblo R, Russell RB, Bourne PE, Bork P, Preissner R.
SuperTarget and Matador: resources for exploring drug-target relationships.
Nucleic Acids Res. 2008 Jan;36(Database issue):D919-22. Epub 2007 Oct 16. 

For the other three datasets, please cite:
Yamanishi,Y., Araki,M., Gutteridge,A., Honda,W. and Kanehisa,M.
Prediction of drug-target interaction networks from the integration of chemical and genomic spaces 
 
**We include all the data here just for the sake of completeness.  

**To run the code: 
main.m   	(to predict the interactions)


**Using the Augmented Lagrange Multiplier (ALM) algorithm, please cite the original paper:
 
Lin, Zhouchen, Minming Chen, and Yi Ma. "The augmented lagrange multiplier method for exact recovery of
corrupted low-rank matrices." arXiv preprint arXiv:1009.5055 (2010).

