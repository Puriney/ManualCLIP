Model and Predict RBP Binding Preference
==================================================

Hidden Markov Model
-------------------------

Toy model: HMM and Splicing Site Recognition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section I would like to introduce a toy application of HMM: splicing site recognition which is a paper from Sean R Eddy so that you will not get confused by the literal definitions. Moreover, you will actually learn how to hack around HMM since I would provide `R` code to reproduce his paper's results. 

Our mission is to locate the most possible 5'-end splicing site given a DNA string. What we have clearly known so far are two things. One is the raw sequence you have observed, `CTTCATGTGAAAGCAGACGTAAGTCA` for example. The other is the fact that sequence comprises of three categories of nucleotides: exon, intron, splicing site. 

We also know the distribution of the `A`, `C`, `G` and `T` is different among exons, introns and splicing site. Though the composition would be estimated from the genome annotations, let's just set them as followed: 

- Exon: uniform base composition on average, 0.25 each base; 
- Intron: A and T are 0.4 each, G and C are 0.1 each; 
- 5'ss : A is 0.05, G is 0.95, others are zero. 






.. image:: _static/HMMandSplicingSite.png
	:scale: 50 %
	:align: center

(img credit: Sean R Eddy)



SVM
-------------------------

Graph Kernel
-------------------------

Random Forest
-------------------------