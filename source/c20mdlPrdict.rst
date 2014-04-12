Model and Predict RBP Binding Preference
==================================================

With methods mentioned in previous chapter we could hopefully distinguish the informative signals from the sequencing noises. Results achieved at this step could be enough for most exploratory researches. However, what if we take a further step? What could we learn from the confident results? Could we learn about the binding preference of RNA binding proteins and construct models to unveil why RBP binds here but not there? Could we even predict the global binding sites of RBP? It does not necessarily means making CLIP-seq *in silico* possible, but machine learning methods could yield promising supplement to the real signals which could be left out by the empirical thresholding. 


Hidden Markov Model
-------------------------

In the beginning of this section I would like to introduce HMM by a toy application: 5' end splicing site recognition, which was originated from NBT paper [Eddy2004]_. Admittedly this toy HMM is naive, it hopefully could help you take first step knowing how to infer the most possible splicing site, simply starting from raw DNA/RNA sequence together with some prior knowledges. For now we begin with the sequence: *CTTCATGTGAAAGCAGACGTAAGTCA*. Series of ideas and concepts would be illustrated step by step so that you could see how the biological question was translated to mathematical questions. Hopefully after you have known HMM I would introduce how HMM was employed in RBP science researches. 

Sequence Path v.s. State Path
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sequence could be read by sequencing and thus observed. Each nucleotide residing on the DNA string belongs to some bio-functional regions. They are exons to be transcribed in the mature RNA, introns to be spliced out and splicing site, the boundary of exon and intron. These regions of different functions could be described as different states. 

Four basic types of nucleotides (A, T, C, G) come together to form a sequence path. Likewise, the state path is composed of "E", "I", and "5", denoting "Exon", "Intron" and "5'-end Splicing Site" respectively. For convenience we disregard 3' end splicing sites here. 

For example:: 

	Sequence Path: CTTCATGTGAAAGCAGACGTAAGTCA
	State    Path: EEEEEEEEEEEEEEEEEE5IIIIIII

Therefore the ultimate goal of locating the splicing site is to figure out the state of each nucleotide and find the one under splicing site state. Generally speaking, given the sequence path we try to infer the state path. 

Why Hidden? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At first sight whether GC-rich or AT-rich of a DNA string could be quickly determined given the sequence path, however you could barely say the state path immediately. This is where the "hidden" comes from. 

The observed sequence path alone cannot help us to reveal the unobserved state path. With additional prior knowledge, for example, the distribution of A/T/G/C in each state, we could unveil the most possible E/I/5 state path in a **top-bottom** fashion. 

Think in HMM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose we have understood the unknown dark matter in the bottom and all the prior knowledge. Let's think in HMM to generate a DNA sequence in **bottom-up** fashion. Oh, calm down, just playing roulette wheel. 

Emission Probs
""""""""""""""""""""""""""""""""""""""""""""""""""""

Let's play with intron state and under intron state we set the bases distribution as: ::

	    A   C   G   T
	I 0.4 0.1 0.1 0.4

.. figure:: _static/c20HMM_Emission_Prob1.png
	:scale: 40 %
	:align: center

	Distribution of A, C, G, T in toy intron. *(See R scripts to generate pie chart in final section in case of your curiosity)*

The **emission probabilities** under intron state could be described as the above roulette wheel. To generate a DNA string given the emission probabilities of intron, it is straightforward to implement sampling intron sequence of 26nt using ``R`` language:  

.. code-block:: r 
	
	> set.seed(1026)
	> dna.bases <- c("A", "C", "G", "T")
	> intron.emission.prob <- c(0.40, 0.10, 0.10, 0.40)
	> seq.len <- 26
	> sample.seq <- sample(x = dna.bases, prob = intron.emission.prob, 
							size = seq.len, replace = T)
	> cat(sample.seq, "\n")
	> table(sample.seq)
	A C A G A A A A A A A T G A T A A A C A G G G A A T 

You could even double-check whether this string follow your prior set distribution. I give you result of 1000nt sequence by setting ``seq.len <- 1000`` without changing other codes:  

.. code-block:: r

	> table(sample.seq)
	sample.seq
	  A   C   G   T 
	396  81 101 422

Emission probabilities of each state is dedicated to generate observed base. Don't forget we assume we are HMM now. If we clearly know state status of each base, we could immediately sample its nucleotide given the emission probabilities. Previously we set ``seq.len`` as ``26`` or ``1000`` to generate sequence of corresponding length. In the same way, we could give the individual base a nucleotide by setting ``seq.len <- 1``.  

This story under intron state could be exactly applied to exon, splicing sites states as well. Thus we have three wheels of emission probabilities under three different states. 

Transition Probs
""""""""""""""""""""""""""""""""""""""""""""""""""""

HMM has a feature: the next state is determined by ongoing state, or ongoing state is determined by previous one. For example, now we know we are under "Exon" state. When we considering generating next one, we are hesitating between whether continuing "Exon", or changing to "Splicing Site" state. Moving from one state to another state is exactly determined by **transition probabilities**. 

Though emission and transition probabilities have different biological missions: generating sequence and state path respectively. Their mathematical principle are exactly same - sampling something following multinomial model. 

Let's assume exon status has the following transition probabilities: 

.. code-block:: r

	> state.bases <- c("E", "5", "I")
	> exon.transition.prob <- c(0.9, 0.1, 0)
	> exon.transition.matrix <- matrix(exon.transition.prob, nrow = 1, 
										byrow = T)
	> colnames(exon.transition.matrix) <- state.bases
	> row.names(exon.transition.matrix) <- "E"
	> print(exon.transition.matrix)
	    E   5 I
	E 0.9 0.1 0

It means the probability of continuing exon state is 0.9 whereas changing to 5' end splicing site state is 0.1. Note that in our toy model exon cannot directly jump to intron state. 

You cannot be unfamiliar with the following codes to choose state for next base: 

.. code-block:: r

	> set.seed(1026)
	> sample.state <- sample(x = state.bases, prob = exon.transition.prob, 
							size = 1, replace = T)
	> cat(sample.state)
	E

Generate sequence
""""""""""""""""""""""""""""""""""""""""""""""""""""

So far hopefully you clear have learned: 
	
- Generating sequence path purely under one state, i.e. intron
- Determine state of next base when current base is under certain state, i.e. exon

Next before generating sequences and of course state path we need make the emission and transition probabilities to the full story: 


.. figure:: _static/HMMandSplicingSite_Mx.png
	:scale: 50 %
	:align: center

	Diagram showing emission and transition probabilities matrix [Eddy2004]_. 

In particular there are additional "Start" and "End" states for model purpose. These two states have no base assignment and they serve as the ``^`` and ``$`` in the regular expression in programming. In the figure the transition probability from "Start" to "Exon" state is ``1.0`` which means in our toy example here we are going to generate sequence starting with "Exon". Likewise the probability of continuing "Intron" is ``0.9``. 

The full emission probabilities matrix would be: 

.. code-block:: bash

		A    C    G    T
	S   0.00 0.00 0.00 0.00
	E   0.25 0.25 0.25 0.25
	SS5 0.05 0.00 0.95 0.00
	I   0.40 0.10 0.10 0.40
	N   0.00 0.00 0.00 0.00

The full transition probabilities matrix would be: 

.. code-block:: bash

	    S   E SS5   I   N
	S   0 1.0 0.0 0.0 0.0
	E   0 0.9 0.1 0.0 0.0
	SS5 0 0.0 0.0 1.0 0.0
	I   0 0.0 0.0 0.9 0.1
	N   0 0.0 0.0 0.0 1.0

Check the source R script to find ``GenerateHMMSeq`` function. After you type the following in the R console, TA-DA, a sequence with HMM in mind is generated. 

.. code-block:: r

	> set.seed(1026)
	> GenerateHMMSeq(emmisionMx=myEmisMtx, transitionMx=myTransMtx, len=20)
	Start
	Base 1 is A under E 
	Base 2 is A under E 
	Base 3 is C under E 
	Base 4 is G under E 
	Base 5 is G under E 
	Base 6 is T under E 
	Base 7 is C under E 
	Base 8 is C under E 
	Base 9 is C under E 
	Base 10 is G under SS5 
	Base 11 is G under I 
	Base 12 is A under I 
	Base 13 is T under I 
	Base 14 is T under I 
	Base 15 is T under I 
	End
	Seq Path: A A C G G T C C C G G A T T T 
	StatPath: E E E E E E E E E SS5 I I I I I N 
	[1] 1

You may noticed that in this case the length of output sequence is 15 though 20 is pre-defined. Is it wrong? No, because the probability of changing from "Intron" to "End" state is ``0.1``, not ``0``, Ongoing "Intron" state is possible to choose "End" state for next base. 

Judge Possible State Path
^^^^^^^^^^^^^^^^^^^^^^^^^^^

So far hopefully you have known how to think in HMM to generate sequence path and state path given the transition and emission probabilities matrix available. Now back to our original goal: locating 5'-end splicing site. Now out mind goes in **top-bottom** direction and our mission is to find the state path with highest probabilities given the observed sequence path. 

Compute Probability of State Path
"""""""""""""""""""""""""""""""""

Recall our example: 

.. code-block:: bash

	Sequence Path: CTTCATGTGAAAGCAGACGTAAGTCA
	State    Path: EEEEEEEEEEEEEEEEEE5IIIIIII

Given the observed sequence of 26nt, how likely do you think is the state path among the hundreds of combinations? In general :math:`{\mathbf P}\left(\mathbf{s},\pi \mid \text{HMM},\theta\right)` of an Hidden Marcov Model with parameters ( :math:`\theta` ) generates a state path ( :math:`\pi` ) and an observed sequence ( :math:`\mathbf{s}` ) , is the product of all the emission probabilities and transition probabilities that were used: 

.. math::
	
	\left( 0.25^{18} * 0.95 * 0.4^3 * 0.1 * 0.1 * 0.4^2 \right) * \left(1 * 0.9^{17} * 0.1 * 1 * 0.9^6 * 0.1 \right)

26nt sequence has 26 emissions and 27 transitions (note the final end state).

See ``ComputeStatePathProb`` function for details and it could yield: 

.. code-block:: r
	> demo.seq.path <- c("CTTCATGTGAAAGCAGACGTAAGTCA")
	> demo.state.path <- c("EEEEEEEEEEEEEEEEEE5IIIIIII")
	> ComputeStatePathProb(seqPath = demo.seq.path, statePath = demo.state.path, 
	                     transMtx = myTransMtx, emisMtx = myEmisMtx, log = TRUE)
	[1] -41.21968

Enumerate State Paths 
""""""""""""""""""""""""

.. figure:: _static/HMMandSplicingSite.png
	:scale: 50 %
	:align: center

	Most possible state path. 

Among all the possible combinations of state path, how to pick up the most possible state path? There is one naive way: enumerate all the candidates, compute their probability, and pick up the highest one. 

.. code-block:: r

	> InspectorHMM(demo.seq.path, myTransMtx, myEmisMtx)
	[1] "All State Path with non-zero Probability: "
	[1] "EEEE5IIIIIIIIIIIIIIIIIIIII"
	[1] 2.904208e-21
	[1] "EEEEEE5IIIIIIIIIIIIIIIIIII"
	[1] 8.621866e-20
	[1] "EEEEEEEE5IIIIIIIIIIIIIIIII"
	[1] 1.347167e-19
	[1] "EEEEEEEEE5IIIIIIIIIIIIIIII"
	[1] 4.431469e-21
	[1] "EEEEEEEEEE5IIIIIIIIIIIIIII"
	[1] 2.769668e-21
	[1] "EEEEEEEEEEE5IIIIIIIIIIIIII"
	[1] 1.731043e-21
	[1] "EEEEEEEEEEEE5IIIIIIIIIIIII"
	[1] 8.222452e-20
	[1] "EEEEEEEEEEEEEE5IIIIIIIIIII"
	[1] 6.761885e-21
	[1] "EEEEEEEEEEEEEEE5IIIIIIIIII"
	[1] 3.211895e-19
	[1] "EEEEEEEEEEEEEEEE5IIIIIIIII"
	[1] 1.056545e-20
	[1] "EEEEEEEEEEEEEEEEEE5IIIIIII"
	[1] 1.254647e-18
	[1] "EEEEEEEEEEEEEEEEEEEE5IIIII"
	[1] 2.579454e-20
	[1] "EEEEEEEEEEEEEEEEEEEEE5IIII"
	[1] 1.612159e-20
	[1] "EEEEEEEEEEEEEEEEEEEEEE5III"
	[1] 7.657755e-19
	Best State Path and original Sequence 
	EEEEEEEEEEEEEEEEEE5IIIIIII 
	CTTCATGTGAAAGCAGACGTAAGTCA 
	Maximum Probability: 1.25464666093437e-18 
	Maximum Prob (log):  -41.2196776860225 
	Posterior Decoding:  0.461971754327234 
	[1] "EEEEEEEEEEEEEEEEEE5IIIIIII"

So far hopefully you have understood the principles of HMM. Obviously we cannot use the enumerate methods because when the length of sequence is increasing we could barely afford the time and computing resources to do that. We do have smarter algorithms. Keep reading. 

Viterbi algorithm
""""""""""""""""""""

SVM
-------------------------

Graph Kernel
-------------------------

Random Forest
-------------------------

References
-------------------------
.. [Eddy2004] Eddy, S. R. What is a hidden Markov model? Nat. Biotechnol. 22, 1315–1316 (2004).

Source Scripts
-------------------------

.. literalinclude:: _scripts/c20HMM_Wheel.R
    :language: r