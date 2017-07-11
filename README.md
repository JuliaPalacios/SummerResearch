# Undergraduate Summer Research in Statistics 2017
June 26 - September 1

*Instructor*: Julia Palacios, Sequoia Hall 141

*Group meetings*: Thursday, 11:00 am - 12:00 pm

*Individual meetings*: Monday

*Location*: Sequoia Hall 207

*Logistics*: We will have group meetings every Thursday and individual meetings every Monday. You are expected to show progress in every meeting.

*Seminars of interest*: Statistics seminar, Tuesdays 4:30PM, Sequoia Hall room 200 (refreshments provided in room 100 at 4PM).

*General instructions*: 
The main components will be reading papers, using in R with many graphics and simulations, and the analysis of real data. You can work on adapting a statistical method to analyze some data in order to answer a scientific question, develop statistical methodology to answer a scientific question or explore the statistical properties of existent methods.

You are expected to start working on your ﬁnal report from the start, so that you have a complete ﬁnal document when your research time is over (and not afterwards), so keep completing your report every week.

Keep an R-Markdown document or Jupyter notebook with your code and comment your code. Keep an updated git repository with your information. We will talk about this on June 29.

*Expectations*:
You will present your project proposal on June 29 in 10 minutes. You will describe the questions you want to answer and how you plan to answer them. You will also have a list of references you will read during the summer program.

Give two oral presentations, the first presentation of 20 minutes will be on Thursday July 20 about your progress and your final presentation of 40 minutes will be on Thursday August 31. 

Your first complete draft is due on August 25 and your final publishable paper is due on September 1st.



# Journal reading

Part of our summer research will involve reading and discussing some scientific papers. I want to emphasize the importance of being a scholar. I ask you to read at least two papers a week, some of them will be discussed during group meetings and some will be uniquely relevant for your projects. When you finish reading your paper, write a short summary of your findings and annotate them here under Description. Please edit this list as needed.


| Date | Link to paper | Reader | Description
|------|--------------|-----------|-----------|
| June 26 | [Methods and models for unravelling human evolutionary history ](https://www.nature.com/nrg/journal/v16/n12/full/nrg4005.html)| All | It will be discussed on Thursday June 29. It presents a high level overview of methods in population genetics. It does not cover all methods out there.|
| June 27 | [A metric on phylogenetic tree shapes](http://biorxiv.org/content/early/2017/02/09/054544)| Julia |Given a tree shape, the authors propose a labeling scheme of internal nodes traversing the tree from tips to root such that the label of the root uniquely identifies the tree shape. With this single "number" they construct different metrics and define addition and multiplication of tree shapes with different mappings. This labeling scheme can be useful for counting ranked tree shapes that are compatible with the data. Not clear how to "adjust" for the ranking or if ranking can come as a second step in the counting |
| June 27 | [Frequency spectrum neutrality tests: one for all and all for one](https://doi.org/10.1534/genetics.109.104042) | Alan | Generalizes SFS-based statistical estimators of mutation rates, and shows how classical tests of neutrality fit within a more general family of test statistics. Constructs an example of a new test statistic that yielded higher power than all previous statistics against evolutionary hypothesis for the lactase gene in European human populations |
| June 28 | [Chapter 1. Fast Sequential Monte Carlo Methods by Rubinstein et al.]| Julia | Just introductory chapter |
| June 28 | [LAMARC 2.0: maximum likelihood and Bayesian estimation of population parameters](https://academic.oup.com/bioinformatics/article-abstract/22/6/768/296494/LAMARC-2-0-maximum-likelihood-and-Bayesian) |Julia|It is a very short paper that describes the type of models incorporated in the software. It has limited functionals of Ne and assumes finite sites mutation model (isochronous sampling). It has Bayesian and frequentist implementations|
| June 28 |[nextflu: real-time tracking of seasonal influenza virus evolution in humans](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv381) | Rayne | Overview of the NextFlu project outlining the data pipeline used to load and clean sequence data and build a phylogenetic tree|
| June 28 |[Gaussian Process-Based Bayesian Nonparametric Inference of Population Size Trajectories from Gene Genealogies](http://onlinelibrary.wiley.com/doi/10.1111/biom.12003/full) | Rayne | Given a geneology, the paper introduces a Bayesian nonparametric estimation of the effective population size|
| June 28 | [Integrated Nested Laplace Approximation for Bayesian Nonparametric Phylodynamics](https://arxiv.org/abs/1210.4908) | Rayne | Improves upon the previous method of estimating the effective population size trajectory  using  INLA instead of MCMC. This improves both accuracy and computational efficiency of the previous method|
| June 29 | [Importance Sampling for the Infinite Sites Model](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2832804/pdf/nihms164971.pdf)|Julia|They used the Griffiths and Tavare mtDNA for testing another importance sampling method for inferring theta under the infinite sites mutation model.|
| June 29 | [Empirical Bayesian analysis of simultaneous changepoints in multiple data sequences](https://arxiv.org/abs/1508.01280)|Julia|An example of MCEM for Empirical Bayes estimation of hyperparameters in the multiple changepoints setting. The Metropolis-Hastings algorithm is very original and tailored to this problem. It exploits a dynamic programming way of computing and sampling from the full posterior distribution|
| July 3 | [Stochastic enumeration method for counting trees](https://people.smp.uq.edu.au/DirkKroese/ps/tree-counting.pdf)| Julia | |
| July 4 | [Footprints of ancient-balanced polymorphisms in genetic variation data from closely related species](http://onlinelibrary.wiley.com/doi/10.1111/evo.12567/epdf)|Julia|They derive approximations to the distribution of some statistics that are indicators of trans-species polymorphism (and signature of balancing selection)--reading incomplete. |
| July 6 | [Optimal neutrality tests based on the frequency spectrum](https://doi.org/10.1534/genetics.110.118570) | All | Building on (Achaz, 2009), paper proposes choosing weights based on fractional frequencies instead of usual frequencies, and also shows how a most optimal test can be to test against an alternative evolutionary hypothesis |
| July 6 |[Phylogenetic MCMC Algorithms Are Misleading on Mixtures of Trees](https://pdfs.semanticscholar.org/bcaf/fce47bdda89be7501c796db8ae00f8e7aed3.pdf)|Julia|They show that the Mixing time is exponential when sampling from the "wrong" posterior distribution. The result is not surprising since it is the wrong model (ignoring the fact that data comes from a mixture) but the mathematical approach is a good example for showing mixing time calculations|
| July 7 | [Non-parametric estimation of population size changes from the site frequency spectrum](http://www.biorxiv.org/content/biorxiv/early/2017/05/02/125351.full.pdf)|Julia||
| July 10 | [Compressive fluorescence microscopy for biological and hyperspectral imaging](http://www.pnas.org/content/109/26/E1679.full.pdf)| Julia | They implement compressive sensing in fluorescence microscopy|
| July 13  | [FastTree](https://academic.oup.com/mbe/article/26/7/1641/1128976/FastTree-Computing-Large-Minimum-Evolution-Trees)| All | |
| TBD | [An expanded view of complex traits](http://www.cell.com/cell/abstract/S0092-8674(17)30629-3)| All | |
| TBD | [Impacts of Neanderthal-Introgressed Sequences on the Landscape of Human Gene Expression](http://www.cell.com/cell/fulltext/S0092-8674(17)30128-9)| All | |
| TBD   | [Topological metrizations of trees, and new quartet methods of tree inference](https://arxiv.org/abs/1704.02004) | Alan | Using different metrizations to design new phylogenetic inference methods |
| TBD   | [A generalized Watterson estimator for next-generation sequencing: From trios to autopolyploids](https://doi.org/10.1016/j.tpb.2015.01.001) | Alan | Uniting new NGS-based estimators proposed previously (see first section) via a common mathematical framework |
| TBD   | [Population genetic analysis of shotgun assemblies of genomic sequences from multiple individuals](https://doi.org/10.1101/gr.074187.107) | Alan | Introducing method for obtaining estimators of nucleotide diversity \theta from shotgun sequencing data |
| TBD   | [Population genomics from pool sequencing](https://doi.org/10.1111/mec.12522) | Alan | Formulating estimators of summary statistics for data arising from NGS of pooled samples |
| TBD   | [Neutrality tests for sequences with missing data](https://doi.org/10.1534/genetics.112.139949) | Alan | SFS-based statistical estimators for incomplete genetic data |
| TBD   | [Fundamental limits on the accuracy of demographic inference based on the sample frequency spectrum](https://doi.org/10.1073/pnas.1503717112) | Alan | Information-theoretic error bound on distance between population growth function and its estimator, and slow convergence of estimators |
| TBD | [The ratio of human X chromosome to autosome diversity..](http://www.nature.com/ng/journal/v42/n10/full/ng.651.html)|Julia||
| TBD | [Analyses of X-linked and autosomal..](http://www.nature.com/ng/journal/v43/n8/full/ng.877.html)| Julia | |
| TBD |[Ascertainment bias in studies of human genome-wide polymorphism](http://genome.cshlp.org/content/15/11/1496.full.pdf+html)|Alan||


Put yourself on the calendar like this:
```
| May 16   | [paper](https://...) | Julia | A very interesting paper |
