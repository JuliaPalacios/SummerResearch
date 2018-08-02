# Undergraduate Summer Research in Statistics 2018
June 25 - August 17, 40hrs per week

*Instructor*: Julia Palacios, Sequoia Hall 141

*Group meetings*: Thursday, 11:00 am - 12:00 pm

*Individual meetings*: Monday

*Location*: Sequoia Hall 207

*Logistics*: We will have group meetings every Thursday and individual meetings every Monday. You are expected to show progress in every meeting.

*Seminars of interest*: Statistics seminar, Tuesdays 4:30PM, Sequoia Hall room 200 (refreshments provided in room 100 at 4PM).

*General instructions*: 
The main components will be reading papers, use R for graphics and simulations and the analysis of real data. You can work on adapting a statistical method to analyze some data in order to answer a scientific question, develop statistical methodology to answer a scientific question or explore the statistical properties of existent methods.

You are expected to start working on your ﬁnal report from the start, so that you have a complete ﬁnal document when your research time is over (and not afterwards), so keep completing your report every week.

Keep an R-Markdown document or Jupyter notebook with your code and comment your code. Keep an updated git repository with your information. We will talk about this on June 29.

*Expectations*:
You will present your project proposal on June 29 in 10 minutes. You will describe the questions you want to answer and how you plan to answer them. You will also have a list of references you will read during the summer program.

Give two oral presentations, the first presentation of 20 minutes will be on Thursday July 27 about your progress and your final presentation of 40 minutes will be on Thursday August 16. 

Your first complete draft is due on August 10 and your final publishable paper is due on August 16.



# Journal reading

Part of our summer research will involve reading and discussing some scientific papers. I want to emphasize the importance of being a scholar. I ask you to read at least two papers a week, some of them will be discussed during group meetings and some will be uniquely relevant for your projects. When you finish reading your paper, write a short summary of your findings and annotate them here under Description. Please edit this list as needed.


| Date | Link to paper | Reader | Description
|------|--------------|-----------|-----------|
| June 27 | [Stochastic models and descriptive statistics for phylogenetic trees, from Yule to today](https://projecteuclid.org/euclid.ss/998929474)| Kyle | This paper was a great introduction to trees and stochastic processes, but I found it also quite complicated and often confusing. Yule's 1924 paper proposes the idea of Yule processes and observes that phylogenetic trees are unbalanced, with a long-tailed distribution. Aldous also spoke of the branching process, based on a tree with shape (m,i), where m = size of parent clade and s = size of smaller daughter clade. Additionally, Aldous touches on selection and systematic biases as well as adaptive radiation. Aldous claims that more elaborate models are necessary, as numbers of species don't tend to increase or decrease very quickly, except in cases of mass extinction, which tend to be completely random. In summary, Aldous concludes that phylogenetic trees have observed imbalance. If macroevolutionary processes are to be amenable to mathematical modeling, we will have to continue refining our methods and metrics to model such processes. The Yule process is used to model macroevolutionary events . Confusions/questions include the snakes table on page 24, the clades on pg. 27, the idea of adaptive radiation as well as the beta-splitting model (Appendix). |
|July 5| [Ranked Tree Shapes, Nonrandom Extinctions, and the Loss of Phylogenetic Diversity](https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syy030/4972608)| Kyle | Phylogenetic diversity (PD) is a measure of the evolutionary legacy of a group of species, which can be used to define conservation priorities. It has been shown that an important loss of species diversity can sometimes lead to a much less important loss of PD, depending on the topology of the species tree and on the distribution of its branch lengths. Amaury et al. (2017) described the concept of ranked tree shapes, an extension of the Beta-splitting model described by Aldous (2001), except with two additional parameters: parameter a, which quantifies the correlation between the richness of a clade and its relative depth, and a parameter η, which quantifies the correlation between the clade and its frequency. The loss of PD is dependent on i) ranked shape of the clade tree, and ii) non-random extinctions. Further research needs to be done on metrics of these ranked tree shapes to reach further conclusions w.r.t. population genomics as well as evolutionary dynamics. |
| June 12   | [A Metric on Phylogenetic Tree Shapes](https://academic.oup.com/sysbio/article/67/1/113/3788885) | All (Kyle will present) | Kyle: A Metric on Phylogenetic Tree Shapes is written by Coljin and focuses on designing metrics, particularly of distance, on phylogenetic tree shapes. Coljin starts from first principles, explaining the labeling scheme and definitions, before presenting several distance functions, over the space of rooted unlabeled shapes. Coljin first builds functions for trees with same size, and then later incorporates tree size, branch lengths and other properties. Data is simulated on seasonal vs. tropical influenza, as the resulting trees have significant difference due to seasonal effects. As with many other papers I've read, issues in detecting difference ensue upon looking at larger sample sizes. |
| June 12   | [Statistics for Phylogenetic Trees](http://pages.stat.wisc.edu/~larget/phylogeny/holmes-thPopBio-2003.pdf) | All | Kyle: In this paper, Susan Holmes writes about the challenges inherent in solving the problems of inference on tree space, and addresses some of the possible ways to deal with these problems. Estimations for trees have been contested for a great amount of time, as trees provide an excellent unit of summary information. Phylogenetic trees have both a discreteness and a complexity that justify recourse to a series of tools from statistics and geometry. Basic information is coded into distances, namely the expected number of mutations. Distance-based methods are semiparametric, on the spectrum of nonparametric to parametric models. |
| June 12   | [Computational Tools for Evaluating Phylogenetic and Hierarchical Clustering Trees](https://arxiv.org/pdf/1006.1015.pdf) | All (Kuhan will present) | Kyle: This paper presents an implementation of the geodesic distance between trees developed by Billera, Holmes, and Vogtmann (2001). It also describes tools for estimating phylogenetic trees from Data. Holmes evaluates trees by comparing two hierarchical clusterings through creating a weighted matching score from the matrix of matchings and by implementing the polynomial time geodesic path algorithm. The goal is to get the geodesic distance between two trees. Very interesting read.  |



Put yourself on the calendar like this:

 ```
 | May 16   | [paper](https://...) | Julia | A very interesting paper |


