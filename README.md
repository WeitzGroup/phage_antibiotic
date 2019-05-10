Code for: Quantitative Models of Phage-Antibiotic Combination Therapy
=============================

This manuscript develops and analyzes a novel model of phage-antibiotic combination therapy, specifically adapted to an *in vivo* context. The objective is to explore the Turner lab's development and clinical application of combination therapy utilizing phage that target antibiotic efflux pumps in *Pseudomonas aeruginosa* (Chan et al., Sci Reports 2016).  The mathematical model we develop leverages our findings of "immunophage synergy" in P.a. (Leung & Weitz, JTB 2017; Roach, Leung, et al., Cell Host Microbe 2017).  We adapted this model specifically to the case of phage-antibiotic treatment of pathogens that have co-linked variation in their susceptibility to phage and antibiotics. In doing so, the manuscript addresses three key questions: 

(i)	How robust is combination therapy to variation in the resistance profiles of pathogens?
(ii)	What is the role of immune responses in shaping therapeutic outcomes?
(iii)	What levels of phage and antibiotics are necessary for curative success?

As we show, combination therapy outperforms either phage or antibiotic alone and therapeutic effectiveness is enhanced given interaction with innate immune responses. Notably, therapeutic success can be achieved even at sub-inhibitory concentrations of antibiotic such as ciprofloxacin. These in-silico findings provide further support to the nascent application of combination therapy to treat MDR bacterial infections, while highlighting the role of innate immunity in shaping therapeutic outcomes.

Code Usage
=========================
All the scripts used to generate article's figures are written in MATLAB (R2018b). Furthermore, the code was reviewed internally using MATLAB (R2019a).

1. Run the script 'Gen_all_figs.m' outside the 'code_comb_therapy' directory to generate the figures from the article as well as supplementary figures. Moreover, the figures will be automatically saved in 'fig' and 'epsc' format inside a 'figures' directory.

'code_comb_therapy' directory contains all the scripts (.m) and functions necessary to generate the article's figures.
