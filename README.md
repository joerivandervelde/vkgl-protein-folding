# Protein folding on VKGL data
![Protein structure example of the CFTR gene](img/pymol-cftr-small.png)

Deep learning has greatly increased our coverage of human protein structures<sup>[1](https://www.nature.com/articles/s41586-021-03819-2) </sup> while sophisticated protein stability prediction algorithms have become accessible to novice users on commodity hardware<sup>[2](https://doi.org/10.1093/bioinformatics/btz184) </sup>.
Decreased stability and consequent protein misfolding by DNA mutations is important mechanism for pathogenicity<sup>[3](https://pubs.acs.org/doi/10.1021/jacs.5b03743) </sup>
The potential for genome diagnostics has been shown<sup>[4](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01082-2) </sup> but ΔΔG thresholds differ greatly between proteins<sup>[4](https://doi.org/10.1186/1471-2105-10-S8-S8) </sup>.



Perform folding of GRCh37 reference proteins and VKGL variant protein.
and calculate the difference in Gibbs free energy change (ΔΔG)
increase in ΔΔG less favourable folding, pathogenic

Gene selection criteria: interesting for re-interpretation by number of VUS variants
top 1000 meeste VUS (until 12 VUS variants), 160 genes with full structure, 71 partial structure
select for monomers, enough initial benign/pathogenic variants
supplemented with a number of hand-picked highly relevant genes such as BRCA1

Example

![CFTR folding on VKGL variants](img/ddg_vkgl_CFTR.png)

Youden's J statistic placed the threshold at 1.39.
According to the data, the chance that a variant is correctly labeled 'pathogenic' above this threshold is 90% (PPV) and 68% of all pathogenic variants can be found this way (sensitivity).

Results

In total 3299 variants were folded for 55 genes
Based current limited set of 35 genes, mean DDG gene threshold: 1.48 (95%CI: 0.97-2)
which is similar to previously found 1.58 and 1.50 ([link](https://www.nature.com/articles/s41598-020-72404-w))

which genes significantly deviate from this?
select genes with enough samples withs (n >= 10 benign and n >=10 pathogenic) and PPV >= 90% (found 10: CFTR, MLH1, ATP7B, LDLR, SCN1A, FGFR2, F8, SLC12A3, MSH2, NPC1)
select for genes above or below the 95%CI of threshold (> 2 or < 0.97)
found two:
SLC12A3 (folded 11 benign and 47 pathogenic variants, threshold 0.41 with PPV 93%, NPV 47%, sens 81%  and spec 73%)
MSH2 (folded 46 benign and 18 pathogenic variants, threshold 3.66 with PPV 100%, NPV 87% NPV, sens 61%  and spec 100%)

Discussion
results seem to indicate that 1 in 5 genes has deviating DDG threshold

## Data used

### VKGL public consensus release April 2023
* Fokkema, IFAC, van der Velde, KJ, Slofstra, MK, et al. Dutch genome diagnostic laboratories accelerated and improved variant interpretation and increased accuracy by sharing data. Human Mutation. 2019; 40: 2230–2238. https://doi.org/10.1002/humu.23896
* [Direct download link](https://downloads.molgeniscloud.org/downloads/VKGL/VKGL_public_consensus_202304.tsv)

### AlphaFold2 human proteome v4
* Tunyasuvunakool, K., Adler, J., Wu, Z. et al. Highly accurate protein structure prediction for the human proteome. Nature 596, 590–596 (2021). https://doi.org/10.1038/s41586-021-03828-1
* [Direct download link](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar)

## Software

### FoldX 5.0
* Javier Delgado and others, FoldX 5.0: working with RNA, small molecules and a new graphical interface, Bioinformatics, Volume 35, Issue 20, October 2019, Pages 4168–4169. https://doi.org/10.1093/bioinformatics/btz184
* [Download under license](https://foldxsuite.crg.eu)

### R
* R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
  Copyright (C) 2023 The R Foundation for Statistical Computing
  Platform: aarch64-apple-darwin20 (64-bit)
* [Download via R-project](https://www.r-project.org)

### Java
* Oracle OpenJDK version 18.0.1
* [Download](https://www.oracle.com/java/technologies/downloads)

