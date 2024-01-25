# Protein folding on VKGL data
![Protein structure example of the CFTR gene](img/pymol-cftr-small.png)

Coverage of human protein structures has recently greatly increased by deep learning<sup>[1](https://www.nature.com/articles/s41586-021-03819-2) </sup> while sophisticated algorithms to predict protein stability have become accessible to mainstream users on commodity hardware<sup>[2](https://doi.org/10.1093/bioinformatics/btz184) </sup>.
Protein stability may be decreased by coding DNA variation and even lead to protein misfolding, representing an important mechanism for pathogenicity<sup>[3](https://pubs.acs.org/doi/10.1021/jacs.5b03743) </sup>.
The potential for genome diagnostics to recognize and report such variation has been shown<sup>[4](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01082-2) </sup>, but the ΔΔG thresholds used for interpretion in this context differ greatly between proteins<sup>[4](https://doi.org/10.1186/1471-2105-10-S8-S8) </sup>.

Here, we perform protein folding on DNA variation from shared by Dutch genome diagnostic laboratories in the VKGL Data Sharing working group.
Essentially, we calculate the difference in Gibbs free energy change (ΔΔG) between wild-type protein sequences and variant sequences.
An increase in ΔΔG indicates that more energy is required for folding, making it less favourable and prone to pathogenic misfolding.
The amino acid changes of the DNA variation is based on GRCh37 and introduced in the AlphaFold2 human proteome.

We selected genes for which many Variants of Unknown Significance (VUS) have been reported for potential re-interpretation or otherwise having a high clinical interest.
In addition, we require a substantial amount of initial benign and pathogenic variants to increase chances of success.
Lastly, the selected genes had protein products consist of single-fragment monomers.

We calculate the ΔΔG for benign and pathogenic variants and use Youden's J statistic to estimate an optimal threshold between these two groups.
For example, for CFTR variants the threshold is placed at 1.39.
According to these data, the chance that a new variant is correctly labeled as 'pathogenic' above this threshold is 90% (positive predictive value, PPV) and 68% of all pathogenic variants can be found this way (sensitivity).

![CFTR folding on VKGL variants](img/ddg_vkgl_CFTR.png)

In total we have folded [5869 mutant proteins](out/all_folded_variants.txt) for DNA variation in [55 genes](out/ddg_vkgl_gene_results.txt).
These variants were classified as 1621 likely benign or benign (LB/B), 1678 likely pathogenic or pathogenic (LP/P), 2482 VUS and 88 conflicting, i.e. multiple classifications on the same protein change that are not identical.
The mean ΔΔG threshold across these genes was 1.48 (95%CI: 0.97-2), which is comparable to the previously found threholds of 1.58 and 1.50<sup>[5](https://www.nature.com/articles/s41598-020-72404-w) </sup>.

We defined genes to have a trustworthy threshold if enough samples were used for estimation (LB/B _n_ >= 10 and LP/P _n_ >=10) and a PPV of >= 90%.
This resulted in 10 genes: CFTR, MLH1, ATP7B, LDLR, SCN1A, FGFR2, F8, SLC12A3, MSH2, and NPC1.
We investigated the potential for re-classification of 456 VUS present in the 10 genes.
Applying the respective gene thresholds resulted in [166 VUS](out/vus_lp_candidates.txt) that might be considered candidates for re-classification as likely pathogenic:

| Gene    | Candidates for re-classification |
|---------|----------------------------------|
| ATP7B   | 17                               |
| CFTR    | 22                               |
| F8      | 14                               |
| FGFR2   | 6                                |
| LDLR    | 17                               |
| MLH1    | 26                               |
| MSH2    | 23                               |
| NPC1    | 4                                |
| SCN1A   | 21                               |
| SLC12A3 | 16                               |

Interestingly, two genes were found to have a ΔΔG above or below the 95%CI of the mean threshold (> 2 or < 0.97).
For [SLC12A3](out/ddg_vkgl_SLC12A3.pdf) a ΔΔG threshold of 0.41 based on folded structures for 11 benign and 47 pathogenic variants, having a PPV of 93%, an NPV of 47%, a sensitivity of 81% and a specificity of 73%, and
for [MSH2](out/ddg_vkgl_MSH2.pdf) a ΔΔG threshold of 3.66 based on folded structures for 46 benign and 18 pathogenic variants, having a PPV 100%, an NPV of 87%, a sensitivity of 61% and a specificity of 100%.

While these results seem to confirm the potential of protein folding for genome diagnostics, the used sample size is relatively small and AlphaFold2 structures might not accurately represent actual biological protein structures nor take into account the highly complex context in which these proteins perform their function.

## Acknowledgements
* __Rene Mulder__
* Jan D.H. Jongbloed
* Kristin M. Abbott
* Birgit Sikkema-Raddatz
* Helga Westers

## Data used

### VKGL public consensus release April 2023
* Fokkema, IFAC, van der Velde, KJ, Slofstra, MK, et al. Dutch genome diagnostic laboratories accelerated and improved variant interpretation and increased accuracy by sharing data. Human Mutation. 2019; 40: 2230–2238. https://doi.org/10.1002/humu.23896
* [Direct download link](https://downloads.molgeniscloud.org/downloads/VKGL/VKGL_public_consensus_202304.tsv)

### AlphaFold2 human proteome v4
* Tunyasuvunakool, K., Adler, J., Wu, Z. et al. Highly accurate protein structure prediction for the human proteome. Nature 596, 590–596 (2021). https://doi.org/10.1038/s41586-021-03828-1
* [Direct download link](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar)

## Software used

### FoldX 5.0
* Javier Delgado and others, FoldX 5.0: working with RNA, small molecules and a new graphical interface, Bioinformatics, Volume 35, Issue 20, October 2019, Pages 4168–4169. https://doi.org/10.1093/bioinformatics/btz184
* [Download under license](https://foldxsuite.crg.eu)

### R 4.2.3
* R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"
  Copyright (C) 2023 The R Foundation for Statistical Computing
  Platform: aarch64-apple-darwin20 (64-bit)
* [Download via R-project](https://www.r-project.org)

### Java 18.0.1
* Oracle OpenJDK version 18.0.1
* [Download](https://www.oracle.com/java/technologies/downloads)
