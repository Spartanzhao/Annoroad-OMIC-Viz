---
title: 'Annoroad Browser: A biological multi-omic visualization tool'
tags:
- Python
- bioinformatics
- multi omic
- Hi-C
- visualization
authors:
- name: Zhao Yue [co-first author] 
  orcid: 0000-0001-5293-820X
  affiliation: 1
- name: Yuan Zan [co-first author] 
  affiliation: 1
bibliography: paper.bib 
date: 18 October 2021
affiliations:
- name: Zhao Yue, Senior bioinformatic engineer , Beijing Annoroad Genomics
  index: 1
---



# Summary
Different types of high throughput sequencing such as Hi-C, ChIP-seq, ATAC-seq, Cut&Tag, transcriptome, BS-seq, WGS have revolutionized multi-omic science. In turn, these technologies also generated some challenges, such as how to deal with a large amount of data and how to do the integration analysis of different omics data, and after all, find the correlation between different omic data. Big data visualization thus is very important for researchers to discover useful clues. However, drawing a multi-omic track plot that reaches the demand for publication requirement is a big challenge. Current published tools we can use are WashU browsers,pyGenomeTrack,IGV tools and so on. So even though these tools can help us visualize omic datas, however, these tools have various disadvantages. WashU can only run a few species, IGV can't draw the heatmap, pyGenomeTrack is unfriendly to use. All of these tools can not compare differences between two sample's Hi-C data and epigenome or RNA data. Further more, these tools can only draw one chromosome area at the same time. If we have so many interesting areas on different chromosomes, then it will be a huge work to draw these pictures one by one.  Thus, we developed Annoroad Browser, a command line based tool to solve this issue. Annoroad Browser is a powerful and flexible tool that perfectly meets the demand of researchers who wants to do integration, visualization, and analysis of epigenetic,3D genomic, methylome WGS and transcriptome data. Besides, our tools can easily combine different tracks in any order according to user needs and generate a vector graph that fits the requirement for publication. Moreover, the picture we generated is highly reproducible. Except for heat map, all the other data tracks can change the color. Users don’t need to install our package. It only needs python3.6 or above. Besides, the required packages are pandas, numpy, scipy and matplotlib. 
Among 3D genomics, there are different kinds of data we can perform such as contact matrix, A/B compartment, delta z score matrix between two samples, insulating score, TAD boundary, loops, 4C display. In addition, we can display different kinds of epigenome and transcriptome data, and calculate difference between two data then show on the track. We support to show deletion, duplication, inversion and insertion from WGS on the track. Based on this handy and powerful tool, we can face the challenge in multi-omic data analysis. For example, if we draw a track plot, someplace on the chromosome has a deletion, then the heat map support there is a TAD fusion, another sample has no such deletion and thus has no TAD fusion. On the track of H3K27ac, we find there is a peak in one side of deletion, and after TAD fusion, the H3K27ac peak has a strong interaction with a gene on the other side of deletion, which formed enhancer hijacking. But the other sample has on such phenomenon. We can further identify it by using delta zscore heatmap in our tool. The transcriptome data showed gene expression is dramatically enhanced after enhancer hijacking. Thus based on our tool, we can easily find such interesting phenomena very quickly. 

# Statement of need
 Annoroad Browser is a python3 based tool. It can perform multi-omic visualization.  It can display, Hi-C heatmaps. delta-HiC of two samples, loops, AB compartment, TAD boundary, structure variation, epigenomic data, gene track, gene expression and log2 fold change of two samples. Besides, users can also flexibly determine the display sequence of their omic data. 

# Features
Data that Annoroad Browser can show:

- `Contact matrix` : 3D genome, especially Hi-C data is a N by N matrix, this matrix can be raw matrix or ice normalized matrix

- `Zscore Delta Matrix` : Some times we want to compare the difference of two matrix, this funciont can first do the zscore of two matrices and then do the subtraction

- `Gene track`: This function can draw the genes on its location of the chromosome which include their exons and introns

- `Bedgraph` : (Epigenomic data or transcriptome data) epigenome or transcriptome data can be displayed by this function

- `Loop`: many softwares such as HiCCUPs can find out the loop from chromomse interaction matrix. We can show the loops on the chromosome.

- `TAD boundary`: chromsome was composed by many different TADs(topological associated domain), thus we want to use this function to perform the TAD boundary

- `Structure Variation`: structurue variation include deletion, duplication, inservion, inversion. Sometimes structure variation can trigger the 3D alteration of chromosome. This function can show each type's structure varation's postion and thus, we can find relationship between 3D structure and structure variation.

- `A/B compartment`: A/B compartment is a very import concept of 3D genome. Chromosme usually have two different kinds of compartment, A compartment usually have higher gene expresion level. B compartment have lower gene expression level.

- `4C show`: 4C means one to many. If you only interest in one postion and want to know how this position interact with other areas, then you can use this function to find out.

- `Signal compare`: calculate the log2 fold change of two sample's epigenome data or transcriptome data.

Annoroad Browser can display different omic data in one picture, thus we can easier to figure out the association between different omic data. Detailed instructions can found in https://github.com/Spartanzhao/Annoroad-Browser/.

As the shown below, FigureA shows three sample's chromosome interaction heatmap first, top two are from cancer cell line, We can clearly see deletions compared with normal cell line at the bottom. Below the heatmap, is the 4C display of the MIR31HG. PANC1 has much stronger interaction on the other side of deletion compared with HPDE6C7. H3K27ac 


<p align="left">
<img src="./result/MIR31HG_2.png?raw=true" width="300" height="800"><img src="./result/BxPC3.SV.png?raw=true" width="300" height="800"><img src="./result/BxPC3.diff_gene.png?raw=true" width="300" height="800">
</p>

# Acknowledgements
I want to Acknowledge my colleagues in Annoroad Genomics, Fan Xuning helped me to improve the speed of heatmap drawing, Yuan Zan made my heatmap's color scheme more beautiful than before. I also want appreciate the dedicated people of these python packages which include pandas, matplotlib, numpy, scipy.

# References
