---
title: 'Annoroad Browser: A biological multi-omic visualization tool'
tags:
- Python
- bioinformatics
- multi omic
- Hi-C
- visualization
authors:
- name: Zhao Yue^[co-first author] 
    orcid: 0000-0001-5293-820X
    affiliation: 1
- name: Yuan Zan^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 1
affiliations:
- name: Zhao Yue, Senior bioinformatic engineer , Beijing Annoroad Genomics
   index: 1
---
date: 18 October 2021
bibliography: paper.bib

# Summary
Different types of high throughput sequencing such as Hi-C, ChIP-seq, ATAC-seq, Cut&Tag, transcriptome, BS-seq, WGS have revolutionized multi-omic science. In turn, these technologies also generated some challenges, such as how to deal with a large amount of data and how to do the integration analysis of different omics data, and after all, find the correlation between different omic data. Big data visualization thus is very important for researchers to discover useful clues. However, drawing a multi-omic track plot that reaches the demand for publication requirement is a big challenge. Thus, we developed Annoroad Browser, a command line based tool to solve this issue. Annoroad Browser is a powerful and flexible tool that perfectly meets the demand of researchers who wants to do integration, visualization, and analysis of epigenetic,3D genomic, methylome WGS and transcriptome data. Besides, our tools can easily combine different tracks in any order according to user needs and generate a vector graph that fits the requirement for publication. Moreover, the picture we generated is highly reproducible. Except for heat map, all the other data tracks can change the color. Users don’t need to install our package. It only needs python3.6 or above. Besides, the required packages are pandas, numpy, scipy and matplotlib. 
Among 3D genomics, there are different kinds of data we can perform such as contact matrix, A/B compartment, delta z score matrix between two samples, insulating score, TAD boundary, loops, 4C display. In addition, we can display different kinds of epigenome and transcriptome data, and calculate difference between two data then show on the track. We support to show deletion, duplication, inversion and insertion from WGS on the track. Based on this handy and powerful tool, we can face the challenge in multi-omic data analysis. For example, if we draw a track plot, someplace on the chromosome has a deletion, then the heat map support there is a TAD fusion, another sample has no such deletion and thus has no TAD fusion. On the track of H3K27ac, we find there is a peak in one side of deletion, and after TAD fusion, the H3K27ac peak has a strong interaction with a gene on the other side of deletion, which formed enhancer hijacking. But the other sample has on such phenomenon. We can further identify it by using delta zscore heatmap in our tool. The transcriptome data showed gene expression is dramatically enhanced after enhancer hijacking. Thus based on our tool, we can easily find such interesting phenomena very quickly. 

# Statement of need
 Annoroad Browser is a python3 based tool. It can perform multi-omic visualization.  It can display, Hi-C heatmaps. delta-HiC of two samples, loops, AB compartment, TAD boundary, structure variation, epigenomic data, gene track, gene expression and log2 fold change of two samples. Besides, users can also flexibly determine the display sequence of their omic data. 
