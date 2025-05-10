# OpenTn5: Open-Source Resource for Robust and Scalable Tn5 Transposase Purification and Characterization
**Jan Soroczynski<sup>1</sup>, Lauren J. Anderson<sup>&ast;1</sup>, Joanna L. Yeung<sup>&ast;1</sup>, Justin Rendleman<sup>1</sup>, Deena A. Oren<sup>2</sup>, Hide A. Konishi<sup>3</sup>, and Viviana I. Risca<sup>‡1</sup>**  
1. Laboratory of Genome Architecture and Dynamics, The Rockefeller University, New York, NY
2. Structural Biology Resource Center, The Rockefeller University, New York, NY
3. Laboratory of Chromosome and Cell Biology, The Rockefeller University, New York, NY  
&ast; These authors contributed equally.  
‡ To whom correspondence should be addressed. Email: vrisca@rockefeller.edu
## Abstract: 
Tagmentation streamlines sequencing library preparation by leveraging the bacterial cut-and-paste Tn5 transposase for DNA fragmentation and sequencing adapter addition. Here we present an open-source protocol for the generation of multi-purpose hyperactive Tn5 transposase, pG-Tn5E54K, L372P, yielding multi-milligram quantities per liter of E. coli culture, sufficient for thousands of tagmentation reactions. Benchmarking shows effectiveness for CUT&Tag, as well as bulk and single-cell ATAC-seq, with the enzyme retaining activity in storage for more than a year.
#### Link to preprint: https://www.biorxiv.org/content/10.1101/2024.07.11.602973v1
## Description of Code:
- **OpenTn5_ATACseq_MCF7_analysis.R**
   
  ⤷ Comparison of ATAC-seq quality control metrics between in-house pG-Tn5 & Illumina's TDE1. Downstream analysis after preprocessing was done with the ATACseq pipeline: https://github.com/riscalab/pipeSeq/blob/master/ATACseq.sh 
  1. `count fragments under peaks`
  2. `Spearman rank correlation on reads in peaks & different sized tiles`
  3. `Amount of E.Coli reads`
  4. `Spearman rank correlation of K562 cells using different volumes of pG-Tn5 transposome, as compared to standard OmniATAC-seq protocol 2.50 µL volume of Illumina TDE1, calculated over 1 kb bins genome-wide`
- **OpenTn5_CUTnTag_Analysis.Rmd**
  
  ⤷ Comparison of CUT&Tag quality control metrics between in-house pG-Tn5 & Illumina's TDE1.
  1. `peakcalling with SEACR with stringent threshold of 1%`
  2. `calculate peak width`
  3. `calculate percent of reproducible peaks between replicates`
  4. `calculate FRIP`
  5. `count fragments under peaks`
  6.  `Spearman rank correlation on reads in peaks & different sized tiles`
  7.  `plot fragment length distribution from sequenced samples & tapestation traces` 
- **OpenTn5_CUTnTag_MCF7_AdditionalQC.Rmd**
  
    ⤷ Additional QC metrics from SUMMARY output file generated from CUTnTag pipeline: https://github.com/riscalab/pipeSeq/blob/master/CUTnTag.sh
