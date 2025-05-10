# OpenTn5: Open-Source Resource for Robust and Scalable Tn5 Transposase Purification and Characterization
### Abstract: 
Tagmentation streamlines sequencing library preparation by leveraging the bacterial cut-and-paste Tn5 transposase for DNA fragmentation and sequencing adapter addition. Here we present an open-source protocol for the generation of multi-purpose hyperactive Tn5 transposase, pG-Tn5E54K, L372P, yielding multi-milligram quantities per liter of E. coli culture, sufficient for thousands of tagmentation reactions. Benchmarking shows effectiveness for CUT&Tag, as well as bulk and single-cell ATAC-seq, with the enzyme retaining activity in storage for more than a year.
### Description of Code:
- `OpenTn5_ATACseq_MCF7_analysis.R`  
  ⤷ *Comparison of ATAC-seq quality control metrics between in-house pG-Tn5 & Illumina's TDE1
  1. `count fragments under peaks`
  2. `Spearman rank correlation on reads in peaks & different sized tiles`
  3. `Amount of E.Coli reads`
  4. `Spearman rank correlation of  K562 cells using different volumes of pG-Tn5 transposome, as compared to standard OmniATAC-seq protocol 2.50 µL volume of Illumina TDE1, calculated over 1 kb bins genome-wide`
- `OpenTn5_CUTnTag_Analysis.Rmd`
  ⤷ *Comparison of CUT&Tag quality control metrics between in-house pG-Tn5 & Illumina's TDE1
  1. `peakcalling with SEACR with stringent threshold of 1%`
  2. `calculate peak width`
  3. `calculate percent of reproducible peaks between replicates`
  4. `calculate FRIP`
  5. `count fragments under peaks`
  6.  `Spearman rank correlation on reads in peaks & different sized tiles`
  7.  `plot fragment length distribution from sequenced samples & tapestation traces`
  8.  `plot coverage around TSS`
