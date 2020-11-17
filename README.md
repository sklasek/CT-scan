## CT-scan

How X-ray CT scanning affects microbial communities in sediments

Authors: Erica Ewton, Scott Klasek, Erin Peck, Jason Wiest, and Frederick Colwell   
Publication: Not published yet.   

**Reproducible workflow:**
This repo shows codes for analyzing sediment microbial community responses to X-ray CT scanning.   
Sequence data are 16S rRNA V4 hypervariable regions are available [on NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA533633)   
For code showing the processing of ASVs from raw fastq sequences and the generation of phyloseq objects, see [sequence processing](https://github.com/sklasek/CT-scan/blob/master/markdowns/01_sequence_processing.md)      
For code showing analysis of microbial community data starting from the phyloseq object ps4, see [sequence analysis](https://github.com/sklasek/CT-scan/blob/master/markdowns/02_sequence_analysis.md)    
If you're not interested in processing these sequences and want to skip to the analysis, you can begin with the phyloseq object ps4, which has been removed of unwanted taxa, decontaminated, and pruned of samples with very poor read depth.   
