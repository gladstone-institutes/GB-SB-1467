# GB-SB-1467
Mouse single-nucleus RNA-seq data analysis (GEM-X Flex Reagent Kit)  
  
## Citation
If you find this code or our paper useful in your research, please cite:

**Vitamin B2 metabolism promotes FSP1 stability to prevent ferroptosis**  
_Deol, K.K., Harris, C.A., Tomlinson, S.J. et al. Nat Struct Mol Biol 33, 525–536 (2026)._  
https://doi.org/10.1038/s41594-026-01759-x

## Researchers
- Jain lab: Skyler Blume, Ankur Garg, Isha Jain
- Corces lab: Adam Turner, Ryan Corces
- Bioinformatics core: Ayushi Agrawal, Michela Traglia, Reuben Thomas, Alex Pico

## Experimental details
- WT and KO mice treated with either PBS or Vitamin B3
- Collection of brains was done in batches
- Some mice are from the same litter
- Whole brain of mice was extracted and one half was used for transciptomics
- The new GEM-X Flex Reagent Kit was used, which in theory will obtain more nuclei per sample. Corces lab processed all 16 samples together.
- Used the Mouse GEM-X Flex Mouse Transcriptome Probe Kit (16 samples)
- 16 mouse samples, aimed to capture 20,000+ nuclei per mouse sample
- 4-plex workflow for 16 samples was used. 4 pools (4 GEM wells) were made, with 4 samples per pool. The samples going into each GEM well were randomized.
- After library preparation, the final Flex snRNA-seq libraries for the 4 GEM wells were pooled together. This pool was then sequenced across 3 lanes of a 25B flow cell on the NovaSeq X.
 
## Relevant folders

1. Hive: /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
 
## References
1. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
2. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html


</br>


