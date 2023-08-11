# Cell-specific_mQTLs

## Useful links

- Project document: https://docs.google.com/document/d/1nfo-oH2fhhFgFX3g32FkXf0iXA1ICXOms1ylqHheYb8/edit?pli=1
- RDSF location: /projects/MRC-IEU/research/projects/ieu3/p4/003
- BC4 location: /user/uy23281/scratch
- epifranklin location: /home/uy23281


## Background

Research question: 

Aims: 

- Compare mQTL associations across contexts
- Can we use predicted cell counts in blood to faithfully reproduce mQTLs in cell-sorted data



## Analysis

### Normalisation methylation data

#### Sorted cell types

To normalise use meffil

```
Rscript scripts/meffil.R
```

#### Bulk tissue

Already normalised, details - functional normalisation using meffil R package 

### Generate matrixeqtl inputs

Use steps 1-3 in GoDMC pipeline. Separate clone for each cell type.

### Perform mQTL analysis 

Use specific scripts in /resources/ folder for each cloned pipeline repository

### Comparison
- Are the GoDMC SNP-CpG pairs replicating in bulk+int?
- Are the GoDMC SNP-CpG pairs replicating in cell-specific? 
- Discovery in bulk+int replicating in cell-specific?
- Discovery in cell-specific replicating in bulk?
  
```
Rscript Compare_mQTLmodels.Rmd
```


