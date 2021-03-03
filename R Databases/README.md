Each file on this folder contains a database ready to use with the R package TFEA.ChIP.

To use these databases with TFEA.ChIP, load the file in your R session and use the function *set_user_data()*:
```
load("path/to/file/database.RData")
set_user_data( MetaData, Mat01 )
```


Currently available databases are:

### Human ChIP-Seq databases:

ReMap+GH_doubleElite.Rdata: links genes to ChIP-Seq experiments with peaks overlapping regulatory regions stored in GeneHancer Double Elite [1]. ChIP-seq source: ReMap 2020 collection [3]

ReMap+GH_doubleElite.Rdata: links genes to ChIP-Seq experiments with peaks overlapping regulatory regions stored in GeneHancer Double Elite [1]. ChIP-seq source: ReMap 2018 collection [2]

TC_TSS1kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 1Kb from a gene.

TC_TSS5kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 5Kb from a gene.

TC_TSS5kb+distantSites.Rdata: includes binding to enhancers linked to known genes[1] 

TC_TSS10kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 10Kb from a gene.

TC_TSS10kb+distantSites.Rdata: includes binding to enhancers linked to known genes[1]


RM_TSS1kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 1Kb from a gene. ChIP-seq source: ReMap 2018 collection [2]

RM_TSS1kb+distantSites.Rdata: includes binding to enhancers linked to known genes[1]. ChIP-seq source: ReMap 2018 collection [2]

RM_TSS5kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 5Kb from a gene. ChIP-seq source: ReMap 2018 collection [2]

RM_TSS5kb+distantSites.Rdata: includes binding to enhancers linked to known genes[1]. ChIP-seq source: ReMap 2018 collection [2]

RM_TSS10kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 10Kb from a gene. ChIP-seq source: ReMap 2018 collection [2]

RM_TSS10kb+distantSites.Rdata: includes binding to enhancers linked to known genes[1]. ChIP-seq source: ReMap 2018 collection [2]

### Mouse ChIP-seq databases:

mm_TSS1kb:links genes to ChIP-Seq experiments with peaks up to 1Kb from a gene.

mm_TSS5kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 5Kb from a gene.

mm_TSS10kb.Rdata: links genes to ChIP-Seq experiments with peaks up to 10Kb from a gene.

### References

[1] Distant binding sites are taken from GeneHancer database:
GeneHancer: genome-wide integration of enhancers and target genes in GeneCards
Simon Fishilevich Ron Nudel Noa Rappaport Rotem Hadar Inbar Plaschkes Tsippi Iny Stein Naomi Rosen Asher Kohn Michal Twik Marilyn Safran
Database, Volume 2017, 1 January 2017, bax028, https://doi.org/10.1093/database/bax028

[2] ReMap 2018: an updated atlas of regulatory regions from an integrative analysis of DNA-binding ChIP-seq experiments
Jeanne Chèneby Marius Gheorghe Marie Artufel Anthony Mathelier Benoit Ballester
Nucleic Acids Research, Volume 46, Issue D1, 4 January 2018, Pages D267–D275
http://tagc.univ-mrs.fr/remap/

[3] ReMap 2020: a database of regulatory regions from an integrative analysis of Human and Arabidopsis DNA-binding sequencing experiments
Jeanne Chèneby, Zacharie Ménétrier, Martin Mestdagh, Thomas Rosnet, Allyssa Douida, Wassim Rhalloussi, Aurélie Bergon, Fabrice Lopez, Benoit Ballester
Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D180–D188, https://doi.org/10.1093/nar/gkz945
