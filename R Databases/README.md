Each file on this folder contains a database ready to use with the R package TFEA.ChIP.

To use these databases with TFEA.ChIP, load the file in your R session and use the function *set_user_data()*:
```
load("path/to/file/database.RData")
set_user_data( MetaData, ChIPDB )
```
# TFEA.ChIP 1.15.2 Update

As time went by, the ChIP-seq collections used to build our database have grown considerably. In order to minimize the resources needed to use TFEA.ChIP, we are changing the format of our ChIP-Gene databases, from the original matrices:

```
> Mat01
       ChIPseq 1 ChIPseq 2 ChIPseq 3 ChIPseq 4
Gene 1         1         1         0         0
Gene 2         1         0         1         0
Gene 3         0         0         0         0
Gene 4         0         1         0         1
```
To a list-based format containing two elements:
* Gene Keys: vector of gene IDs
* ChIP Targets: list of vectors, one per ChIP-seq experiment in the, database, containing the putative targets assigned. Each target is coded as its position in the vector 'Gene Keys':
```
> ChIPDB
$`Gene Keys`
[1] "Gene 1" "Gene 2" "Gene 3" "Gene 4"

$`ChIP Targets`
$`ChIP Targets`$`ChIPseq 1`
[1] 1 2

$`ChIP Targets`$`ChIPseq 2`
[1] 1 4

$`ChIP Targets`$`ChIPseq 3`
[1] 2

$`ChIP Targets`$`ChIPseq 4`
[1] 4
```

Old databases are still compatible with the current version of the package.

Currently available databases are:

### Human ChIP-Seq databases:

**ReMap2022+EnsTSS+CellTypeEnh.Rdata**: links genes to ChIP-Seq experiments with peaks overlapping Transcription Starting Sites annotated in Ensemble Release 105 and cell type specific regulatory regions generated using ABC-Enhancer-Gene-Prediction [5]. ChIP-seq source: ReMap 2022 collection [4]

**ReMap2022+EnsTSS+GH.Rdata**: links genes to ChIP-Seq experiments with peaks overlapping regulatory regions stored in GeneHancer Double Elite [1] and Transcription Starting Sites annotated in Ensemble Release 105. ChIP-seq source: ReMap 2022 collection [4]

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

[4] ReMap 2022: a database of Human, Mouse, Drosophila and Arabidopsis regulatory regions from an integrative analysis of DNA-binding sequencing experiments
Fayrouz Hammal, Pierre De Langen, Aurélie Bergon, Fabrice Lopez, Benoit Ballester. Nucleic Acids Research, 2021 Nov 9. https://doi.org/10.1093/nar/gkab996

[5] Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat Genet 51, 1664–1669 (2019). Fulco, C.P., Nasser, J., Jones, T.R. et al. https://doi.org/10.1038/s41588-019-0538-0
