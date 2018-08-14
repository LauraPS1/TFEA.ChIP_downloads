# TFEA.ChIP
TFEA.ChIP is an R package in developement. Its purpose is to analyze transcription factor enrichment in a set of differentially expressed genes.

A web implementation of TFEA.ChIP is available at https://www.iib.uam.es/TFEA.ChIP/

Specifically TFEA.ChIP, uses information derived from the hundreds of ChIP-Seq experiments from the 
ENCODE Consortium<sup>[1]</sup>  expanded to include additional
datasets contributed to GEO database<sup>[2][3]</sup> by individual laboratories
representing the binding sites of factors not assayed by ENCODE. The package includes a set of tools
to customize the ChIP data, perform enrichment analysis and visualize the results. The package implements
two enrichment analysis methods:

* Analysis of the association of TFBS and differential expression from 2x2 tables recording the presence
of binding sites for a given TF in DE and control genes. The statistical significance of the association
for each factor determined by a Fisher’s exact test.

* GSEA analysis, based on the core function of the GSEA algorithm for R<sup>[4][5]</sup>, *GSEA.EnrichmentScore*.
<br>
1: ENCODE Project Consortium (2012) Nature 489, 57-74<br>
2: Edgar, R et al. (2002) Nucleic Acids Res. 30:207-10<br>
3: Barrett, T et al. (2013) Nucleic Acids Res. 41(Database issue):D991-5<br>
4: Subramanian, Tamayo, et al. (2005) PNAS 102, 15545-15550<br>
5: Mootha, Lindgren, et al. (2003) Nat Genet 34, 267-273<br>
