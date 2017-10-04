################### IMPORTS ####
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges GRanges distanceToNearest
#' @importFrom IRanges IRanges
#' @importFrom biomaRt select
#' @importFrom dplyr "%>%" arrange
#' @importFrom grDevices colorRamp
#' @importFrom stats fisher.test p.adjust
#' @importFrom utils data
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import org.Hs.eg.db
################### FUNCTIONS ####

txt2GR<-function(fileTable,format,fileMetaData,alpha=NULL){

    #' @title Function to filter a ChIP-Seq input.
    #' @description Function to filter a ChIP-Seq output (in .narrowpeak or
    #' MACS's peaks.bed formats) and then store the peak coordinates in a
    #' GenomicRanges object, associated to its metadata.
    #' @param fileTable data frame from a txt/tsv/bed file
    #' @param format "narrowPeak" or "macs".
    #' narrowPeak fields:
    #' 'chrom','chromStart','chromEnd','name','score','strand','signalValue',
    #' 'pValue','qValue','peak'
    #' macs fields:
    #' 'chrom','chromStart','chromEnd','name','qValue'
    #' @param fileMetaData Data frame/matrix/array contaning the following
    #' fields: 'Name','Accession','Cell','Cell Type','Treatment','Antibody',
    #' 'TF'.
    #' @param alpha max p-value to consider ChIPseq peaks as significant and
    #' include them in the database. By default alpha is 0.05 for narrow peak
    #' files and 1e-05 for MACS files
    #' @return The function returns a GR object generated from the ChIP-Seq
    #' dataset input.
    #' @export txt2GR
    #' @examples
    #' data("ARNT.peaks.bed","ARNT.metadata",package = "TFEA.ChIP")
    #' ARNT.gr<-txt2GR(ARNT.peaks.bed,"macs",ARNT.metadata)

    requireNamespace("GenomicRanges")
    requireNamespace("IRanges")
    requireNamespace("dplyr")

    if (is.data.frame(fileMetaData)==FALSE){
        if(is.matrix(fileMetaData)){
            if(length(fileMetaData[1,])==7){
                fileMetaData<-as.data.frame(
                    fileMetaData, stringsAsFactors=FALSE)
            }else{
                warning("fileMetaData format error: 'fileMetaData' must be a",
                        " data frame/matrix/array with 7 atributes: 'Name',",
                        "'Accession', 'Cell', 'Cell Type','Treatment',",
                        "'Antibody','TF'")
                break
            }
        }else if(is.array(fileMetaData)){
            if(length(fileMetaData)==7){
                fileMetaData<-as.data.frame(
                    fileMetaData,stringsAsFactors=FALSE)
            }else{
                warning("fileMetaData format error: 'fileMetaData' must be a",
                        " data frame/matrix/array with 7 atributes: 'Name',",
                        "'Accession', 'Cell', 'Cell Type','Treatment',",
                        "'Antibody','TF'")
                break
            }
        }else{
            warning("fileMetaData format error: 'fileMetaData' must be a data",
                "frame/matrix/array with 7 atributes: 'Name', 'Accession',",
                " 'Cell', 'Cell Type','Treatment','Antibody','TF'")
            break
        }
        }
    format<-tolower(format)

    if(format=="narrowpeak"){

        if(fileTable[1,8]==-1 & fileTable[1,9]==-1){
            warning ("The ChIP-Seq input file ",fileMetaData$Name,
                " does not include p-value or Q-value for each peak.",
                " Please, make sure the peaks in the input file have been",
                " previously filtered according to their significance")
            fileTable<-fileTable[,1:3]
            Stat<-"no score"
            fileTable$score=rep(NA,length(fileTable[,1]))
            colnames(fileTable)[1:3]<-c("chr","start","end")
        }else if(fileTable[1,8]==-1){
            fileTable<-fileTable[,c(1,2,3,9)]
            colnames(fileTable)<-c("chr","start","end","score")
            Stat<-"log10(p-Value)"
            if(is.null(alpha)){valLimit<-1.3}else{valLimit<-(-log10(alpha))}
            fileTable<-fileTable[fileTable$score>valLimit,]
        }else if(fileTable[1,9]==-1){
            fileTable<-fileTable[,c(1,2,3,8)]
            colnames(fileTable)<-c("chr","start","end","score")
            Stat<-"p-Value"
            if(is.null(alpha)){valLimit<-0.05}else{valLimit<-alpha}
            fileTable<-fileTable[fileTable$score<valLimit,]
        }

        fileMetaData<-c(fileMetaData,Stat)
        MDframe<-as.data.frame(lapply(
            fileMetaData, rep,length(fileTable[,1])))
        colnames(MDframe)<-c("Name","Accession","Cell","Cell Type",
            "Treatment", "Antibody","TF","Score Type")

        gr<-GenomicRanges::GRanges(
            seqnames=fileTable$chr,
            ranges=IRanges::IRanges(fileTable$start,end=fileTable$end),
            score=fileTable$score,
            mcols=MDframe
        )
        return(gr)

    }else if(format=="macs"){
        if(is.null(alpha)){valLimit<-50}else{valLimit<-(-10*log10(alpha))}
        if(length(fileTable[1,])==5){
            fileTable<-fileTable[,c(1,2,3,5)]
            colnames(fileTable)<-c("chr","start","end","score")
            fileTable<-fileTable[fileTable$score>valLimit,]
            Stat<-"-10*log.Pvalue"
        }else if (length(fileTable[1,])==4 & is.character(fileTable[1,4])){
            # if the 4th column consists of peak names
            warning ("The ChIP-Seq input file does not include p-value or",
                " Q-value for each peak. Please, make sure the peaks",
                " in the input file have been previously filtered",
                " according to their significance")
            fileTable<-fileTable[,1:3]
            Stat<-"no score"
            fileTable$score=rep(NA,length(fileTable[,1]))
            colnames(fileTable)[1:3]<-c("chr","start","end")
        }else if (length(fileTable[1,])==4 & !is.character(fileTable[1,4])){
            # if the 4th column consists of adjusted p-values
            colnames(fileTable)<-c("chr","start","end","score")
            fileTable<-fileTable[fileTable$score>valLimit,]
            Stat<-"-10*log.Pvalue"
        }

        fileMetaData<-c(fileMetaData,Stat)
        MDframe<-as.data.frame(lapply(fileMetaData, rep,length(fileTable[,1])))
        colnames(MDframe)<-c("Name","Accession","Cell","Cell Type","Treatment",
            "Antibody","TF","Score Type")

        gr<-GenomicRanges::GRanges(
            seqnames=fileTable$chr,
            ranges=IRanges::IRanges(fileTable$start,end=fileTable$end),
            score=fileTable$score,
            mcols=MDframe
        )
        return(gr)

    }else{
        warning("Wrong file format. Only narrowPeak or MACS output",
            " ('_peaks.bed') are supported.\nPlease choose either",
            " 'narrowPeak' or 'Macs'.")
        break
    }
        }

GR2tfbs_db<-function(Ref.db,gr.list,distanceMargin=10){

    #' @title Makes a TFBS-gene binding database
    #' @description GR2tfbs_db generates a TFBS-gene binding database
    #' through the association of ChIP-Seq peak coordinates (provided
    #' as a GenomicRange object) to overlapping genes or gene-associated
    #' Dnase regions (Ref.db).
    #' @param Ref.db GenomicRanges object containing a database of reference
    #' elements (either Genes or gene-associate Dnase regions) including
    #' a gene_id metacolumn
    #' @param gr.list List of GR objects containing ChIP-seq peak coordinates
    #' (output of txt2GR).
    #' @param distanceMargin Maximum distance allowed between a gene or DHS to
    #' assign a gene to a ChIP-seq peak. Set to 10 bases by default.
    #' @return List of vectors, one for every ChIP-Seq, storing the IDs of the
    #' genes to which the TF bound in the ChIP-Seq.
    #' @export GR2tfbs_db
    #' @examples
    #' data("DnaseHS_db","gr.list", package="TFEA.ChIP")
    #' GR2tfbs_db(DnaseHS_db, gr.list)

    if(!requireNamespace("S4Vectors", quietly = TRUE)){
        stop("S4Vectors package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }
    requireNamespace("S4Vectors")

    m<-0

    for(i in 1:length(gr.list)){

        gr<-gr.list[[i]]

        nearest_index<-suppressWarnings(
            GenomicRanges::distanceToNearest(gr,Ref.db,select="all"))
        nearest_index<-nearest_index[
            !is.na(nearest_index@elementMetadata@listData$distance)]
        nearest_index<-nearest_index[
            nearest_index@elementMetadata@listData$distance<=distanceMargin]
        inSubject<-S4Vectors::subjectHits(nearest_index)

        if(length(inSubject)==0){
            # in case any ChIP-Seq dataset does not have any genes to be
            # assigned.
            m<-m+1
            next
        }
        assigned_genes<-Ref.db[inSubject]$gene_id

        if ((i-m)==1){
            TFgenes_list<-list(assigned_genes)
            names(TFgenes_list)[i-m]<-as.character(
                gr.list[[i]]$mcols.Accession[1])
        }else{
            TFgenes_list<-c(TFgenes_list,list(assigned_genes))
            names(TFgenes_list)[i-m]<-as.character(
                gr.list[[i]]$mcols.Accession[1])
        }
    }

    return(TFgenes_list)
}

SearchID<-function(GeneID,id_db){

    #' @title Searchs for an Entrez gene ID in a TF-gene database
    #' @description Searchs for an Entrez gene ID in a db that contains
    #' the IDs of the genes a TF binds to.
    #' @param GeneID gene Entrez ID.
    #' @param id_db TF - gene binding database.
    #' @return 1/0 row. Each element represents a TF ChIP-Seq experiment.

    # examples
    # SearchID(GeneID,id_db)

    TF_row<-rep(0,length=length(id_db))

    for (i in 1:length(id_db)){
        if (GeneID %in% id_db[[i]]==TRUE){TF_row[i]<-1}
    }
    return(TF_row)
}

makeTFBSmatrix<-function(GeneList,id_db){

    #' @title Function to search for a list of entrez gene IDs.
    #' @description Function to search for a list of entrez gene IDs
    #' in a  TF-gene binding data base.
    #' @param GeneList Array of gene Entrez IDs
    #' @param id_db TF - gene binding database.
    #' @return 1/0 matrix. Each row represents a gene, each column,
    #' a ChIP-Seq file.
    #' @export makeTFBSmatrix
    #' @examples
    #' data("tfbs.database","Entrez.gene.IDs",package = "TFEA.ChIP")
    #' makeTFBSmatrix(Entrez.gene.IDs,tfbs.database)

    for (i in 1:length(GeneList)){

        TFrow<-SearchID(GeneList[i],id_db)

        if (i==1){
            TF_matrix<-matrix(ncol=length(id_db))
            TF_matrix<-TFrow
            rm(TFrow)
        }else{
            TF_matrix<-rbind(TF_matrix,TFrow)
            rm(TFrow)
        }
    }

    colnames(TF_matrix)<-names(id_db)
    rownames(TF_matrix)<-GeneList
    return(TF_matrix)
}

set_user_data<-function(metadata,binary_matrix){
    #' @title Sets the data objects as default.
    #' @description Function to set the data objects provided by the user
    #' as default to the rest of the functions.
    #' @param metadata Data frame/matrix/array contaning the following fields:
    #' 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'.
    #' @param binary_matrix Matrix[n,m] which rows correspond to all the human
    #' genes that have been assigned an Entrez ID, and its columns, to every
    #' ChIP-Seq experiment in the database. The values are 1 – if the ChIP-Seq
    #' has a peak assigned to that gene – or 0 – if it hasn’t –.
    #' @return sets the user's metadata table and TFBS matrix as the variables
    #' "MetaData" and "Mat01", used by the rest of the package.
    #' @export set_user_data
    #' @examples
    #' data("MetaData","Mat01",package="TFEA.ChIP")
    #' # For this example, we will usethe variables already included in the
    #' # package.
    #' set_user_data(MetaData,Mat01)

    assign("MetaData", metadata, envir = .GlobalEnv)
    assign("Mat01",binary_matrix,envir = .GlobalEnv)
}

GeneID2entrez<-function(gene.IDs,return.Matrix = FALSE){

    #' @title Translates gene IDs from Gene Symbol or Ensemble ID to Entrez ID.
    #' @description Translates gene IDs from Gene Symbol or Ensemble Gene ID
    #' to Entrez Gene ID using the IDs approved by HGNC. When translating from
    #' Gene Symbol, keep in mind that many genes have been given more than one
    #' symbol through the years. This function will return the Entrez ID
    #' corresponding to the currently approved symbols if they exist, otherwise
    #' NA is returned. In addition some genes might map to more than one Entrez
    #' ID, in this case gene is assigned to the first match and a warning
    #' is displayed.
    #' @param gene.IDs Array of Gene Symbols or Ensemble Gene IDs.
    #' @param return.Matrix T/F. When TRUE, the function returns a matrix[n,2],
    #' one column with the gene symbols or Ensemble IDs, another with their
    #' respective Entrez IDs.
    #' @return Vector or matrix containing the Entrez IDs(or NA) corresponding
    #' to every element of the input.
    #' @export GeneID2entrez
    #' @examples
    #' GeneID2entrez(c("TNMD","DPM1","SCYL3","FGR","CFH","FUCA2","GCLC"))

    requireNamespace("biomaRt")
    requireNamespace("GenomicFeatures")

    Genes<-GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    suppressMessages(GeneNames<-biomaRt::select(
        org.Hs.eg.db, Genes$gene_id, c("SYMBOL", "ENSEMBL")))
    # suppressWarnings added to avoid 'select()' returned 1:many mapping
    # between keys and columns

    if(grepl("ENSG0",gene.IDs[1])==TRUE){ID.type<-"ENSEMBL"
    }else{ID.type<-"SYMBOL"}

    tmp<-match(as.character(gene.IDs),GeneNames[,ID.type])
    tmp2<-match(GeneNames[,ID.type],as.character(gene.IDs))

    if(sum(duplicated(tmp2[!is.na(tmp2)]))>0){
        warning("Some genes returned 1:many mapping to ENTREZ ID. ",
                "Genes were assigned the first ENTREZ ID match found.\n",
                call. =FALSE)
        }
    cat("Done! ",length(tmp[!is.na(tmp)])," genes of ",length(tmp),
        " successfully translated.\n")

    if(return.Matrix==TRUE){
        if (length(tmp[is.na(tmp)])>0){cat("Couldn't find Entrez IDs for ",
            length(tmp[is.na(tmp)])," genes (NAs returned instead).\n")}
        return(data.frame(
            GENE.ID=gene.IDs,ENTREZ.ID=GeneNames[tmp,"ENTREZID"]))
    }else{
        if (length(tmp[is.na(tmp)])>0){cat("Couldn't find Entrez IDs for ",
            length(tmp[is.na(tmp)]),"genes.\n")}
        return(GeneNames[tmp[!is.na(tmp)],"ENTREZID"])
    }
}

get_chip_index<-function(database = "g",TFfilter = NULL){

    #' @title Creates df containing accessions of ChIP-Seq datasets and TF.
    #' @description Function to create a data frame containing the ChIP-Seq
    #' dataset accession IDs and the transcription factor tested in each ChIP.
    #' This index is used in functions like “contingency_matrix” and “GSEA_run”
    #' as a filter to select specific ChIPs or transcription factors to run an
    #' analysis.
    #' @param database Name of the database used: "encode"/"e" or "general"/"g".
    #' @param TFfilter (Optional) Transcription factors of interest.
    #' @return Data frame containig the accession ID and TF for every ChIP-Seq
    #' experiment included in the metadata files.
    #' @export get_chip_index
    #' @examples
    #' get_chip_index(database="e")
    #' get_chip_index(database="general",TFfilter=c("SMAD2","SMAD4"))

    requireNamespace("dplyr")
    requireNamespace("utils")
    if(!exists("MetaData")){
        data("MetaData",package = "TFEA.ChIP",envir = environment())
        }
    if(is.null(TFfilter)){
        Index<-dplyr::select(MetaData,Accession,TF)
        database<-tolower(database)
        if(database=="encode"|database=="e"){
            Index<-Index[grepl("^wg.*",Index$Accession),]
        }
        return(Index)

    }else{
        Index<-dplyr::select(MetaData,Accession,TF)
        Index<-Index[Index$TF %in% TFfilter,]
        database<-tolower(database)
        if(database=="encode"|database=="e"){
            Index<-Index[grepl("^wg.*",Index$Accession),]
        }
        return(Index)
    }
}

contingency_matrix<-function(test_list,control_list,chip_index=get_chip_index()){

    #' @title Computes 2x2 contingency matrices
    #' @description Function to compute contingency 2x2 matrix by the partition
    #' of the two gene ID lists according to the presence or absence of the
    #' terms in these list in a ChIP-Seq binding database.
    #' @param test_list List of gene Entrez IDs
    #' @param control_list If not provided, all human genes not present in
    #' test_list will be used as control.
    #' @param chip_index Output of the function “get_chip_index”, a data frame
    #' containing accession IDs of ChIPs on the database and the TF each one
    #' tests. If not provided, the whole internal database will be used
    #' @return List of contingency matrices, one CM per element in chip_index
    #' (i.e. per ChIP-seq dataset).
    #' @export contingency_matrix
    #' @examples
    #' data("Genes.Upreg",package = "TFEA.ChIP")
    #' CM_list_UP <- contingency_matrix(Genes.Upreg)

    requireNamespace("utils")

    if (missing(control_list)){
        # Generating control gene list in case is not provided.
        Genes<-GenomicFeatures::genes(
            TxDb.Hsapiens.UCSC.hg19.knownGene)$gene_id
        control_list<-Genes[!(Genes %in% test_list)]
        rm(Genes)
    }else{
        control_list<-control_list[!(control_list %in% test_list)]
    }
    if (exists("Mat01")==FALSE){
        data("Mat01",package = "TFEA.ChIP",envir = environment())
        }
    Matrix1<-Mat01[rownames(Mat01)%in%test_list,
        colnames(Mat01)%in%chip_index$Accession]
    Matrix2<-Mat01[rownames(Mat01)%in%control_list,
        colnames(Mat01)%in%chip_index$Accession]

    for (i in 1:length(chip_index[,1])){
        chip.vector1<-Matrix1[,chip_index$Accession[i]]
        chip.vector2<-Matrix2[,chip_index$Accession[i]]

        pos1<-length(chip.vector1[chip.vector1==1])
        pos2<-length(chip.vector2[chip.vector2==1])
        neg1<-length(chip.vector1[chip.vector1==0])
        neg2<-length(chip.vector2[chip.vector2==0])

        contMatrix<-cbind(c(pos1,pos2),c(neg1,neg2))
        rownames(contMatrix)<-c("Test","Control")
        colnames(contMatrix)<-c("Positive","Negative")

        if (i==1){
            contMatrix_list<-list(contMatrix)
            names(contMatrix_list)[i]<-as.character(chip_index$Accession[i])
            rm(contMatrix)
        }else{
            contMatrix_list<-c(contMatrix_list,list(contMatrix))
            names(contMatrix_list)[i]<-as.character(chip_index$Accession[i])
            rm(contMatrix)
        }
    }
    return(contMatrix_list)
}

getCMstats<-function(contMatrix_list,chip_index=get_chip_index()){

    #' @title Generate statistical parameters from a contingency_matrix output
    #' @description From a list of contingency matrices, such as the output
    #' from “contingency_matrix”, this function computes a fisher's exact test
    #' for each matrix and generates a data frame that stores accession ID of a
    #' ChIP-Seq experiment, the TF tested in that experiment, the p-value and
    #' the odds ratio resulting from the test.
    #' @param contMatrix_list Output of “contingency_matrix”, a list of
    #' contingency matrix.
    #' @param chip_index Output of the function “get_chip_index”, a data frame
    #' containing accession IDs of ChIPs on the database and the TF each one
    #' tests. If not provided, the whole internal database will be used
    #' @return Data frame containing accession ID of a ChIP-Seq experiment, the
    #' TF tested in that experiment, raw p-value (-10*log10 pvalue), odds-ratio
    #' and FDR-adjusted p-values (-10*log10 adj.pvalue).
    #' @export getCMstats
    #' @examples
    #' data("CM_list",package = "TFEA.ChIP")
    #' stats_mat_UP <- getCMstats(CM_list)

    requireNamespace("stats")

    for (i in 1:length(contMatrix_list)){
        pval<-stats::fisher.test(x=contMatrix_list[[i]])
        if (i==1){
            pval_list<-list(pval)
            names(pval_list)[i]<-names(contMatrix_list)[i]
        }else{
            pval_list<-c(pval_list,list(pval))
            names(pval_list)[i]<-names(contMatrix_list)[i]
        }
    }

    statMat<-data.frame(
        Accession=chip_index$Accession,
        TF=chip_index$TF,
        p.value=NA,
        OR=NA)
    for (idx in names(contMatrix_list)){
        FTres<-try({stats::fisher.test(x=contMatrix_list[[idx]])},
            silent = TRUE)
        if(class(FTres)=="htest"){
            statMat$p.value[which(statMat$Accession==idx)]<-FTres$p.value
            statMat$OR[which(statMat$Accession==idx)]<-FTres$estimate
        }
    }
    statMat$adj.p.value<-stats::p.adjust(statMat$p.value,"fdr")
    statMat$log.adj.pVal<-(-1*(log10(statMat$adj.p.value)))

    return(statMat)
}

GSEA_EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 0, correl.vector = NULL) {

    # Computes the weighted GSEA score of gene.set in gene.list.
    # Developed by The Broad Institute

    #' @title Computes the weighted GSEA score of gene.set in gene.list.
    #' @description Computes the weighted GSEA score of gene.set in gene.list.
    #' @param gene.list The ordered gene list
    #' @param gene.set A gene set, e.g. gene IDs corresponding to a ChIP-Seq
    #' experiment's peaks.
    #' @param weighted.score.type Type of score: weight:
    #' 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
    #' @param correl.vector A vector with the coorelations (such as signal to
    #' noise scores) corresponding to the genes in the gene list
    #' @return list of:
    #' ES: Enrichment score (real number between -1 and +1)
    #' arg.ES: Location in gene.list where the peak running enrichment occurs
    #' (peak of the "mountain")
    #' RES: Numerical vector containing the running enrichment score for all
    #' locations in the gene list
    #' tag.indicator: Binary vector indicating the location of the gene sets
    #' (1's) in the gene list
    #' @export GSEA_EnrichmentScore
    #' @examples
    #' GSEA_EnrichmentScore(gene.list=c("3091","2034","405","55818"),
    #' gene.set=c("2034","112399","405"))

    tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <-  N - Nh
    if (weighted.score.type == 0) {
        correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector**alpha)
    sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
    if(sum.correl.tag>0){
        norm.tag    <- 1.0/sum.correl.tag
        norm.no.tag <- 1.0/(N-sum.correl.tag)
        RES <- cumsum(tag.indicator * correl.vector * norm.tag -
            no.tag.indicator * norm.no.tag)
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (abs(max.ES) > abs(min.ES)) {
            #      ES <- max.ES
            ES <- signif(max.ES, digits = 5)
            arg.ES <- which.max(RES)
        } else {
            #      ES <- min.ES
            ES <- signif(min.ES, digits=5)
            arg.ES <- which.min(RES)
        }
        return(list(ES = ES, arg.ES = arg.ES, RES = RES,
            indicator = tag.indicator))

    }else{
        RES<-rep(NaN,length(gene.list))
        indicator<-rep(0,length(gene.list))
        ES<-NA
        arg.ES<-NA
        return(list(ES = ES, arg.ES = arg.ES, RES = RES,
            indicator = tag.indicator))
    }
}

GSEA_Shuffling<-function(gene.list,permutations){

    #' @title Function to create shuffled gene lists to run GSEA.
    #' @description Function to create shuffled gene lists to run GSEA.
    #' @param gene.list Vector of gene Entrez IDs.
    #' @param permutations Number of shuffled gene vector IDs required.
    #' @return Vector of randomly arranged Entrez IDs.
    # examples
    # GSEA_Shuffling(gene.list,1000)

    shuffledGL<-list()
    for(i in 1:permutations){
        tmp<-sample(gene.list)
        shuffledGL<-c(shuffledGL,list(tmp))
    }
    return(shuffledGL)
}

GSEA_run<-function(gene.list,chip_index=get_chip_index(),get.RES = FALSE,RES.filter = NULL){

    #' @title Function to run a GSEA analysis
    #' @description Function to run a GSEA to analyze the distribution of TFBS
    #' across a sorted list of genes.
    #' @param gene.list List of Entrez IDs ordered by their fold change.
    #' @param chip_index Output of the function “get_chip_index”, a data frame
    #' containing accession IDs of ChIPs on the database and the TF each one
    #' tests. If not provided, the whole internal database will be used
    #' @param get.RES (Optional) boolean. If TRUE, the function stores RES of
    #' all/some TF.
    #' @param RES.filter (Optional) chr vector. When get.RES==TRUE, allows to
    #' choose which TF's RES to store.
    #' @return a list of:
    #' Enrichment.table: data frame containing accession ID, TF name,
    #' enrichment score, p-value, and argument of every ChIP-Seq experiment.
    #' RES (optional): list of running sums of every ChIP-Seq
    #' indicators (optional): list of 0/1 vectors that stores the matches (1)
    #' and mismatches (0) between the gene list and the gene set.
    #' @export GSEA_run
    #' @examples
    #' data("Entrez.gene.IDs",package = "TFEA.ChIP")
    #' chip_index<-get_chip_index(TFfilter = c("HIF1A","EPAS1","ARNT"))
    #' GSEA.result <- GSEA_run(Entrez.gene.IDs,chip_index,get.RES = TRUE)

    requireNamespace("stats")
    requireNamespace("utils")

    if (!exists("Mat01")){
        data("Mat01",package = "TFEA.ChIP",envir = environment())
    }
    Mat01<-Mat01[,colnames(Mat01)%in%chip_index$Accession]

    shuffledGL<-GSEA_Shuffling(gene.list,1000)
    # Generate random gene lists to get a p-value for ESs.
    enrichmentScore<-vector()
    pval<-vector()
    enrichmentArg<-vector()

    if(get.RES==TRUE){
        res<-list()
        ind<-list()
    }
    for (i in 1:length(chip_index$Accession)){

        chip.genes<-Mat01[,colnames(Mat01)==chip_index$Accession[i]]
        chip.genes<-names(chip.genes[chip.genes==1])

        if(length(chip.genes)>10){
            result<-GSEA_EnrichmentScore(gene.list, chip.genes)
            shuffled.ES<-vector()
            for (j in 1:length(shuffledGL)){
                # Get ES for the random gene list to compute pval.
                shuffled.ES[j]<-GSEA_EnrichmentScore(
                    shuffledGL[[j]], chip.genes)$ES
            }

            enrichmentScore[i]<-result$ES
            pval[i]<-sum(abs(shuffled.ES) >= abs(result$ES)) / 1000
            enrichmentArg[i]<-result$arg.ES

            if(get.RES==TRUE & missing(RES.filter)){
                # Store running sums of selected TFs.
                res<-c(res,list(result$RES))
                names(res)[length(res)]<-chip_index$Accession[i]
                ind<-c(ind,list(result$indicator))
                names(ind)[length(ind)]<-chip_index$Accession[i]
            }else if(get.RES==TRUE & missing(RES.filter)==FALSE){
                if(chip_index$TF[i]%in%RES.filter){
                    res<-c(res,list(result$RES))
                    names(res)[length(res)]<-chip_index$Accession[i]
                    ind<-c(ind,list(result$indicator))
                    names(ind)[length(ind)]<-chip_index$Accession[i]
                }
            }
        }else{chip_index<-chip_index[-i,]}

        if(i==abs(length(chip_index[,1])*0.25)){
            cat("|||||||||| 25%\n")
        }else if (i==abs(length(chip_index[,1])*0.5)){
            cat("|||||||||||||||||||| 50%\n")
        }else if (i==abs(length(chip_index[,1])*0.75)){
            cat("|||||||||||||||||||||||||||||| 75%\n")
        }else if (i==length(chip_index[,1])){
            cat("|||||||||||||||||||||||||||||||||||||||| Done!\n")}
    }
    pval.adj<-stats::p.adjust(pval,"fdr") # Adjust pvalues

    enrichmentTable<-cbind(chip_index$Accession,chip_index$TF,
        as.numeric(enrichmentScore),as.numeric(pval.adj),
        as.numeric(enrichmentArg))

    enrichmentTable<-as.data.frame(enrichmentTable,stringsAsFactors=FALSE)
    colnames(enrichmentTable)<-c("Accession","TF","ES","pval.ES","Arg.ES")
    enrichmentTable$ES<-as.numeric(enrichmentTable$ES)
    enrichmentTable$pval.ES<-as.numeric(enrichmentTable$pval.ES)
    enrichmentTable$Arg.ES<-as.numeric(enrichmentTable$Arg.ES)
    enrichmentTable<-enrichmentTable[!is.na(enrichmentTable$pval.ES),]

    if(get.RES==TRUE){
        GSEA_results<-list(enrichmentTable,res,ind)
        names(GSEA_results)<-c("Enrichment.table","RES","indicators")
        return(GSEA_results)
    }else{
        return(enrichmentTable)
    }
}

plot_CM<-function(CM.statMatrix,plot_title = NULL,specialTF = NULL,TF_colors = NULL){

    #' @title Makes an interactive html plot from an enrichment table.
    #' @description Function to generate an interactive html plot from a
    #' transcription factor enrichment table, output of the function
    #' "getCMstats".
    #' @param CM.statMatrix Output of the function "getCMstats".
    #' A data frame storing: Accession ID of every ChIP-Seq tested,
    #' Transcription Factor,Odds Ratio, p-value and adjusted p-value.
    #' @param plot_title The title for the plot.
    #' @param specialTF (Optional) Named vector of TF symbols -as written in
    #' the enrichment table- to be highlighted in the plot. The name of each
    #' element of the vector specifies its color group, i.e.: naming elements
    #' HIF1A and HIF1B as "HIF" to represent them with the same color.
    #' @param TF_colors (Optional) Nolors to highlight TFs chosen in specialTF.
    #' @return plotly scatter plot.
    #' @export plot_CM
    #' @examples
    #' data("stat_mat",package = "TFEA.ChIP")
    #' plot_CM(stat_mat)

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }
    requireNamespace("plotly")

    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}
    if (is.null(specialTF)){
        CM.statMatrix$highlight<-rep("Other",length(CM.statMatrix[,1]))
        markerColors<-c("azure4")
        names(markerColors)<-c("Other")
    }
    if (!is.null(specialTF) & is.null(TF_colors)){
        TF_colors<-c("red","blue","green","hotpink","cyan","greenyellow",
            "gold","darkorchid","chocolate1","black","lightpink","seagreen")
        TF_colors<-TF_colors[1:length(unique(names(specialTF)))]
        highlight_list<-highlight_TF(CM.statMatrix,2,specialTF,TF_colors)
        CM.statMatrix$highlight<-highlight_list[[1]]
        markerColors<-highlight_list[[2]]
    }
    if(!is.null(specialTF)&!is.null(TF_colors)){
        highlight_list<-highlight_TF(CM.statMatrix,2,specialTF,TF_colors)
        CM.statMatrix$highlight<-highlight_list[[1]]
        markerColors<-highlight_list[[2]]
    }
    if (!exists("MetaData")){
        data("MetaData",package = "TFEA.ChIP",envir = environment())
    }
    MetaData<-MetaData[MetaData$Accession%in%CM.statMatrix$Accession,]
    MetaData<-dplyr::arrange(MetaData,Accession)
    CM.statMatrix<-dplyr::arrange(CM.statMatrix,Accession)
    CM.statMatrix$Treatment<-MetaData$Treatment
    CM.statMatrix$Cell<-MetaData$Cell
    rm(MetaData)

    if(length(CM.statMatrix[CM.statMatrix$OR==Inf,1])>0){
        warn_number<-length(CM.statMatrix[CM.statMatrix$OR==Inf,1])
        CM.statMatrix[CM.statMatrix$OR==Inf,]$OR<-rep(
            max(CM.statMatrix[CM.statMatrix$OR!=Inf,]$OR),
            length(CM.statMatrix[CM.statMatrix$OR==Inf,1]))
        warning(warn_number," elements have an Odds Ratio of Inf.",
                " Maximum value for OR introduced instead.")
    }
    if(length(CM.statMatrix[CM.statMatrix$OR==-Inf,1])>0){
        warn_number<-length(CM.statMatrix[CM.statMatrix$OR==-Inf,1])
        CM.statMatrix[CM.statMatrix$OR==-Inf,]$OR<-rep(
            min(CM.statMatrix[CM.statMatrix$OR!=-Inf,]$OR),
            length(CM.statMatrix[CM.statMatrix$OR==-Inf,1]))
        warning(warn_number," elements have an Odds Ratio of -Inf. Minimum",
            " value for OR introduced instead.")
    }
    if(length(CM.statMatrix[CM.statMatrix$adj.p.value==0,1])>0){
        warn_number<-length(CM.statMatrix[CM.statMatrix$adj.p.value==0,1])
        CM.statMatrix[CM.statMatrix$p.value==0,]$log.adj.pVal<-rep(
            max(CM.statMatrix[CM.statMatrix$adj.p.value!=0,]$log.adj.pVal),
            length(CM.statMatrix[CM.statMatrix$adj.p.value==0,1]))
        warning(warn_number," elements have a -log(p-Value) of Inf. ",
            "Maximum value for -log(p-Val) introduced instead.")
    }

    if (length(markerColors)>1){
        CM.statMatrix_highlighted<-CM.statMatrix[
            CM.statMatrix$highlight!="Other",]
        CM.statMatrix_other<-CM.statMatrix[
            CM.statMatrix$highlight=="Other",]

        p<-plotly::plot_ly(CM.statMatrix_other, x=~log.adj.pVal,
            y=~OR,type="scatter", mode="markers",
            text=paste0(
                CM.statMatrix_other$Accession, ": ", CM.statMatrix_other$TF,
                '<br>Treatment: ', CM.statMatrix_other$Treatment,
                '<br>Cell: ', CM.statMatrix_other$Cell),
            color = ~highlight, colors=markerColors)

        p<-plotly::add_markers(p,x=CM.statMatrix_highlighted$log.adj.pVal,
            y=CM.statMatrix_highlighted$OR,type="scatter", mode="markers",
            text=paste0(
                CM.statMatrix_highlighted$Accession, ": ",
                CM.statMatrix_highlighted$TF,
                '<br>Treatment: ',CM.statMatrix_highlighted$Treatment,
                '<br>Cell: ',CM.statMatrix_highlighted$Cell),
            color = CM.statMatrix_highlighted$highlight, colors=markerColors
            )%>%
            plotly::layout(title=plot_title)

    }else if (length(markerColors)==1){
        p<-plotly::plot_ly(CM.statMatrix, x=~log.adj.pVal,y=~OR,type="scatter",
            mode="markers",
            text=paste0(
                CM.statMatrix$Accession,": ",CM.statMatrix$TF,
                '<br>Treatment: ',CM.statMatrix$Treatment,
                '<br>Cell: ',CM.statMatrix$Cell),
            color = ~highlight, colors=markerColors)%>%
        plotly::layout(title=plot_title)
    }
    p
    return(p)
}

plot_ES<-function(GSEA_result,LFC,plot_title = NULL,specialTF = NULL,TF_colors = NULL,Accession=NULL,TF=NULL){

    #' @title Plots Enrichment Score from the output of GSEA.run.
    #' @description Function to plot the Enrichment Score of every member of
    #' the ChIPseq binding database.
    #' @param GSEA_result Returned by GSEA_run
    #' @param LFC Vector with log2(Fold Change) of every gene that has an
    #' Entrez ID. Arranged from higher to lower.
    #' @param plot_title (Optional) Title for the plot
    #' @param specialTF (Optional) Named vector of TF symbols -as written in
    #' the enrichment table- to be highlighted in the plot. The name of each
    #' element specifies its color group, i.e.: naming elements HIF1A and HIF1B
    #' as "HIF" to represent them with the same color.
    #' @param TF_colors (Optional) Colors to highlight TFs chosen in specialTF.
    #' @param Accession (Optional) restricts plot to the indicated list dataset
    #' IDs.
    #' @param TF (Optional) restricts plot to the indicated list transcription
    #' factor names.
    #' @return Plotly object with a scatter plot -Enrichment scores- and a
    #' heatmap -log2(fold change) bar-.
    #' @export plot_ES
    #' @examples
    #' data("GSEA.result","log2.FC",package = "TFEA.ChIP")
    #' TF.hightlight<-c("EPAS1")
    #' names(TF.hightlight)<-c("EPAS1")
    #' col<- c("red")
    #' plot_ES(GSEA.result,log2.FC,specialTF = TF.hightlight,TF_colors = col)

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }
    requireNamespace("dplyr")
    requireNamespace("plotly")

    if(is.list(GSEA_result)==TRUE){
        enrichmentTable<-GSEA_result$Enrichment
    }else if(is.data.frame(GSEA_result)==TRUE){
        enrichmentTable<-GSEA_result
    }

    if (!is.null(Accession) | !is.null(TF)){
        if(is.null(Accession)){Accession<-enrichmentTable$Accession}
        if(is.null(TF)){TF<-enrichmentTable$TF}
        SS<-((enrichmentTable$Accession %in% Accession) &
            (enrichmentTable$TF %in% TF))
        enrichmentTable<-enrichmentTable[which(SS),]
    }

    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}

    if (is.null(specialTF)){
        highlight_list<-rep("Other",length(enrichmentTable[,1]))
        enrichmentTable$highlight<-highlight_list
        markerColors<-c("azure4")
        names(markerColors)<-c("Other")
    }

    if (!is.null(specialTF) & is.null(TF_colors)){
        TF_colors<-c("red","blue","green","hotpink","cyan","greenyellow",
            "gold","darkorchid","chocolate1","black","lightpink","seagreen")
        TF_colors<-TF_colors[1:length(unique(names(specialTF)))]
        highlight_list<-highlight_TF(enrichmentTable,2,specialTF,TF_colors)
        enrichmentTable$highlight<-highlight_list[[1]]
        markerColors<-highlight_list[[2]]
    }
    if(!is.null(specialTF)&!is.null(TF_colors)){
        highlight_list<-highlight_TF(enrichmentTable,2,specialTF,TF_colors)
        enrichmentTable$highlight<-highlight_list[[1]]
        markerColors<-highlight_list[[2]]
    }

    simbolo<-rep("pVal>0.05",times=length(enrichmentTable[,1]))
    for (i in 1:length(enrichmentTable[,1])){
        if(!is.na(enrichmentTable[i,4])){
            if (enrichmentTable[i,4]<=0.05){simbolo[i]<-"pVal<0.05"}
        }
    }
    enrichmentTable$symbol<-simbolo
    if (!exists("MetaData")){
        data("MetaData",package = "TFEA.ChIP",envir = environment())
    }
    MetaData<-MetaData[MetaData$Accession%in%enrichmentTable$Accession,]
    MetaData<-dplyr::arrange(MetaData,Accession)
    enrichmentTable<-dplyr::arrange(enrichmentTable,Accession)
    enrichmentTable$Treatment<-MetaData$Treatment
    enrichmentTable$Cell<-MetaData$Cell
    rm(MetaData)

    if(length(markerColors>1)){
        enrichmentTable_highlighted<-enrichmentTable[
            enrichmentTable$highlight!="Other",]
        enrichmentTable_other<-enrichmentTable[
            enrichmentTable$highlight=="Other",]

        p<-plotly::plot_ly(enrichmentTable_other,
            x=enrichmentTable_other$Arg.ES,
            y=enrichmentTable_other$ES, type="scatter", mode="markers",
            text=paste0(
                enrichmentTable_other$Accession,": ",enrichmentTable_other$TF,
                '<br>Pval: ',round(enrichmentTable_other$pval.ES,3),
                '<br>Treatment: ',enrichmentTable_other$Treatment,
                '<br>Cell: ',enrichmentTable_other$Cell),
            color=enrichmentTable_other$highlight, colors=markerColors,
            symbol=enrichmentTable_other$symbol, symbols=c("x","circle"))

        p<-plotly::add_markers(p,x=enrichmentTable_highlighted$Arg.ES,
            y=enrichmentTable_highlighted$ES,type="scatter", mode="markers",
            text=paste0(
                enrichmentTable_highlighted$Accession,": ",
                enrichmentTable_highlighted$TF,
                '<br>Pval: ',round(enrichmentTable_highlighted$pval.ES,3),
                '<br>Treatment: ',enrichmentTable_highlighted$Treatment,
                '<br>Cell: ',enrichmentTable_highlighted$Cell),
            color=enrichmentTable_highlighted$highlight, colors=markerColors,
            symbol=enrichmentTable_highlighted$symbol, symbols=c("x","circle")
        )%>%
        plotly::layout(title=plot_title,
            xaxis = list(title = "Argument"),
            yaxis = list (title = "ES"))
    }else if (length(markerColors)==1){
        p<-plotly::plot_ly(enrichmentTable, x=enrichmentTable$Arg.ES,
            y=enrichmentTable$ES, type="scatter", mode="markers",
            text=paste0(
                enrichmentTable$Accession,": ",enrichmentTable$TF,
                '<br>Pval: ',round(enrichmentTable$pval.ES,3),
                '<br>Treatment: ',enrichmentTable$Treatment,
                '<br>Cell: ',enrichmentTable$Cell),
            color=enrichmentTable$highlight, colors=markerColors,
            symbol=enrichmentTable$symbol, symbols=c("x","circle"))%>%
        plotly::layout(title=plot_title,
            xaxis = list(title = "Argument"),
            yaxis = list (title = "ES"))
    }
    LFC.bar<-get_LFC_bar(LFC)

    graf<-plotly::subplot(p, LFC.bar, shareX=TRUE, nrows=2,
        heights=c(0.95, 0.05), titleY=TRUE)
    graf
    return(graf)
}

plot_RES<-function(GSEA_result,LFC,plot_title = NULL,line.colors = NULL,line.styles = NULL,Accession=NULL,TF=NULL){

    #' @title Plots all the RES stored in a GSEA_run output.
    #' @description Function to plot all the RES stored in a GSEA_run output.
    #' @param GSEA_result Returned by GSEA_run
    #' @param LFC Vector with log2(Fold Change) of every gene that has an
    #' Entrez ID. Arranged from higher to lower.
    #' @param plot_title (Optional) Title for the plot.
    #' @param line.colors (Optional) Vector of colors for each line.
    #' @param line.styles (Optional) Vector of line styles for each line
    #' ("solid"/"dash"/"longdash").
    #' @param Accession (Optional) restricts plot to the indicated list dataset
    #' IDs.
    #' @param TF (Optional) restricts plot to the indicated list transcription
    #' factor names.
    #' @return Plotly object with a line plot -running sums- and a
    #' heatmap -log2(fold change) bar-.
    #' @export plot_RES
    #' @examples
    #' data("GSEA.result","log2.FC",package = "TFEA.ChIP")
    #' plot_RES(GSEA.result,log2.FC,TF=c("EPAS1"),
    #'     Accession=c("GSM2390642","GSM2390643"))

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work.",
            " Please install it.", call. = FALSE)
    }
    requireNamespace("utils")
    requireNamespace("dplyr")
    requireNamespace("plotly")

    if (!is.null(Accession) | !is.null(TF)){
        if(is.null(Accession)){
            Accession<-GSEA_result$Enrichment.table$Accession}
        if(is.null(TF)){TF<-GSEA_result$Enrichment.table$TF}
        GSEA_result$Enrichment.table<-GSEA_result$Enrichment.table[
            GSEA_result$Enrichment.table$Accession %in% Accession &
            GSEA_result$Enrichment.table$TF %in% TF,]
        GSEA_result$RES<-GSEA_result$RES[names(GSEA_result$RES) %in% Accession]
        GSEA_result$indicators<-GSEA_result$indicators[
            names(GSEA_result$indicators %in% Accession)]
    }else{
        Accession<-GSEA_result$Enrichment.table$Accession
    }

    if(is.null(line.colors)){
        line.colors<-c("red","blue","green","hotpink","cyan","greenyellow",
            "gold","darkorchid","chocolate1","black","lightpink","seagreen")
        line.colors<-line.colors[1:length(names(GSEA_result$RES))]
    }
    if(is.null(line.styles)){
        line.styles<-rep("solid",length(names(GSEA_result$RES)))
    }
    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}

    GSEA.runningSum<-GSEA_result$RES

    chip_index<-get_chip_index("g")
    tf<-chip_index[chip_index[,1]%in%Accession,]
    rm(chip_index)
    if (!exists("MetaData")){
        data("MetaData",package = "TFEA.ChIP",envir = environment())
    }
    MetaData<-MetaData[MetaData$Accession%in%Accession,]
    Cell<-vector()
    Treatment<-rep("no treatment",length(Accession))
    TF<-vector()
    RES<-list()
    for(i in 1:length(Accession)){
        Cell[i]<-MetaData[MetaData$Accession==Accession[i],]$Cell
        if(MetaData[MetaData$Accession==Accession[i],]$Treatment!="none"){
            Treatment[i]<-MetaData[MetaData$Accession==Accession[i],5]
        }
        TF[i]<-tf[tf[,1]==Accession[i],2]
        RES<-c(RES,GSEA.runningSum[names(GSEA.runningSum)%in%Accession[i]][1])
    }
    tabla<-data.frame(Accession,Cell,Treatment,TF,stringsAsFactors = FALSE)
    tabla$RES<-RES

    rm(Cell,Treatment,TF,RES)

    if(length(tabla[,1])>1){
        for(i in 1:length(Accession)){
            if (i==1){
                grafica<-plotly::plot_ly(tabla, x=c(1:length(tabla$RES[[1]])),
                    y=tabla$RES[[Accession[1]]], type="scatter", mode="lines",
                    line=list(color=line.colors[1], dash=line.styles[1]),
                    name=paste0(tabla$Accession[1]," - ", tabla$TF[1]),
                    text=paste0(
                        tabla$Accession[1]," - ", tabla$TF[1],
                        '<br>Cell: ',tabla$Cell[1],
                        '<br>Treatment: ', tabla$Treatment[1]))
            }else if (i>1 & i<length(Accession)){
                grafica<-plotly::add_trace(p = grafica,
                    y=tabla$RES[[Accession[i]]],
                    type="scatter", mode="lines",
                    line=list(color=line.colors[i], dash=line.styles[i]),
                    name=paste0(tabla$Accession[i]," - ", tabla$TF[i]),
                    text=paste0(
                        tabla$Accession[i]," - ", tabla$TF[i],
                        '<br>Cell: ', tabla$Cell[i],
                        '<br>Treatment: ', tabla$Treatment[i]))
            }else if (i==length(Accession)){
                grafica<-plotly::add_trace(p = grafica,
                    y=tabla$RES[[Accession[i]]], type="scatter", mode="lines",
                    line=list(color=line.colors[i], dash=line.styles[i]),
                    name=paste0(tabla$Accession[i]," - ", tabla$TF[i]),
                    text=paste0(
                        tabla$Accession[i]," - ", tabla$TF[i],
                        '<br>Cell: ', tabla$Cell[i],
                        '<br>Treatment: ', tabla$Treatment[i]))%>%
                plotly::layout(title=plot_title,
                    xaxis = list(title = "Argument"),
                    yaxis = list (title = "ES"))
            }
        }
    }else{
        grafica<-plotly::plot_ly(tabla,x=c(1:length(tabla$RES[[1]])),
            y=tabla$RES[[Accession[1]]], type="scatter", mode="lines",
            line=list(color=line.colors[1], dash=line.styles[1]),
            name=paste0(tabla$Accession[1]," - ", tabla$TF[1]),
            text=paste0(
                tabla$Accession[i]," - ", tabla$TF[i],
                '<br>Cell: ', tabla$Cell[i],
                '<br>Treatment: ', tabla$Treatment[i]))%>%
        plotly::layout(title=plot_title,
            xaxis = list(title = "Argument"),
            yaxis = list (title = "ES"))
    }
    LFC.bar<-get_LFC_bar(LFC)

    graf<-plotly::subplot(grafica, LFC.bar, shareX = TRUE,
        nrows = 2, heights = c(0.95, 0.05), titleY = TRUE)
    graf
    return(graf)
}

highlight_TF<-function(table,column,specialTF,markerColors){

    #' @title Highlight certain transcription factors in a plotly graph.
    #' @description Function to highlight certain transcription factors using
    #' different colors in a plotly graph.
    #' @param table Enrichment matrix/data.frame.
    #' @param column Column # that stores the TF name in the matrix/df.
    #' @param specialTF Named vector containing TF names as they appear in the
    #' enrichment matrix/df and nicknames for their color group.
    #' Example:
    #'           specialTF<-c("HIF1A","EPAS1","ARNT","SIN3A")
    #'           names(specialTF)<-c("HIF","HIF","HIF","SIN3A")
    #' @param markerColors Vector specifying the shade for every color group.
    #' @return List of two objects:
    #' A vector to attach to the enrichment matrix/df pointing out the color
    #' group of every row.
    #' A named vector connecting each color group to the chosen color.
    # examples
    # highlight_TF(CM.statMatrix_UP,4,specialTF,colors)


    highlight<-rep("Other",length(table[,1]))
    for (i in 1:length(specialTF)){
        for (j in 1:length(table[,1])){
            if(!is.na(table[j,column])){
                if (table[j,column]==specialTF[i]){
                    highlight[j]<-names(specialTF)[i]
                }
            }
        }
    }
    markerColors<-c("azure4",markerColors)
    names(markerColors)<-c("Other",unique(names(specialTF)))
    return(list(highlight,markerColors))
}

get_LFC_bar<-function(LFC){

    #' @title Plots a color bar from log2(Fold Change) values.
    #' @description Function to plot a color bar from log2(Fold Change)
    #' values from an expression experiment.
    #' @param LFC Vector of log2(fold change) values arranged from higher
    #' to lower. Use ony the values of genes that have an Entrez ID.
    #' @return Plotly heatmap plot -log2(fold change) bar-.
    # examples
    # get_LFC_bar(arranged.log2FC.array)

    if(!requireNamespace("scales", quietly = TRUE)){
        stop("scales package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }
    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. ",
            "Please install it.", call. = FALSE)
    }
    requireNamespace("grDevices")
    requireNamespace("dplyr")
    requireNamespace("plotly")
    requireNamespace("scales")

    vals <- scales::rescale(LFC)
    o <- order(vals, decreasing = FALSE)
    cols1 <- scales::col_numeric(grDevices::colorRamp(c("mistyrose","red3")),
        domain = NULL)(vals[1:length(LFC[LFC>0])])
    cols2 <- scales::col_numeric(grDevices::colorRamp(c("navy","lightcyan")),
        domain = NULL)(vals[length(LFC[LFC>0])+1:length(LFC)])
    cols<-c(cols1,cols2)

    colz <-data.frame(vals[o], cols[o])

    LFC.bar<-plotly::plot_ly(x=c(1:length(LFC)), y=rep(1,length(LFC)),
        z = LFC, type = "heatmap",colorscale=colz,showscale = FALSE)%>%
    plotly::layout(yaxis=list(visible=FALSE))

    return(LFC.bar)
}
