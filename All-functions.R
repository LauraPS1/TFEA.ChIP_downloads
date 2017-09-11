################### FUNCTIONS ####

txt2GR<-function(fileTable,format,GRfolder,fileMetaData){

    #' @title Function to filter a ChIP-Seq output.
    #' @description Function to filter a ChIP-Seq output (in .narrowpeak or MACS's peaks.bed formats) and
    #' then store the peak coordinates in a GenomicRanges object, associated to its metadata.
    #' @param fileTable data frame from a txt/tsv/bed file
    #' @param format "narrowPeak" or "macs".
    #' narrowPeak fields:'chrom','chromStart','chromEnd','name','score','strand','signalValue','pValue','qValue','peak'
    #' macs fields: 'chrom','chromStart','chromEnd','name','qValue'
    #' @param GRfolder Path to the folder to save the GR objects
    #' @param MetaData Data frame/matrix/array contaning the following fields: 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'.
    #' @return The function saves the GR object generated in a .Rdata file in GRfolder.
    #' @export txt2GR
    #' @examples
    #' txt2GR(smad4_peaks.bed,"macs","~/folder/GR",SMAD4_metaData)

    requireNamespace("GenomicRanges")
    requireNamespace("IRanges")
    requireNamespace("dplyr")

    if (is.data.frame(fileMetaData)==F){
        if(is.matrix(fileMetaData)){
            if(length(fileMetaData[1,])==7){
                fileMetaData<-as.data.frame(fileMetaData,stringsAsFactors=F)
            }else{
                warning("fileMetaData format error: 'fileMetaData' must be a data frame/matrix/array
                        with 7 atributes: 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'")
                break
            }
        }else if(is.array(fileMetaData)){
            if(length(fileMetaData)==7){
                fileMetaData<-as.data.frame(fileMetaData,stringsAsFactors=F)
            }else{
                warning("fileMetaData format error: 'fileMetaData' must be a data frame/matrix/array
                        with 7 atributes: 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'")
                break
            }
        }else{
            warning("fileMetaData format error: 'fileMetaData' must be a data frame/matrix/array
                    with 7 atributes: 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'")
            break
        }
        }
    format<-tolower(format)

    if(format=="narrowpeak"){

        if(fileTable[1,8]==-1 & fileTable[1,9]==-1){
            warning ("The ChIP-Seq input file ",fileMetaData$Name," does not include p-value or Q-value for each peak. Please, make sure the peaks in the input file have been previously filtered according to their significance")
            fileTable<-dplyr::select(fileTable,V1,V2,V3)
            Stat<-"no score"
            fileTable$score=rep(NA,length(fileTable[,1]))
            colnames(fileTable)[1:3]<-c("chr","start","end")
        }else if(fileTable[1,8]==-1){
            fileTable<-dplyr::select(fileTable,V1,V2,V3,V9)
            colnames(fileTable)<-c("chr","start","end","score")
            Stat<-"log10(p-Value)"
            valLimit<-1.3
            fileTable<-fileTable[fileTable$score>valLimit,]
        }else if(fileTable[1,9]==-1){
            fileTable<-dplyr::select(fileTable,V1,V2,V3,V8)
            colnames(fileTable)<-c("chr","start","end","score")
            Stat<-"p-Value"
            valLimit<-0.05
            fileTable<-fileTable[fileTable$score<valLimit,]
        }

        fileMetaData<-c(fileMetaData,Stat)

        MDframe<-as.data.frame(lapply(fileMetaData, rep,length(fileTable[,1])))
        colnames(MDframe)<-c("Name","Accession","Cell","Cell Type","Treatment","Antibody","TF","Score Type")

        gr<-GenomicRanges::GRanges(

            seqnames=fileTable$chr,
            ranges=IRanges::IRanges(fileTable$start,end=fileTable$end),
            score=fileTable$score,
            mcols=MDframe
        )

        save(gr,file = paste0(GRfolder,"/",fileMetaData$Accession,".Rdata"))

    }else if(format=="macs"){

        if(length(fileTable[1,])==5){
            fileTable<-dplyr::select(fileTable,V1,V2,V3,V5)
            colnames(fileTable)<-c("chr","start","end","score")
            fileTable<-fileTable[fileTable$score>50,]
            Stat<-"10*log10(p-Value)"
        }else if (length(fileTable[1,])==4 & is.character(fileTable[1,4])){ # if the 4th column consists of peak names
            warning ("The ChIP-Seq input file ",fileMetaData$Name," does not include p-value or Q-value for each peak. Please, make sure the peaks in the input file have been previously filtered according to their significance")
            fileTable<-dplyr::select(fileTable,V1,V2,V3)
            Stat<-"no score"
            fileTable$score=rep(NA,length(fileTable[,1]))
            colnames(fileTable)[1:3]<-c("chr","start","end")
        }else if (length(fileTable[1,])==4 & !is.character(fileTable[1,4])){ # if the 4th column consists of adjusted p-values
            colnames(fileTable)<-c("chr","start","end","score")
            fileTable<-fileTable[fileTable$score>50,]
            Stat<-"10*log10(p-Value)"
        }
        
        fileMetaData<-c(fileMetaData,Stat)
        
        MDframe<-as.data.frame(lapply(fileMetaData, rep,length(fileTable[,1])))
        colnames(MDframe)<-c("Name","Accession","Cell","Cell Type","Treatment","Antibody","TF","Score Type")

        gr<-GenomicRanges::GRanges(

            seqnames=fileTable$chr,
            ranges=IRanges::IRanges(fileTable$start,end=fileTable$end),
            score=fileTable$score,
            mcols=MDframe
        )

        save(gr,file = paste0(GRfolder,"/",fileMetaData$Accession,".Rdata"))

    }else{
        warning("Wrong file format. Only narrowPeak or MACS output ('_peaks.bed') are supported.
              Please choose either 'narrowPeak' or 'Macs'.")
        break
    }
}

GR2id_db<-function(Dnase.db,gr.list,GRfolder){

    #' @title Makes a TF-gene binding database
    #' @description GR2id_db generates a TF-gene binding database from ChIP-Seq peak coordinates
    #' in GenomicRange objects.
    #' @param Dnase.db GenomicRanges object containing a database of Dnase hipersensitive sites
    #' @param gr.list Vector of paths to .Rdata files containing each one a GenomicRanges object (output of txt2GR)
    #' @param GRfolder Path to the folder storing the .Rdata files.
    #' @return List of vectors, one for every ChIP-Seq, storing the IDs of the genes to which the TF bound in the ChIP-Seq.
    #' @export GR2id_db
    #' @examples
    #' GR2id_db(Dnase.db=DnaseHS.db, gr.list=dir("~/folder/GR"),GRfolder="~/folder/GR")

    if(!requireNamespace("S4Vectors", quietly = TRUE)){
        stop("S4Vectors package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    requireNamespace("S4Vectors")

    m<-0

    for(i in 1:length(gr.list)){

        load(paste0(GRfolder,"/",gr.list[i]))

        index_cercanos<-suppressWarnings(GenomicRanges::distanceToNearest(gr,Dnase.db,select="all"))
        index_cercanos<-index_cercanos[!is.na(index_cercanos@elementMetadata@listData$distance)]
        index_cercanos<-index_cercanos[index_cercanos@elementMetadata@listData$distance<10]
        inSubject<-S4Vectors::subjectHits(index_cercanos)

        if(length(inSubject)==0){
            m<-m+1
            next
        }

        listaID<-Dnase.db[inSubject]$gene_id


        if (i==1){
            ID_list<-list(listaID)
            names(ID_list)[i]<-substr(gr.list[i],1,nchar(gr.list[i])-6)
            rm(index_cercanos,inSubject,listaID)
        }else{
            ID_list<-c(ID_list,list(listaID))
            names(ID_list)[i-m]<-substr(gr.list[i],1,nchar(gr.list[i])-6)
            rm(index_cercanos,inSubject,listaID)
        }
    }

    return(ID_list)
}

SearchID<-function(GeneID,id_db){

    #' @title Searchs for an Entrez gene ID in a TF-gene database
    #' @description Searchs for an Entrez gene ID in a db that contains the IDs of the genes a TF binds to.
    #' @param GeneID gene Entrez ID.
    #' @param id_db TF - gene binding database.
    #' @return 1/0 row. Each element represents a TF ChIP-Seq experiment.

    # examples
    # SearchID(GeneID,id_db)

    TF_row<-rep(0,length=length(id_db))

    for (i in 1:length(id_db)){
        if (GeneID %in% id_db[[i]]==T){TF_row[i]<-1}
    }
    return(TF_row)
}

SearchIDlist<-function(GeneList,id_db){

    #' @title Function to search for a list of entrez gene IDs.
    #' @description Function to search for a list of entrez gene IDs in a  TF-gene binding data base.
    #' @param GeneList Array of gene Entrez IDs
    #' @param id_db TF - gene binding database.
    #' @return 1/0 matrix. Each row represents a gene, each column, a ChIP-Seq file.
    #' @export SearchIDlist
    #' @examples
    #' SearchIDlist(GeneList,id_db)

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
    #' @description Function to set the data objects provided by the user as default to the rest of the functions.
    #' @param metadata Data frame/matrix/array contaning the following fields: 'Name','Accession','Cell','Cell Type','Treatment','Antibody','TF'.
    #' @param binary_matrix Matrix[n,m] which rows correspond to all the human genes that have been assigned an Entrez ID, and its columns, to every ChIP-Seq experiment in the database.
    #' The values are 1 – if the ChIP-Seq has a peak assigned to that gene – or 0 – if it hasn’t –.
    #' @export set_user_data
    #' @examples
    #' set_user_data(MetaData,Mat01)

    assign("MetaData", metadata, envir = .GlobalEnv)
    assign("Mat01",binary_matrix,envir = .GlobalEnv)
}

get_chip_index<-function(database = "g",TFfilter = NULL){

    #' @title Creates a data frame containing accession IDs of ChIP-Seq experiments and TF tested
    #' @description Function to create a data frame containing accession IDs of ChIP-Seq experiments and
    #' the transcription factor tested in each ChIP. This index is required for functions
    #' like “contingency_matrix” and “GSEA.run” and can also be used as a filter to select
    #' specific ChIPs or transcription factors to run an analysis.
    #' @param database Name of the database used: "encode"/"e" or "general"/"g".
    #' @param TFfilter (Optional) Transcription factors of interest.
    #' @return Data frame containig the accession # and TF for every ChIP-Seq experiment included
    #' in the metadata files.
    #' @export get_chip_index
    #' @examples
    #' get_chip_index(database="e")
    #' get_chip_index(database="general",TFfilter=c("SMAD2","SMAD4"))

    requireNamespace("dplyr")
    requireNamespace("utils")

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
    #' @description Function to compute contingency 2x2 matrix based on the matches between two gene ID lists
    #' and a ChIP-Seq binding database.
    #' @param test_list List of gene Entrez IDs
    #' @param control_list If not provided, all human genes not present in test_list will be used as control.
    #' @param chip_index Output of the function “get_chip_index”, a data frame containing
    #' accession IDs of ChIPs on the database and the TF each one tests. If not provided, the whole internal database will be used
    #' @return List of contingency matrices, one CM per element in chip_index (i.e. per ChIP-seq dataset).
    #' @export contingency_matrix
    #' @examples
    #' contingency_matrix(Genes.Upreg)

    requireNamespace("utils")

    if (missing(control_list)){ # Generating control gene list in case is not provided.
        suppressMessages(require("TxDb.Hsapiens.UCSC.hg19.knownGene",quietly = T))
        Genes<-GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)$gene_id
        control_list<-Genes[!(Genes %in% test_list)]
        rm(Genes)
    }else{
        control_list<-control_list[!(control_list %in% test_list)]
    }

    Matrix1<-Mat01[rownames(Mat01)%in%test_list,colnames(Mat01)%in%chip_index$Accession]
    Matrix2<-Mat01[rownames(Mat01)%in%control_list,colnames(Mat01)%in%chip_index$Accession]


    for (i in 1:length(chip_index[,1])){

        chip.vector1<-Matrix1[,chip_index$Accession[i]]
        chip.vector2<-Matrix2[,chip_index$Accession[i]]

        pos1<-length(chip.vector1[chip.vector1==1])
        pos2<-length(chip.vector2[chip.vector2==1])
        neg1<-length(chip.vector1[chip.vector1==0])
        neg2<-length(chip.vector2[chip.vector2==0])
        CM_gr<-cbind(c(pos1,pos2),c(neg1,neg2))
        rownames(CM_gr)<-c("Test","Control")
        colnames(CM_gr)<-c("Positive","Negative")

        if (i==1){
            CM_list<-list(CM_gr)
            names(CM_list)[i]<-as.character(chip_index$Accession[i])
            rm(CM_gr)
        }else{
            CM_list<-c(CM_list,list(CM_gr))
            names(CM_list)[i]<-as.character(chip_index$Accession[i])
            rm(CM_gr)
        }
    }

    return(CM_list)
}

getCMstats<-function(CM_list,chip_index=get_chip_index()){

    #' @title Generates a data frame with Accession, TF, OR and p-val from a contingency_matrix output
    #' @description From a list of contingency matrices, such as the output from “contingency_matrix”, this function computes a fisher's exact test for each matrix and generates a data frame that stores
    #' accession ID of a ChIP-Seq experiment, the TF tested in that experiment, and the p-value and the odds ratio resulting
    #' from the test.
    #' @param CM_list Output of “contingency_matrix”, a list of contingency matrix.
    #' @param chip_index Output of the function “get_chip_index”, a data frame containing
    #' accession IDs of ChIPs on the database and the TF each one tests. If not provided, the whole internal database will be used
    #' @return Data frame containing accession ID of a ChIP-Seq experiment, the TF tested
    #' in that experiment, raw p-value (-10*log10 pvalue), odds-ratio and FDR-adjusted p-values (-10*log10 adj.pvalue).
    #' @export getCMstats
    #' @examples
    #' getCMstats(CM_list_UP)

    requireNamespace("stats")

    for (i in 1:length(CM_list)){
        pval<-stats::fisher.test(x=CM_list[[i]])
        if (i==1){
            pval_list<-list(pval)
            names(pval_list)[i]<-names(CM_list)[i]
        }else{
            pval_list<-c(pval_list,list(pval))
            names(pval_list)[i]<-names(CM_list)[i]
        }
    }

    statMat<-data.frame(Accession=chip_index$Accession,TF=chip_index$TF,p.value=NA,OR=NA)
    for (idx in names(CM_list)){
        FTres<-try({stats::fisher.test(x=CM_list[[idx]])},silent = T)
        if(class(FTres)=="htest"){
            statMat$p.value[which(statMat$Accession==idx)]<-FTres$p.value
            statMat$OR[which(statMat$Accession==idx)]<-FTres$estimate
        }
    }
    statMat$adj.p.value<-stats::p.adjust(statMat$p.value,"fdr")
    statMat$log.adj.pVal<-(-1*(log(statMat$adj.p.value)))
    return(statMat)

}

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 0, correl.vector = NULL) {

    # Computes the weighted GSEA score of gene.set in gene.list.
    # The weighted score type is the exponent of the correlation
    # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is
    # necessary to input the correlation vector with the values in the same order as in the gene list.
    #
    # Inputs:
    #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)
    #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset)
    #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
    #   correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
    #
    # Outputs:
    #   ES: Enrichment score (real number between -1 and +1)
    #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
    #   RES: Numerical vector containing the running enrichment score for all locations in the gene list
    #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
    #
    # The Broad Institute
    # SOFTWARE COPYRIGHT NOTICE AGREEMENT
    # This software and its documentation are copyright 2003 by the
    # Broad Institute/Massachusetts Institute of Technology.
    # All rights are reserved.
    #
    # This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for
    # its use, misuse, or functionality.

    #' @title Computes the weighted GSEA score of gene.set in gene.list.
    #' @description Computes the weighted GSEA score of gene.set in gene.list.
    #' @param gene.list The ordered gene list
    #' @param gene.set A gene set, e.g. gene IDs corresponding to a ChIP-Seq experiment's peaks.
    #' @param weighted.score.type Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)
    #' @param correl.vector A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list
    #' @return list of:
    #' ES: Enrichment score (real number between -1 and +1)
    #' arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain")
    #' RES: Numerical vector containing the running enrichment score for all locations in the gene list
    #' tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list
    #' @export GSEA.EnrichmentScore
    #' @examples
    #' GSEA.EnrichmentScore(gene.list, ARNT.ChIP.genes, weighted.score.type = 0, correl.vector = NULL)

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
        RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
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
        return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))

    }else{
        RES<-rep(NaN,length(gene.list))
        indicator<-rep(0,length(gene.list))
        ES<-NA
        arg.ES<-NA
        return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
    }
}

GSEA.Shuffling<-function(gene.list,permutations){

    #' @title Function to create shuffled gene lists to run GSEA.
    #' @description Function to create shuffled gene lists to run GSEA.
    #' @param gene.list Vector of gene Entrez IDs.
    #' @param permutations Number of shuffled gene vector IDs required.
    #' @return Vector of randomly arranged Entrez IDs.
    # examples
    # GSEA.Shuffling(gene.list,1000)

    shuffledGL<-list()
    for(i in 1:permutations){
        tmp<-sample(gene.list)
        shuffledGL<-c(shuffledGL,list(tmp))
    }
    return(shuffledGL)
}

GSEA.run<-function(gene.list,chip_index=get_chip_index(),get.RES = F,RES.filter = NULL){

    #' @title Function to run a GSEA analysis with an ID database.
    #' @description Function to run a GSEA analysis with an TF-gene binding database.
    #' @param gene.list List of Entrez IDs ordered by their fold change.
    #' @param chip_index Data.frame indicating the TF corresponding to each ChIP-Seq experiment in the database.
    #' @param get.RES (Optional) boolean. If TRUE, the function stores RES of all/some TF.
    #' @param RES.filter (Optional) chr vector. When get.RES==TRUE, allows to choose which TF's RES to store.
    #' @return a list of:
    #' Enrichment.table: data frame containing accession ID, TF name, enrichment score, p-value, and argument of every ChIP-Seq experiment.
    #' RES (optional): list of running sums of every ChIP-Seq
    #' indicators (optional): list of 0/1 vectors that stores the matches (1) and mismatches (0) between the gene list and the gene set.
    #' @export GSEA.run
    #' @examples
    #' GSEA.run(Arranged.Gene.List, chip.index)
    #' GSEA.run(Arranged.Gene.List, chip.index, get.RES = T, RES.filter =c("SMAD2","SMAD4") )

    requireNamespace("stats")
    requireNamespace("utils")

    Mat01<-Mat01[,colnames(Mat01)%in%chip_index$Accession]

    shuffledGL<-GSEA.Shuffling(gene.list,1000)  # Generate random gene lists to
                                                # get a p-value for ESs.
    EnrichmentScore<-vector()
    pval<-vector()
    EnrichmentArg<-vector()

    if(get.RES==T){
        res<-list()
        ind<-list()
    }

    for (i in 1:length(chip_index$Accession)){

        chip.genes<-Mat01[,colnames(Mat01)==chip_index$Accession[i]]
        chip.genes<-names(chip.genes[chip.genes==1])

        if(length(chip.genes)>10){
            resultado<-GSEA.EnrichmentScore(gene.list, chip.genes)

            shuffled.ES<-vector()
            for (j in 1:length(shuffledGL)){ # Get ES for the random gene list to compute pval.
                shuffled.ES[j]<-GSEA.EnrichmentScore(shuffledGL[[j]], chip.genes)$ES
            }

            EnrichmentScore[i]<-resultado$ES
            pval[i]<-sum(abs(shuffled.ES) >= abs(resultado$ES)) / 1000
            EnrichmentArg[i]<-resultado$arg.ES

            if(get.RES==T & missing(RES.filter)){ # Store running sums of selected TFs.
                res<-c(res,list(resultado$RES))
                names(res)[length(res)]<-chip_index$Accession[i]
                ind<-c(ind,list(resultado$indicator))
                names(ind)[length(ind)]<-chip_index$Accession[i]
            }else if(get.RES==T & missing(RES.filter)==F){
                if(chip_index$TF[i]%in%RES.filter){
                    res<-c(res,list(resultado$RES))
                    names(res)[length(res)]<-chip_index$Accession[i]
                    ind<-c(ind,list(resultado$indicator))
                    names(ind)[length(ind)]<-chip_index$Accession[i]
                }
            }
        }else{chip_index<-chip_index[-i,]}
    }
    pval.adj<-stats::p.adjust(pval,"fdr") # Adjust pvalues

    tablaEnriquecimiento<-cbind(chip_index$Accession,chip_index$TF,as.numeric(EnrichmentScore),
                                as.numeric(pval.adj),as.numeric(EnrichmentArg))

    tablaEnriquecimiento<-as.data.frame(tablaEnriquecimiento,stringsAsFactors=F)
    colnames(tablaEnriquecimiento)<-c("Accession","TF","ES","pval.ES","Arg.ES")
    tablaEnriquecimiento$ES<-as.numeric(tablaEnriquecimiento$ES)
    tablaEnriquecimiento$pval.ES<-as.numeric(tablaEnriquecimiento$pval.ES)
    tablaEnriquecimiento$Arg.ES<-as.numeric(tablaEnriquecimiento$Arg.ES)

    tablaEnriquecimiento<-tablaEnriquecimiento[!is.na(tablaEnriquecimiento$pval.ES),]

    if(get.RES==T){
        resultados<-list(tablaEnriquecimiento,res,ind)
        names(resultados)<-c("Enrichment.table","RES","indicators")
        return(resultados)
    }else{
        return(tablaEnriquecimiento)
    }
}

GeneID2entrez<-function(gene.IDs,return.Matrix = F){

    #' @title Translates gene IDs from Gene Symbol or Ensemble Gene ID to Entrez Gene ID.
    #' @description Translates gene IDs from Gene Symbol or Ensemble Gene ID to Entrez Gene ID using the IDs approved by HGNC.
    #' When translating from Gene Symbol, keep in mind that many genes have been given more than one symbol through the years.
    #' This function will return the Entrez ID corresponding to the currently approved symbols if they exist, otherwise NA is returned
    #' In addition some genes might map to more than one Entrez ID, in this case gene is assigned to the first match and a warning is displayed.
    #' @param gene.IDs Array of Gene Symbols or Ensemble Gene IDs.
    #' @param return.Matrix T/F. When TRUE, the function returns a matrix[n,2], one column with the gene symbols or Ensemble IDs, another with their respective Entrez IDs.
    #' @return Vector or matrix containing the Entrez IDs(or NA) corresponding to every element of the input.
    #' @export GeneID2entrez
    #' @examples
    #' GeneID2entrez(c("KLHL13","PROM1","ZNF663P","GYS1","LPHN3","ERO1L","EGLN3"))
    #' GeneID2entrez(c("ENSG00000121410","ENSG00000156006","ENSG00000196839"),return.Matrix=T)

    requireNamespace("biomaRt")
    requireNamespace("GenomicFeatures")

    Genes<-GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    suppressMessages(GeneNames<-biomaRt::select(org.Hs.eg.db, Genes$gene_id, c("SYMBOL", "ENSEMBL")))#suppressWarnings added to avoid 'select()' returned 1:many mapping between keys and columns

    if(grepl("ENSG0",gene.IDs[1])==T){ID.type<-"ENSEMBL"
    }else{ID.type<-"SYMBOL"}

    tmp<-match(as.character(gene.IDs),GeneNames[,ID.type])
    tmp2<-match(GeneNames[,ID.type],as.character(gene.IDs))

    if(sum(duplicated(tmp2[!is.na(tmp2)]))>0){warning("Some genes returned 1:many mapping to ENTREZ ID. Genes were assigned the first ENTREZ ID match found.",call. =F)}

    cat("Done! ",length(tmp[!is.na(tmp)])," genes of ",length(tmp)," successfully translated.\n")
    if (length(tmp[is.na(tmp)])>0){cat("Couldn't find Entrez IDs for ",length(tmp[is.na(tmp)])," genes (NAs returned instead).\n")}

    if(return.Matrix==T){
        return(data.frame(GENE.ID=gene.IDs,ENTREZ.ID=GeneNames[tmp,"ENTREZID"]))
    }else{
        return(GeneNames[tmp[!is.na(tmp)],"ENTREZID"])
    }
}

highlight.TF<-function(table,column,specialTF,colores){

    #' @title Highlight certain transcription factors in a plotly graph.
    #' @description Function to highlight certain transcription factors using different colors in a plotly graph.
    #' @param table Enrichment matrix/data.frame.
    #' @param column Column # that stores the TF name in the matrix/df.
    #' @param specialTF Named vector containing TF names as they appear in the enrichment matrix/df and nicknames for their color group.
    #' Example:
    #'           specialTF<-c("HIF1A","HIF1A-hx","EPAS1","EPAS1-hx","ARNT","ARNT-hx","SIN3A","SAP30","MXI1")
    #'           names(specialTF)<-c("HIF","HIF","HIF","HIF","HIF","HIF","SIN3A","SAP30","MXI1")
    #' @param colores Vector specifying the shade for every color group.
    #' @return List of two objects:
    #' A vector to attach to the enrichment matrix/df pointing out the color group of every row.
    #' A named vector connecting each color group to the chosen color.
    # examples
    # highlight.TF(pval_mat_UP,4,specialTF,colors)


    highlight<-rep("Other",length(table[,1]))
    for (i in 1:length(specialTF)){

        for (j in 1:length(table[,1])){

            if(!is.na(table[j,column])){
                if (table[j,column]==specialTF[i]){highlight[j]<-names(specialTF)[i]}
            }
        }
    }
    colores<-c("azure4",colores)
    names(colores)<-c("Other",unique(names(specialTF)))
    return(list(highlight,colores))
}

plot_CM<-function(CM.statMatrix,plot_title = NULL,specialTF = NULL,TF_colors = NULL){

    #' @title Makes an interactive html plot from an enrichment table.
    #' @description Function to make an interactive html plot from a transcription
    #' factor enrichment table, output of the function "getPvalMat".
    #' @param CM.statMatrix Output of the function "getCMstats".
    #' A data frame storing: Accession ID of every ChIP-Seq tested, Transcription Factor,Odds Ratio, p-value and adjusted p-value.
    #' @param plot_title The title for the plot.
    #' @param specialTF (Optional) Named vector of TF symbols -as written in the enrichment table- to be highlighted in the plot.
    #' The name of each element of the vector specifies its color group, i.e.: naming elements HIF1A and HIF1B as "HIF" to represent them with the same color.
    #' @param TF_colors (Optional) Nolors to highlight TFs chosen in specialTF.
    #' @return plotly scatter plot.
    #' @export plot_CM
    #' @examples
    #' plot_CM(pval_mat_UP, "Enriquecimiento TF en GSE69303 (upreg)", specialTF,colors)
    #' plot_CM(pval_mat_UP)

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    requireNamespace("plotly")

    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}
    if (is.null(specialTF)){
        CM.statMatrix$highlight<-rep("Other",length(CM.statMatrix[,1]))
        colores<-c("azure4")
        names(colores)<-c("Other")
    }
    if (!is.null(specialTF) & is.null(TF_colors)){
        TF_colors<-c("red","blue","green","hotpink","cyan","greenyellow","gold",
                     "darkorchid","chocolate1","black","lightpink","seagreen")
        TF_colors<-TF_colors[1:length(unique(names(specialTF)))]
        highlight_list<-highlight.TF(CM.statMatrix,2,specialTF,TF_colors)
        CM.statMatrix$highlight<-highlight_list[[1]]
        colores<-highlight_list[[2]]
    }
    if(!is.null(specialTF)&!is.null(TF_colors)){
        highlight_list<-highlight.TF(CM.statMatrix,2,specialTF,TF_colors)
        CM.statMatrix$highlight<-highlight_list[[1]]
        colores<-highlight_list[[2]]
    }

    MetaData<-MetaData[MetaData$Accession%in%CM.statMatrix$Accession,]
    MetaData<-dplyr::arrange(MetaData,Accession)
    CM.statMatrix<-dplyr::arrange(CM.statMatrix,Accession)
    CM.statMatrix$Treatment<-MetaData$Treatment
    CM.statMatrix$Cell<-MetaData$Cell
    rm(MetaData)
    
    if(length(CM.statMatrix[CM.statMatrix$OR==Inf,1])>0){
        warn<-length(CM.statMatrix[CM.statMatrix$OR==Inf,1])
        CM.statMatrix[CM.statMatrix$OR==Inf,]$OR<-rep(max(CM.statMatrix[CM.statMatrix$OR!=Inf,]$OR),
                                                      length(CM.statMatrix[CM.statMatrix$OR==Inf,1]))
        warning(warn," elements have an Odds Ratio of Inf. Maximum value for OR introduced instead.")
    }
    if(length(CM.statMatrix[CM.statMatrix$OR==-Inf,1])>0){
        warn<-length(CM.statMatrix[CM.statMatrix$OR==-Inf,1])
        CM.statMatrix[CM.statMatrix$OR==-Inf,]$OR<-rep(min(CM.statMatrix[CM.statMatrix$OR!=-Inf,]$OR),
                                                       length(CM.statMatrix[CM.statMatrix$OR==-Inf,1]))
        warning(warn," elements have an Odds Ratio of -Inf. Minimum value for OR introduced instead.")
    }
    if(length(CM.statMatrix[CM.statMatrix$adj.p.value==0,1])>0){
        warn<-length(CM.statMatrix[CM.statMatrix$adj.p.value==0,1])
        CM.statMatrix[CM.statMatrix$adj.p.value==0,]$log.adj.pVal<-rep(max(CM.statMatrix[CM.statMatrix$adj.p.value!=0,]$log.adj.pVal),
                                                                  length(CM.statMatrix[CM.statMatrix$adj.p.value==0,1]))
        warning(warn," elements have a -log(p-Value) of Inf. Maximum value for -log(p-Val) introduced instead.")
    }
    
    if (length(colores)>1){
        CM.statMatrix_highlighted<-CM.statMatrix[CM.statMatrix$highlight!="Other",]
        CM.statMatrix_other<-CM.statMatrix[CM.statMatrix$highlight=="Other",]
    
        p<-plotly::plot_ly(CM.statMatrix_other, x=~log.adj.pVal,y=~OR,type="scatter",mode="markers",
                           text=paste0(CM.statMatrix_other$Accession,": ",CM.statMatrix_other$TF,
                                       '<br>Treatment: ',CM.statMatrix_other$Treatment,
                                       '<br>Cell: ',CM.statMatrix_other$Cell),
                           color = ~highlight, colors=colores)
        p<-plotly::add_markers(p,x=CM.statMatrix_highlighted$log.adj.pVal, y=CM.statMatrix_highlighted$OR,type="scatter", mode="markers",
                               text=paste0(CM.statMatrix_highlighted$Accession,": ",CM.statMatrix_highlighted$TF,
                                           '<br>Treatment: ',CM.statMatrix_highlighted$Treatment,
                                           '<br>Cell: ',CM.statMatrix_highlighted$Cell),
                               color = CM.statMatrix_highlighted$highlight, colors=colores)%>%
            plotly::layout(title=plot_title)
    }else if (length(colores)==1){
        p<-plotly::plot_ly(CM.statMatrix, x=~log.adj.pVal,y=~OR,type="scatter",mode="markers",
                           text=paste0(CM.statMatrix$Accession,": ",CM.statMatrix$TF,
                                       '<br>Treatment: ',CM.statMatrix$Treatment,
                                       '<br>Cell: ',CM.statMatrix$Cell),
                           color = ~highlight, colors=colores)
    }
    p
    return(p)
}

plot_GSEA_ES<-function(GSEA_result,LFC,plot_title = NULL,specialTF = NULL,TF_colors = NULL){

    #' @title Plots Enrichment Score from the output of GSEA.run.
    #' @description Function to plot the Enrichment Score of every member of the ChIPseq binding database.
    #' @param GSEA_result Returned by GSEA.run
    #' @param LFC Vector with log2(Fold Change) of every gene that has an Entrez ID. Arranged from higher to lower.
    #' @param plot_title (Optional) Title for the plot
    #' @param specialTF (Optional) Named vector of TF symbols -as written in the enrichment table- to be highlighted in the plot.
    #' The name of each element specifies its color group, i.e.: naming elements HIF1A and HIF1B as "HIF" to represent them with the same color.
    #' @param TF_colors (Optional) Colors to highlight TFs chosen in specialTF.
    #' @return Plotly object with a scatter plot -Enrichment scores- and a heatmap -log2(fold change) bar-.
    #' @export plot_GSEA_ES
    #' @examples
    #' plot_GSEA_ES(GSEA_result,LFC,"Transcription Factor Enrichment",specialTF,colors)
    #' plot_GSEA_ES(GSEA_result=GSEA_result,LFC=LFC)

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    requireNamespace("dplyr")
    requireNamespace("plotly")

   if(is.list(GSEA_result)==T){
        tabla.Enr<-GSEA_result$Enrichment
    }else if(is.data.frame(GSEA_result)==T){
        tabla.Enr<-GSEA_result
    }

    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}

    if (is.null(specialTF)){
        highlight_list<-rep("Other",length(tabla.Enr[,1]))
        tabla.Enr$highlight<-highlight_list
        colores<-c("azure4")
        names(colores)<-c("Other")
    }

    if (!is.null(specialTF) & is.null(TF_colors)){
        TF_colors<-c("red","blue","green","hotpink","cyan","greenyellow","gold",
                     "darkorchid","chocolate1","black","lightpink","seagreen")
        TF_colors<-TF_colors[1:length(unique(names(specialTF)))]
        highlight_list<-highlight.TF(tabla.Enr,2,specialTF,TF_colors)
        tabla.Enr$highlight<-highlight_list[[1]]
        colores<-highlight_list[[2]]
    }
    if(!is.null(specialTF)&!is.null(TF_colors)){
        highlight_list<-highlight.TF(tabla.Enr,2,specialTF,TF_colors)
        tabla.Enr$highlight<-highlight_list[[1]]
        colores<-highlight_list[[2]]
    }

    simbolo<-rep("o",times=length(tabla.Enr[,1]))
    for (i in 1:length(tabla.Enr[,1])){
        if(!is.na(tabla.Enr[i,4])){
            if (tabla.Enr[i,4]<=0.05){simbolo[i]<-"x"}
        }
    }
    tabla.Enr$symbol<-simbolo

    MetaData<-MetaData[MetaData$Accession%in%tabla.Enr$Accession,]
    MetaData<-dplyr::arrange(MetaData,Accession)
    tabla.Enr<-dplyr::arrange(tabla.Enr,Accession)
    tabla.Enr$Treatment<-MetaData$Treatment
    rm(MetaData)

        if(length(colores>1)){
        tabla.Enr_highlighted<-tabla.Enr[tabla.Enr$highlight!="Other",]
        tabla.Enr_other<-rbind(tabla.Enr[tabla.Enr$highlight=="Other",],tabla.Enr_highlighted)
    
        p<-plotly::plot_ly(tabla.Enr_other, x=tabla.Enr_other$Arg.ES,y=tabla.Enr_other$ES, type="scatter", mode="markers",
                   text=paste0(tabla.Enr_other$Accession,": ",tabla.Enr_other$TF,
                              '<br>Pval: ',round(tabla.Enr_other$pval.ES,3),
                              '<br>Treatment: ',tabla.Enr_other$Treatment),
                   color=tabla.Enr_other$highlight, colors=colores, symbol=tabla.Enr_other$symbol, symbols=c("circle","x"))
    
        p<-plotly::add_markers(p,x=tabla.Enr_highlighted$Arg.ES, y=tabla.Enr_highlighted$ES,type="scatter", mode="markers",
                               text=paste0(tabla.Enr_highlighted$Accession,": ",tabla.Enr_highlighted$TF,
                                           '<br>Pval: ',round(tabla.Enr_highlighted$pval.ES,3),
                                           '<br>Treatment: ',tabla.Enr_highlighted$Treatment),
                               color=tabla.Enr_highlighted$highlight, colors=colores,
                               symbol=tabla.Enr_highlighted$symbol, symbols=c("circle","x"))%>%
            plotly::layout(title=plot_title,
                   xaxis = list(title = "Argument"),
                   yaxis = list (title = "ES"))
    }else if (length(colores)==1){
        p<-plotly::plot_ly(tabla.Enr, x=tabla.Enr$Arg.ES,y=tabla.Enr$ES, type="scatter", mode="markers",
                           text=paste0(tabla.Enr$Accession,": ",tabla.Enr$TF,
                                       '<br>Pval: ',round(tabla.Enr$pval.ES,3),
                                       '<br>Treatment: ',tabla.Enr$Treatment),
                           color=tabla.Enr$highlight, colors=colores, symbol=tabla.Enr$symbol, symbols=c("circle","x"))
    }

    LFC.bar<-get_LFC_bar(LFC)

    graf<-plotly::subplot(p,LFC.bar,shareX = TRUE,nrows = 2,heights = c(0.95, 0.05),titleY = TRUE)
    graf
    return(graf)
}

plot_GSEA_RES<-function(GSEA_result,LFC,plot_title = NULL,line.colors = NULL,line.styles = NULL){

    #' @title Plots all the RES stored in a GSEA.run output.
    #' @description Function to plot all the RES stored in a GSEA.run output.
    #' @param GSEA_result Returned by GSEA.run
    #' @param LFC Vector with log2(Fold Change) of every gene that has an Entrez ID. Arranged from higher to lower.
    #' @param plot_title (Optional) Title for the plot.
    #' @param line.colors (Optional) Vector of colors for each line.
    #' @param line.styles (Optional) Vector of line styles for each line ("solid"/"dash"/"longdash").
    #' @return Plotly object with a line plot -running sums- and a heatmap -log2(fold change) bar-.
    #' @export plot_GSEA_RES
    #' @examples
    #' plot_GSEA_RES(GSEA_result,LFC,"Transcription Factor Enrichment",colors.RES,lines.RES)
    #' plot_GSEA_RES(GSEA_result=GSEA_result,LFC=LFC)

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    requireNamespace("utils")
    requireNamespace("dplyr")
    requireNamespace("plotly")

    if(is.null(line.colors)){
        line.colors<-c("red","blue","green","hotpink","cyan","greenyellow","gold",
                       "darkorchid","chocolate1","black","lightpink","seagreen")
        line.colors<-line.colors[1:length(names(GSEA_result$RES))]
    }
    if(is.null(line.styles)){line.styles<-rep("solid",length(names(GSEA_result$RES)))}
    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}

    accessions<-names(GSEA_result$RES)
    GSEA.runningSum<-GSEA_result$RES

    chip_index<-get_chip_index("g")
    tf<-chip_index[chip_index[,1]%in%accessions,]
    rm(chip_index)

    MetaData<-MetaData[MetaData$Accession%in%accessions,]
    Accessions<-vector()
    Cell<-vector()
    Treatment<-rep("no treatment",length(accessions))
    TF<-vector()
    RES<-list()
    for(i in 1:length(accessions)){
        Accessions[i]<-accessions[i]
        Cell[i]<-MetaData[MetaData$Accession==accessions[i],3]
        if(MetaData[MetaData$Accession==accessions[i],5]!="none"){Treatment[i]<-MetaData[MetaData$Accession==accessions[i],5]}
        TF[i]<-tf[tf[,1]==accessions[i],2]
        RES<-c(RES,GSEA.runningSum[names(GSEA.runningSum)%in%accessions[i]][1])
    }
    tabla<-data.frame(Accessions,Cell,Treatment,TF,stringsAsFactors = F)
    tabla$RES<-RES

    rm(Accessions,Cell,Treatment,TF,RES)

    if(length(tabla[,1])>1){
        for(i in 1:length(accessions)){
            if (i==1){
                grafica<-plotly::plot_ly(tabla,x=c(1:length(tabla$RES[[1]])),
                                 y=tabla$RES[[accessions[1]]],
                                 type="scatter", mode="lines", line=list(color=line.colors[1],dash=line.styles[1]),
                                 name=paste0(tabla$Accessions[1]," - ",tabla$TF[1]),
                                 text=paste0(tabla$Accessions[1]," - ",tabla$TF[1],'<br>Cell: ',tabla$Cell[1],", ",tabla$Treatment[1]))
            }else if (i>1 & i<length(accessions)){
                grafica<-plotly::add_trace(p = grafica,y=tabla$RES[[accessions[i]]],
                                   type="scatter", mode="lines",line=list(color=line.colors[i],dash=line.styles[i]),
                                   name=paste0(tabla$Accessions[i]," - ",tabla$TF[i]),
                                   text=paste0(tabla$Accessions[i]," - ",tabla$TF[i],'<br>Cell: ',tabla$Cell[i],", ",tabla$Treatment[i]))
            }else if (i==length(accessions)){
                grafica<-plotly::add_trace(p = grafica,y=tabla$RES[[accessions[i]]],
                                   type="scatter", mode="lines",line=list(color=line.colors[i],dash=line.styles[i]),
                                   name=paste0(tabla$Accessions[i]," - ",tabla$TF[i]),
                                   text=paste0(tabla$Accessions[i]," - ",tabla$TF[i],'<br>Cell: ',tabla$Cell[i],", ",tabla$Treatment[i]))%>%
                    plotly::layout(title=plot_title,
                           xaxis = list(title = "Argument"),
                           yaxis = list (title = "ES"))
            }
        }
    }else{
        grafica<-plotly::plot_ly(tabla,x=c(1:length(tabla$RES[[1]])),
                         y=tabla$RES[[accessions[1]]],
                         type="scatter", mode="lines", line=list(color=line.colors[1],dash=line.styles[1]),
                         name=paste0(tabla$Accessions[1]," - ",tabla$TF[1]))%>%
            plotly::layout(title=plot_title,
                   xaxis = list(title = "Argument"),
                   yaxis = list (title = "ES"))
    }

    LFC.bar<-get_LFC_bar(LFC)

    graf<-plotly::subplot(grafica,LFC.bar,shareX = TRUE,nrows = 2,heights = c(0.95, 0.05),titleY = TRUE)
    graf
    return(graf)
}

get_LFC_bar<-function(LFC){

    #' @title Plots a color bar from log2(Fold Change) values.
    #' @description Function to plot a color bar from log2(Fold Change) values from an expression experiment.
    #' @param LFC Vector of log2(fold change) values arranged from higher to lower. Use ony the values of genes that have an Entrez ID.
    #' @return Plotly heatmap plot -log2(fold change) bar-.
    # examples
    # get_LFC_bar(arranged.log2FC.array)

    if(!requireNamespace("scales", quietly = TRUE)){
        stop("scales package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    requireNamespace("grDevices")
    requireNamespace("dplyr")
    requireNamespace("plotly")
    requireNamespace("scales")

    vals <- scales::rescale(LFC)
    o <- order(vals, decreasing = F)
    cols1 <- scales::col_numeric(grDevices::colorRamp(c("mistyrose","red3")), domain = NULL)(vals[1:length(LFC[LFC>0])])
    cols2 <- scales::col_numeric(grDevices::colorRamp(c("navy","lightcyan")), domain = NULL)(vals[length(LFC[LFC>0])+1:length(LFC)])
    cols<-c(cols1,cols2)

    colz <-data.frame(vals[o], cols[o])

    LFC.bar<-plotly::plot_ly(x=c(1:length(LFC)), y=rep(1,length(LFC)),
                     z = LFC, type = "heatmap",colorscale=colz,showscale = FALSE)%>%
        plotly::layout(yaxis=list(visible=F))

    return(LFC.bar)
}

plot_RES<-function(GSEA.runningSum,LFC,plot_title = NULL,line.colors = NULL,line.styles = NULL){

    #' @title Plots selected RES from the output of the function GSEA.run.
    #' @description Function to plot selected RES from the output of the function GSEA.run.
    #' @param GSEA.runningSum List of Running Enrichment Score vectors.
    #' @param plot_title Title for the plot.
    #' @param LFC Vector with log2(Fold Change) of every gene that has an Entrez ID. Arranged from higher to lower
    #' @param line.colors (Optional) Vector of colors for each line.
    #' @param line.styles (Optional) Vector of line styles for each line (solid/dash/longdash).
    #' @return Plotly object with a line plot -running sums- and a heatmap -log2(fold change) bar-.
    #' @export plot_RES
    #' @examples
    #' plot_RES(EPAS1.runningSums,LFC,"Transcription Factor Enrichment",colors,lines)
    #' plot_RES(GSEA.runningSum=EPAS1.runningSums,LFC=LFC)

    if(!requireNamespace("plotly", quietly = TRUE)){
        stop("plotly package needed for this function to work. Please install it.",
             call. = FALSE)
    }
    requireNamespace("utils")
    requireNamespace("dplyr")
    requireNamespace("plotly")

    if(is.null(line.colors)){
        line.colors<-c("red","blue","green","hotpink","cyan","greenyellow","gold",
                       "darkorchid","chocolate1","black","lightpink","seagreen")
        line.colors<-line.colors[1:length(names(GSEA.runningSum))]
    }
    if(is.null(line.styles)){line.styles<-rep("solid",length(names(GSEA.runningSum)))}
    if (is.null(plot_title)){plot_title<-"Transcription Factor Enrichment"}

    accessions<-names(GSEA.runningSum)

    chip_index<-get_chip_index("g")
    tf<-chip_index[chip_index[,1]%in%accessions,]
    rm(chip_index)

    MetaData<-MetaData[MetaData$Accession%in%accessions,]

    Accessions<-vector()
    Cell<-vector()
    Treatment<-rep("no treatment",length(accessions))
    TF<-vector()
    RES<-list()
    for(i in 1:length(accessions)){
        Accessions[i]<-accessions[i]
        Cell[i]<-MetaData[MetaData$Accession==accessions[i],3]
        if(MetaData[MetaData$Accession==accessions[i],5]!="none"){Treatment[i]<-MetaData[MetaData$Accession==accessions[i],5]}
        TF[i]<-tf[tf[,1]==accessions[i],2]
        RES<-c(RES,GSEA.runningSum[names(GSEA.runningSum)%in%accessions[i]][1])
    }
    tabla<-data.frame(Accessions,Cell,Treatment,TF,stringsAsFactors = F)
    tabla$RES<-RES

    rm(Accessions,Cell,Treatment,TF,RES)

    if(length(tabla[,1])>1){
        for(i in 1:length(accessions)){
            if (i==1){
                grafica<-plotly::plot_ly(tabla,x=c(1:length(tabla$RES[[1]])),
                                 y=tabla$RES[[accessions[1]]],
                                 type="scatter", mode="lines", line=list(color=line.colors[1],dash=line.styles[1]),
                                 name=paste0(tabla$Accessions[1]," - ",tabla$TF[1]),
                                 text=paste0(tabla$Accessions[1]," - ",tabla$TF[1],'<br>Cell: ',tabla$Cell[1],", ",tabla$Treatment[1]))
            }else if (i>1 & i<length(accessions)){
                grafica<-plotly::add_trace(p = grafica,y=tabla$RES[[accessions[i]]],
                                   type="scatter", mode="lines",line=list(color=line.colors[i],dash=line.styles[i]),
                                   name=paste0(tabla$Accessions[i]," - ",tabla$TF[i]),
                                   text=paste0(tabla$Accessions[i]," - ",tabla$TF[i],'<br>Cell: ',tabla$Cell[i],", ",tabla$Treatment[i]))
            }else if (i==length(accessions)){
                grafica<-plotly::add_trace(p = grafica,y=tabla$RES[[accessions[i]]],
                                   type="scatter", mode="lines",line=list(color=line.colors[i],dash=line.styles[i]),
                                   name=paste0(tabla$Accessions[i]," - ",tabla$TF[i]),
                                   text=paste0(tabla$Accessions[i]," - ",tabla$TF[i],'<br>Cell: ',tabla$Cell[i],", ",tabla$Treatment[i]))%>%
                    plotly::layout(title=plot_title,
                           xaxis = list(title = "Argument"),
                           yaxis = list (title = "ES"))
            }
        }
    }else{
        grafica<-plotly::plot_ly(tabla,x=c(1:length(tabla$RES[[1]])),
                         y=tabla$RES[[accessions[1]]],
                         type="scatter", mode="lines", line=list(color=line.colors[1],dash=line.styles[1]),
                         name=paste0(tabla$Accessions[1]," - ",tabla$TF[1]))%>%
            plotly::layout(title=plot_title,
                   xaxis = list(title = "Argument"),
                   yaxis = list (title = "ES"))
    }

    LFC.bar<-get_LFC_bar(LFC)

    graf<-plotly::subplot(grafica,LFC.bar,shareX = TRUE,nrows = 2,heights = c(0.95, 0.05),titleY = TRUE)
    graf
    return(graf)
}

