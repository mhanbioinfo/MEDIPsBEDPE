##########################################################################
##Function calculates density of C's, G's, and CpGs in the given sequences as well as in the reference genome
##Returns CpG enrichment values for the given sequences.
##########################################################################
##Input:	bam file or tab (|) separated txt file "chr | start | stop  | strand"
##Param:	BSgenome, pattern (character), chromosomes, shift, extend
##Output:	CpGenrich results object
##Requires:	BSgenome, Biostrings, GenomicRanges, MEDIPS.Bam2GRanges, MEDIPS.txt2Granges, MEDIPS.GenomicCoordinates
##Modified:	05/21/2020
##Author:	Joern Dietrich, Lukas Chavez

## new version by SPPearce
## https://github.com/chavez-lab/MEDIPS/issues/2

## edited by: Ming Han (2022_07_22)
## do no re-read BAM or BEDPE input if 'GRange.Reads' already exists

MEDIPS.CpGenrichNew <-
  function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){

    ## Proof correctness....
    if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}

    ## Read region file
    fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
    path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/")
    if(path==""){path=getwd()}
    if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}

    dataset = get(ls(paste("package:", BSgenome, sep = ""))[1])

    ## do no re-read BAM or BEDPE input if 'GRange.Reads' already exists
    if (exists("GRangeObj")){
      GRange.Reads = GRangeObj
    } else {
      if(!paired){
        GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)
        ## saves 'GRanges.Reads' to global env
        ## so don't have to run getGRange() again if other functions need it
        assign("GRangeObj", GRange.Reads, envir = .GlobalEnv)
      }	else{
        GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)
        ## saves 'GRanges.Reads' to global env
        ## so don't have to run getPairedGRange() again if other functions need it
        assign("GRangeObj", GRange.Reads, envir = .GlobalEnv)
      }
    }

    ## Sort chromosomes
    if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
    if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}

    ## Get chromosome lengths for all chromosomes within data set.
    cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))

    chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])

    ranges(GRange.Reads) <- restrict(ranges(GRange.Reads),+1)

    ##Calculate CpG density for regions
    total=length(chromosomes)
    cat("Calculating CpG density for given regions...\n")

    readsChars <- unlist(getSeq(dataset, GRange.Reads, as.character=TRUE))

    regions.CG = sum(vcountPattern("CG",readsChars))
    regions.C  = sum(vcountPattern("C",readsChars))
    regions.G  = sum(vcountPattern("G",readsChars))
    all.genomic= sum(width(readsChars))

    nReads <- length(readsChars)

    regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100
    regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))

    CG <- DNAStringSet("CG")
    pdict0 <- PDict(CG)
    params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
    genome.CG=sum(bsapply(params, pdict = pdict0))
    params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
    alphabet=bsapply(params)
    genome.l=sum(as.numeric(alphabet))
    genome.C=as.numeric(sum(alphabet[2,]))
    genome.G=as.numeric(sum(alphabet[3,]))
    genome.relH=genome.CG/genome.l*100
    genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);

    ##Calculate CpG density for reference genome

    enrichment.score.relH=regions.relH/genome.relH
    enrichment.score.GoGe=regions.GoGe/genome.GoGe

    gc()
    return(list(genome=BSgenome,
                regions.CG=regions.CG, regions.C=regions.C, regions.G=regions.G,
                regions.relH=regions.relH, regions.GoGe=regions.GoGe,
                genome.C=genome.C, genome.G=genome.G, genome.CG=genome.CG,
                genome.relH=genome.relH, genome.GoGe=genome.GoGe,
                enrichment.score.relH=enrichment.score.relH,
                enrichment.score.GoGe=enrichment.score.GoGe))

  }