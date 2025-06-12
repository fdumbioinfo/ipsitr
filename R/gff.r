#' @title Ensembl parser for gene, transcript, exon and prot.
#' @description ggf file parser
#' @param gff character path for gff ensembl file
#' @param pepfa character path for protein ensembl file
#' @return data.frame
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach foreach
#' @importFrom stringr str_to_lower
#' @importFrom utils head
#' @importFrom rlang .data
#' @export
gff <- function( gff = NULL , pepfa = NULL )
{
  i=j=1
  #
  # gff gene, transcript, exon
  #
  gff %>% moal:::nbline(compression = "gz") -> NbLine
  gff %>% gzfile -> f
  f %>% readLines(NbLine) -> rl0
  close(f)
  gff %>% basename %>% sub("^(.).*_(.).*\\..*.gz$","\\1\\2",.) %>% stringr::str_to_lower(.) -> Species
  gff %>% basename %>% gsub("\\.","_",.) %>% gsub("_gff3_gz","",.) -> DirName0
  DirName0 %>% dir.create
  # log genome info
  Log0 <- list()
  rl0 %>% grep("^#!genome-",.,value = T) %>% paste0(collapse="\n") -> Log0[[1]]
  #
  # table 9 : annotations
  #
  # rl0 %>% strsplit("\t") %>% lapply(length) %>% unlist %>% table
  rl0 %>% strsplit("\t") %>% lapply(length) %>% unlist %>% "=="(9) -> sel
  rl0[sel] -> rl1
  rl1 %>% strsplit("\t") %>% lapply("[",3) %>% unlist -> FeatureType0
  c( FeatureType0 %>% table %>% names,FeatureType0 %>% table %>% as.numeric ) %>%
    matrix(ncol=2,byrow = F) %>% apply(1,paste0,collapse=": ") %>% paste0(collapse = "\n") -> Log0[[2]]
  #
  # gene
  #
  rl1 %>% grep("ID=gene:",.) -> sel
  rl1[sel] -> rl2
  rl2 %>% strsplit("\t") %>% lapply("[",3) %>% unlist -> FeatureType0
  c( FeatureType0 %>% table %>% names,FeatureType0 %>% table %>% as.numeric ) %>%
    matrix(ncol=2,byrow = F) %>% apply(1,paste0,collapse=": ") %>% paste0(collapse = "\n") -> Log0[[3]]
  rl2 %>% strsplit("\t") %>% unlist %>% matrix(ncol=9,byrow=T) %>% data.frame -> rl3
  rl3$X9 %>% strsplit(";") %>% lapply(length) %>% unlist %>% table -> t
  a0 <- foreach(i=1:length(t)) %do%
    {
      rl3$X9 %>% strsplit(";") %>% lapply(length) %>% unlist %>% "=="(t %>% names %>% "["(i) ) %>% which -> sel
      rl3[sel,] -> rl4
      rl4$X9  %>% strsplit(";") %>% unlist %>% sub("^(.*)=.*","\\1",.) %>% sub("^(transcript:).*","\\1",.) %>% sub("^(gene:).*","\\1",.) %>%
        matrix(ncol=t %>% names %>% "["(i) %>% as.numeric ,byrow = T) %>% data.frame %>% apply(MARGIN = 1,paste0,collapse="") -> HeaderAll0
      HeaderAll0 %>% table -> HeaderAll1
      foreach( j=1:length(HeaderAll1) ) %do%
        {
          HeaderAll0 %>% "=="(HeaderAll1[j] %>% names ) -> sel
          rl4[sel,] -> rl5
          rl5[,-c(9,6,8)] %>% setNames( c("chr","dbtype","genetype","start","end","strand") ) %>% data.frame -> t0
          rl5$X9[1] %>% strsplit(";") %>% unlist %>% sub("(.*)=.*","\\1",.) %>% sub("(.*):.*","\\1",.) -> Header0
          rl5$X9 %>% strsplit(";") %>% unlist %>% sub("^.*=(.*)","\\1",.) %>% sub("^transcript:(.*)","\\1",.) %>% sub("^gene:(.*)","\\1",.) %>% 
            matrix(ncol=t %>% names %>% "["(i) %>% as.numeric ,byrow = T) %>% data.frame %>% setNames(Header0) -> tt0
          data.frame(tt0,t0)
        }
    }
  a0 %>% unlist(recursive = F) -> a1
  a2 <- foreach(i=1:(length(a1)-1) ) %do% { left_join(a1[[i]],a1[[length(a1)]]) %>% data.frame(check.names = F) }
  a3 <- foreach(i=1:length(a2) ) %do% { a1[[length(a1)]] %>% colnames %>% match(a2[[i]] %>% colnames) %>% a2[[i]][,.] %>% data.frame(check.names = F) }
  a4 <- foreach(i=1:length(a3),.combine = "rbind" ) %do% { a3[[i]] }
  a4 %>% rbind(a1[[length(a1)]],.) -> a5
  colnames(a5)[1] <- "ENSGID"
  colnames(a5)[2] <- "Symbol"
  a5[,c(1:4,8,11:13)] -> gene0
  #
  # transcript
  #
  rl1 %>% grep("ID=transcript",.) -> sel
  rl1 %>% "["(sel) -> rl2
  rl2 %>% strsplit("\t") %>% lapply("[",3) %>% unlist -> FeatureType0
  c( FeatureType0 %>% table %>% names,FeatureType0 %>% table %>% as.numeric ) %>%
    matrix(ncol=2,byrow = F) %>% apply(1,paste0,collapse=": ") %>% paste0(collapse = "\n") -> Log0[[4]]
  rl2 %>% strsplit("\t") %>% unlist %>% matrix(ncol=9,byrow = T) %>% data.frame -> rl3
  rl3$X9 %>% strsplit(";") %>% lapply(length) %>% unlist %>% table -> t
  a0 <- foreach(i=1:length(t)) %do%
    {
      rl3$X9 %>% strsplit(";") %>% lapply(length) %>% unlist %>% "=="(t %>% names %>% "["(i) ) %>% which -> sel
      rl3[sel,] -> rl4
      rl4$X9  %>% strsplit(";") %>% unlist %>% sub("^(.*)=.*","\\1",.) %>% sub("^(transcript:).*","\\1",.) %>% sub("^(gene:).*","\\1",.) %>%
        matrix(ncol=t %>% names %>% "["(i) %>% as.numeric ,byrow = T) %>% data.frame %>% apply(MARGIN = 1,paste0,collapse="") -> HeaderAll0
      HeaderAll0 %>% table -> HeaderAll1
      foreach( j=1:length(HeaderAll1) ) %do%
        {
          HeaderAll0 %>% "=="(HeaderAll1[j] %>% names ) -> sel
          rl4[sel,] -> rl5
          rl5[,-c(9,6,8)] %>% setNames( c("chr","dbtype","genetype","start","end","strand") ) %>% data.frame -> t0
          rl5$X9[1] %>% strsplit(";") %>% unlist %>% sub("(.*)=.*","\\1",.) %>% sub("(.*):.*","\\1",.) -> Header0
          rl5$X9 %>% strsplit(";") %>% unlist %>% sub("^.*=(.*)","\\1",.) %>% sub("^transcript:(.*)","\\1",.) %>% sub("^gene:(.*)","\\1",.) %>% 
            matrix(ncol=t %>% names %>% "["(i) %>% as.numeric ,byrow = T) %>% data.frame %>% setNames(Header0) -> tt0
          data.frame(tt0,t0)
        }
    }
  a0 %>% unlist(recursive = F) -> a1
  a2 <- foreach(i=1:(length(a1)-1) ) %do% { left_join(a1[[i]],a1[[length(a1)]]) %>% data.frame(check.names = F) }
  a3 <- foreach(i=1:length(a2) ) %do% { a1[[length(a1)]] %>% colnames %>% match(a2[[i]] %>% colnames) %>% a2[[i]][,.] %>% data.frame(check.names = F) }
  a4 <- foreach(i=1:length(a3),.combine = "rbind" ) %do% { a3[[i]] }
  a4 %>% rbind(a1[[length(a1)]],.) -> a5
  colnames(a5)[1] <- "ENSTID"
  colnames(a5)[2] <- "ENSGID"
  colnames(a5)[3] <- "TranscriptSymbol"
  a5 %>% colnames
  if( a5 %>% colnames %>% grepl("ccds",.) %>% any ){ a5[,c(1:4,10,13:15)] -> transcript0 }else{ a5[,c(1:4,8,11:13)] -> transcript0 }
  #
  # exon
  #
  rl1 %>% strsplit("\t") %>% lapply("[",3) %>% "=="("exon") %>% which -> sel
  rl1 %>% "["(sel) -> rl2
  rl2 %>% strsplit("\t") %>% lapply("[",3) %>% unlist -> FeatureType0
  c( FeatureType0 %>% table %>% names,FeatureType0 %>% table %>% as.numeric ) %>%
    matrix(ncol=2,byrow = F) %>% apply(1,paste0,collapse=": ") %>% paste0(collapse = "\n") -> Log0[[5]]
  rl2 %>% strsplit("\t") %>% unlist %>% matrix(ncol=9,byrow = T) %>% data.frame -> rl3
  rl3$X9 %>% strsplit(";") %>% lapply(length) %>% unlist %>% table -> t
  a0 <- foreach(i=1:length(t)) %do%
    {
      rl3$X9 %>% strsplit(";") %>% lapply(length) %>% unlist %>% "=="(t %>% names %>% "["(i) ) %>% which -> sel
      rl3[sel,] -> rl4
      rl4$X9  %>% strsplit(";") %>% unlist %>% sub("^(.*)=.*","\\1",.) %>% sub("^(transcript:).*","\\1",.) %>% sub("^(gene:).*","\\1",.) %>%
        matrix(ncol=t %>% names %>% "["(i) %>% as.numeric ,byrow = T) %>% data.frame %>% apply(MARGIN = 1,paste0,collapse="") -> HeaderAll0
      HeaderAll0 %>% table -> HeaderAll1
      foreach( j=1:length(HeaderAll1) ) %do%
        {
          HeaderAll0 %>% "=="(HeaderAll1[j] %>% names ) -> sel
          rl4[sel,] -> rl5
          rl5[,-c(9,6,8)] %>% setNames( c("chr","dbtype","genetype","start","end","strand") ) %>% data.frame -> t0
          rl5$X9[1] %>% strsplit(";") %>% unlist %>% sub("(.*)=.*","\\1",.) %>% sub("(.*):.*","\\1",.) -> Header0
          rl5$X9 %>% strsplit(";") %>% unlist %>% sub("^.*=(.*)","\\1",.) %>% sub("^transcript:(.*)","\\1",.) %>% sub("^gene:(.*)","\\1",.) %>% 
            matrix(ncol=t %>% names %>% "["(i) %>% as.numeric ,byrow = T) %>% data.frame %>% setNames(Header0) -> tt0
          data.frame(tt0,t0)
        }
    }
  a0 %>% unlist(recursive = F) %>% data.frame -> a1
  colnames(a1)[1] <- "ENSTID"
  colnames(a1)[2] <- "ENSEID"
  a1 %>% dplyr::select(c(2,1,12:14)) -> exon0
  #
  # prot
  #
  pepfa %>% moal:::nbline(compression = "gz") -> NbLine
  pepfa %>% gzfile -> f
  f %>% readLines(NbLine) -> rl0
  close(f)
  rl0 %>% grep("^>",.) -> sel
  rl0[sel] -> rl1
  # rl1 %>% strsplit(" ") %>% lapply(length) %>% unlist %>% table
  # no description
  rl1 %>% grep("description",.,invert = T) -> sel
  rl1[sel] -> rl2
  rl2 %>% strsplit(" ") %>% lapply(length) %>% unlist %>% table -> t
  a0 <- foreach( i=1:length(t) ) %do%
    {
      rl2 %>% strsplit(" ") %>% lapply(length) %>% unlist %>% "=="(t[i] %>% names) %>% which -> sel
      rl2[sel] -> rl3
      rl3 %>% strsplit(" ") %>% lapply("[",1) %>% unlist %>% gsub(">","",.) %>% sub("(.*)\\..$","\\1",.) -> ENSPID
      rl3 %>% strsplit(" ") %>% lapply("[",3) %>% unlist %>% sub("^scaffold:(.*)","\\1",.) -> chr
      rl3[1] %>% strsplit(" ") %>% lapply("[",-c(1,2,3)) %>% 
        unlist %>% sub("^(.*):.*","\\1",.) %>% sub("^(scaffold):.*","\\1",.) -> Header0
      rl3 %>% strsplit(" ") %>% lapply("[",-c(1,2,3)) %>% unlist %>% sub("^scaffold:(.*)","\\1",.) %>% 
        sub("^.*:(.*)","\\1",.) %>% sub("(.*)\\..*","\\1",.) %>% matrix(ncol= t[i] %>% names %>% as.numeric %>% "-"(3) ,byrow = T) %>% data.frame %>% setNames(Header0) -> rl3
      data.frame(ENSPID,chr,rl3)
    }
  if( length(a0)>1 )
  {
    a1 <- foreach(i=1:(length(a0)-1) ) %do% { left_join(a0[[i]],a0[[length(a0)]]) %>% data.frame(check.names = F) }
    a2 <- foreach(i=1:length(a1) ) %do% { a0[[length(a0)]] %>% colnames %>% match(a1[[i]] %>% colnames) %>% a1[[i]][,.] %>% data.frame(check.names = F) }
    a3 <- foreach(i=1:length(a2),.combine = "rbind" ) %do% { a2[[i]] }
    a3 %>% rbind(a0[[length(a0)]],.) -> d0
  }else
  {
    a0[[1]] -> d0
  }
  # description
  rl1 %>% grep("description",.) -> sel
  rl1[sel] -> rl2
  rl2 %>% sub("^(.*)description:.*","\\1",.) %>% strsplit(" ") %>% 
    lapply("[",-c(1,2,3)) %>% lapply(length) %>% unlist %>% table -> t
  a0 <- foreach( i=1:length(t) ) %do%
    {
      rl2 %>% sub("^(.*)description:.*","\\1",.) %>% strsplit(" ") %>% 
        lapply("[",-c(1,2,3)) %>% lapply(length) %>% unlist %>% "=="(t[i] %>% names) %>% which -> sel
      rl2[sel] -> rl3
      rl3 %>% strsplit(" ") %>% lapply("[",1) %>% unlist %>% gsub(">","",.) %>% sub("(.*)\\..*$","\\1",.) -> ENSPID
      rl3 %>% strsplit(" ") %>% lapply("[",3) %>% unlist -> chr
      rl3[1] %>% sub("^(.*)description:.*","\\1",.) %>% strsplit(" ") %>% lapply("[",-c(1,2,3)) %>% unlist %>% 
        sub("^(.*):.*","\\1",.) -> Header0
      rl3 %>% sub("^(.*)description:.*","\\1",.) %>% strsplit(" ") %>% lapply("[",-c(1,2,3)) %>% unlist %>%
        sub("^.*:(.*)","\\1",.) %>% sub("(.*)\\..*$","\\1",.) %>% matrix(ncol= t[i] %>% names %>% as.numeric ,byrow = T) %>%
        data.frame %>% setNames(Header0) -> rl4
      data.frame(ENSPID,chr,rl4)
    }
  a1 <- foreach(i=1:(length(a0)-1) ) %do% { left_join(a0[[i]],a0[[length(a0)]]) %>% data.frame(check.names = F) }
  a2 <- foreach(i=1:length(a1) ) %do% { a0[[length(a0)]] %>% colnames %>% match(a1[[i]] %>% colnames) %>% a1[[i]][,.] %>% data.frame(check.names = F) }
  a3 <- foreach(i=1:length(a2),.combine = "rbind" ) %do% { a2[[i]] }
  a3 %>% rbind(a0[[length(a0)]],.) -> a4
  a4 %>% head
  d0 %>% head
  if( d0 %>% colnames %>% grepl("^gene_symbol$",.) %>% any ){ data.frame(d0 ) %>% rbind(a0[[2]],a1[[1]],. ) -> prot0 }else
  { data.frame(d0, gene_symbol = rep(NA,nrow(d0))  ) %>% rbind(a0[[2]],a1[[1]],. ) -> prot0 }
  prot0[,c(1,3,4,6,7,2)] %>% setNames(c("ENSPID","ENSGID","ENSTID","biotype","Symbol","chr"))-> prot1
  #
  # output
  #
  # gene
  transcript0[,c(1,2)] %>% group_by(.data$ENSGID) %>% 
    mutate( ENSTIDs = paste0(.data$ENSTID %>% unique ,collapse = "|") ) %>%
    dplyr::slice(1) %>% data.frame -> t0
  exon0[,c(1,2)] %>% group_by(.data$ENSTID) %>% 
    mutate( ENSEIDs = paste0(.data$ENSEID %>% unique ,collapse = "|") ) %>%
    dplyr::slice(1) %>% dplyr::select(-1) %>% data.frame -> tt0
  t0 %>% inner_join(tt0) %>% dplyr::select(-1) -> t1
  prot1[,c(1,2)] %>% group_by(.data$ENSGID) %>% 
    mutate( ENSPIDs = paste0(.data$ENSPID %>% unique ,collapse = "|") ) %>%
    dplyr::slice(1) %>% dplyr::select(-1) %>% data.frame -> ttt0
  t1 %>% dplyr::left_join(ttt0, by = "ENSGID") -> t2
  t2 %>% inner_join(gene0,.) -> gene1
  paste("annotations_gtf_",Species,"_gene_",gene1 %>% nrow,".tsv",sep="") -> FileName0
  DirName0 %>% file.path(FileName0) -> FileName1
  gene1 %>% output(FileName1)
  # transcripts
  transcript0 %>% left_join(t1) -> transcript1
  gene1 %>% colnames
  transcript1 %>% left_join(gene1[,c(1,2,4)], by = "ENSGID") %>% data.frame -> transcript2
  paste("annotations_gtf_",Species,"_transcript_",transcript2 %>% nrow,".tsv",sep="") -> FileName0
  DirName0 %>% file.path(FileName0) -> FileName1
  transcript2 %>% output(FileName1)
  # protein
  paste("annotations_gtf_",Species,"_protein_",prot1 %>% nrow,".tsv",sep="") -> FileName0
  DirName0 %>% file.path(FileName0) -> FileName1
  prot1 %>% output(FileName1)
  # exon
  exon0 %>% dplyr::left_join(transcript2[,c(1,11)], by = "ENSTID") -> exon1
  paste("annotations_gtf_",Species,"_exon_",exon1 %>% nrow,".tsv",sep="") -> FileName0
  DirName0 %>% file.path(FileName0) -> FileName1
  exon1 %>% output(FileName1)
  # log
  file.path(DirName0,"log.txt") %>% file("w") -> f
  foreach(i=1:length(Log0)) %do% { Log0[[i]] %>% paste(.,"\n",sep="") %>% writeLines(f) }
  close(f)
}
