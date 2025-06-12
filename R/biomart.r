#' @title Biomart wrapper.
#' @description Biomart wrapper.
#' @param species character 
#' @return data.frame
#' @examples
#' # not run
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach foreach
#' @importFrom stringr str_to_lower
#' @importFrom utils head
#' @importFrom rlang .data
#' @import moal
#' @import biomaRt
#' @export
biomart <- function( species = NULL )
{
  if( !is.null(species) )
  { 
    moal:::orthoinfo -> OrthoInfo0
    OrthoInfo0 %>% lapply("[",1) %>% "=="(species) %>% which %>% OrthoInfo0[.] %>% unlist %>%
      "["(3) %>% gsub("(.).* (.*)","\\1\\2",.) %>% stringr::str_to_lower(.) %>% paste("_gene_ensembl",sep="") -> Species0
  }else{ species = "hs" ; Species0 <- "hsapiens_gene_ensembl" }
  #
  # db
  #
  biomaRt::useMart(biomart="ensembl",dataset=Species0) -> mart
  c("ensembl_gene_id","ensembl_transcript_id","ensembl_peptide_id",
    "external_gene_name","description","chromosome_name","transcript_biotype","transcript_count","external_synonym") -> attr
  getBM( mart = mart, attributes = attr, uniqueRows = T ) -> a1
  c("ensembl_gene_id","ensembl_transcript_id","ensembl_peptide_id","entrezgene_id","uniprotswissprot","uniprotsptrembl") -> attr
  getBM( mart = mart, attributes = attr, uniqueRows = T ) -> aa1
  a1 %>% dplyr::left_join(aa1, relationship = "many-to-many") %>% unique -> a2
  c("ENSGID","ENSTID","ENSPID","Symbol","Description","Chromosome",
    "TranscriptType","TranscriptCount","Syn","NCBIGeneID","UniProtSP","UniProtTr") -> Header0
  a2 %>% setNames(Header0) -> a2
  a2 %>% head
  a2 %>% group_by(.data$ENSGID) %>%
    mutate(
      ENSTIDs = paste0(.data$ENSTID %>% unique ,collapse = "|"),
      ENSPIDs = paste0(.data$ENSPID %>% unique ,collapse = "|"),
      Synonyms = paste0(.data$Syn %>% unique ,collapse = "|")) %>% dplyr::select(.data$ENSTIDs,.data$ENSPIDs,.data$Synonyms) %>% data.frame -> Collapse0
  a2 %>% cbind(Collapse0[,-1]) -> All0
  # gene
  All0 %>% dplyr::select(-.data$ENSTID,-.data$ENSPID,-.data$Syn) %>% dplyr::group_by(.data$ENSGID) %>% dplyr::slice(1) %>% data.frame  -> gene0
  # transcript
  All0 %>% dplyr::select(c(2,1,3:ncol(All0))) %>% dplyr::group_by(.data$ENSTID) %>% dplyr::slice(1) %>% data.frame  -> transcript0
  # peptide
  All0$ENSPID %>% grep("^$",.,invert = T) %>% All0[.,] -> All1
  All1 %>% dplyr::select(c(3,1,2,4:ncol(All0))) %>% dplyr::group_by(.data$ENSPID) %>% dplyr::slice(1) %>% data.frame -> pep0
  list(gene0,transcript0,pep0) %>% return()
}
