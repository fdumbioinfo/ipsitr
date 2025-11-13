#' @title extract rawdata for agilent array file
#' @details Extraction values from rawdata file
#' @param files character rawdata file path
#' @param outputdata logical
#' @param arraytype character agilent arraytype "rna" by default or "mirna"
#' @param path character output direcotry
#' @return lsit of 5 elements:
#'  - all data
#'  - rawdata matrix with control probes
#'  - rawdata matrix without control probes
#'  - matrix rawdata wihout control dimmension
#'  - annotations
#' @examples
#' # not run
#' # rawdatafilepath %>% rawdataextraction
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom moal output
#' @importFrom magrittr %>%
#' @importFrom limma read.maimages
#' @importFrom dplyr select slice group_by summarise_all
#' @importFrom rlang .data
#' @export
agilentrawdataextraction <- function( files = NULL , outputdata = TRUE , arraytype = "rna" , path = "." )
{
  #
  # rawdata extraction
  #
  list(G="gProcessedSignal") -> c0
  files %>% limma::read.maimages(source="agilent",columns=c0) -> a0
  data.frame(ProbeName=a0$genes$ProbeName,a0$E,stringsAsFactors=F) -> m0
  c("ProbeName",paste("s",1:length(files),sep="")) -> colnames(m0)
  # remove control probes
  m0 %>% dplyr::slice(which(a0$genes$ControlType == 0)) -> m1
  # summarization for duplicated probes with median
  m1 %>% dplyr::group_by(.data$ProbeName) %>% dplyr::summarise_all(mean) %>% data.frame -> m2
  if(outputdata)
  {
    paste("rawdata_",dim(m2)[2]-1,"_",dim(m2)[1],".tsv",sep="") -> FileName
    path %>% file.path(FileName) -> FileName
    m2 %>% moal::output(FileName)
    # annotations
    if(arraytype == "rna"){ a0$gene %>% dplyr::select(c(7,1,2,3,4,5,6,8,9,10)) -> annot0 }
    if(arraytype == "mir"){ a0$gene %>% dplyr::select(c(4,1,2,3,5)) -> annot0 }
    annot0 %>% dplyr::slice(which(a0$genes$ControlType == 0)) -> annot1
    annot1 %>% dplyr::group_by(.data$ProbeName) %>% dplyr::slice(1) -> annot2
    paste("annotations_",dim(annot2)[1],".tsv",sep="") -> FileName
    path %>% file.path(FileName) -> FileName
    annot2 %>% moal::output( FileName )
  }
  list(a0=a0,m0=m0,m1=m2,dimpdf=dim(m2)[1],annotations= annot0) %>% return
}




