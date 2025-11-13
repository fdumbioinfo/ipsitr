#' @title agilent workflow for rawdata extraction , QCs and normalization
#' @description QC analysis for Agilent rawdata and normalisation for agilent chip rawdata
#' @param files character list of rawdata files
#' @param sif data.frame sample meta data
#' @param bgfilter numeric threshold for background
#' @param dirname character
#' @param path character
#' @details
#' Results directory contain boxplot, acp, hierachical clustering for rawdata and dotplot of QC probes
#' A sample information file (sif) is created to complete with experimental factor variable.
#' Re run QCAgilent with competed sif for generated plot with that correspond to factor levels.
#' @return
#' Create directory QC results
#' @examples
#' # not run
#' # pathfile %>% qcagilent
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom grDevices pdf graphics.off
#' @importFrom graphics legend title
#' @importFrom stats sd
#' @export
agilentwfmir <- function( files , sif = NULL , bgfilter = 16 , dirname = NULL , path = "." )
{
  "agilentwf" -> DirName
  if(!is.null(dirname)){ paste(DirName,"_",dirname,sep="") -> DirName }
  path %>% file.path(DirName) -> Path
  Path %>% dir.create
  ### QC rawdata
  files %>% qcagilent(sif=sif,path=Path)
  ### filter bg
  file.path(Path,"QCAgilent") %>% list.files(full.names=T) -> l0
  l0 %>% basename %>% grep("^rawdata",.) %>% l0[.] %>% input -> m0
  m0 %>% dplyr::select(-1) %>% apply(1,mean) %>% ">"(.,bgfilter) %>% which -> sel
  m0[sel,] -> m1
  ### normalization
  m1 %>% dplyr::select(-1) %>% norm(.) -> m2
  m2 %>% apply(1,sd) %>% "=="(.,0) %>% which -> selsd0
  if(length(selsd0) > 0)
  { data.frame("rowID"=m1$ProbeName[-selsd0],m2[-selsd0,],stringsAsFactors=F) -> m3 }else
    { data.frame("rowID"=m1$ProbeName,m2,stringsAsFactors=F) -> m3 }
  paste("normdata_",ncol(m3)-1,"_", nrow(m3),".tsv",sep="") -> FileName
  Path %>% file.path(FileName) -> FileName
  m3 %>% output( FileName )
  ### output annotations
  l0 %>% basename %>% grep("^annotations" , . ) %>% l0[.] %>% input -> a0
  a0 %>% colnames %>% grep("GeneName" , . ) -> num
  colnames(a0)[num] <- "Symbol"
  a0[sel,] -> a1
  if( length( selsd0 ) > 0 ){ a1[-selsd0,] -> a2 }else{ a1 -> a2 }
  paste("annotations_", nrow(a2),".tsv" , sep = "" ) -> FileName
  Path %>% file.path( FileName ) -> FileName
  a2 %>% dplyr::select( .data$ProbeName , .data$Symbol , .data$SystematicName , .data$Description ) %>% output( FileName )
  ### QC normdata
  m3 %>% qc( sif = sif , dirname = "QCnormdata" , path = Path , inputdata = F)
}
