#' @title To extract data from spectrometry mzxml file
#' @param .files character vector of mzxml file
#' @param .path character
#' @param .dirname character
#' @return data.frame
#' @examples
#' # not run
#' # xcmsmzxml( listfiles )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach  %do%
#' @importFrom xcms xcmsSet group retcor fillPeaks peakTable
#' @export
orbitrapwf1 <- function(
  files , .path = "." , .dirname = "orbitrapwf1" )
{
  .path %>% file.path( .dirname ) -> Path
  Path %>% dir.create
  # sif
  Path %>% file.path(paste( "metadata_",length(.files),".tsv",sep="") ) -> FileName
  data.frame( "SampleID" = paste("s",1:length(.files),sep=""),"FileName" = .files %>% basename ) -> s0
  s0 %>% output( FileName )
  # extraction
  .files %>% xcmsSet -> x0
  x0 %>% group -> x1
  x1 %>% retcor( family = "symmetric" , plottype = "mdevden" ) -> x2
  x2 %>% group( bw = 10 ) -> x3
  x3 %>% fillPeaks -> x4
  x4 %>% peakTable -> x5
  # output
  data.frame("rowID" = paste("row",1:nrow(x5) , sep = ""),x5,stringsAsFactors = F) -> x6
  colnames(x6)[ 10:( length( .files ) + 9 ) ] <- paste("s",1:length(.files),sep="")
  # rawdata
  x6[ , c(1,10:(length(.files)+9) ) ] -> Dat0
  Path %>% file.path( paste("rawdata_",length( .files ),"_",nrow( x5 ),".tsv",sep="") ) -> FileName
  Dat0 %>% output( FileName )
  # annotation
  x6[ , c(1:9) ] -> annot0
  Path %>% file.path( paste("annotations_",nrow(x5),".tsv",sep="") ) -> FileName
  annot0 %>% output( FileName )
}






