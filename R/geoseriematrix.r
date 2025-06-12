#' @title Principal Component Analysis
#' @description Principal Component Analysis
#' @param file character path
#' @return files
#' @examples
#' # not run
#' # geoseriematrix( file )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @export
geoseriematrix <- function( file )
{
  file %>% moal:::nbline(compression = "gz") -> NbLine
  file %>% gzfile -> f
  f %>% readLines(NbLine) -> rl0
  close(f)
  rl0 %>% grep("^!Series_geo_accession",.,value=T ) %>% strsplit("\t") %>% unlist %>% "["(2) %>% gsub('\\"',"",.) -> geoid
  #
  # sample data
  #
  rl0 %>% grep("^!Sample_title",.,value=T) %>% strsplit("\t") %>% unlist %>% "["(-1) %>% gsub('\\"',"",.) -> GEOTitle
  rl0 %>% grep("^!Sample_geo_accession",.,value=T) %>% strsplit("\t") %>% unlist %>% "["(-1) %>% gsub('\\"',"",.) -> GEOID
  paste(geoid,"_sampledata_",length(GEOTitle),".tsv",sep="") -> FileName0
  paste("s",1:length(GEOTitle),sep="") -> SampleID
  data.frame(SampleID,GEOID,GEOTitle) -> s0
  s0 %>% output(FileName0)
  #
  # normdata
  #
  rl0[ ( rl0 %>% grep( '"ID_REF"', . ) ):length(rl0)  ] -> normdata
  rl0 %>% grep('"ID_REF"',.,value=T) %>% strsplit("\t") %>% unlist %>% gsub('\\"',"",.) %>% "["(-1) -> header
  if( header %>% "=="(s0$GEOID) %>% all )
  {
    rl0[( rl0 %>% grep('"ID_REF"',.)+1):length(rl0) ] -> m0
    m0 %>% strsplit("\t") %>% unlist %>% gsub('\\"',"",.) %>% matrix(nrow = length(m0) , byrow = T ) %>% data.frame -> m1
    m1 %>% setNames(c("rowID",header)) -> m2
    paste(geoid,"_normdata_",length(GEOTitle),"_",nrow(m2),".tsv",sep="") -> FileName0
    m2 %>% output(FileName0)
  }else{ paste("no match sample and data" %>% cat) }
}
