#' @title install custom cdf
#' @description
#' Function to install custom cdf from brainarray site.
#' @param arraytype character came from affynorm function
#' @param version character custom cdf version to install
#' @param annotationsource character
#' @details install custom cdf from brainarray for affymetrix .cel file normalization
#' @return character package name
#' @examples
#' # not run
#' # celfiles %>% customcdfInstall
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom  xml2 read_html
#' @importFrom rvest html_nodes html_attrs
#' @importFrom utils download.file install.packages
#' @importFrom dplyr select filter
#' @importFrom utils installed.packages
#' @importFrom rlang .data
#' @export

customcdfInstall <- function( arraytype = NULL , version = NULL , annotationsource = "entrezg" )
{
  INSTALL <- TRUE
  # package name
  customcdfList %>% dplyr::filter( .data$arrayType == arraytype ) %>% dplyr::select( .data$brainarrayName ) %>% as.character -> PackageName0
  customcdfList %>% dplyr::filter( .data$arrayType == arraytype ) %>% dplyr::select( .data$organism ) %>% as.character -> Species
  paste( "pd.",PackageName0,".",Species,".",annotationsource,sep="") -> PackageName
  # package list
  installed.packages() %>% as.data.frame %>% "["(.,1) %>% unlist -> PackageList
  read_html("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp") %>%
    html_nodes("a") %>% html_attrs() %>% grep("#v" , . , value = T ) %>% '['(1) %>%
    sub("#v" , "" , . ) %>% paste( . , ".0.0" , sep = ""  ) -> Version
  if( !(grepl( PackageName , PackageList ) %>% any) )
  {
    tmpDir = tempdir()
    # url
    paste( "http://mbni.org/customcdf/" , Version, "/", annotationsource, ".download/",
           PackageName, "_" , Version , ".tar.gz", sep = "") -> url
    download.file( url , paste( tmpDir , "\\", PackageName , sep = "") )
    install.packages( paste( tmpDir , "\\", PackageName , sep = ""), repos = NULL, type = "source")
    print( paste( "customcdf used: ", PackageName , "version: " , Version ) )
    INSTALL <- FALSE
  }
  if( INSTALL & !is.null(version) )
  {
    tmpDir = tempdir()
    # url
    paste( "http://mbni.org/customcdf/" , version, "/", annotationsource, ".download/",
           PackageName, "_" , version , ".tar.gz", sep = "") -> url
    download.file( url , paste( tmpDir , "\\", PackageName , sep = "") )
    install.packages( paste( tmpDir , "\\", PackageName , sep = ""), repos = NULL, type = "source")
    print( paste( "customcdf used: ", PackageName , "version: " , version) )
    Version <- version
  }
  list( PackageName , Version ) %>% return()
}







