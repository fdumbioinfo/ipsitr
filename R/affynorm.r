#' @title affynorm
#' @description normalize affymetrix rawdata cel files using rma
#' and brainarry custom cdf
#' @param celfiles character vector contening cel files path
#' @param customcdf boolean TRUE by default use the last version of custom cdf. see details.
#' @param versioncdf character if you want to use old version of custom cdf. see details
#' @param norm logical apply quantile normalization
#' @param bg logical apply rma background correction
#' @details The function by default use customcdf = TRUE meaning that the it will install the last
#' version of custom cdf. You can choose an older version using for exemple versioncdf = 23.0.0.
#' customcdf = FALSE will use affymetrix probset design.
#' custom cdf version use ncbi entrez gene id.
#' @return data.frame of normalized data with features id as first column.
#' @examples
#' # not run
#' # affynorm(celfiles)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom affyio read.celfile.header
#' @importFrom oligo read.celfiles rma
#' @importFrom Biobase exprs
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export
#' 
affynorm <- function( 
    celfiles = NULL , customcdf = TRUE , versioncdf = NULL,
    norm = TRUE , bg = TRUE)
{
  # Normalization with custom cdf
  if(customcdf)
  {
    # array type
    affyio::read.celfile.header(celfiles[1], info="full") -> ArrayInfos
    ArrayInfos[[1]] -> ArrayType
    # custom cdf install
    customcdfInstall(arraytype=ArrayType, version=versioncdf) -> PackageInfo
    PackageInfo[[1]] -> PackageName
    PackageInfo[[2]] -> PackageVersion
    paste("cdf used: ",PackageName," version: ",PackageVersion,sep="") %>% cat
    # rma
    oligo::read.celfiles(filenames=celfiles, pkgname=PackageName) %>%
      oligo::rma(normalize=norm, background=bg) %>% Biobase::exprs(.) -> m0
    data.frame(rownames(m0) %>% sub("_at","",.), m0, stringsAsFactors=F) %>%
      setNames(c("rowID",paste("s",1:dim(m0)[2],sep=""))) %>% return()
    }else
      {
        # affymetrix design
        oligo::read.celfiles(filenames=celfiles) %>% 
          oligo::rma(normalize=norm, background=bg) %>% Biobase::exprs(.) -> m0
        data.frame(rownames(m0) %>% sub("_at","",.), m0, stringsAsFactors=F) %>%
          setNames(c("rowID", paste("s",1:dim(m0)[2], sep=""))) %>% return()
      }
}













