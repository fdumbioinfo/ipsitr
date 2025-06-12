#' @title replace value by group median
#' @param dat character vector of mzxml file
#' @param backgroundvalue numeric value to substitute
#' @param method character
#' @param background logical
#' @param nareplacevalue numeric replace value for NA 0 by default
#' @return data.frame
#' @examples
#' # not run
#' # replacegroupmed(dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom foreach  %do% foreach
#' @importFrom stats median setNames
#' @importFrom edgeR calcNormFactors
#' @export
removerownoise <- function( dat , backgroundvalue = 10, method = "median", nareplacevalue = 0 )
{
  dat[,-1] -> Dat0
  # replace NA nareplacevalue
  Dat0 %>% unlist %>% is.na %>% which -> sel
  if( length(sel) > 0 ){
    Dat0 %>% unlist -> t
    t %>% replace(sel,nareplacevalue) -> tt
    tt %>% matrix( nrow = nrow(Dat0) ) -> Dat1 }else{ Dat0 -> Dat1 }
  # replace "NA" nareplacevalue
  dat %>% unlist %>% grep("^NA$",.) -> sel 
  if( length(sel) > 0 ){
    Dat1 %>% unlist -> t
    t %>% replace(sel,nareplacevalue) -> tt
    tt %>% matrix( nrow = nrow(Dat0) ) -> Dat2 }else{ Dat1 -> Dat2 }
  # remove background
  Dat2 %>% apply(1,method) -> Mean0
  Mean0 %>% ">"(backgroundvalue) %>% which -> sel
  if( length(sel) > 0 ){ 
    Dat2[sel,] -> t
    t %>% data.frame("rowID"=dat[sel,1],.) %>% setNames(dat %>% colnames) -> Dat3
    }else{ Dat2 %>% data.frame("rowID"=dat[sel,1],.) %>% setNames(dat %>% colnames) -> Dat3 }
  Dat3 %>% head
  Dat3 %>% return()
}
