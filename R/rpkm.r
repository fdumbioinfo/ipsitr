#' @title read per kilobase per millions
#' @description Read per millions calculation from RNASeq raw data count.
#' @param dat matrix numeric
#' @param length numeric vector of base pair length
#' @return data.frame
#' @examples
#' # not run
#' # rpkm( mat , length )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select 
#' @export

rpkm <- function( dat , length = NULL )
{
  dat %>% apply(2,sum) -> Sum0
  dat %>% "*"(10^6) -> Dat1
  Dat1 %>% "*"(10^3) -> Dat2
  Dat2 %>% "/"(Sum0) -> Dat3
  Dat3 %>% "/"(length) -> Dat4
  Dat4 %>% return()
}
