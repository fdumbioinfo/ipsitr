#' @title Transcripts per millions
#' @description Read per millions calculation from RNASeq raw data count.
#' @param dat matrix numeric
#' @param length numeric vector of base pair length
#' @return data.frame
#' @examples
#' # not run
#' # tpm( mat , length )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select 
#' @export

tpm <- function( dat , length = NULL )
{
  dat %>% "*"(10^3) -> Dat1
  Dat1 %>% "/"(length) -> Dat2
  Dat2 %>% apply(2,sum) -> Sum0
  Dat2 %>% "/"(Sum0) -> Dat3
  Dat3 %>% "*"(10^6) -> Dat4
  Dat4 %>% return()
}
