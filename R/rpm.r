#' @title read per millions
#' @description Read per millions calculation from RNASeq raw data count.
#' @param dat data.frame raw data count with row ID
#' @return data.frame 
#' @examples
#' # not run
#' # rpm( dat )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr select 
#' @export

rpm <- function( dat )
{
  dat %>% apply(2,sum) -> Sum0
  dat %>% "*"(10^6) -> Dat1
  Dat1 %>% "/"(Sum0) -> Dat2
  Dat2 %>% return()
}
