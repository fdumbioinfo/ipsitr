#' @title grep in a data.frame
#' @description Function that grep in a data.frame for each column
#' @param pattern character to search in the data.frame
#' @param dt data.frame
#' @param value logical if FALSE return the data.frame of matching raws
#' @details Function that make grep for a data.frame.
#' @return numeric vector of matching raws or values
#' @examples
#' # not run
#' # grepdt("a", dt , .value = TRUE )
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @export

grepdt <- function( pattern , dt , value = FALSE )
{
  list() -> l0
  for( i in 1:dim( dt )[2] )
  {
    dt[,i] %>% as.character -> x
    grep( pattern = pattern , x = x  ) -> l0[[i]] }
    if( value )
    {
      l0 %>% unlist %>% unique -> sel
      dt[ sel , ] %>% return()
    }else { l0 %>% unlist %>% unique %>% return() }
}






