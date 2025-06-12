#' @title replace value in data.frame
#' @param dat data.frame
#' @param comptype character for comparison operator
#' @param compvalue numeric or character value
#' @param replacevalue numeric or character value for replacement
#' @param dopar logical to use parallel computing
#' @return data.frame
#' @examples
#' # not run
#' # replacemat(.dat)
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @export

replacemat <- function( dat , comptype = "<=" , compvalue = 0 , replacevalue = 1 , dopar = TRUE )
{
  dat %>% unlist %>% sapply( comptype , compvalue ) %>% which -> sel
  dat %>% unlist %>% replace(  list = sel , values = replacevalue ) -> Dat0
  Dat0 %>% unlist %>% matrix( nrow = nrow(dat) ) %>% as.data.frame -> Dat1
  Dat1 %>% return()
}










