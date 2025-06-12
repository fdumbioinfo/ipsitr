#' @title replace value by group median
#' @param dat character vector of mzxml file
#' @param value numeric value to substitute
#' @param factor character
#' @details
#' replace value by column median if group size is one and replace by row group median if group size > 1
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
replacegroupmed <- function( dat , value = 0 , factor = NULL )
{
  i=j=1
  dat %>% "=="(value) %>% which(arr.ind=T) %>% data.frame -> Coord
  Coord %>% head
  Coord %>% dim
  Coord$col %>% table
  # replace with median group
  if( nrow(Coord) > 0 & !is.null(factor) )
  {
    edgeR::calcNormFactors(dat) -> NormFactors
    factor %>% table %>% "=="(1) %>% which %>% names -> Levels
    Levels %>% paste("^",.,"$",sep="") %>% paste0(collapse = "|") -> Grep0
    factor %>% as.character %>% grep(Grep0,.,invert = T) -> Selcol
    Selcol %>% paste("^",.,"$",sep="") %>% paste0(collapse = "|") -> Grep1
    Coord$col %>% grep(Grep1,.) -> Selrow
    Coord[Selrow,] -> Coord1
    Coord1 %>% head
    Coord1 %>% dim
    Coord1$col %>% table
    foreach(i=1:nrow(Coord1)) %do%
      {
        Coord1[i,]
        Coord[i,2]
        factor %>% table
        factor[Coord1[i,2]] %>% as.character %>% paste("^",.,"$",sep="") -> Grep0
        factor %>% grep(Grep0,.) -> sel
        dat[i,sel] %>% as.numeric %>% median -> Median0
        Median0 * NormFactors[Coord[i,2]] -> Median1
        dat[Coord1[i,1],Coord1[i,2]] <- Median1
      }
    
    Dat1 <- foreach( i=1:length(factor %>% levels) , .combine = "cbind") %do%
      {
        factor %>% levels %>% "["(.,i) %>% "=="(.,factor) %>% which -> sel
        dat[,sel] %>% as.data.frame -> Dat0
        Dat0 %>% unlist %>% median -> MedGroupCol
        Dat0 %>% "=="( ., value ) %>% which( arr.ind=T ) -> selvalue
        if( nrow(selvalue) > 0 )
        {
          foreach( j=1:nrow(selvalue) ) %do%
            {
              Dat0[selvalue[j,1],] %>% as.numeric %>% median -> MedGroupRow
              if( MedGroupRow %>% "=="(.,0) ){ Dat0[ selvalue[j,1] , selvalue[j,2] ] <- MedGroupCol }else
              { Dat0[ selvalue[j,1] , selvalue[j,2] ] <- MedGroupRow }
            }
        }
        Dat0 %>% as.data.frame
      }
  }
  # 
  Coord$col %>% table %>% names %>% as.numeric -> Ncol
  if( nrow(Coord) > 0 )
  {
    foreach( i=1:length(Ncol) ) %do%
      { 
        Coord %>% dplyr::filter(col == Ncol[i] ) %>% dplyr::select(row) %>% unlist -> sel
        dat[sel,Ncol[i]] <- NormFactors[Ncol[i]] 
      }
  }
  
  
  Dat1 %>% setNames( colnames(dat) ) %>% return()
}