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
#' @importFrom graphics hist
#' @importFrom foreach  %do% foreach
#' @importFrom stats median setNames quantile
#' @importFrom edgeR calcNormFactors
#' @export
replacegroupmed2 <- function( dat , value = 0 , factor = NULL )
{
  i=j=1
  dat %>% "=="(value) %>% which(arr.ind=T) %>% data.frame -> Coord
  Coord %>% head
  Coord %>% dim
  Coord$col %>% table
  factor %>% table -> FactorTable
  factor %>% levels %>% length -> NGroup
  dat %>% unlist %>% stats::quantile(.) -> QuantAll
  FactorTableSel <- foreach(i=1:length(FactorTable)) %do%
    {
      FactorTable[i] %>% names %>% paste("^",.,"$",sep="") -> Grep0
      factor %>% grep(Grep0,.)
    }
  dat %>% apply(1,quantile) -> QuantRow
  dat %>% apply(1,quantile) -> QuantCol
  dat %>% unlist %>% graphics::hist(breaks = 10) -> h0
  h0 %>% class
  
  dat %>% apply(1,hist,breaks=100) -> QuantCol10
  QuantCol10[[1]][] %>% head
  dat %>% ncol %>% "/"(3) %>% "*"(2)
  dat %>% apply(2,"==",0) %>% apply(1,which) %>% lapply(length) %>% unlist -> n0
  n0 %>% length
  n0 %>% head
  n0 %>% ">"(0) %>% which -> sel
  sel %>% length
  dat %>% dplyr::slice(sel) %>% data.frame -> dat1
  dat1 %>% dim
  dat -> dato
  foreach(i=1:nrow(Coord)) %do%
    {
      factor[Coord[i,2]] %>% paste("^",.,"$",sep="") -> Grep0
      factor %>% grep(Grep0,.) -> selt
      dat[Coord[i,1],selt] %>% "=="(value) %>% which -> seltt
      if(seltt %>% length %>% "/"(length(selt)) %>% ">"(0.6))
      {
        QuantAll %>% "+"(QuantCol[Coord[i,2]]) %>% "/"(2)
      }
      
    }
  
  
  foreach(i=1:nrow(dat)) %do%
    {
      foreach(j=1:NGroup) %do%
        {
          dat[i,FactorTableSel[j] %>% unlist]
          
          
        }
      dat1  
      
      
      
    }
  
  dat %>% apply(1,"==",0) %>% lapply(which) %>% length
  dat %>% apply(2,quantile)
  # 1 remove row with to much 0 2/3
  # replace with median group
  # si 1 replicat -> remplace par une valeur globale de bruit pondere -> moyenne( 1er quart de la matrice * 1 quartile de la colonne * 1er quatile de ligne)
  # si replicat > 1 -> remplace par la mediane du groupe
  # si replicat > 1 et que des 0 -> remplace par une valeur globale de bruit pondere -> moyenne( 1er quart de la matrice * 1 quartile de la colonne * 1er quatile de ligne)
  #
  # 1
  #
  # dat %>% apply(==,0)
  
  
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