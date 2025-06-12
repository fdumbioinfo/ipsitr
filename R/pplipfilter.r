#' @title filter for phospholipid class
#' @description filter for phospholipid class
#' @param dat matrix numeric
#' @param factor factor
#' @param samplename character
#' @return results directory
#' @examples
#' # not run
#' # pplipfilter(dat)
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @export
pplipfilter <- function( dat , rtinfmin = NULL , rtsupmin = NULL , dirname = NULL, path = "." )
{  
  dat -> Dat0
  i=1
  #
  ifelse( is.null(dirname), "pplipfilter" -> DirName0, paste("pplipfilter_",dirname,sep="") -> DirName0 )
  path %>% file.path(DirName0) -> Path0
  if(!dir.exists(Path0)){ Path0 %>% dir.create }
  #
  if(is.null(rtinfmin) & is.null(rtsupmin))
  {
    RtInfMin <- c("9.6","11.0","11.5","13.1","14.8","16.0","17.9","19.8","22.0","12.5","8.0")
    RtSupMin <- c("10.1","11.5","12.6","14.5","16.3","17.0","19.8","22.2","23.5","13.0","9.0")
  }else{
    RtInfMin <- RtInfMin
    RtSupMin <- RtSupMin }
  # création de la matrice contenant les valeurs de filtre rt et mz
  data.frame(
    Famille = c("PG","CL","PI","PE","PS","LPE","PC","SM","LPC","MLCL","PA"),
    RtInfMin = RtInfMin,
    RtSupMin = RtSupMin,
    stringsAsFactors = F ) -> FilterO
  # pour convertir les minutes en seconde
  FilterO %>% 
    mutate( 
      RtInfsec = as.numeric(sub("(.*)\\.(.*)","\\1",RtInfMin,perl=T))*60+as.numeric(sub("(.*)\\.(.*)","\\2",RtInfMin,perl=T))*6,
      RtSupsec = as.numeric(sub("(.*)\\.(.*)","\\1",RtSupMin,perl=T))*60+as.numeric(sub("(.*)\\.(.*)","\\2",RtSupMin,perl=T))*6) %>%
    mutate(
      MzInf = c(650,1200,800,650,700,400,750,650,400,1100,650),
      MzSup = c(850,1600,950,900,900,700,950,900,700,1310,800) ) %>%
    mutate( PAIR = c("IMPAIR","CL","IMPAIR","PAIR","PAIR","PAIR","PAIR","IMPAIR","PAIR","PAIR","PAIR")) -> Filter1
  # output
  Filter1 %>% output("pplipfilter_rtmz.tsv")
  # add decimal column and mz interger for pair/impair filtering
  Dat0 %>% mutate(mzEnt=sub("(.*)\\..*","\\1",mz) %>% as.numeric,
                  mzDec=substr(sub( ".*\\.(.*)","\\1",mz),1,1) %>% as.numeric) -> Dat1
  #
  foreach(i=1:nrow(Filter1)) %do%
  {
    # Filter1 %>% "[["(i,"Famille") -> famlip0
    # Filter1 %>% "[["(i,"PAIR") -> pair
    # # filtre rt
    # Filter1 %>% "[["(i,"RtInfsec") -> RtInf
    # Filter1 %>% "[["(i,"RtSupsec") -> RtSup
    # Dat1$rt
    # Filter1$RtInfMin[i]
    # Dat1 %>% dplyr::filter(.data$rt > Filter1$RtInfsec[i] & .data$rt < Filter1$RtSupsec[i]) -> Dat2
    Dat1 %>% dplyr::filter(.data$rt > Filter1$RtInfMin[i] %>% as.numeric & .data$rt < Filter1$RtSupMin[i] %>% as.numeric) -> Dat2
    # # filtre mz
    # Filter1$Famille[i]
    # Filter1$RtInfMin[i]
    # Filter1 %>% dplyr::filter(.data$Famille == Filter1$Famille[i]) %>% dplyr::select(Filter1$RtInfMin[i]) %>% as.numeric -> MzInf
    # Filter1 %>% dplyr::filter(.data$Famille == famlip0) %>% dplyr::select(MzSup) %>% as.numeric -> MzSup
    #
    if(Filter1$Famille[i] == "CL")
    {
      Dat2 %>% dplyr::filter(mz > Filter1$RtInfMin[i] & mz < Filter1$RtSupMin[i]) -> t
      which(
        ( (Filter1$mzEnt[i] %% 2) == 0 & Filter1$mzDec[i] <= 1 ) |
        ( (Filter1$mzEnt[i] %% 2) != 0 & Filter1$mzDec[i] >= 9 ) ) -> sel
      t[sel,] -> Dat3
    }else{
      if(Filter1$PAIR[i] == "PAIR")
      {
        Dat2 %>% dplyr::filter(mz > Filter1$MzInf[i] & mz < Filter1$MzSup[i]) %>%
          dplyr::slice(which(Filter1$mzEnt[i] %% 2 ==0 )) -> Dat3
      }else{
        Dat2 %>% dplyr::slice( which(!c( Filter1$mzEnt[i] %% 2 == 0 ))) %>%
          dplyr::filter( mz > Filter1$MzInf[i] & mz < Filter1$MzSup[i] ) -> Dat3 }
    }
    #
    Dat3 %>% dplyr::select( -c(mzEnt,mzDec) ) -> Dat4
    # output
    paste(Filter1$Famille[i],sep="") -> DirName0
    Path0 %>% file.path(DirName0) %>% dir.create
    paste("peaktablexcms_",Filter1$Famille[i],"_",dim(Dat3)[1],".txt",sep="") -> FileName0
    Path0 %>% file.path(DirName0,FileName0) -> FileName1
    Dat4 %>% output(FileName1)
  }
}