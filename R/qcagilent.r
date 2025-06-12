#'  @title QC analysis for Agilent rawdata
#' @description QC analysis for Agilent rawdata expression array.
#' @param files character list of rawdata files
#' @param sif data.frame with sample information file
#' @param path character
#' @details
#' Results directory contain boxplot, acp, hierachical clustering for rawdata and dotplot of QC probes
#' A sample information file (sif) is created to complete with experimental factor variable.
#' Re run QCAgilent with competed sif for generated plot with that correspond to factor levels.
#' @return
#' Create directory QC results
#' @examples
#' # not run
#' # pathfile %>% qcagilent
#' @author Florent Dumont <florent.dumont@universite-paris-saclay.fr>
#' @importFrom magrittr %>%
#' @import dplyr
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom utils stack
#' @importFrom grDevices pdf graphics.off
#' @importFrom graphics legend title
#' @importFrom magrittr %>%
#' @export

qcagilent <- function( files, sif = NULL, path = "." )
{
   path %>% file.path( "QCAgilent" ) -> Path
   Path %>% dir.create
   files %>% agilentrawdataextraction( output = T , path = Path ) -> re0
   re0$m0 -> m0
   re0$m1 -> m1
   ### agilent QC control
   # E1A
   m0 %>% slice( grep( "E1A" , m0$ProbeName ) ) -> e1a0
   data.frame(
      E1AType  = c("3","a104","a107","a135","a20","a22","a97","n11","n9","1" ),
      E1AConcentrationpgL = c(0.04,0.4,4,40,133,400,1333,4000,13333,40000),
      stringsAsFactors = F) -> e1ainfo
   # ea1 type
   e1a0 %>%
      data.frame(
         ProbeName = e1a0$ProbeName,
         E1AType = e1a0$ProbeName %>% sub( ".*_.*_(.*)" , "\\1" , . ),
         . , stringsAsFactors = F ) -> e1a1
   # joining e1a values with concentration
   e1a1 %>%
      inner_join( e1ainfo , by = "E1AType" ) %>%
      select( c(1,2,(3+length(files)+1),4:(3+length( files ) ) ) ) %>%
      arrange( .data$E1AConcentrationpgL )  -> e1a2

   # recording rawdata in tsv

   paste("QCpos_E1A_rawdata_",dim(e1a2)[2]-3 , "_" , dim(e1a2)[1],".tsv" , sep = "" ) -> FileName
   Path %>% file.path( FileName ) -> FileName
   e1a2 %>% output( FileName  )

   # creating sample information file

   if( is.null( sif ) )
   {
      files %>% basename -> FileName

      data.frame(
         "SampleID" = colnames(m1)[-1],
         FileName,
         stringsAsFactors = F ) -> Sif

      paste( "metadata_", length(files),".tsv", sep = "" ) -> FileName
      Path %>% file.path( FileName ) -> FileName
      Sif %>% output( FileName  )

   }else
      {
         paste("metadata_",length( files ),".tsv",sep="") -> FileName
         Path %>% file.path( FileName ) -> FileName
         sif %>% output( FileName )
      }

   # E1A plot

   rep( as.character(e1a2$E1AType), dim(e1a2)[2]-3  ) -> color

   # joining dot in the plot

   rep(
      paste( "e1a" , 1:dim(e1a2)[1] , sep = ""),
      dim( e1a2 )[2]-3 ) -> group

   e1a2 %>%
      select( -c(1,2,3) ) %>%
      log2 %>%
      data.frame %>%
      stack %>%
      ggplot( aes(x = .data$ind , y = .data$values , color = color , group = group ) ) +
      geom_point( fill = "white" , size = 0.2) +
      geom_line( size = 0.2) +
      xlab("Samples") +
      theme( axis.text.x = element_text( face = "plain" , color="black" , size = 7, angle = 90 ) ) +
      theme( axis.text.y = element_text( face = "plain" , color = "black" , size = 8, angle = 0 ) ) +
      theme( legend.position = "right" ) +
      ggtitle( paste( "QCpos_rawdata_log2_E1A_", dim(e1a2)[2]-3 , "_10x35" , sep = "" ) ) +
      guides( color = guide_legend( "E1A Type" )  )

   # recording plot

   paste( "QCpos_rawdata_log2_E1A_", dim(m0)[2]-1,"_10x35.pdf" , sep = "") -> FileName
   Path %>% file.path( FileName ) -> FileName
   ggsave( FileName )

   # 3xSLv1

   m0 %>%
      slice( grep("3xSLv1" , m0$ProbeName) ) -> slv10

   paste( "QCneg_3xSLv1_rawdata_", dim(slv10)[2]-1 , "_" , dim(slv10)[1], ".tsv" , sep = "" ) -> FileName
   Path %>% file.path( FileName ) -> FileName
   slv10 %>% output( FileName )

   # 3xSLv1 plot

   slv10 %>%
      select(-1 ) %>%
      data.frame %>%
      stack %>%
      ggplot( aes(x = .data$ind , y = .data$values  ) ) +
      geom_point( fill = "white") +
      xlab("Samples") +
      theme( axis.text.x = element_text( face = "plain" , color="black" , size = 7, angle = 90 ) ) +
      theme(axis.text.y = element_text( face = "plain" , color = "black" , size = 8, angle = 0 ) ) +
      theme(legend.position = "none") +
      ggtitle( paste("QCneg_rawdata_3xSLv1_",dim(m0)[2]-1,"_1x308" , sep = "") )


   paste( "QCneg_rawdata_3xSLv1_", dim(m0)[2]-1,"_1x308.pdf" , sep = "") -> FileName
   Path %>% file.path( FileName ) -> FileName
   ggsave( FileName )

   ### QC

   if( is.null( sif ) )
   {
      qc( dat = m1 , path = Path , dirname = "QCrawdata" )
   }else
      { qc( dat = m1 , sif = sif, path = Path, inputdata = F, dirname = "QCrawdata" ) }

}






