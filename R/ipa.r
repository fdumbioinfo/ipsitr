#' @title ipa
#' @description parse and create plot from ipa all enrichment results
#' @param file character path to the file
#' @param top numeric top feature to display
#' @param labsize numeric feature size 
#' @param dpi character jpeg resolution between retina print screen 
#' @param path character for relative path of output directory
#' @param dirname character name for output
#' @return no values
#' @examples
#' # not run
#' # ipa( file , path "" )
#' @author Florent Dumont <florent.dumont@univresite-paris-saclay.fr>
#' @import ggplot2
#' @importFrom rlang .data
#' @export
ipa <- function( 
    file = "", top = 80, labsize = 8 , dpi = "retina",
    overlapmin = 2 , enaScoremin = 1, background = 25000,
    path = "." , dirname = NULL )
{
  i=1
  # output directory
  ifelse( is.null(dirname), "ena_ipa" -> DirName, paste("ena_IPA_",dirname,sep="") -> DirName )
  path %>% file.path( DirName ) -> Path ; Path %>% dir.create
  #
  # read file
  #
  file %>% file("r") -> f
  f %>% readLines() -> rl0
  close(f)
  #
  rl0 %>% grep("^$",.) -> sel
  c("Analysis Details","Ingenuity Canonical Pathways" , "Upstream Regulators",
    "Diseases and Bio Functions","Tox Functions","Regulator Effects","Analysis Ready Molecules") -> Grep0
  #
  # Analysis setting
  #
  rl0 %>% grep(Grep0[1],.) -> Grep1
  rl0[Grep1:(sel[3]-1)] -> rl1
  "IPACoreAnalysisSettings.txt" %>% file.path(Path , . ) %>% file("w") -> f
  rl1 %>% writeLines(f)
  close(f)
  #
  # input List
  #
  rl0 %>% grep(Grep0[7],.) -> Grep1
  rl0[(Grep1+1):(sel[length(sel)]-1)] -> rl1
  rl1[-1] %>% strsplit("\t") %>% unlist -> rl2
  length(rl1) -> listsize
  rl1[1] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  header[1] <- "rowID"
  rl2 %>% matrix( ncol = header %>% length , byrow = T ) %>% data.frame %>% setNames(header) -> rl3
  rl3$rowID %>% strsplit("\\/") %>% lapply("[",1) %>% unlist -> rl3$rowID
  rl3$Symbol %>% strsplit("\\/") %>% lapply("[",1) %>% unlist -> rl3$Symbol
  paste("InputList_",listsize,".tsv",sep="") %>% file.path(Path,.) -> FileName0
  rl3 %>% output(FileName0)
  #
  # Canonical Pathways
  #
  # create results table
  rl0 %>% grep(Grep0[2],.) -> Grep2
  rl0[Grep2:(sel[4]-1)] -> rl1
  rl1[-1] %>% strsplit("\t") %>% unlist -> rl2
  rl1[1] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  header[1] <- "Pathways"
  header[2] <- "log10pval"
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  fO <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  fO -> rl3$Molecules
  rl3$Ratio %>% as.numeric -> rl3$Ratio
  rl3$log10pval %>% as.numeric -> rl3$log10pval
  # compute enrichment score
  f0 <- foreach( i=1:nrow(rl3), .combine = "rbind" ) %do% 
    {
      rl3[i,5] %>% strsplit("\\|") %>% unlist %>% length %>% as.numeric -> OverlapSize
      OverlapSize %>% "/"( rl3[i,4] %>% as.numeric ) %>% round -> GenesetSize
      OverlapSize / GenesetSize -> OverlapRatio
      listsize * ( GenesetSize / background ) -> exp
      OverlapSize / exp -> ENAScore
      c(OverlapSize,listsize-OverlapSize,GenesetSize-OverlapSize,background-OverlapSize-GenesetSize) %>%
        matrix(ncol = 2) %>% fisher.test %>% unlist %>% "["(1) %>% as.numeric -> pval
      c(OverlapSize, GenesetSize, OverlapRatio, ENAScore, pval)
    }
  f0[,5] %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10 %>% '*'(.,-1) %>% round(.,4) -> log10pvalFDR
  c("Name","SymbolList","OverlapSize","GenesetSize","OverlapRatio","ENAScore","pval","pvalFDR","log10pvalFDR") -> Header0
  rl3 %>% data.frame(f0,pvalFDR,log10pvalFDR) %>% dplyr::select( c(1,5,6:10,11,12) ) %>% setNames(Header0) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enaScoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("CanonicalPathways_",nrow(rl5),".tsv",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  rl5 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(c(1,9,6)) -> rl6
  ipabarplotcp( dat = rl6, top = top , labsize = labsize ) -> p
  # output
  paste("CanonicalPathways_barplot_top",top,".jpeg",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  ggsave( filename=FileName1 , plot = p , width = 12 , height = 10 , dpi = dpi )
  #
  # upstream regulators
  #
  # create results table
  rl0 %>% grep(Grep0[3],.) -> Grep3
  rl0[Grep3:(sel[5]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$UpstreamRegulator -> Name
  rl3$MoleculeType -> MoleculeType
  rl3$`p-valueofoverlap` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  rl3$Targetmoleculesindataset %>% lapply(FUN=gsub,pattern=",",replacement="\\|") %>% unlist -> SymbolList
  f0 <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Targetmoleculesindataset[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  f0 -> SymbolList
  SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist -> OverlapSize
  OverlapSize -> ENAScore
  data.frame(Name,MoleculeType,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enaScoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("UpstreamRegulators_",nrow(rl5),".tsv",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  rl5 %>% output(FileName1)
  # barplot
  rl5 %>% colnames
  rl5 %>% dplyr::select( Name, log10pvalFDR, OverlapSize ) -> rl6
  ipabarplotup( dat=rl6, top=top, labsize=labsize ) -> p
  # output
  paste("UpstreamRegulators_barplot_",top,".jpeg",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  ggsave( filename=FileName1 , plot = p , width = 12 , height = 10 , dpi = "retina" )
  #
  # Functions
  #
  # create results table
  rl0 %>% grep(Grep0[4],.) -> Grep4
  rl0[Grep4:(sel[7]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$DiseasesorFunctionsAnnotation -> Name
  rl3$Categories -> IPACategories
  rl3$Functions -> IPAFunctions
  rl3$Molecules %>% lapply(FUN=gsub,pattern=",",replacement="\\|") %>% unlist -> SymbolList
  SymbolList <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist %>% as.numeric -> OverlapSize
  OverlapSize -> ENAScore
  rl3$`p-Value` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  data.frame( Name,IPACategories,IPAFunctions,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR ) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enaScoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("Functions_",nrow(rl5),".tsv",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  rl5 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(Name,log10pvalFDR,OverlapSize) -> rl6
  ipabarplotfun( dat = rl6, top = top , labsize = labsize ) -> p
  # output
  paste("Functions_barplot_",top,".jpeg",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  ggsave( filename=FileName1 , plot = p , width = 12 , height = 10 , dpi = "retina" )
  #
  # Tox function
  #
  # create results table
  rl0 %>% grep(Grep0[5],.) -> Grep5
  rl0[Grep5:(sel[8]-1)] -> rl1
  rl1[-c(1:2)] %>% strsplit("\t") %>% unlist -> rl2
  rl1[2] %>% strsplit("\t") %>% unlist %>% gsub(" ","",.) -> header
  rl2 %>% matrix(ncol = header %>% length , byrow = T) %>% data.frame %>% setNames(header) -> rl3
  rl3$DiseasesorFunctionsAnnotation -> Name
  rl3$Categories -> IPACategories
  rl3$Functions -> IPAFunctions
  SymbolList <- foreach(i=1:nrow(rl3),.combine="c") %do% 
    { 
      rl3$Molecules[i] %>% gsub(" \\(includes others\\)","",.) %>% strsplit(",") %>% unlist %>%
        strsplit("\\/") %>% lapply("[",1) %>% unlist %>% paste0(collapse = "|")
    }
  SymbolList %>% strsplit("\\|") %>% lapply(length) %>% unlist %>% as.numeric -> OverlapSize
  OverlapSize -> ENAScore
  rl3$`p-Value` %>% as.numeric -> pval
  pval %>% p.adjust(method = "fdr") -> pvalFDR
  pvalFDR %>% log10(.) %>% "*"(-1) -> log10pvalFDR
  data.frame( Name,IPACategories,IPAFunctions,SymbolList,OverlapSize,ENAScore,pval,pvalFDR,log10pvalFDR ) -> rl4
  rl4 %>% dplyr::filter(OverlapSize > overlapmin & ENAScore > enaScoremin) %>% dplyr::arrange( desc(log10pvalFDR ) ) -> rl5
  # output
  paste("ToxFunctions_",nrow(rl5),".tsv",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  rl5 %>% output(FileName1)
  # barplot
  rl5 %>% dplyr::select(Name,log10pvalFDR,OverlapSize) -> rl6
  ipabarplottox( dat=rl6, top=top, labsize=labsize ) -> p
  # output
  paste("ToxFunctions_barplot_",top,".jpeg",sep="") -> FileName0
  FileName0 %>% file.path(Path,.) -> FileName1
  ggsave( filename=FileName1 , plot = p , width = 12 , height = 10 , dpi = "retina" )
}

