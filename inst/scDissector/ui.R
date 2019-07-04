library(ggvis)
library(shinyTree)
# For dropdown menu

actionLink <- function(inputId, ...) {
  tags$a(href='javascript:void',
         id=inputId,
         class='action-button',
         ...)
}

#vers=read.csv("data/versions.csv",stringsAsFactors = F)
#annnot_fn=paste(strsplit(vers$path[1],"\\.")[[1]][1],".annots",sep="")
#if (file.exists(annnot_fn)){
#  clustnames=read.table(annnot_fn)
#}else {
#  clustnames=1:vers$k[1]
#}
mainPanel(width=12,  
  fluidPage(
   
  #    h4("scDissector v0.99"),
    titlePanel("scDissector v1.00"),
      navbarPage("",id = "inMain",footer = fluidPage(
        img(src = 'sinai_logo.png',height = '70px', width = '70px'),
        p(em("(c) 2019 Ephraim (Effi) Kenigsberg. Department of Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai"))),
  #   tabsetPanel(id="inMain",
                  tabPanel("Data", fluidRow(column(12,
                                  br(),
                                  textInput("inDatapath",width="75%", "Data path:", value =""), 
                                  wellPanel(
                                           selectInput("inModelVer", "Model Version:",choices = c(),width=300),
                                          selectInput("inSampleToAdd", "Samples:",choices = c(),width=300),
                                          actionButton("inAddSample","Add"),
                                          textInput("inSamples", width="100%", "Samples To Load:"),
                                          textInput("inMinUmis", width=100, "Minimum #UMIs/cell:",value = 800),
                                          textInput("inMaxUmis", width=100, "Maximum #UMIs/cell:",value = 25000),
                                          actionButton("inLoad","Load")
                                  ),
                                  br(),br(),
                  #                h5("Import Data and Clustering"),
                  #                wellPanel(
                  #                  textInput("inImportUmitabFile",width=1000, "UMI Table File:", value =""),
                  #                  textInput("inImportCellToClusterFile",width=1000, "Cell To Cluster File:", value =""),
                  #                  textInput("inImportVersionName", "Imported Version Name:", value =""),
                  #                  actionButton("inImport","Import")
                  #                ),
                  br(),br(),br(),br(),br(),br(),br(),br())
                                  )),
          tabPanel("MetaData",
                   fluidRow(
                     DT::dataTableOutput("mytable",height = 600)),fluidRow(
                       actionButton("selectAllSamples","Select All"),
                       actionButton("selectSamples","Copy Selected Ids"),
                       textInput("inSamplesToShow", width=2000, "Samples:")
                     )),
           tabPanel("Gating",
                           fluidRow(
                             selectInput("inGatingSample", "Sample:",choices = c()),
                             selectInput("inGatingShowClusters", "Cluster to highlight:",choices = c()),
                             uiOutput("gaiting_plots_dynamic")
                           )
                  ),
                  tabPanel("Basics", 
                           fluidRow(
                             plotOutput("correlation_betwen_clusters",width=600,height=600),
                             plotOutput("ncells_barplot",width="100%",height=200),
                             plotOutput("UMI_boxplot",width="100%",height=200)
                          )
                   ),
                  tabPanel("Clusters", 
                    fluidRow(
                          selectInput("inModelOrAverage",label = NULL,choices=c("Model parameters","Batch-corrected Average","Average"),selected = "Model parameters"),
                          column(1,offset=0, style='padding:0px;',uiOutput("samples_enrichment")),column(10,uiOutput("avg_profile_plot")),column(1,offset=0, style='padding:0px;',uiOutput("cluster_sizes")),
                          column(4,plotOutput("modelSampleLegend",height = 200)
                                 ),
                          column(4,
                          sliderInput("inModelColorScale","Log2(expression/mean)",min = -8,max = 8,step = 1,value = c(-4,4)),
                          selectInput("inAbsOrRel","Expression Values",choices=c("Relative","Absolute"))),
                          column(4,
                                 plotOutput("colorLegendModel",height = 70,width=200)
                                 ),
                          column(12,
                                  hr()),
                          column(6,
                            textInput("inGenes", width=2000, "Genes:"),
                            wellPanel(
                            selectInput("inGeneSets", "Gene Sets:",choices = c()),
                            actionButton("inClusterGenes", "Reorder Genes"),
                            selectInput("inReorderingMethod","Reordering method",choices=c("Hierarchical clustering","Diagonal"))),
                            wellPanel(
                              h4("Gene Selection"),
                              selectInput("inSelectGenesFrom","Select From:",choices=c("All genes","Without RPs","TFs","Surface markers")),
                              selectInput("inNgenes", "N=",choices =c(50,100,200,300,400,500,1000),width=100),
                              actionButton("inBlindChisqSelectedClusters", "Chi sq. screen on selected clusters"),
                              actionButton("inVarMeanScreenSelectedClusters", "VarMean screen on selected clusters"),
                              sliderInput("inMinExprForScreen","Min expression (log10)",min = -7,max = -1,step = .1,value=c(-5,-1)),
                              br(),
                               actionButton("inFC", "Fold Change screen (FG/BG)"),
                            textInput("inFC_fgClusts","Fold Change FG clusters"),
                            textInput("inFC_bgClusts","Fold Change BG clusters")
                            ),
                            selectInput("categorizeSamplesBy",label = "Categorize Samples By:",choices=c())
                           ),
                          column(6,
                            textInput("inClusters","Clusters:",width=2000),
                            wellPanel( 
                            actionButton("inOrderClusters", "Reorder Clusters"),
                            selectInput("inReorderingClustersMethod","Reordering method",choices=c("Hierarchical clustering","Diagonal","Cluster-sets","Formula")),
                            textInput("inOrderClusetersByGenes","Order by gene formula:",width=2000),
                            actionButton("inResetClusters", "Reset Clusters"),
                            actionButton("inSaveClusterOrder","Save Order")
                            ),
                            wellPanel(
                              h4("Cluster-sets"),
                              fluidRow(column(6,textInput("inAddClusterSet", "Type in name",value="")),column(6,h2(""),actionButton("inAddClusterSetButton","Add"))),
                              fluidRow(column(6,selectInput("removeClusterSetSelectInput","Select to remove",choices = c())),column(6,h2(""),actionButton("inRemoveClusterSetButton","Remove"))),
                              actionButton("saveClusterSetButtion","Save Cluster-sets"),
                              actionButton("reloadClusterSetButtion","Reload Cluster-sets")
                            ),
                            shinyTree("clusters_sets_shinytree", theme="proton", themeIcons = FALSE, themeDots = FALSE,dragAndDrop=T)
          #                  wellPanel(
          #                    selectInput("inAnnotateClusterNum","Cluster",choices=c()),
          #                  textInput("inClustAnnot", "New Annotation:"),
          #                  actionButton("inAnnotateCluster","Annotate"),
          #                  actionButton("inSaveAnnot", "Save Annotations"))
                            ),
                    column(12,
                     hr(),
                     hr()) 
            )), 
                tabPanel("Cells", fluidRow(
                      plotOutput("truthplot",width="100%",height = "900"),
                      plotOutput("colorLegendTruth",height = 70,width=200),
                      plotOutput("truthSampleLegend",width=200,height=200),
                      sliderInput("inTruthColorScale","log2(1+#UMIs)",min = 0,max = 8,step = .1,value = c(0,2)),
                      textInput("inGenesScore", width=2000, "Order cells within clusters by the following genes:"),
                      checkboxInput("inScoreReverseOrder",label = "Reverse Order",value = F),
                      selectInput("inTruthDownSamplingVersion",label="Down-sampling version:",choices = c(""),width=200),
                      checkboxInput("inTruthShowSeparatorBars",label = "Show Separator Bars",value = T),
                      wellPanel(radioButtons("inTruthEqualRepresent", "Equally represent:",
                                             c("Samples" = "equally_represent_samples",
                                               "Clusters" = "equally_represent_clusters")),
                        selectInput("inTruthNcellsPer",label="#Cells sampled",choices=params$nrandom_cells_per_sample_choices,selected = 250,width=200)),
                      downloadButton('downloadTruthPlot', 'Download Plot')
                  )),
                tabPanel("Samples",
                  column(12,wellPanel(
                    selectInput("inReorderSamplesBy", "ClusterSet:",choices=c()),
                    actionButton("inSamplesAddClusterSet","Add clusterSet"),
                    textAreaInput("inSamplesClusterSets",label = "ClusterSets to view:",value = ""),
                    actionButton("inReorderSamples","Reorder samples")
                  )),
                  column(12,
                  uiOutput("subtype_freqs"),
                  hr()),
                  column(12,uiOutput("sample_avg_profile_plot")),
                  column(12,sliderInput("inSamplesColorScale","Log2(expression/mean)",min = -8,max = 8,step = 1,value = c(-4,4)))),
                tabPanel("Clustering QC",fluidRow(
                      column(4,
                        wellPanel(
                          selectInput("inQCClust", "Cluster:",choices=c()),
                          selectInput("inQCDownSamplingVersion",label="Down-sampling version:",choices = c(""),width=200),
                          selectInput("inQCNcellsPerSample",label="#Cells Per Sample=",choices=params$nrandom_cells_per_sample_choices,selected = 250,width=200)
                          
                        )
                      ),
                      column(8,
                        plotOutput("cellcor")
                      ),
                      column(3,
                        wellPanel(sliderInput("mean", "Log10(mean #UMIs)", -3, 2, value = c(-1, 1), step=.1),
                        sliderInput("varmean", "Log2(variance/mean)", -1, 10, value = c(0, 8), step = .5),
                        h4("two genes:"),
                        textInput("inGene1",  "Gene1:", value = "IRF8"),
                        textInput("inGene2",  "Gene2:", value = "ACTB"))
                      ),
                      column(9,
                        ggvisOutput("varmeanplot"),
                        h4("Table (log2(1+#UMIs)"),
                        textOutput("gene2"),
                        textOutput("gene1"),
                        tableOutput("table")
                       )
                      )
                  ),
                tabPanel("Gene Modules",fluidRow(
                wellPanel(
                  plotOutput("varMeanThreshPlot"),
                  selectInput("inModulesDownSamplingVersion",label="Down-sampling version:",choices = c(""),width=200),
                  fluidRow(column(3,sliderInput("inVarMeanXlim",min=-6,max=2,step=.1,value = c(-2,2),label = "X-axis range",width(100)))),
                  textInput("inVarMean_MeanThresh",  "Min Log10(mean):", value = "-2"),
                  textInput("inVarMean_varmeanThresh",  "Min Log2(var/mean):", value = "0.50")
                 ),
                wellPanel(actionButton("inModulesGetCormap","Get Correlation Map"),
                uiOutput("cor_gene_module_plot")),
                wellPanel(
                  textInput("inNUmberOfGeneModules","Number of Modules:",value = 50),
                  uiOutput("cor_module_plot"),
                  uiOutput("avg_module_plot")),
                sliderInput("inAvgModuleColorScale","Log2(expression/(Module_mean))",min = -4,max = 4,step = 1,value = c(-2,2)),
                       textInput("inModules", width=2000, "Modules:"),
                       wellPanel(
                         actionButton("inClusterModules", "Reorder Modules"),
                         selectInput("inReorderingModulesMethod","Reordering method",choices=c("Hierarchical clustering","Diagonal")),
                         downloadButton('inDownloadModuleList', 'Download Modules')
                         ),
                      wellPanel(
                        selectInput("inModuleSelect","Show Module:",choices=c()),
                        textOutput("textModuleGenes")
                      )
              ))
            )
         )
          
  
  
)

          
  

