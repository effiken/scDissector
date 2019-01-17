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
    titlePanel("scDissector v0.99"),
      navbarPage("",id = "inMain",footer = fluidPage(
        img(src = 'sinai_logo.png',height = '70px', width = '70px'),
        p(em("(c) 2018 Ephraim (Effi) Kenigsberg. Department of Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai"))),
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
                  tabPanel("Gating",
                           fluidRow(
                             selectInput("inGatingSample", "Sample:",choices = c()),
                             selectInput("inGatingShowClusters", "Cluster to highlight:",choices = c()),
                             uiOutput("gaiting_plots_dynamic")
                           )
                  ),
                  tabPanel("Basics", 
                           fluidRow(
                             plotOutput("ncells_barplot",width="100%",height=200),
                             plotOutput("UMI_boxplot",width="100%",height=200)
                          )
                   ),
                  tabPanel("Clusters", 
                    fluidRow(
                          selectInput("inModelOrAverage",label = NULL,choices=c("Model parameters","Batch-corrected Average","Average")),
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
                              br(),
                               actionButton("inFC", "Fold Change screen (FG/BG)"),
                            textInput("inFC_fgClusts","Fold Change FG clusters"),
                            textInput("inFC_bgClusts","Fold Change BG clusters")
                            )
                           ),
                          column(6,
                            textInput("inClusters","Clusters:",width=2000),
                            wellPanel( 
                            actionButton("inOrderClusters", "Reorder Clusters"),
                            selectInput("inReorderingClustersMethod","Reordering method",choices=c("Hierarchical clustering","Diagonal","Cluster-sets","Formula")),
                            textInput("inOrderClusetersByGenes","Order by gene formula:",width=2000),
                            actionButton("inResetClusters", "Reset Clusters"),
                            actionButton("inSaveClusterOrder","Save Order"),
                            actionButton("inUndoClusterChanges","Undo")
                            ),
                            wellPanel(
                              actionButton("inAddClusterSetButton","Add a cluster-set"),
                              textInput("inAddClusterSet", "Name:",value=""),
                              selectInput("removeClusterSetSelectInput","Remove",choices = c()),
                              actionButton("inRemoveClusterSetButton","Remove a cluster-set"),
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
                     textInput("inSamplesToShow", width=2000, "Samples:"),
                     textInput("inSampleColors", width=2000, "Colors:"),
                     actionButton("inResetSamples","Reset"),
                     hr()) 
            )), 
                tabPanel("Cells", fluidRow(
                      plotOutput("truthplot",width="100%",height = "900"),
                      plotOutput("colorLegendTruth",height = 70,width=200),
                      plotOutput("truthSampleLegend",width=200,height=200),
                      sliderInput("inTruthColorScale","log2(1+#UMIs)",min = 0,max = 8,step = .1,value = c(0,2)),
                      selectInput("inTruthDownSamplingVersion",label="Down-sampling version:",choices = c(""),width=200),
                      checkboxInput("inTruthShowSeparatorBars",label = "Show Separator Bars",value = T),
                      selectInput("inTruthNcellsPerSample",label="#Cells Per Sample=",choices=params$nrandom_cells_per_sample_choices,selected = 1000,width=200),
                      downloadButton('downloadTruthPlot', 'Download Plot')
                  )),
                tabPanel("ClusterSets",
                        fluidRow(
                        plotOutput("correlation_betwen_clusters",width=600,height=600)
                        
                        )),
                tabPanel("Clustering QC",fluidRow(
                      column(4,
                        wellPanel(
                          selectInput("inQCClust", "Cluster:",choices=c()),
                          selectInput("inQCDownSamplingVersion",label="Down-sampling version:",choices = c(""),width=200),
                          selectInput("inQCNcellsPerSample",label="#Cells Per Sample=",choices=params$nrandom_cells_per_sample_choices,selected = 1000,width=200)
                          
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
                  textInput("inVarMean_MeanThresh",  "Min Log10(mean):", value = "-2"),
                  textInput("inVarMean_varmeanThresh",  "Min Log2(var/mean):", value = "0.50"),
                  sliderInput("inVarMeanXlim",min=-6,max=2,step=.1,value = c(-1.5,2),label = "X-axis range",width(100))
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
              )),
                tabPanel("Samples",fluidRow(
                         uiOutput("subtype_freqs"),
                         column(6,wellPanel(
                           selectInput("inReorderSamplesBy", "ClusterSet:",choices=c()),
                           actionButton("inSamplesAddClusterSet","Add clusterSet"),
                           textAreaInput("inSamplesClusterSets",label = "ClusterSets to view:",value = ""),
                           actionButton("inReorderSamples","Reorder samples")
                         )),
                          column(12,uiOutput("sample_avg_profile_plot")),
                          column(12,sliderInput("inSamplesColorScale","Log2(expression/mean)",min = -8,max = 8,step = 1,value = c(-4,4))),
                          column(6,textInput("inProjectSampleGroup1","Samples Group 1:")),
                          column(6,textInput("inProjectSampleGroup2","Samples Group 2:")),
                          column(12,selectInput("inProjectPlotType","Plot Type",choices=c("Side by Side","Tile")),
                          uiOutput("projection_barplot_plot"),
                   selectInput("inClustForDiffGeneExprsProjVsRef","Cluster for GE analysis:",choices=c()),
                   wellPanel(textInput("Gene1ForExprsTableRefVsProj","Gene:"),
                   tableOutput("Gene1ExprsTableRefVsProj")),
                   plotOutput("DiffGeneExprsProjVsRef",width="100%",height=500))))
              )
         )
          
  
  
)

          
  

