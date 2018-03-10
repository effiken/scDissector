library(ggvis)

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
   
      h4("scDissector v0.62"),
     tabsetPanel(id="inMain",
                  tabPanel("Data", fluidRow(column(12,
                                  br(),
                                  textInput("inDatapath",width="75%", "Data path:", value =""), 
                                  wellPanel(
                                           selectInput("inModelVer", "Model Version:",choices = c(),width=300),
                                          selectInput("inSampleToAdd", "Samples:",choices = c(),width=300),
                                          actionButton("inAddSample","Add"),
                                          textInput("inSamples", width="100%", "Samples To Load:"),
                                          textInput("inMinUmis", width=100, "Minimum #UMIs/cell:",value = 250),
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
                             #   plotOutput("BatchHeatmap",width="150%",height=600)
                          
                          )
                   ),
                  
                  tabPanel("Clusters", 
                    fluidRow(
                          selectInput("inModelOrAverage","View",choices=c("Model","Batch-corrected Average","Average")),
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
                              selectInput("inSelectGenesFrom","Select From:",choices=c("All genes","TFs","Surface markers")),
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
                            selectInput("inReorderingClustersMethod","Reordering method",choices=c("Hierarchical clustering","Diagonal","Formula")),
                            textInput("inOrderClusetersByGenes","Order by gene formula:",width=2000),
                            actionButton("inResetClusters", "Reset Clusters"),
                            actionButton("inSaveClusterOrder","Save Order"),
                            actionButton("inUndoClusterChanges","Undo")
                            ),
                            wellPanel(
                              selectInput("inAnnotateClusterNum","Cluster",choices=c()),
                            textInput("inClustAnnot", "New Annotation:"),
                            actionButton("inAnnotateCluster","Annotate"),
                            actionButton("inSaveAnnot", "Save Annotations"))
                            ),
                    column(12,
                     hr(),
                     textInput("inSamplesToShow", width=2000, "Samples:"),
                     actionButton("inResetSamples","Reset"),
                     hr()) 
                         )#,
#                    fluidRow(column(12,
#                            wellPanel(
#                            downloadButton("downloadExprsTable","Download Full Expression Table"),
#                            downloadButton("downloadCellCountsTable","Download Cell Counts Table")),
#                            wellPanel(
#                            actionButton("detectDiffExrs","screen Diff. expressed Genes"),
#                            downloadButton("downloadSignifExprsTable","Download Diffrential Expression Table"),
#                            selectInput("inputFDR",label="FDR=",choices=c(1e-1,5e-2,1e-2,1e-3,1e-4,1e-5),selected = 1e-2,width=100),
#                            selectInput("inputMinFC",label="Fold Change >",choices=c(1.5,2,4,8),selected = 2,width=100),
#                            selectInput("inMinAvgExprs",label="Avg. Exprs >",choices=c(1e-3,5e-3,1e-2,5e-2,1e-1),selected = 2,width=100))
#                            ))
                        ), 
                tabPanel("Truth", fluidRow(
                      plotOutput("truthplot",width="100%",height = "700"),
                      plotOutput("colorLegendTruth",height = 70,width=200),
                      plotOutput("truthSampleLegend",width=200,height=200),
                      sliderInput("inTruthColorScale","log2(1+#UMIs)",min = 0,max = 8,step = .1,value = c(0,2)),
                      selectInput("inTruthDownSamplingVersion",label="Down-sampling version:",choices = c(""),width=200),
                      checkboxInput("inTruthShowSeparatorBars",label = "Show Separator Bars",value = T),
                      selectInput("inTruthNcellsPerSample",label="#Cells Per Sample=",choices=params$nrandom_cells_per_sample_choices,selected = 1000,width=200)
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
                  textInput("inVarMean_varmeanThresh",  "Min Log2(var/mean):", value = "0.15"),
                  sliderInput("inVarMeanXlim",min=-6,max=2,step=.1,value = c(-1.5,2),label = "X-axis range")
                ),
                wellPanel(
                  selectInput("inNUmberOfGeneModules","Number of Modules:",c(10,20,50,100,200)),
                  actionButton("inGetModules","Get Modules")
                ),
                uiOutput("avg_module_plot"),
                       textInput("inModules", width=2000, "Modules:"),
                       wellPanel(
                         actionButton("inClusterModules", "Reorder Modules"),
                         selectInput("inReorderingModulesMethod","Reordering method",choices=c("Hierarchical clustering","Diagonal"))
                         ),
                      wellPanel(
                        selectInput("inModuleSelect","Show Module:",choices=c()),
                        textOutput("textModuleGenes")
                      )
              )
              ),
                tabPanel("Samples",fluidRow(
                          column(12,uiOutput("sample_avg_profile_plot")),
                          column(12,sliderInput("inSamplesColorScale","Log2(expression/mean)",min = -8,max = 8,step = 1,value = c(-4,4))),
                          column(6,textInput("inProjectSampleGroup1","Samples Group 1:")),
                          column(6,textInput("inProjectSampleGroup2","Samples Group 2:")),
                          column(12,selectInput("inProjectPlotType","Plot Type",choices=c("Side by Side","Tile")),
                          uiOutput("projection_barplot_plot"),
                   selectInput("inClustForDiffGeneExprsProjVsRef","Cluster for GE analysis:",choices=c()),
                   wellPanel(textInput("Gene1ForExprsTableRefVsProj","Gene:"),
                   tableOutput("Gene1ExprsTableRefVsProj")),
                   plotOutput("DiffGeneExprsProjVsRef",width="100%",height=500)))),
                     
                   fluidRow(column(12,br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br(),br())
              ),
          br()),
          img(src = 'sinai_logo.png',height = '70px', width = '70px'),
          p(em("(c) 2018 Ephraim (Effi) Kenigsberg. Department of Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai"))
      )
  
)

          
  

