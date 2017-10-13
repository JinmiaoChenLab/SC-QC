library(shiny)
library(shinyFiles)
library(gplots)
library(ggplot2)
library(reshape2)
library(fastqcr)

shinyUI(fluidPage(
    #titlePanel("Quality examination of RNAseq Data"),
    titlePanel("scSeqQC"),
    hr(),
    
    fluidRow(
        ## =========================================================================================================##
        ## Input
        ## =========================================================================================================##
        
        column(3,
               h2('fastqc'),
               wellPanel(
                   h4("Step 1: run fastqc (.fastq/ .fastq.gz)"),
                   shinyDirButton('directory', 'Folder select', 'Please select a folder'),
                   verbatimTextOutput('directorypath'),
                   
                   fluidRow(
                       column(6,
                              actionButton("runFastqc", "run", icon = icon("hand-o-right"))
                       )
                   ),
                   hr(),
                   h4("Step 2: summarize fastqc output (.txt)"),
                   fileInput(inputId = 'FastqcStat_files',
                             label = "",
                             multiple = TRUE,
                             accept = c(".txt")),
                   actionButton("summaryFastqcStat", "summarize", icon = icon("hand-o-right"))
                   
               ),
               hr(),
               h2("sample meta data"),
               wellPanel(
                   fileInput(inputId = 'sample_file',
                             label = "Sample meta data file(txt)",
                             multiple = FALSE,
                             accept = c(".txt")),
                   uiOutput("comparison_group"),
                   actionButton('combine_groups',"Create new group"),
                   uiOutput("comparison_group_order"),
                   actionButton('update_order',"Update"),
                   
                   h6("User input order"),
                   verbatimTextOutput("placeholder", placeholder = TRUE)
               ),
               hr(),
               h2("expression matrix: gene detection rate"),
               wellPanel(
                   fileInput(inputId = 'expr_file',
                             label = "exprs matrix gene-by-sample (txt)",
                             multiple = FALSE,
                             accept = c(".txt")),
                   
                   textInput("exprs_detectionLimit", "Expression detection limit", "1"),
                   #verbatimTextOutput("exprs_detectionLimit")
                   actionButton('computeButton',"Compute")
                   
                   
               ),
               hr(),
               h2("Statistic results saving"),
               wellPanel(
                   actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),
                   actionButton("OpenDir_save", "Open results folder", icon = icon("folder"))
                   ),
               
               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ",
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        
        ## =========================================================================================================##
        ## Visualization
        ## =========================================================================================================##
        column(9,
               tabsetPanel(type = "pills",
                           tabPanel("Total reads", fluidPage(
                               hr(),
                               tabsetPanel(id="D_totalReads",
                                           tabPanel(title="Barplot", value="B_panel1",
                                                    
                                                    hr(),
                                                    actionButton("updateBarplot_totalReads", "Update Plot", icon = icon("hand-pointer-o")),
                                                    
                                                    fluidRow(
                                                        column(12,plotOutput("TotalReadsBarPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFBarplot_totalReads", "Download Barplot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="R_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="R_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirR", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           ),
                                           tabPanel(title="Violinplot", value="V_panel1",
                                                    hr(),
                                                    actionButton("updateVlnplot_totalReads", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("totalReadsVlnPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFVlnplot_totalReads", "Download Violin plot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="R_vln_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="R_vln_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirR_vln", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           )
                                           # tabPanel(title="fdaNorm", value="N_panel2")
                               )
                           )),
                           tabPanel("Gene detection", fluidPage(
                               hr(),
                               tabsetPanel(id="D_geneDetection",
                                           
                                           tabPanel(title="Barplot", value="B_panel2",
                                                    # br(),
                                                    # uiOutput("markerFilter"),
                                                    hr(),
                                                    actionButton("updateBarplot_nGene", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("GeneNoBarPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFBarplot_nGene", "Download Barplot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="G_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="G_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirG", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           ),
                                           tabPanel(title="Violinplot", value="V_panel2",
                                                    hr(),
                                                    actionButton("updateVlnplot_nGene", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("GeneNoVlnPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFVlnplot_nGene", "Download Violin plot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="G_vln_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="G_vln_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirG_vln", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           )
                               )
                               
                           )),
                           
                           tabPanel("Deduplicated %", fluidPage(
                               hr(),
                               tabsetPanel(id="D_dedupPerc",
                                           
                                           tabPanel(title="Barplot", value="B_panel3",
                                                    # br(),
                                                    # uiOutput("markerFilter"),
                                                    hr(),
                                                    actionButton("updateBarplot_dedup", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("DedupBarPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFBarplot_dedup", "Download Barplot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="D_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="D_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirD", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           ),
                                           tabPanel(title="Violinplot", value="V_panel3",
                                                    hr(),
                                                    actionButton("updateVlnplot_dedup", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("DedupVlnPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFVlnplot_dedup", "Download Violin plot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="D_vln_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="D_vln_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirD_vln", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           )
                               )
                               
                           )),
                           tabPanel("Per base Quality", fluidPage(
                               hr(),
                               actionButton("updateHeatmap_perBaseQ", "Update Plot", icon = icon("hand-pointer-o")),
                               
                               fluidRow(
                                   column(12,plotOutput("perBaseQPlot", height = 800)    ),
                                   column(12, 
                                          h3("Export to PDF:"),
                                          wellPanel(
                                              actionButton("PDF_perBaseQ", "Download Barplot in PDF", icon = icon("download")),
                                              
                                              hr(),
                                              fluidRow(
                                                  column(6,
                                                         sliderInput(inputId="Q_tab1_w", label = "PDF width(in):", 
                                                                     min=3, max=30, value=20, width=100, ticks=F)
                                                  ),
                                                  column(6, 
                                                         sliderInput(inputId="Q_tab1_h", label = "PDF height(in):", 
                                                                     min=3, max=20, value=8, width=100, ticks=F)
                                                  )),
                                              
                                              actionButton("OpenDirQ", "Open download folder", icon = icon("folder"))
                                          )
                                   )
                                   
                               )
                               
                           )),
                           tabPanel("fastqc results summary", fluidPage(
                               hr(),
                               actionButton("updateHeatmap_summary", "Update Plot", icon = icon("hand-pointer-o")),
                               
                               fluidRow(
                                   column(12,plotOutput("fastqcSummaryPlot", height = 800)    ),
                                   column(12, 
                                          h3("Export to PDF:"),
                                          wellPanel(
                                              actionButton("PDF_fastqcSummary", "Download Barplot in PDF", icon = icon("download")),
                                              
                                              hr(),
                                              fluidRow(
                                                  column(6,
                                                         sliderInput(inputId="S_tab1_w", label = "PDF width(in):", 
                                                                     min=3, max=30, value=20, width=100, ticks=F)
                                                  ),
                                                  column(6, 
                                                         sliderInput(inputId="S_tab1_h", label = "PDF height(in):", 
                                                                     min=3, max=20, value=8, width=100, ticks=F)
                                                  )),
                                              
                                              actionButton("OpenDirS", "Open download folder", icon = icon("folder"))
                                          )
                                   )
                                   
                               )
                               
                           )),

                           tabPanel("Total fastqc errors", fluidPage(
                               hr(),
                               tabsetPanel(id="D_totalErrors",
                                           tabPanel(title="Barplot", value="B_panel4",
                                                    
                                                    hr(),
                                                    actionButton("updateBarplot_totalErrors", "Update Plot", icon = icon("hand-pointer-o")),
                                                    
                                                    fluidRow(
                                                        column(12,plotOutput("TotalErrorBarPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFBarplot_totalError", "Download Barplot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="E_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="E_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirE", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           ),
                                           tabPanel(title="Violinplot", value="V_panel4",
                                                    hr(),
                                                    actionButton("updateVlnplot_totalErrors", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("totalErrorsVlnPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFVlnplot_totalErrors", "Download Violin plot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="E_vln_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="E_vln_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirE_vln", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           )
                               )
                           )),
                           tabPanel("Total fastqc warnings", fluidPage(
                               hr(),
                               tabsetPanel(id="D_totalWarnings",
                                           tabPanel(title="Barplot", value="B_panel5",
                                                    
                                                    hr(),
                                                    actionButton("updateBarplot_totalWarnings", "Update Plot", icon = icon("hand-pointer-o")),
                                                    
                                                    fluidRow(
                                                        column(12,plotOutput("TotalWarningBarPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFBarplot_TotalWarning", "Download Barplot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="W_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="W_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirW", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           ),
                                           tabPanel(title="Violinplot", value="V_panel5",
                                                    hr(),
                                                    actionButton("updateVlnplot_totalWarnings", "Update Plot", icon = icon("hand-pointer-o")),
                                                    fluidRow(
                                                        column(12,plotOutput("totalWarningsVlnPlot", height = 800)    ),
                                                        column(12, 
                                                               h3("Export to PDF:"),
                                                               wellPanel(
                                                                   actionButton("PDFVlnplot_totalWarnings", "Download Violin plot in PDF", icon = icon("download")),
                                                                   
                                                                   hr(),
                                                                   fluidRow(
                                                                       column(6,
                                                                              sliderInput(inputId="W_vln_tab1_w", label = "PDF width(in):", 
                                                                                          min=3, max=30, value=20, width=100, ticks=F)
                                                                       ),
                                                                       column(6, 
                                                                              sliderInput(inputId="W_vln_tab1_h", label = "PDF height(in):", 
                                                                                          min=3, max=20, value=8, width=100, ticks=F)
                                                                       )),
                                                                   
                                                                   actionButton("OpenDirW_vln", "Open download folder", icon = icon("folder"))
                                                               )
                                                        )
                                                        
                                                    )
                                           )
                               )
                           )),
                           tabPanel("Quality features comparison", fluidPage(
                               br(),
                               uiOutput("featureFilter"),
                               checkboxInput("featureSelection", "features chosen", TRUE),
                               hr(),
                               
                               conditionalPanel(
                                   condition = "input.featureSelection",
                                   uiOutput("ordering_feature")
                               ),
                               hr(),
                               
                            actionButton("updateLinePlot_featureOverlay", "Update Plot", icon = icon("hand-pointer-o")),
                               
                               fluidRow(
                                   column(9,plotOutput("featureOverlay_linePlot", height = 800)    ),
                                   column(3, 
                                          h3("Export to PDF:"),
                                          wellPanel(
                                              actionButton("PDF_featureOverlay", "Download line plot in PDF", icon = icon("download")),
                                              
                                              hr(),
                                              fluidRow(
                                                  column(6,
                                                         sliderInput(inputId="O_tab1_w", label = "PDF width(in):", 
                                                                     min=3, max=30, value=20, width=100, ticks=F)
                                                  ),
                                                  column(6, 
                                                         sliderInput(inputId="O_tab1_h", label = "PDF height(in):", 
                                                                     min=3, max=20, value=8, width=100, ticks=F)
                                                  )),
                                              
                                              actionButton("OpenDirO", "Open download folder", icon = icon("folder"))
                                          )
                                   )
                                   
                               )
                               
                           ))
               )
        )
    )
))