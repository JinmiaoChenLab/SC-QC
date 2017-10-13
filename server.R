## max data size
options(shiny.maxRequestSize=1024^3)
options(shiny.launch.browser = T)
#options(shiny.host = "0.0.0.0")

shinyServer(function(input, output, session) {
    v <- reactiveValues(exprsData = NULL, nGene = NULL, spl_info = NULL, stats = NULL, meanQ = NULL,
                        user_group = NULL, user_groupby_level = NULL, selectedFeatures = NULL, summary_table = NULL)
    volumes <- getVolumes()
    shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))
    output$directorypath <- renderPrint({parseDirPath(volumes,input$directory)})
    output$directorypath_fastqcStats <- renderPrint({parseDirPath(volumes,input$directory)})
    
    
    ##-------------------Side Panel-------------------
    
    observeEvent(input$runFastqc, {
        path1 <- parseDirPath(volumes,input$directory)
        withProgress(message="Running fastqc...", value=0, {
            if(length(list.files(path = file.path(path1, "FASTQC"), pattern = '_fastqc.txt')>0)){
                cat('fastqc results are found!\n')
            }else{
                #system(paste0("/data1/home/mclau/pp/RNAseq_QCandNormalization/inhouse_QCwrapper/fastqc.sh -e2 ", path1))
                fastqc(fq.dir = path1)
                ## extract txt files with renaming
                qc_unzip(file.path(path1, "FASTQC"))
                all_folders <- list.dirs(path = file.path(path1,"FASTQC"), recursive=FALSE)
                #all_folders <- list.files(path = file.path(path1,"FASTQC"))
                #all_folders <- all_folders[-grep(x = all_folders, pattern = 'html')]
                for (i in all_folders){
                    #new_fn <- gsub(x=i,pattern = '.zip', replacement = '.txt')
                    new_fn <- paste0(i,'.txt')
                    file.rename(from = file.path(i,fsep=.Platform$file.sep,'fastqc_data.txt'), to =new_fn )
                }
            }
            
        })
    })
    
    
    observeEvent(input$summaryFastqcStat, {
        fastqcStatFiles <- input$FastqcStat_files
        print(fastqcStatFiles$datapath)
        
        withProgress(message="Summarizing fastqc...", value=0, {
            if(length(fastqcStatFiles$datapath)==0){
                cat('fastqc results are NOT found!\n')
            }else{
                res <- summarize_fastqc_stats(fastqc_files = fastqcStatFiles$datapath)
                stats <- res$summary
                stats$SampleID <- stats$Filename
                stats$Filename <- NULL
                meanQ <- res$meanQ
                meanQ$SampleID <- meanQ$Filename
                meanQ$Filename <- NULL
                v$stats <- stats
                v$meanQ <- meanQ
            }
        })
    })
    
    observeEvent(input$expr_file, {
        
        expr_file <- input$expr_file
        print(expr_file$datapath)
        
        exprsData <- read.table(file = expr_file$datapath, row.names=1, header = TRUE, as.is = TRUE, sep = '\t', check.names = FALSE)
        v$exprsData <- exprsData
        print(dim(v$exprsData))
        
        
    })
    
    observeEvent(input$sample_file, {
        sampleFile <- input$sample_file
        print(sampleFile$datapath)
        
        spl_info <- read.table(file = sampleFile$datapath, header = TRUE, as.is = TRUE, sep = '\t', check.names = FALSE)
        v$spl_info <- spl_info
        #print(v$spl_info)
        
    })
    
    
    output$comparison_group <- renderUI({
        
        gpNames <- colnames(v$spl_info)
        checkboxGroupInput('comparison_group', strong('Select comparison group(s):'),
                           gpNames, inline = TRUE)
    })
    
    observeEvent(input$combine_groups,{
        
        if(length(input$comparison_group)>1){
            
            comparison_group <- (input$comparison_group)
            spl <- (v$spl_info)
            spl[,paste(comparison_group, collapse = '_')] <- apply(spl[,comparison_group], MARGIN = 1, function(y){paste(y,collapse = '_')})
            # if(length(input$comparison_group)>1){
            #     v$user_group <- paste(comparison_group, collapse = '_')
            # }else{
            #     v$user_group <- input$comparison_group
            # }
            # 
            v$spl_info <- spl
        }else{
            # v$user_group <- input$comparison_group
        }
    })
    
    output$comparison_group_order <- renderUI({
        
        gp_levels <- unique(v$spl_info[,input$comparison_group])
        selectizeInput('comparison_group_order',
                       'Select order',
                       choices = gp_levels,
                       multiple = TRUE)
    })
    
    observeEvent(input$update_order,{
        
        user_groupby_level <- as.factor(input$comparison_group_order)
        v$spl_info[,input$comparison_group] <- factor(v$spl_info[,input$comparison_group], levels = user_groupby_level)
        output$placeholder <- renderText({ as.character(user_groupby_level) })
        v$user_groupby_level <- factor(user_groupby_level, levels = user_groupby_level) 
        v$user_group <- input$comparison_group
    })
    
    observeEvent(input$computeButton,{
        withProgress(message="Computing gene detection... ", value=0, {
            if(is.null(v$exprsData)){cat('exprs data is missing!\n'); return;}
            exprsData <- isolate(v$exprsData)
            if(is.null(v$spl_info)) {stop('Missing sample info!\n')}
            spl <- v$spl_info
            
            pairedEnd_samples <- spl$SampleID[grep(x = spl$SampleID, pattern = '_[1-2]')]
            if(length(pairedEnd_samples) >0){
                ss <- unlist(strsplit(x = pairedEnd_samples, split = '_'))
                sss <- ss[seq(1, length(ss), 2)]
                unique_pairedEnd_spl <- unique(sss)
                unique_id <- which(colnames(exprsData) %in% unique_pairedEnd_spl)
                
                end2_TPM <- exprsData[,unique_id]
                colnames(end2_TPM) <- paste0(colnames(end2_TPM),'_2')
                colnames(exprsData)[unique_id] <- paste0(colnames(exprsData)[unique_id],'_1')
                
                exprsData <- cbind(exprsData,end2_TPM)
                v$exprsData <- exprsData 
            }
            
            
            exprsData <- exprsData[colSums(exprsData)>0,]
            nGene <- colSums(exprsData > input$exprs_detectionLimit)
            print(head(nGene))
            v$nGene <- nGene
        })
    })   
    
    
    observeEvent(input$OpenDir_save, {
        resDir <- paste0(getwd(), .Platform$file.sep, "results_", Sys.Date())
        if(dir.exists(resDir)){
            opendir(resDir)
        }else{
            dir.create(resDir)
            opendir(resDir)
        }
    })
    
    observeEvent(input$saveButton, {
        if(!is.null(v$summary_table)){
            withProgress(message="Saving statistic summary table...", value=0, {
                print(getwd())
                
                #dir.create(paste0("results_", Sys.Date()))
                resDir <- paste0(getwd(), .Platform$file.sep, "results_", Sys.Date())
                if(dir.exists(resDir)){
                    opendir(resDir)
                }else{
                    dir.create(resDir)
                    opendir(resDir)
                }
                filename1 <- paste0(resDir, .Platform$file.sep, "stat_summary_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(resDir, .Platform$file.sep,
                                        "stat_summary_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".txt");
                    i = i + 1;
                }
                write.table(v$summary_table, file = filename1, row.names = FALSE, col.names = TRUE, sep = '\t')
            })
        }else{
            warning('no summary table found\n')
        }
    })
    
    ##---------------Total Reads tabset-------------------
    ####################
    ## Bar plot
    ####################    
    observeEvent(input$updateBarplot_totalReads, {
        
        V_TotalReadsBarPlotInput <- function(){
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            
            ## order by no. of reads
            df$Number_of_Reads_factor <- as.factor(df$Number_of_Reads)
            df <- df[order(df[,v$user_group],df$Number_of_Reads_factor),]
            
            
            withProgress(message="Generating Bar Plot...", value=0, {
                
                df <- within(df, SampleID <- factor(SampleID, levels=(SampleID)))
                gp <- ggplot(data=df, aes_string(x="SampleID", y="Number_of_Reads", fill=v$user_group)) +
                    labs(y = "Total no. of reads", x = 'Cells') +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    geom_bar(stat="identity")+ theme(#axis.text.x = element_text(size=15,angle=90,hjust=.5,vjust=.5,face="plain"),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                        axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                        axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                
                plot(gp)
            })
            
        }
        output$TotalReadsBarPlot <- renderPlot({
            V_TotalReadsBarPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFBarplot_totalReads, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                #dir.create(paste0("PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "totalReads_barplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "totalReads_barplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$R_tab1_w),
                    height=as.numeric(input$R_tab1_h))
                V_TotalReadsBarPlotInput()
                dev.off()
            })
            
        })
        observeEvent(input$OpenDirR, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
                #stop("PDF not created yet!")
            }
        })
        
    })    
    ####################
    ## Violin plot
    ####################
    observeEvent(input$updateVlnplot_totalReads, {
        V_totalReadsVlnPlotInput <- function(){
            
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            
            withProgress(message="Generating Violin Plot...", value=0, {
                dodge <- position_dodge(width = 0.4)
                gp <- ggplot(data = df, aes_string(x = v$user_group, y = "Number_of_Reads", fill = v$user_group)) +
                    geom_violin(position = dodge) +
                    geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
                    geom_point(position = position_jitter(width = 0.2)) +
                    labs(y = "Total no. of reads") +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    theme(axis.text.x = element_text(size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
                          axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                          axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                          axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                plot(gp)
            })
            
            
        }
        output$totalReadsVlnPlot <- renderPlot({
            V_totalReadsVlnPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFVlnplot_totalReads, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                #dir.create(paste0("PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "totalReads_vlnplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "totalReads_vlnplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$R_vln_tab1_w),
                    height=as.numeric(input$R_vln_tab1_h))
                V_totalReadsVlnPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirR_vln, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
                #stop("PDF not created yet!")
            }
        })
    })
    
    ##---------------Gene detection tabset-------------------
    ####################
    ## Bar plot
    ####################
    observeEvent(input$updateBarplot_nGene, {
        V_GeneNoBarPlotInput <- function(){
            g <- v$nGene
            spl_info <- v$spl_info
            data <- data.frame('SampleID'=names(g), 'nGene'=g)
            if(sum(!data$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in exprs file are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% data$SampleID) > 0){stop('Some cells in sample info are missing from exprs file\n')}
            
            df <- merge(data,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            
            ## order by no. of reads
            df$nGene_factor <- as.factor(df$nGene)
            df <- df[order(df[,v$user_group],df$nGene_factor),]
            
            withProgress(message="Generating Bar Plot...", value=0, {
                
                df <- within(df, SampleID <- factor(SampleID, levels=(SampleID)))
                gp <- ggplot(data=df, aes_string(x="SampleID", y="nGene", fill=v$user_group)) +
                    labs(y = "Total no. of genes", x = 'Cells') +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    geom_bar(stat="identity")+ theme(#axis.text.x = element_text(size=15,angle=90,hjust=.5,vjust=.5,face="plain"),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                        axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                        axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                
                plot(gp)
            })
            
            
        }
        output$GeneNoBarPlot <- renderPlot({
            V_GeneNoBarPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFBarplot_nGene, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                #dir.create(paste0("PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "nGene_barplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "nGene_barplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$G_tab1_w),
                    height=as.numeric(input$G_tab1_h))
                V_GeneNoBarPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirG, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
                #stop("PDF not created yet!")
            }
        })
    })
    ####################
    ## Violin plot
    ####################
    observeEvent(input$updateVlnplot_nGene, {
        V_nGeneVlnPlotInput <- function(){
            
            g <- v$nGene
            spl_info <- v$spl_info
            data <- data.frame('SampleID'=names(g), 'nGene'=g)
            if(sum(!data$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in exprs file are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% data$SampleID) > 0){stop('Some cells in sample info are missing from exprs file\n')}
            df <- merge(data,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            withProgress(message="Generating Violin Plot...", value=0, {
                dodge <- position_dodge(width = 0.4)
                
                gp <- ggplot(data = df, aes_string(x = v$user_group, y = "nGene", fill = v$user_group)) +
                    geom_violin(position = dodge) +
                    geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
                    geom_point(position = position_jitter(width = 0.2)) +
                    labs(y = "Total no. of genes") +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    theme(axis.text.x = element_text(size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
                          axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                          axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                          axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                plot(gp)
            })
            
        }
        output$GeneNoVlnPlot <- renderPlot({
            V_nGeneVlnPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFVlnplot_nGene, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                #dir.create(paste0("PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "nGene_vlnplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "nGene_vlnplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$G_vln_tab1_w),
                    height=as.numeric(input$G_vln_tab1_h))
                V_nGeneVlnPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirG_vln, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
    })
    
    ##---------------Deduplicated Percentage tabset-------------------
    
    ####################
    ## Bar plot
    ####################
    observeEvent(input$updateBarplot_dedup, {
        V_dedupBarPlotInput <- function(){
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            
            ## order by no. of reads
            df$dedup_perc_factor <- as.factor(df$dedup_perc)
            df <- df[order(df[,v$user_group],df$dedup_perc_factor),]
            
            withProgress(message="Generating Bar Plot...", value=0, {
                
                df <- within(df, SampleID <- factor(SampleID, levels=(SampleID)))
                gp <- ggplot(data=df, aes_string(x="SampleID", y="dedup_perc", fill=v$user_group)) +
                    labs(y = "Deduplicated  % (unique/total)", x = 'Cells') +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    geom_bar(stat="identity")+ theme(#axis.text.x = element_text(size=15,angle=90,hjust=.5,vjust=.5,face="plain"),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                        axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                        axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                plot(gp)
            })
            
            
        }
        output$DedupBarPlot <- renderPlot({
            V_dedupBarPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFBarplot_dedup, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                #dir.create(paste0("PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "deduplicatedPerc_barplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "deduplicatedPerc_barplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$D_tab1_w),
                    height=as.numeric(input$D_tab1_h))
                V_dedupBarPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirD, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
    })
    
    ####################
    ## Violin plot
    ####################
    observeEvent(input$updateVlnplot_dedup, {
        V_dedupVlnPlotInput <- function(){
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            withProgress(message="Generating Violin Plot...", value=0, {
                dodge <- position_dodge(width = 0.4)
                
                gp <- ggplot(data = df, aes_string(x = v$user_group, y = "dedup_perc", fill = v$user_group)) +
                    geom_violin(position = dodge) +
                    geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
                    geom_point(position = position_jitter(width = 0.2)) +
                    labs(y = "Deduplicated % (unique/total)") +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    theme(axis.text.x = element_text(size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
                          axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                          axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                          axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                
                plot(gp)
            })
            
            
        }
        output$DedupVlnPlot <- renderPlot({
            V_dedupVlnPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFVlnplot_dedup, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "deduplicatedPerc_vlnplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "deduplicatedPerc_vlnplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$D_vln_tab1_w),
                    height=as.numeric(input$D_vln_tab1_h))
                V_dedupVlnPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirD_vln, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
                #stop("PDF not created yet!")
            }
        })
    })
    
    
    ##---------------Quality features comparison tabset-------------------
    
    ####################
    ## line plot
    ####################
    output$featureFilter <- renderUI({
        stats <- v$stats
        g <- v$nGene; nGene <- data.frame('SampleID'=names(g), 'nGene'=g)
        spl_info <- v$spl_info
        
        # checks
        if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
        if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
        
        if(sum(!nGene$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in exprs file are missing from sample info\n')}
        if(sum(!spl_info$SampleID %in% nGene$SampleID) > 0){stop('Some cells in sample info are missing from exprs file\n')}
        
        
        df <- merge(stats,spl_info,by="SampleID")
        df <- merge(nGene,df,by="SampleID")
        df <- df[order(df[,v$user_group]),]
        new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
        v$summary_table <- new_summary
        features <- colnames(df)
        features <- features[order(features)]
        selectizeInput('features', 'Choose features:',
                       choices = features, selected = features,
                       multiple = TRUE, width = "100%")
        
        
    })
    output$ordering_feature <- renderUI({
        
        
        selectizeInput('key_feature', 'Choose ONE feature for ordering:',
                       choices = input$features, selected = input$features,
                       multiple = TRUE, width = "100%")
        
    })
    
    observeEvent(input$updateLinePlot_featureOverlay, {
        V_featureOverlayPlotInput <- function(){
            stats <- v$stats
            g <- v$nGene; nGene <- data.frame('SampleID'=names(g), 'nGene'=g)
            spl_info <- v$spl_info
            
            # checks
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            
            if(sum(!nGene$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in exprs file are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% nGene$SampleID) > 0){stop('Some cells in sample info are missing from exprs file\n')}
            
            
            df <- merge(stats,spl_info,by="SampleID")
            df <- merge(nGene,df,by="SampleID")
            df <- df[order(df[,input$key_feature]),]
            # df <- df[,c('SampleID',v$user_group,input$features)]
            
            ## ordering of features
            feature_defaultOrder <- c('Number_of_Reads','nGene','dedup_perc','Number_of_errors','Number_of_warnings','Length')
            feature_defaultOrder <- feature_defaultOrder[-which(feature_defaultOrder == input$key_feature)]
            feature_defaultOrder <- feature_defaultOrder[feature_defaultOrder %in% input$features]
            df <- df[,c('SampleID',v$user_group,input$key_feature,feature_defaultOrder)]
            
            
            # if(input$key_feature == 'nGene'){
            #     title_feature <- 'no of gene detected'    
            # }else if(input$key_feature == 'dedup_perc'){
            #     title_feature <- 'deduplicated percentage'    
            # }else {
            #     title_feature <- input$key_feature
            # }
            
            withProgress(message="Generating Line Plot...", value=0, {
                
                df <- within(df, SampleID <- factor(SampleID, levels=(SampleID)))
                df_long <- melt(df,id=c(v$user_group,"SampleID")) # convert to long format
                colnames(df_long)[which(colnames(df_long) == 'variable')] <- 'Quality_feature'
                #gp <- ggplot(df_long, aes_string(x="SampleID",y="value", colour="Quality_feature", shape= v$user_group)) +
                gp <- ggplot(df_long, aes_string(x="SampleID",y="value", colour= v$user_group)) +
                    geom_line(aes(group=Quality_feature))+ geom_point(size=5)+
                    labs(y = "Feature values", x = 'Cells') +
                    theme(#axis.text.x = element_text(size=15,angle=90,hjust=.5,vjust=.5,face="plain"),
                        axis.text.x = element_blank(),
                        axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                        axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                        axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
                        plot.title = element_text(color="blue", size=30, face="bold.italic"))+
                    facet_grid(Quality_feature ~ ., scales="free") +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) 
                #+ ggtitle(paste0("Ordered by ", title_feature))
                
                plot(gp)
            })
            
            
        }
        output$featureOverlay_linePlot <- renderPlot({
            V_featureOverlayPlotInput()
        }, height = 800)
        
        observeEvent(input$PDF_featureOverlay, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "featureOverlay_lineplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "featureOverlay_lineplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$O_tab1_w),
                    height=as.numeric(input$O_tab1_h))
                V_featureOverlayPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirO, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
    })
    
    
    ##---------------Per Base Quality tabset-------------------
    observeEvent(input$updateHeatmap_perBaseQ, {
        stats <- v$stats
        spl_info <- v$spl_info
        if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
        if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
        
        if(length(unique(stats$Length))==1){
            V_perBaseQPlotInput <- function(){
                
                meanQ <- v$meanQ
                df <- spl_info
                df <- df[order(df[,v$user_group]),]
                meanQ <- meanQ[match(df$SampleID, meanQ$SampleID),]
                
                if(!identical(as.character(df$SampleID), as.character(meanQ$SampleID)))stop("ERROR!")
                rownames(meanQ) <- meanQ$SampleID
                meanQ$SampleID <- NULL
                
                ## plot setups
                group_by_color <-rainbow(length(v$user_groupby_level))[as.integer(df[,v$user_group])]
                if(nrow(meanQ)<30) {sepcolor <- "black"
                }else{
                    sepcolor <- "white"
                }
                sepwidth <- c(0,min(0.5/nrow(meanQ),0.01))
                if(nrow(meanQ)<60){rowsep <- 1:nrow(meanQ)
                }else{
                    rowsep <- 0
                }
                ## legend
                group_by_legend <- paste(v$user_group,levels(df[,v$user_group]))
                
                dd <- t(apply(meanQ, MARGIN = 1, as.numeric))
                colnames(dd) <- colnames(meanQ)
                withProgress(message="Generating Bar Plot...", value=0, {
                    ## plot
                    h<-heatmap.2(as.matrix(dd), trace = 'none', scale = 'none',
                                 dendrogram = 'none', Rowv=FALSE,Colv=FALSE, cexCol = 1.2,
                                 margin=c(12,15),labRow="",
                                 RowSideColors = group_by_color,
                                 rowsep=rowsep, sepcolor = sepcolor, sepwidth=sepwidth)
                    legend(#"bottomright",
                            "bottomleft", 
                           legend = group_by_legend,
                           col = rainbow(length(v$user_groupby_level)),
                           lty= 1, horiz = FALSE,           
                           lwd = 8, cex=0.8, y.intersp=1.4	
                    )
                    h
                })
                
                
                
                
            }
            output$perBaseQPlot <- renderPlot({
                V_perBaseQPlotInput()
            }, height = 800)
            
            observeEvent(input$PDF_perBaseQ, {
                
                withProgress(message="Downloading plot PDF files...", value=0, {
                    print(getwd())
                    pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                    if(dir.exists(pdfDir)){
                        opendir(pdfDir)
                    }else{
                        dir.create(pdfDir)
                        opendir(pdfDir)
                    }
                    filename1 <- paste0(pdfDir, .Platform$file.sep, "perBaseQPlot_", Sys.Date(), ".pdf")
                    i = 0
                    while(file.exists(filename1)){
                        filename1 <- paste0(pdfDir, .Platform$file.sep,
                                            "perBaseQPlot_",
                                            Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                        i = i + 1;
                    }
                    pdf(filename1,
                        width=as.numeric(input$Q_tab1_w),
                        height=as.numeric(input$Q_tab1_h))
                    V_perBaseQPlotInput()
                    dev.off()
                })
                
            })
        }else{
            print('Different read lengths!!\n')
        }
        
        observeEvent(input$OpenDirQ, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
        
    })
    
    ##---------------FASTQC Results Summary tabset-------------------
    observeEvent(input$updateHeatmap_summary, {
        
        V_fastqcSummaryPlotInput <- function(){
            fastqc_terms <- read.table(file = 'fastqc_terms.txt', header = TRUE, sep = '\t', as.is = TRUE)
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            
            ## Error terms
            err_mat <- matrix(data = 0, nrow = nrow(stats), ncol = nrow(fastqc_terms))
            rownames(err_mat) <- stats$SampleID
            colnames(err_mat) <- fastqc_terms[,"names"]
            err_found <- strsplit(x=as.character(stats$Errors), split=', ')
            for(i in 1:nrow(err_mat)){
                err_mat[i,which(colnames(err_mat) %in% err_found[[i]])] <- 2
            }
            ## Warning terms
            warning_mat <- matrix(data = 0, nrow = nrow(stats), ncol = nrow(fastqc_terms))
            rownames(warning_mat) <- stats$SampleID
            colnames(warning_mat) <- fastqc_terms[,"names"]
            warning_found <- strsplit(x=as.character(stats$Warnings), split=', ')
            for(i in 1:nrow(warning_mat)){
                warning_mat[i,which(colnames(warning_mat) %in% warning_found[[i]])] <- 1
            }
            
            ## combine
            if(!identical(as.character(rownames(warning_mat)), as.character(rownames(err_mat))))stop("ERROR @error_warming_heatmap\n")
            m <- err_mat+warning_mat
            df <- spl_info
            df <- df[order(df[,v$user_group]),]
            m <- m[match(df$SampleID, rownames(m)),]
            
            ## plot setups
            group_by_color <-rainbow(length(v$user_groupby_level))[as.integer(df[,v$user_group])]
            err_warning_color <- rainbow(ncol(m))
            err_warning_legend_shortnms <- fastqc_terms[,"short_names"]
            
            if(nrow(m)<30) {sepcolor <- "black"
            }else{
                sepcolor <- "white"
            }
            sepwidth <- c(0,min(0.5/nrow(m),0.01))
            if(nrow(m)<60){rowsep <- 1:nrow(m)
            }else{
                rowsep <- 0
            }
            ## legend
            group_by_legend <- paste(v$user_group,levels(df[,v$user_group]))
            
            withProgress(message="Generating Bar Plot...", value=0, {
                
                ## plot
                h<-heatmap.2(as.matrix(m), trace = 'none', scale = 'none',
                             dendrogram = 'none', Rowv=FALSE,Colv=FALSE, cexCol = 1.8,
                             margin=c(16,12),labRow="",
                             RowSideColors = group_by_color,
                             ColSideColors = err_warning_color,  
                             breaks = seq(0,2,length.out = 4), col=c('green','pink', 'red'),
                             key.xlab = 'pass/warning/fail',
                             labCol = err_warning_legend_shortnms,#cexCol = 0.8 + 1/log10(nc),
                             #rowsep=sep_id, 
                             keysize=1.2,
                             rowsep=rowsep, sepcolor = sepcolor, sepwidth=sepwidth)
                legend(#"bottomright", 
                        "bottomleft",
                       legend = group_by_legend,
                       col = rainbow(length(v$user_groupby_level)),
                       lty= 1, horiz = FALSE,           
                       lwd = 8, cex=0.8, y.intersp=1.4	
                )
                
                h
                
            })
            
        }
        output$fastqcSummaryPlot <- renderPlot({
            V_fastqcSummaryPlotInput()
        }, height = 800)
        
        observeEvent(input$PDF_fastqcSummary, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                #dir.create(paste0("PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "fastqcSummaryPlot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "fastqcSummaryPlot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$S_tab1_w),
                    height=as.numeric(input$S_tab1_h))
                V_fastqcSummaryPlotInput()
                dev.off()
            })
            
        })
        observeEvent(input$OpenDirS, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
        
    })  
    
    
    
    
    ##---------------Total fastqc Error tabset-------------------
    
    ####################
    ## Bar plot
    ####################   
    observeEvent(input$updateBarplot_totalErrors, {
        V_TotalErrorBarPlotInput <- function(){
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            
            ## order by no. of reads
            df$Number_of_errors_factor <- as.factor(df$Number_of_errors)
            df <- df[order(df[,v$user_group],df$Number_of_errors_factor),]
            
            withProgress(message="Generating Bar Plot...", value=0, {
                
                df <- within(df, SampleID <- factor(SampleID, levels=(SampleID)))
                gp <- ggplot(data=df, aes_string(x="SampleID", y="Number_of_errors", fill=v$user_group)) +
                    labs(y = "Total no. of errors", x = 'Cells') +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    geom_bar(stat="identity")+ theme(#axis.text.x = element_text(size=15,angle=90,hjust=.5,vjust=.5,face="plain"),
                        axis.text.x=element_blank(),
                        axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                        axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                        axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                
                plot(gp)
            })
            
        }
        output$TotalErrorBarPlot <- renderPlot({
            V_TotalErrorBarPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFBarplot_totalError, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "totalError_barplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "totalError_barplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$E_tab1_w),
                    height=as.numeric(input$E_tab1_h))
                V_TotalErrorBarPlotInput()
                dev.off()
            })
            
        })
        observeEvent(input$OpenDirE, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
        
    })   
    
    ####################
    ## Violin plot
    ####################
    observeEvent(input$updateVlnplot_totalErrors, {
        V_totalErrorsVlnPlotInput <- function(){
            
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            withProgress(message="Generating Violin Plot...", value=0, {
                dodge <- position_dodge(width = 0.4)
                gp <- ggplot(data = df, aes_string(x = v$user_group, y = "Number_of_errors", fill = v$user_group)) +
                    labs(y = "Total no. of errors") +
                    geom_violin(position = dodge) +
                    geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
                    geom_point(position = position_jitter(width = 0.2)) +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    theme(axis.text.x = element_text(size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
                          axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                          axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                          axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                plot(gp)
            })
            
        }
        output$totalErrorsVlnPlot <- renderPlot({
            V_totalErrorsVlnPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFVlnplot_totalErrors, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "totalErrors_vlnplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "totalErrors_vlnplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$E_vln_tab1_w),
                    height=as.numeric(input$E_vln_tab1_h))
                V_totalErrorsVlnPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirE_vln, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
    })
    
    
    ##---------------Total fastqc warning tabset-------------------
    ####################
    ## Bar plot
    ####################    
    observeEvent(input$updateBarplot_totalWarnings, {
        
        V_TotalWarningBarPlotInput <- function(){
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            
            ## order by no. of reads
            df$Number_of_warnings_factor <- as.factor(df$Number_of_warnings)
            df <- df[order(df[,v$user_group],df$Number_of_warnings_factor),]
            
            withProgress(message="Generating Bar Plot...", value=0, {
                
                df <- within(df, SampleID <- factor(SampleID, levels=(SampleID)))
                gp <- ggplot(data=df, aes_string(x="SampleID", y="Number_of_warnings", fill=v$user_group)) +
                    labs(y = "Total no. of warnings", x = 'Cells') +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    geom_bar(stat="identity")+ theme(#axis.text.x = element_text(size=15,angle=90,hjust=.5,vjust=.5,face="plain"),
                        axis.text.x=element_blank(),
                        axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                        axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                        axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                
                plot(gp)
            })
            
        }
        output$TotalWarningBarPlot <- renderPlot({
            V_TotalWarningBarPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFBarplot_TotalWarning, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "totalWarning_barplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "TotalWarning_barplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$W_tab1_w),
                    height=as.numeric(input$W_tab1_h))
                V_TotalWarningBarPlotInput()
                dev.off()
            })
            
        })
        observeEvent(input$OpenDirW, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
        
    })   
    
    ####################
    ## Violin plot
    ####################
    observeEvent(input$updateVlnplot_totalWarnings, {
        V_totalWarningsVlnPlotInput <- function(){
            stats <- v$stats
            spl_info <- v$spl_info
            if(sum(!stats$SampleID %in% spl_info$SampleID) > 0){stop('Some cells in stats are missing from sample info\n')}
            if(sum(!spl_info$SampleID %in% stats$SampleID) > 0){stop('Some cells in sample info are missing from stats\n')}
            df <- merge(stats,spl_info,by="SampleID")
            df <- df[order(df[,v$user_group]),]
            new_summary <- save_statSummary(old=isolate(v$summary_table), new= df)
            v$summary_table <- new_summary
            withProgress(message="Generating Violin Plot...", value=0, {
                dodge <- position_dodge(width = 0.4)
                gp <- ggplot(data = df, aes_string(x = v$user_group, y = "Number_of_warnings", fill = v$user_group)) +
                    labs(y = "Total no. of warnings") +
                    geom_violin(position = dodge) +
                    geom_boxplot(width=.1, outlier.colour=NA, position = dodge) +
                    geom_point(position = position_jitter(width = 0.2)) +
                    theme(legend.text=element_text(size=12), legend.title=element_text(size=18)) +
                    theme(axis.text.x = element_text(size=15,angle=0,hjust=.5,vjust=.5,face="plain"),
                          axis.text.y = element_text(size=15,angle=0,hjust=1,vjust=0,face="plain"),  
                          axis.title.x = element_text(size=20,angle=0,hjust=.5,vjust=0,face="plain"),
                          axis.title.y = element_text(size=20,angle=90,hjust=.5,vjust=.5,face="plain"))
                plot(gp)
            })
            
        }
        output$totalWarningsVlnPlot <- renderPlot({
            V_totalWarningsVlnPlotInput()
        }, height = 800)
        
        observeEvent(input$PDFVlnplot_totalWarnings, {
            
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
                if(dir.exists(pdfDir)){
                    opendir(pdfDir)
                }else{
                    dir.create(pdfDir)
                    opendir(pdfDir)
                }
                filename1 <- paste0(pdfDir, .Platform$file.sep, "totalWarnings_vlnplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                    filename1 <- paste0(pdfDir, .Platform$file.sep,
                                        "totalWarnings_vlnplot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pdf(filename1,
                    width=as.numeric(input$W_vln_tab1_w),
                    height=as.numeric(input$W_vln_tab1_h))
                V_totalWarningsVlnPlotInput()
                dev.off()
            })
            
        })
        
        
        observeEvent(input$OpenDirW_vln, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
                opendir(pdfDir)
            }else{
                dir.create(pdfDir)
                opendir(pdfDir)
            }
        })
    })
})


