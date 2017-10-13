## =====================================================================
## function for updating statistic table
## =====================================================================
save_statSummary <- function(old =  NULL, new = NULL){
    if(is.null(old)) return(new)
    
    if(sum(colnames(old) %in% colnames(new)) < ncol(old)){
        return(merge(old, new))
    }
        
}

## =====================================================================
## function for opening the results directory
## =====================================================================
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}

## =====================================================================
## function for extracting quality results from fastqc txt files
## =====================================================================
summarize_fastqc_stats <- function(fastqc_files = NULL){
  
  f <- fastqc_files
  
  summary_table <-  data.frame('Filename'=rep(NA,length(f)),'Number_of_Reads'=rep(0,length(f)),
                               'Length'= rep(0,length(f)),'Number_of_warnings'=rep(0,length(f)),
                               'Warnings'=rep(NA,length(f)),'Number_of_errors'= rep(0,length(f)),
                               'Errors'= rep(NA,length(f)), 'dedup_perc'=rep(0,length(f)))
  mean_table <- data.frame('Filename'=rep(NA,length(f)))
  diff_read_length <- FALSE
  
  for (i in 1:length(f))
  {
    
    d1 <- readLines(f[i])
    line_END_MODULEs <- grep(x=d1, pattern = '>>END_MODULE')
    
    ################
    ## filename
    ################
    Filename_line <- grep(x=d1, pattern = 'Filename')
    t <- unlist(strsplit(x=d1[Filename_line], split = '\t'))
    Fname <- gsub(x = t[2], pattern = '.fastq.*', replacement = "")
    summary_table[i,'Filename'] <- Fname
    
    
    ################
    ## Summary of pass/fail & basic stats
    ################
    lines_all_measurements <- grep(x=d1, pattern = '>>')
    lines_all_measurements<-d1[setdiff(lines_all_measurements, line_END_MODULEs)]
    t <- unlist(strsplit(x=lines_all_measurements, split = '\t'))
    status <- matrix(t, nrow = length(lines_all_measurements), ncol = 2, byrow = TRUE)
    status <- gsub(x = status, pattern = '>>', replacement = "")
    status <- as.data.frame(status)
    # remove 'Basic Statistics'
    status <- status[-grep(x=status[,1], pattern = 'Basic Statistics'),]
    t <- as.data.frame(table(status[,2])); rownames(t) <- t[,1]
    summary_table[i,'Number_of_errors'] <- t['fail','Freq']
    summary_table[i,'Number_of_warnings'] <- t['warn','Freq']
    summary_table[i,'Errors'] <- paste(status[status[,2] == 'fail',1], collapse = ', ')
    summary_table[i,'Warnings'] <- paste(status[status[,2] == 'warn',1], collapse = ', ')
    
    
    # 'Basic Statistics'
    startline_BasicStats <- grep(x=d1, pattern = '>>Basic Statistics')
    temp_line_END_MODULEs <- line_END_MODULEs[line_END_MODULEs > startline_BasicStats]
    endline_BasicStats <- temp_line_END_MODULEs[which.min(temp_line_END_MODULEs - startline_BasicStats)]
    
    BasicStats <- d1[c((startline_BasicStats+1):(endline_BasicStats-1))]
    t <- unlist(strsplit(x=BasicStats, split = '\t'))
    BasicStats <- matrix(t, nrow = length(t)/2, ncol = 2, byrow = TRUE)
    
    summary_table[i,'Number_of_Reads'] <- BasicStats[grep(x=BasicStats[,1], pattern = 'Total Sequences'),2]
    summary_table[i,'Length'] <- BasicStats[grep(x=BasicStats[,1], pattern = 'Sequence length'),2]
    
    ################
    ## deduplicated percentage
    ################
    dedup_line <- grep(x=d1, pattern = '#Total Deduplicated Percentage')
    t <- unlist(strsplit(x=d1[dedup_line], split = '\t'))
    dedup_perc <- round(as.numeric(t[2]), digits=2)
    summary_table[i,'dedup_perc'] <- dedup_perc
    
    ################
    ## per base quality
    ################
    if(diff_read_length){
      cat('Different read length detected!\n')
      
    }else{
      # get data chunck
      startline_perbaseQ <- grep(x=d1, pattern = '>>Per base sequence quality')
      temp_line_END_MODULEs <- line_END_MODULEs[line_END_MODULEs > startline_perbaseQ]
      endline_perbaseQ <- temp_line_END_MODULEs[which.min(temp_line_END_MODULEs - startline_perbaseQ)]
      
      perbaseQ <- d1[c((startline_perbaseQ+1):(endline_perbaseQ-1))]
      t <- unlist(strsplit(x=perbaseQ, split = '\t'))
      perbaseQ <- matrix(t, nrow = length(t)/7, ncol = 7, byrow = TRUE)
      mean_vec <- perbaseQ[2:nrow(perbaseQ),2]
      
      if(i==1){
        baseCat <- perbaseQ[2:nrow(perbaseQ),1]
        empty_mean <- as.data.frame(matrix(-1, nrow = nrow(mean_table), ncol = length(mean_vec)))
        colnames(empty_mean) <- baseCat
        mean_table <- cbind(mean_table, empty_mean)   
      }
      
      mean_table[i,'Filename'] <- Fname
      if(length(mean_vec) != ncol(mean_table)-1){
          diff_read_length<-TRUE
          mean_table[i, 2:ncol(mean_table)] <- NA  
      }else{
          mean_table[i, 2:ncol(mean_table)] <- mean_vec            
      }
      
    }
    
    
  }
  
  summary_table$Number_of_Reads <- as.numeric(summary_table$Number_of_Reads)
  output<-list('summary'=summary_table, 'meanQ'=mean_table)
  return(output)
}

