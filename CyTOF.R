# Goal: Analysis of CyTOF data

# packages to be installed/loaded
# install.packages("package name")
library(xlsx)

######################### Data Munging #########################
# load data
file_path <- file.path("/","Users","mshruti","Desktop","Stanford","udn627524","CyTOF")
if (getwd() != file_path)
{
  setwd(file_path)
}
cytof_data <- read.xlsx2("Ashley_05022016.xls",1,stringsAsFactors=F)

head(cytof_data)
dim(cytof_data)

# function to subset rows using pattern matching than exact match
subset_rows_with_pattern <- function(data_frame,pattern,col_name){
  pattern_indices <- c()
  for(i in seq_len(nrow(data_frame)))
  {
    if(grepl(pattern,data_frame[i,col_name]))
    {
      pattern_indices <- c(pattern_indices,i)
    }
  }
  return(data_frame[pattern_indices,])
}

# function to extract rows based on cell type, disease condition and stimulant 
getData <- function(data_frame,cell_type,disease_status,stimulant){
  # be careful with reg exporeesion when cell_type contains "+"
  if(disease_status=="healthy")
  {
    motif <- paste("2634_",".*","/",cell_type,"$",sep="")
  }
  else if(disease_status=="disease")
  {
    motif <- paste("UDN_627524_",".*","/",cell_type,"$",sep="")
  }
  else
  {
    print("Error: disease_status can be either healthy or disease")
  }
  
  stimulant_pattern <- paste("^",stimulant,"_",sep="")
  motif_stimulant_pattern <- paste(stimulant_pattern,motif,sep="")
  df <- subset_rows_with_pattern(data_frame,motif_stimulant_pattern,"Name")
  return(df)
}

### data extraction for 4 different blood cells
conditions <- c("healthy","disease")
cells <- c("Basophils","CD16 Hi Monocytes","CD16 Low Monocytes","Lymphocytes")
stimulants <- 1:8
df_name <- cytof_data

# extract data for different conditions
for(condition in conditions){
  for(cell in cells){
    for(stimulant in stimulants){
      assign( make.names(paste(condition,cell,stimulant,sep="_")), getData(df_name,cell,condition,stimulant) ) 
    }
  }
}

# function to extract cell count for each cell
extract_cellcount <- function(data_frame){
  return(as.numeric(data_frame[1,"X.Cells"]))
}

###  data extraction for lymphocytes 

# lymphocytes is a vector containing names of different lymphocytes as named in the given data
lymphocytes_name <- c("CD3+/CD4+/Central Memory","CD3+/CD4+/Effector","CD3+/CD4+/Effector Memory",
                      "CD3+/CD4+/HLA-DR+CD38+","CD3+/CD4+/Naive","CD3+/CD4+/Tregs","CD3+/CD4+CD8+",
                      "CD3+/CD4-CD8-","CD3+/CD8+/Central Memory","CD3+/CD8+/Effector","CD3+/CD8+/Effector Memory",
                      "CD3+/CD8+/HLA-DR+CD38+","CD3+/CD8+/Naive","CD3-/B cells","CD3-/B cells/IgD+ Memory",
                      "CD3-/B cells/IgD-CD27-","CD3-/B cells/Naive","CD3-/B cells/Switched Memory",
                      "CD3-/B cells/Transitional","CD3-/CD19+","CD3-/CD19+/Plasmoblast","CD3-/Non B/DC",
                      "CD3-/Non B/DC/mDC","CD3-/Non B/DC/pDC","CD3-/Non B/NK","CD3-/Non B/NK/CD16 Hi NK",
                      "CD3-/Non B/NK/CD16 Low NK","CD3-/Non B/NK/HLA-DR+NK","NKT")

# substitute multiple different patterns by different characters in a string in a single command.
findReplace <- function(find,replace,string){
  for(i in seq_along(find))
  {
    string <- gsub(find[i], replace[i], string)
  }
  return(string)
}

search_for <- c("\\s+", "[\\/]", "\\+", "\\-")
replace_by <- c("", "_", "positive", "negative")

# extract data for different conditions
for(condition in conditions){
  for(cell in lymphocytes_name){
    for(stimulant in stimulants){
      # how the data should be named
      cell_name <- findReplace(find=search_for, replace=replace_by, cell)
      
      # pattern to search for this cell type in the data frame (cytof_data) 
      cell_pattern <- findReplace(find=c("\\+","\\-","\\/"),replace=c("\\\\+","\\\\-","\\\\/"),cell)
      
      assign( paste(condition, cell_name, stimulant,sep="_"),
              getData(df_name,cell_pattern,condition,stimulant) ) 
    }
  }
}

### data extraction for protein levels in various cell types, under different conditions conditions and stimulants.
proteins <- c("pERK1/2","IkB","pAKT","p-P38","pSTAT3","pS6","pSTAT1","pSTAT5","pPLCg2")

getProteinStimulantData <- function(data_frame,cell_type){
  protein_stimulant_df <- data.frame(matrix(data=NA, ncol=5, dimnames= list(c(), c('cell_type', 'protein', 'stimulant', 'healthy', 'disease'))))
  
  healthy_list <- list()
  for(stimulant in stimulants) 
  {
    healthy_list[[stimulant]] <- getData(data_frame,cell_type,"healthy",stimulant)
  }
  
  disease_list <- list()
  for(stimulant in stimulants) 
  {
    disease_list[[stimulant]] <- getData(data_frame,cell_type,"disease",stimulant)
  }
  
  count <- 1 #row number for final dataframe
  for(protein in proteins) 
  {
    for(stimulant in stimulants) 
    {
      healthy_stimulated <- as.numeric(subset(healthy_list[[stimulant]],Well.ID.==protein,Statistic))
      disease_stimulated <- as.numeric(subset(disease_list[[stimulant]],Well.ID.==protein,Statistic))
      
      protein_stimulant_df[count,'cell_type'] = cell_type
      protein_stimulant_df[count,'protein'] = protein
      protein_stimulant_df[count,'stimulant'] = stimulant
      protein_stimulant_df[count,'healthy'] = healthy_stimulated
      protein_stimulant_df[count,'disease'] = disease_stimulated
      
      count <- count + 1
    }
  }
  return(protein_stimulant_df) 
}

######################### Exploratory Analysis #########################

# Compare cell count Healthy vs Disease in unstimulated condition
cellcounts <- data.frame()
for(cell in cells){
  for(condition in conditions){
    pattern <- make.names( paste(condition,cell,1,sep="_") ) #1 for unstimulated condition
    # get the desired data frame
    df <- eval( parse(text=pattern) )
    # row is cell-type, column is disease condition in cellcounts data frame
    cellcounts[make.names(cell),condition] <- extract_cellcount(df)
  }
}
cellcounts

# plot Change in Cell Count in Disease in Unstimulated condition
par(mfrow=c(2,2))
row_id <- 0
apply(cellcounts,1,function(i) {
  # "<<-" : search to made through parent environments for an existing definition of the variable being assigned.
  (row_id <<- row_id+1); 
  barplot((as.matrix(i)), ylab="cell count", main=rownames(cellcounts[row_id,]), names.arg=colnames(cellcounts), beside=TRUE, col=c("lightblue","brown"), space=c(0.2,0) )})

# plot Change in Cell Count as log2 ratio
logratio_cellcounts <- transform(cellcounts, logratio=log2(disease/healthy))
par(mfrow=c(1,1))
barplot(logratio_cellcounts$logratio, beside=T, col=2:5, legend.text=rownames(logratio_cellcounts), args.legend=list(x="bottomright"), ylab="Log2 Fold Change Disease VS Healthy")

### Compare Lympocyte cell count Healthy vs Disease in unstimulated condition
lymphocytes_cellcounts <- data.frame()

for(cell in lymphocytes_name){
  for(condition in conditions){
    cell_name <- findReplace(search_for, replace_by, cell)
    pattern <- paste(condition,cell_name,1,sep="_") #1 for unstimulated condition
    # get the desired data frame
    df <- eval( parse(text=pattern) )
    # row is cell-type, column is disease condition in cellcounts data frame
    lymphocytes_cellcounts[cell_name,condition] <- extract_cellcount(df)
  }
}

lymphocytes_cellcounts

# plot Change in Lympocyte cells' in Disease in Unstimulated condition
par(mfrow=c(2,2))
row_id <- 0
apply(lymphocytes_cellcounts,1,function(i) {
# "<<-" : search to made through parent environments for an existing definition of the variable being assigned.
(row_id <<- row_id+1);
barplot((as.matrix(i)), ylab="cell count", main=rownames(lymphocytes_cellcounts[row_id,]), names.arg=colnames(lymphocytes_cellcounts), beside=TRUE, col=c("lightblue","brown"), space=c(0.2,0) )})

# plot Change in Cell Count as log2 ratio
logratio_lymphocytes_cellcounts <- transform(lymphocytes_cellcounts, logratio=log2(disease/healthy))
par(mfrow=c(1,1))
barplot(logratio_lymphocytes_cellcounts$logratio, beside=T, col=1:nrow(logratio_lymphocytes_cellcounts), ylab="Log2 Fold Change Disease VS Healthy")
#legend(x="bottomright",legend=rownames(logratio_lymphocytes_cellcounts), col=1:nrow(lymphocytes_cellcounts),fill=TRUE)

# since there are too many lymphocyte cells, select the ones which show atleast 2 fold change in disease
logratio2_lymphocytes_cellcounts <- subset(logratio_lymphocytes_cellcounts,abs(logratio)>1)
barplot(logratio2_lymphocytes_cellcounts$logratio, beside=T, col=1:nrow(logratio2_lymphocytes_cellcounts), ylab="Log2 Fold Change Disease VS Healthy",legend.text=rownames(logratio2_lymphocytes_cellcounts), args.legend=list(x="bottom", cex=0.7) )

######## Change in Phosphorylated Proteins Levels in Disease

# change in phosphorylated protein Levels in Disease, stimulator=1 implies unstimulated
# TODO: rename the function in a better way
proteinStimulantPlot <- function(data_frame,cell_type,stimulator,log=F){
    # log True: implies if you want to plot healthy and disease sideways; False implies to plot log2 fold change in disease/healthy
    df <- getProteinStimulantData(data_frame,cell_type)
    df2 <- subset(df,stimulant==stimulator) 
    df3 <- as.matrix( df2[,c("healthy","disease")] )
    rownames(df3) <- df2$protein
    
    row_id <- 0

    if(log==T)
    {
      logratio_df3 <- transform(df3, logratio=log2(disease/healthy))
      par(mfrow=c(1,1))
      apply(logratio_df3,1,function(i) 
        {
        # "<<-" : search to made through parent environments for an existing definition of the variable being assigned.
        (row_id <<- row_id+1); 
        barplot(logratio_df3$logratio, beside=T, col=2:9, legend.text=rownames(logratio_df3), 
                args.legend=list(x="bottom", cex=0.9, box.col="white"), ylab="Log2 Fold Change Disease VS Healthy", 
                main=c(cell_type,paste("stimulant:",stimulator)) )
        })
    }
    
    else
    {
      par(mfrow=c(2,2))
      apply(df3,1,function(i) 
      {
        # "<<-" : search to made through parent environments for an existing definition of the variable being assigned.
        (row_id <<- row_id+1); 
        barplot((as.matrix(i)), ylab=rownames(df3)[row_id], main=c(cell_type,paste("stimulant:",stimulator)), 
                names.arg=colnames(df3), beside=TRUE, col=c("lightblue","brown"), space=c(0.2,0) )
      })
    }
    #return(df3)
}

# Change in Phosphorylated Proteins Levels in Disease: Basophils (Unstimulated)
# sapply(cells[1:3], function(i) {proteinStimulantPlot(cytof_data,i,1,log=F)})
sapply(cells[1:3], function(i) {proteinStimulantPlot(cytof_data,i,1,log=T)})

# relative change in Phosphorylated Proteins Levels in Disease upon stimulation
proteinDifferentStimulantsPlot <- function(data_frame,cell_type,log=F){
    protein_stimulation_ratios_df <- data.frame(matrix(data=NA, ncol=5, dimnames= list(c(), c('cell_type', 'protein', 'stimulant', 'stimulation_ratio_healthy', 'stimulation_ratio_disease'))))
  
    df <- getProteinStimulantData(data_frame,cell_type)
    count <- 1
    for(prot in proteins)
    {
      healthy_unstimulated <- subset(df,protein==prot & stimulant==1,healthy,drop=T)
      disease_unstimulated <- subset(df,protein==prot & stimulant==1,disease,drop=T)
      # print(healthy_unstimulated)
      # print(disease_unstimulated)
      
      for(stimulator in stimulants[2:length(stimulants)])
      {
        healthy_stimulated <- subset(df,protein==prot & stimulant==stimulator,healthy,drop=T)
        disease_stimulated <- subset(df,protein==prot & stimulant==stimulator,disease,drop=T)

        # difference in protein levels upon stimulation
        healthy_diff = healthy_stimulated - healthy_unstimulated
        disease_diff = disease_stimulated - disease_unstimulated
        # relative change in protein levels upon stimulation
        healthy_ratio = healthy_diff/healthy_unstimulated
        disease_ratio = disease_diff/disease_unstimulated

        protein_stimulation_ratios_df[count,'cell_type'] = cell_type
        protein_stimulation_ratios_df[count,'protein'] = prot
        protein_stimulation_ratios_df[count,'stimulant'] = stimulator
        protein_stimulation_ratios_df[count,'stimulation_ratio_healthy'] = healthy_ratio
        protein_stimulation_ratios_df[count,'stimulation_ratio_disease'] = disease_ratio

        #print(c(prot,stimulator,count,healthy_ratio,disease_ratio))
        count <- count + 1
      }
    }
    return(protein_stimulation_ratios_df)
}



# trouble shoot
protein_stimulation_basophils <- proteinDifferentStimulantsPlot(cytof_data,'Basophils')
# ToDO: check if u got right values from above


x <- protein_stimulation_basophils
x2 <- subset(x,protein=="pSTAT5")
x3 <- as.matrix(x2[,c('stimulation_ratio_healthy','stimulation_ratio_disease')])
rownames(x3) <- x2$stimulant
