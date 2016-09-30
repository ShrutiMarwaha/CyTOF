library(xlsx)

# load data
file_path <- file.path("/","Users","shruti","Desktop","Stanford","udn627524","CyTOF")
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
  # hard coding for unstimulated <- "^1_" 
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
  stimulant_motif_pattern <- paste(stimulant_pattern,motif,sep="")
  df <- subset_rows_with_pattern(data_frame,stimulant_motif_pattern,"Name")
  return(df)
}

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

# Compare cell count Healthy vs Disease in unstimulated condition
cellcounts <- data.frame(matrix(data=NA,nrow=4,ncol=2,dimnames=list(make.names(cells),conditions) ))
cellcounts["Basophils","healthy"] <- extract_cellcount(healthy_Basophils_1)
cellcounts["Basophils","disease"] <- extract_cellcount(disease_Basophils_1)
cellcounts["CD16.Hi.Monocytes","healthy"] <- extract_cellcount(healthy_CD16.Hi.Monocytes_1) 
cellcounts["CD16.Hi.Monocytes","disease"] <- extract_cellcount(disease_CD16.Hi.Monocytes_1)
cellcounts["CD16.Low.Monocytes","healthy"] <- extract_cellcount(healthy_CD16.Low.Monocytes_1) 
cellcounts["CD16.Low.Monocytes","disease"] <- extract_cellcount(disease_CD16.Low.Monocytes_1)
cellcounts["Lymphocytes","healthy"] <- extract_cellcount(healthy_Lymphocytes_1)
cellcounts["Lymphocytes","disease"] <- extract_cellcount(disease_Lymphocytes_1)
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
