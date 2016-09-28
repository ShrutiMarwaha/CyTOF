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
      assign( paste(condition,cell,stimulant,sep="_"), getData(df_name,cell,condition,stimulant) )     
      # NOT WORKING for names with sapce
    }
  }
}

# function to extract cell count for each cell
extract_cellcount <- function(data_frame){
  return(as.numeric(data_frame[1,"X.Cells"]))
}

df_cellcounts <- data.frame(matrix(data=NA,nrow=4,ncol=2,dimnames=list(cells,conditions) ))
df_celltypes["Basophils","Healthy"] <- extract_cellcount(healthy_Basophils_1)
df_celltypes["Basophils","Disease"] <- extract_cellcount(disease_Basophils_1)
# NOT WORKING
df_celltypes["CD16 Hi Monocytes","Healthy"] <- extract_cellcount(healthy_CD16 Hi Monocytes_1) 
df_celltypes["CD16 Hi Monocytes","Disease"] <- as.numeric(disease_HiMonocytes_unstimulated[1,"X.Cells"])
df_celltypes["CD16 Low Monocytes","Healthy"] <- as.numeric(healthy_LowMonocytes_unstimulated[1,"X.Cells"])
df_celltypes["CD16 Low Monocytes","Disease"] <- as.numeric(disease_LowMonocytes_unstimulated[1,"X.Cells"])
df_celltypes["Lymphocytes","Healthy"] <- as.numeric(healthy_Lymphocytes_unstimulated[1,"X.Cells"])
df_celltypes["Lymphocytes","Disease"] <- as.numeric(disease_Lymphocytes_unstimulated[1,"X.Cells"])
df_celltypes
