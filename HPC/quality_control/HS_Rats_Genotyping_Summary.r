################### Read in arguments ###################
args <- commandArgs(TRUE)
args
# args[1] to specific prefix
prefix <- args[1]
# args[2] to specific dir_path
dir_path <- args[2]
# args[3] to specific code path
code_path <- args[3]
# args[4] to specific mode
mode <- args[4]
# other args for sex_outliers_Sample_ID
if(length(args) >= 5){
    sex_outliers_Sample_ID <- args[5: length(args)]
}else{
    sex_outliers_Sample_ID <- c()
}

if(mode == "Part1" || mode == "Part2"){
    options(install.packages.compile.from.source = "always")
    for(pac in c("rmarkdown", "dplyr", "downloadthis", "DT", "crosstalk", "reactable")){
        if (!pac %in% installed.packages()) install.packages(pac, repos='https://cran.us.r-project.org')
    }
    if(mode == "Part1"){
        # Part 1 for this flowcell's demux and alignment
        library(rmarkdown)
        render(paste0(code_path, "/quality_control/HS_Rats_Genotyping_Summary_Part1.Rmd"), 
                params=list(dir_path=dir_path, code_path=code_path, sex_outliers_Sample_ID=sex_outliers_Sample_ID),
                output_dir=paste0(dir_path),
                output_file=paste0("HS_Rats_Genotyping_Summary_", prefix,"_Part1.html"))
    }
    if(mode == "Part2"){
        # Part 2 for genotypes
        library(rmarkdown)
        render(paste0(code_path, "/quality_control/HS_Rats_Genotyping_Summary.Rmd"), 
                params=list(dir_path=dir_path, code_path=code_path, sex_outliers_Sample_ID=sex_outliers_Sample_ID),
                output_dir=paste0(dir_path),
                output_file=paste0("HS_Rats_Genotyping_Summary_", prefix,".html"))
    }
}