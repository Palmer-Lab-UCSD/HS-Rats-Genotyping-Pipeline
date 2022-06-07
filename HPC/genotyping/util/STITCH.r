library(STITCH)

################### Read in arguments ###################
args <- commandArgs(TRUE)
args
# args[1] to specific chromosome
chr <- args[1]
# args[2] to specific K
K <- as.numeric(args[2])
# args[3] to specific number of generation
nGen <- as.numeric(args[3])
# args[4] to specific number of iteration
niterations <- as.numeric(args[4])
# args[5] to specific method selection
method <- args[5]
# args[6] to specific output directory
outputdir <- args[6]
# args[7] to specific temporary directory
tempdir <- args[7]
# args[8] to specific bam files list
bamlist <- args[8]
# args[9] to specific sample name file
sampleNames_file <- args[9]
# args[10] to specific position file
posfile <- args[10]
# args[11] to specific regionStart
regionStart <- as.numeric(args[11])
# args[12] to specific regionEnd
regionEnd <- as.numeric(args[12])
# args[13] to specific nCore
nCore <- as.numeric(args[13])

if(length(args) > 13){
   # args[14] to specific reference panel
   refPnls <- args[14]
   # choose founder haplotype data for STITCH
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   # refPnls should not include the following post-fix
   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   refHap <- paste(refPnls, ".hap.gz", sep="")
   refLgd <- paste(refPnls, ".legend.gz", sep="")
}else{
   refHap <- ""
   refLgd <- ""
}

# Other arguments for STITCH
buffer<- 1e+6
server_environment <- "server"

# create a temporary directory
tempdir <- paste(tempdir, "/stitch", niterations ,"_", chr, "_", regionStart, sep="")
tmd <- paste(tempdir, "/tmd", sep="")
rsd <- paste(tempdir, "/rsd", sep="")
if(file.exists(tempdir))
   system(paste("rm -rf ", tempdir, sep=""))
system(paste0("mkdir -p ", tempdir))
system(paste0("mkdir -p ", tmd))
system(paste0("mkdir -p ", rsd))


if(niterations <= 18){
         st<- STITCH(
                    regionStart = regionStart+1,
                      regionEnd = regionEnd,
                         buffer = buffer,
                         method = method,
                      outputdir = rsd,
                            chr = chr,
                        posfile = posfile,
                        bamlist = bamlist,
               sampleNames_file = sampleNames_file,
       reference_haplotype_file = refHap,
          reference_legend_file = refLgd,
                              K = K,
                    niterations = niterations,
     shuffleHaplotypeIterations = NA,
               refillIterations = NA,
                        tempdir = tmd,
                         nCores = nCore,
                           nGen = nGen,
           inputBundleBlockSize = NA,
       output_haplotype_dosages = TRUE
         )
} else {
         st<- STITCH(
                    regionStart = regionStart+1,
                      regionEnd = regionEnd,
                         buffer = buffer,
                         method = method,
                      outputdir = rsd,
                            chr = chr,
                        posfile = posfile,
                        bamlist = bamlist,
               sampleNames_file = sampleNames_file,
       reference_haplotype_file = refHap,
          reference_legend_file = refLgd,
                              K = K,
                    niterations = niterations,
                        tempdir = tmd,
                         nCores = nCore,
                           nGen = nGen,
           inputBundleBlockSize = NA,
       output_haplotype_dosages = TRUE
         )
}

system(paste0("mv ", rsd, "/stitch.*.vcf.gz ", outputdir))
system(paste0("rm -rf ", tempdir))
