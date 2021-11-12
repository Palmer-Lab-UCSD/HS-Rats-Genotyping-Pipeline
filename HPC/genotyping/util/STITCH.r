###### TODO
###### 1. ${code}/STITCH.r BIG CHANGE HERE
######    a. bcftools path
######    b. conflicts between pos file and reference panel
######    c. currently this has to take a reference panel
######       and the pos file is generated from the reference panel


##########################################################
## Function to split chromosome into pieces to speed up ##
##########################################################
psf<- function(pos, binSize=7e+6, minSnp=1000, quiet=FALSE){
# binSize: bin size
# minSnp: minmum number of SNPs in each bin
   nPt<- floor((max(pos) - min(pos))/binSize)
      #nPt<- round(nPt/nCpu)
      if(nPt < 1) nPt<- 1
      #nPt<- nPt * nCpu
   ps<- seq(min(pos), max(pos), length.out=nPt+1)
      ps<- round(ps)
      ps[1]<- min(pos)-1
      ps[nPt+1]<- max(pos)
   nS<- NULL
   for(n in 1:nPt){
      nS<- c(nS, sum(pos >= ps[n]+1 & pos <= ps[n+1]))
   }
   ms<- min(minSnp, quantile(nS,0.25))
      ms<- ceiling(ms)
   ok<- FALSE; if(nPt < 2) ok<- TRUE
   while(!ok){# merge small ones that are next to each other
      ok<- TRUE
      for(n in 1:(nPt-1)){
         if(nS[n] < minSnp && nS[n+1] < minSnp){
            ps<- ps[-n-1]
            nS[n]<- nS[n] + nS[n+1]
               nS<- nS[-n-1]
            nPt<- nPt - 1

            if(nPt > 1) ok<- FALSE
            break
         }
      }
   }
   ok<- FALSE; if(nPt < 2) ok<- TRUE
   while(!ok){# merge a small one to its neighbor
      ok<- TRUE
      for(n in 1:nPt){
         if(nS[n] < minSnp){
            if(n == 1){
               nS[n+1]<- nS[n] + nS[n+1]
            }else if(n == nPt){
               nS[n-1]<- nS[n] + nS[n-1]
            }else{
               if(nS[n-1] > nS[n+1]){
                  nS[n+1]<- nS[n] + nS[n+1]
               }else{
                  nS[n-1]<- nS[n] + nS[n-1]
               }
            }
            nS<- nS[-n]
            ps<- ps[-n]
            nPt<- nPt - 1

            ok<- FALSE
            break
         }
      }
   }

   if(!quiet){
      cat("Number of jobs: ", nPt, "\n", sep="")
      cat("Bin size (MB):\n"); print(summary(diff(ps)/1e+6))
   }

   list(nPt = nPt, ps = ps)
}

##########################################################
##                     Main process                     ##
##########################################################
##                 with reference panel                 ##
##########################################################
library(parallel)
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
# args[10] to specific reference panel
refPnls <- args[10]

# Other arguments for STITCH
nCore<- 1
buffer<- 1e+6
server_environment <- "server"
inputBundleBlockSize <- NA
genfile <- ""
# reference genome no need for bamfiles
# ref<- "/projects/ps-palmer/reference_genomes/rat/rn6.fa"

# create a temporary directory
tempdir <- paste(tempdir, "/Tmp", chr, sep="")
datadir <- file.path(tempdir, paste(chr, "Data", sep=""))
resultdir <- file.path(tempdir, paste(chr, "Result", sep=""))
if(!file.exists(datadir))
   system(paste0("mkdir -p ", datadir))
if(!file.exists(resultdir))
   system(paste0("mkdir -p ", resultdir))
system(paste0("rm -r ", resultdir, "/*", sep=""), ignore.stderr = TRUE)

# choose founder haplotype data for STITCH
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# In folder refPnls, the file names follow the format phasedVariants_chr
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
refHap<- paste(refPnls, "/phasedVariants_", chr, ".hap.gz", sep="")
refLgd<- paste(refPnls, "/phasedVariants_", chr, ".legend.gz", sep="")

# generate pos file based on reference panel
phyMap<- read.table(refLgd,header=TRUE,stringsAsFactors=FALSE)
phyMap<- cbind(sapply(strsplit(phyMap[,"id"],":"),function(x)x[1]),phyMap[,c(2,1,3,4)])
colnames(phyMap)<- c("chr","pos","id","ref","alt")
summary(phyMap[,"pos"])
posfile <- file.path(datadir, paste("pos", chr, ".txt", sep=""))
write.table(
   phyMap[,-3],
   file = posfile,
   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
)

# split the process into small process based on pos file
# set binSize to be 7 Mb
# and minimum number of SNPs as 10000
ps<- psf(phyMap[,"pos"], binSize=7e+6, minSnp=10000, quiet=FALSE)
   nPt<- ps$nPt
   ps<- ps$ps
nCpu<- min(nPt,5)

# run STITCH in parallel
lst<- list()
nes<- 0
pid<- NULL
cnt<- 0
for(n in 1:nPt){
   cnt<- cnt+1
   pid<- c(pid, mcparallel(
   {###---------------------
      cat("stitch: now process ", n, " out of ", nPt, " blocks. Please wait...\n", sep="")
      tmd<- paste(resultdir, "/tmd", sprintf("%03d", n, sep=""), sep="")
      if(!file.exists(tmd)){
         system(paste0("mkdir -p ", tmd))
      }else{
         system(paste("rm -rf ", tmd, sep=""))
         system(paste0("mkdir -p ", tmd))
      }

      rsd <- file.path(resultdir, paste("rsd", sprintf("%03d",n), sep=""))
      if(!file.exists(rsd)){
         system(paste0("mkdir -p ", rsd))
      }else{
         system(paste("rm -rf ", rsd, sep=""))
         system(paste0("mkdir -p ", rsd))
      }
      if(niterations <= 18){
               st<- STITCH(
                          regionStart = ps[n]+1,
                            regionEnd = ps[n+1],
                               buffer = buffer,
                               method = method,
                            outputdir = rsd,
                                  chr = chr,
                              posfile = posfile,
                              genfile = genfile,
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
                 inputBundleBlockSize = inputBundleBlockSize
               )
      } else {
               st<- STITCH(
                          regionStart = ps[n]+1,
                            regionEnd = ps[n+1],
                               buffer = buffer,
                               method = method,
                            outputdir = rsd,
                                  chr = chr,
                              posfile = posfile,
                              genfile = genfile,
                              bamlist = bamlist,
                     sampleNames_file = sampleNames_file,
             reference_haplotype_file = refHap,
                reference_legend_file = refLgd,
                                    K = K,
                          niterations = niterations,
                              tempdir = tmd,
                               nCores = nCore,
                                 nGen = nGen,
                 inputBundleBlockSize = inputBundleBlockSize
               )
      }
      # Please update bcftools version
      system(paste("/projects/ps-palmer/software/local/src/bcftools-1.13/bcftools index -f -t ", rsd, "/stitch.chr*.vcf.gz", sep=""))
      cat("\n=>", sprintf("%3d",n), "/", nPt, " is done!\n", sep="")
      system(paste("rm -rf ", tmd, sep=""))
      
      paste("Job #", n, ": okay!", sep="")
   },###---------------------
      mc.set.seed = FALSE,
      silent = FALSE
   )$pid)
   Sys.sleep(5)
   if(n == nPt){
      oTmp<- mccollect(pid,wait=TRUE)
      # process...
      lst<- c(lst, oTmp)
      nes<- sum(nes, length(grep("Error",oTmp)))
      gc()
   }else if(cnt %% nCpu == 0){
      while(TRUE){# wait for some job to finish
         Sys.sleep(30)
         idxTmp<- NULL
         for(mTmp in 1:cnt){
            oTmp<- mccollect(pid[mTmp],wait=FALSE)
            if(!is.null(oTmp[[1]])){
               # process...
               lst<- c(lst, oTmp)
               idxTmp<- c(idxTmp, mTmp)
               nes<- sum(nes, length(grep("Error",oTmp)))
            }
         }
         if(length(idxTmp) > 0){
            cnt<- cnt - length(idxTmp)
            pid<- pid[-idxTmp]
            break
         }
         rm(idxTmp, mTmp, oTmp)
      }
      gc()
   }
}

lst

cat("********Number of failed blocks:", nes, "********\n")

fls<- NULL
for(n in 1:nPt){
   f<- list.files(paste(resultdir, "/rsd", sprintf("%03d",n), sep=""), "stitch.*.vcf.gz$")
   if(length(f) != 1){
      cat(n, " ")
      next
   }
   f<- paste(resultdir, "/rsd", sprintf("%03d",n),"/",f,sep="")
   fls<- c(fls, f)
}
cat("...Move on...\n")
if(length(fls) != nPt) stop("File counts wrong")
str<- "/projects/ps-palmer/software/local/src/bcftools-1.13/bcftools concat --no-version -a -d none"
   str<- paste(str, " -O z -o ", outputdir, "/", chr, "_stitch.vcf.gz", sep="")
   str<- paste(str, " ", paste(fls, collapse=" "), sep="")
system(str, ignore.stdout=FALSE)

str<- paste("rm -rf ", tempdir, sep="")
system(str)

q("no")

##################
# the end #
###########
