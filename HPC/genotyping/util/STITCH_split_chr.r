## Function to split chromosome into pieces
psf<- function(pos, binSize=7e+6, minSnp=1000, quiet=FALSE){
# binSize: bin size
# minSnp: minmum number of SNPs in each bin

   # get number of bins
   nPt <- floor((max(pos) - min(pos))/binSize)
   if(nPt < 1) nPt <- 1

   # get bins with rounding
   ps <- seq(min(pos), max(pos), length.out=nPt+1)
   ps <- round(ps)
   ps[1] <- min(pos)-1
   ps[nPt+1] <- max(pos)

   # get number of SNPs for each bin
   nS <- NULL
   for(n in 1:nPt){
      nS <- c(nS, sum(pos >= ps[n]+1 & pos <= ps[n+1]))
   }

   # merge small bins that are next to each other
   ok <- FALSE; if(nPt < 2) ok <- TRUE
   while(!ok){
      ok <- TRUE
      for(n in 1:(nPt-1)){
         if(nS[n] < minSnp && nS[n+1] < minSnp){
            ps <- ps[-n-1]
            nS[n] <- nS[n] + nS[n+1]
            nS <- nS[-n-1]
            nPt <- nPt - 1
            if(nPt > 1) ok<- FALSE
            break
         }
      }
   }

   # merge a small bin to its neighbor
   ok<- FALSE; if(nPt < 2) ok<- TRUE
   while(!ok){
      ok<- TRUE
      for(n in 1:nPt){
         if(nS[n] < minSnp){
            if(n == 1){
               nS[n] <- nS[n] + nS[n+1]
               nS <- nS[-n-1]
               ps <- ps[-n-1]
            }else if(n == nPt){
               nS[n-1] <- nS[n] + nS[n-1]
               nS <- nS[-n]
               ps <- ps[-n]
            }else{
               if(nS[n-1] > nS[n+1]){
                  nS[n]<- nS[n] + nS[n+1]
                  nS <- nS[-n-1]
                  ps <- ps[-n-1]
               }else{
                  nS[n-1]<- nS[n] + nS[n-1]
                  nS <- nS[-n]
                  ps <- ps[-n]
               }
            }
            nPt <- nPt - 1
            ok <- FALSE
            break
         }
      }
   }

   list(nPt = nPt, ps = ps)
}

##########################################################
##                     Main process                     ##
##########################################################

################### Read in arguments ###################
args <- commandArgs(TRUE)
# args[1] to specific position file
posfile <- args[1]
# args[2] to specific outfile
outfile <- args[2]

# read in pos file
positions<- read.table(posfile,header=FALSE,stringsAsFactors=FALSE)
colnames(positions)<- c("chr","pos","ref","alt")

# split the process into small process based on pos file
# set binSize to be 7 Mb
# and minimum number of SNPs as 8000
ps<- psf(positions[,"pos"], binSize=6e+6, minSnp=8000, quiet=FALSE)
ps<- ps$ps

write.table(
   ps,
   file = outfile,
   quote=FALSE, row.names=FALSE, col.names=FALSE
)