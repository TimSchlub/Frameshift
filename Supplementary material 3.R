nuc <- "atggcaaaacaattgcgtgtaactttagatcagaccttaacccagattacaaagagagttgatgaagacaacgacgaaactaatgaaccaatacaaccaaaagaacaagaaaatattcaaagagaaataaattgcatccaagtggaatcattcaatgatacaacggtgaaatgtgaatttcaaagccaatttagacagttaatcacatatccagtgcatatattctacatacacaataatacatcatttacccaagtcgaggcaatcattcaagaaattaacgatagccacgtccacattttaattatttgcattgatacaaatttattcagtaaacatttagctggaattttaaatacacaatcgttgctgatatttgcatttaaaccagtttggatcggtaaagcatttgattttctacttgattcaggagttttaatagaaccagtaacacatgaacagataaattttgacgagataattgaaggtatagaagcacggaaacaaaatgactctattgatgttcataacgtccaagatagtattgtacctatgctacaaagtgggcttgttgtcagtaccatttcgaataatcaagaatttcatcacactacatctatattacgccataacacactcaagaatgctatatcagatatgtcggatgatgttaatattgtcaatgctgtcactaatcgcgatttacgggtacttcttccatcagtcgcttccaacactacagagtctctcagaagtcaagaattaccaaggcattctaagtcgtttcaatgtgacgccagaacagttaaatttgtatcgcaaattcatgacatcagttatagaacaggcagaagaatccaagaacacagactagttgatacaattccagataacgctaacaacgccagctatcacaatattcgatatactgtgttaggtaatgttttactaatggtaccgtatgatgacatgtttcacatgtttgaagagtactcattagatttggcttccgttagtttcgcgttgagtatcaaata"
n <- 2000
n_synon <- 2000

# Defining functions -------------------------------------------------------------------------------------

# This function returns the reverse sequence
strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

# This function returns the complimentary sequence 
complem <- function(x) {
  x <- tolower(x)
  x <- gsub("a","T",x)
  x <- gsub("c","G",x)
  x <- gsub("t","A",x)
  x <- gsub("g","C",x)
  return(tolower(x))
}

findORF <- function(nuc, fs=0) {
# Returns a list of nucleotide open reading frames that are bookended by stop codons
# (TGA, TAG or TAA), after a frame shift of fs is applied. fs must be 0,1, or 2
  if (!(fs %in% c(0,1,2))) stop("fs must take value 0, 1 or 2")
  nuc <- toupper(paste(nuc,collapse=""))
  
  nuc <- substr(nuc,fs+1, nchar(nuc)) # nucleotide sequence after frameshift
  codons <- substring(nuc, seq(1,nchar(nuc),3), seq(3,nchar(nuc),3)) # split by codons
  codons <- gsub("TAA|TAG|TGA","999",codons) # replacing stop codons with 999
  ORFs <- strsplit(paste(codons,collapse=""),"999") # nucleotide sequence of each ORF
  
  return(unlist(ORFs))
}

synonymousMutation <- function(codon, codon_table){
# Returns a random nucleotide sequence that encodes the same codon
  codon <- as.character(codon)
  codon <- toupper(codon)
  aa <- codon_table[codon_table$nucleotide==codon,]$a
  if (length(aa)==1) {
    return(as.character(sample(codon_table[codon_table$a==aa,]$nucleotide ,1)))
  } else return(codon)
}

codon_permute_bootstrap <- function(nuc,n=1000,seed=NULL){
# This function calculates the observed and expected distributions of frameshifted proteins
# This function calculates the observed and expected distributions of the length 
# of ORFs under different reading frames using the codon permutation method.
# nuc is the nucleotide sequence, n is the number of bootstraps (default = 1000), seed is
# used to set the random seed for reproducable results.
  if (n < 1) stop ("n must be 1 or greater")
  if (class(seed) == "numeric") set.seed(seed);
  nuc <- toupper(paste(nuc,collapse=""))
  
  codons <- substring(nuc, seq(1,nchar(nuc),3), seq(3,nchar(nuc),3))
  b.distrib0 <- list(NULL) # Original reading frame
  b.distrib1 <- list(NULL) # +1 reading frame
  b.distrib2 <- list(NULL) # +2 reading frame
  com.b.distrib0 <- list(NULL) # c0 reading frame
  com.b.distrib1 <- list(NULL) # c1 reading frame
  com.b.distrib2 <- list(NULL) # c2 reading frame
  rev.b.distrib0 <- list(NULL) # -0 reading frame
  rev.b.distrib1 <- list(NULL) # -1 reading frame
  rev.b.distrib2 <- list(NULL) # -2 reading frame
  revcom.b.distrib0 <- list(NULL) # -c0 reading frame 
  revcom.b.distrib1 <- list(NULL) # -c1 reading frame
  revcom.b.distrib2 <- list(NULL) # -c2 reading frame
  for (i in 1:n) { # Loops through bootstraps
    
    nuc.permute <- paste(sample(codons,length(codons),replace=FALSE), collapse="") # permutes the codons
    
    com.nuc.permute <- complem(nuc.permute) # complementary sequence (c0, c1, c2)
    rev.nuc.permute <- strReverse(nuc.permute) # reverse sequence (-0, -1, -2)
    revcom.nuc.permute <- complem(rev.nuc.permute) # reverse complementary sequence (-c0, -c1, -c2)
    
    b.ORF0 <- findORF(nuc.permute,0) 
    b.ORF1 <- findORF(nuc.permute,1)
    b.ORF2 <- findORF(nuc.permute,2)
    
    com.b.ORF0 <- findORF(com.nuc.permute,0) 
    com.b.ORF1 <- findORF(com.nuc.permute,1)
    com.b.ORF2 <- findORF(com.nuc.permute,2)
    
    rev.b.ORF0 <- findORF(rev.nuc.permute,0) 
    rev.b.ORF1 <- findORF(rev.nuc.permute,1)
    rev.b.ORF2 <- findORF(rev.nuc.permute,2)
    
    revcom.b.ORF0 <- findORF(revcom.nuc.permute,0) 
    revcom.b.ORF1 <- findORF(revcom.nuc.permute,1)
    revcom.b.ORF2 <- findORF(revcom.nuc.permute,2)
    
    b.distrib0[[i]] <- nchar(b.ORF0)
    b.distrib1[[i]] <- nchar(b.ORF1)
    b.distrib2[[i]] <- nchar(b.ORF2)
    
    com.b.distrib0[[i]] <- nchar(com.b.ORF0)
    com.b.distrib1[[i]] <- nchar(com.b.ORF1)
    com.b.distrib2[[i]] <- nchar(com.b.ORF2)
    
    rev.b.distrib0[[i]] <- nchar(rev.b.ORF0)
    rev.b.distrib1[[i]] <- nchar(rev.b.ORF1)
    rev.b.distrib2[[i]] <- nchar(rev.b.ORF2)
    
    revcom.b.distrib0[[i]] <- nchar(revcom.b.ORF0)
    revcom.b.distrib1[[i]] <- nchar(revcom.b.ORF1)
    revcom.b.distrib2[[i]] <- nchar(revcom.b.ORF2)
    
    rm(nuc.permute, b.ORF0, b.ORF1, b.ORF2
       ,com.nuc.permute, com.b.ORF0, com.b.ORF1, com.b.ORF2
       ,rev.nuc.permute, rev.b.ORF0, rev.b.ORF1, rev.b.ORF2
       ,revcom.nuc.permute, revcom.b.ORF0, revcom.b.ORF1, revcom.b.ORF2)
    #if (i %% round(n/100) == 0){ Sys.sleep(0.01); print(i); flush.console() } # Uncomment this line to track progress
  }
  
  
  returnList <- list( b.distrib0, com.b.distrib0, rev.b.distrib0, revcom.b.distrib0
                     ,b.distrib1, com.b.distrib1, rev.b.distrib1, revcom.b.distrib1
                     ,b.distrib2, com.b.distrib2, rev.b.distrib2, revcom.b.distrib2
                      )
  
  names(returnList) <- c( "b.distrib0", "com.b.distrib0", "rev.b.distrib0", "revcom.b.distrib0"
                         ,"b.distrib1", "com.b.distrib1", "rev.b.distrib1", "revcom.b.distrib1"
                         ,"b.distrib2", "com.b.distrib2", "rev.b.distrib2", "revcom.b.distrib2"
                         )
  return(returnList)
}
  
synonymous_mutation_bootstrap <- function(nuc,n=1000,seed=NULL){
  # This function calculates the observed and expected distributions of frameshifted proteins
  # This function calculates the observed and expected distributions of the length 
  # of ORFs under different reading frames using the synonymous mutation method.
  # nuc is the nucleotide sequence, n is the number of bootstraps (default = 1000), seed is
  # used to set the random seed for reproducable results.
  if (n < 1) stop ("n must be 1 or greater")
  if (class(seed) == "numeric") set.seed(seed);
  nuc <- toupper(paste(nuc,collapse=""))
  
  codons <- substring(nuc, seq(1,nchar(nuc),3), seq(3,nchar(nuc),3))
  b.distrib0 <- list(NULL)
  b.distrib1 <- list(NULL)
  b.distrib2 <- list(NULL)
  com.b.distrib0 <- list(NULL)
  com.b.distrib1 <- list(NULL)
  com.b.distrib2 <- list(NULL)
  rev.b.distrib0 <- list(NULL)
  rev.b.distrib1 <- list(NULL)
  rev.b.distrib2 <- list(NULL)
  revcom.b.distrib0 <- list(NULL)
  revcom.b.distrib1 <- list(NULL)
  revcom.b.distrib2 <- list(NULL)
  
  # Table of codons
  codon_table <- readRDS("codon_table.rds") 
  
  for (i in 1:n) {
    # randomly generates a nucleotide sequence that encodes the same codon sequnce in frame 0
    nuc.permute <- paste(sapply(codons, FUN=synonymousMutation, codon_table),collapse="")
    
    com.nuc.permute <- complem(nuc.permute)
    rev.nuc.permute <- strReverse(nuc.permute)
    revcom.nuc.permute <- complem(rev.nuc.permute)
    
    b.ORF0 <- findORF(nuc.permute,0) 
    b.ORF1 <- findORF(nuc.permute,1)
    b.ORF2 <- findORF(nuc.permute,2)
    
    com.b.ORF0 <- findORF(com.nuc.permute,0) 
    com.b.ORF1 <- findORF(com.nuc.permute,1)
    com.b.ORF2 <- findORF(com.nuc.permute,2)
    
    rev.b.ORF0 <- findORF(rev.nuc.permute,0) 
    rev.b.ORF1 <- findORF(rev.nuc.permute,1)
    rev.b.ORF2 <- findORF(rev.nuc.permute,2)
    
    revcom.b.ORF0 <- findORF(revcom.nuc.permute,0) 
    revcom.b.ORF1 <- findORF(revcom.nuc.permute,1)
    revcom.b.ORF2 <- findORF(revcom.nuc.permute,2)
    
    b.distrib0[[i]] <- nchar(b.ORF0)
    b.distrib1[[i]] <- nchar(b.ORF1)
    b.distrib2[[i]] <- nchar(b.ORF2)
    
    com.b.distrib0[[i]] <- nchar(com.b.ORF0)
    com.b.distrib1[[i]] <- nchar(com.b.ORF1)
    com.b.distrib2[[i]] <- nchar(com.b.ORF2)
    
    rev.b.distrib0[[i]] <- nchar(rev.b.ORF0)
    rev.b.distrib1[[i]] <- nchar(rev.b.ORF1)
    rev.b.distrib2[[i]] <- nchar(rev.b.ORF2)
    
    revcom.b.distrib0[[i]] <- nchar(revcom.b.ORF0)
    revcom.b.distrib1[[i]] <- nchar(revcom.b.ORF1)
    revcom.b.distrib2[[i]] <- nchar(revcom.b.ORF2)
    
    rm(nuc.permute, b.ORF0, b.ORF1, b.ORF2
       ,com.nuc.permute, com.b.ORF0, com.b.ORF1, com.b.ORF2
       ,rev.nuc.permute, rev.b.ORF0, rev.b.ORF1, rev.b.ORF2
       ,revcom.nuc.permute, revcom.b.ORF0, revcom.b.ORF1, revcom.b.ORF2)
  }
  
  returnList <- list(  b.distrib0, com.b.distrib0, rev.b.distrib0, revcom.b.distrib0
                      ,b.distrib1, com.b.distrib1, rev.b.distrib1, revcom.b.distrib1
                      ,b.distrib2, com.b.distrib2, rev.b.distrib2, revcom.b.distrib2
  )
  
  names(returnList) <- c( "b.distrib0", "com.b.distrib0", "rev.b.distrib0", "revcom.b.distrib0"
                          ,"b.distrib1", "com.b.distrib1", "rev.b.distrib1", "revcom.b.distrib1"
                          ,"b.distrib2", "com.b.distrib2", "rev.b.distrib2", "revcom.b.distrib2"
  )
  return(returnList)
  
}

# A function to formats the results
bootResults <- function(nuc.ORF,permute.distrib,synon.distrib,fs){
  x <- data.frame(nuc_ORF = nuc.ORF
                  ,aa_length = round(nchar(nuc.ORF)/3)
                  ,Pval_permute = 1-permute.distrib(nchar(nuc.ORF)-0.1)^length(nuc.ORF)
                  ,Pval_synon = 1-synon.distrib(nchar(nuc.ORF)-0.1)^length(nuc.ORF)
                  ,fs = fs
  )  
  
  nstart <- head(c(1,(cumsum(x$aa_length+1)*3)+1),-1)
  nend <- c(tail(nstart,-1)-1,sum(x$aa_length+1)*3)
  
  
  if (fs %in% c("plus_1","com_plus_1","rev_plus_1","revcom_plus_1")) {
    nstart <- nstart + 1
    nend <- nend + 1
    nend[length(nend)] <- tail(nend,1)-3
  }
  if (fs %in% c("plus_2","com_plus_2","rev_plus_2","revcom_plus_2")) {
    nstart <- nstart + 2
    nend <- nend + 2
    nend[length(nend)] <- tail(nend,1)-3
  }
  x$nstart <- nstart
  x$end <- nend
  return(x)
}


#Runs the analysis ------------------------------------------------------------
com.nuc <- complem(nuc)
rev.nuc <- strReverse(nuc)
revcom.nuc <- complem(rev.nuc)

ORF0 <- findORF(nuc,0) 
ORF1 <- findORF(nuc,1)
ORF2 <- findORF(nuc,2)
com.ORF0 <- findORF(com.nuc,0)
com.ORF1 <- findORF(com.nuc,1)
com.ORF2 <- findORF(com.nuc,2)
rev.ORF0 <- findORF(rev.nuc,0)
rev.ORF1 <- findORF(rev.nuc,1)
rev.ORF2 <- findORF(rev.nuc,2)
revcom.ORF0 <- findORF(revcom.nuc,0)
revcom.ORF1 <- findORF(revcom.nuc,1)
revcom.ORF2 <- findORF(revcom.nuc,2)

distrib.permute <- codon_permute_bootstrap(nuc,n)
distrib.synon <- synonymous_mutation_bootstrap(nuc,n_synon)

Results <- rbind(
  bootResults(ORF0,        ecdf(unlist(distrib.permute$b.distrib0)),        ecdf(unlist(distrib.synon$b.distrib0)),  "plus_0")
  ,bootResults(ORF1,        ecdf(unlist(distrib.permute$b.distrib1)),        ecdf(unlist(distrib.synon$b.distrib1)),  "plus_1")
  ,bootResults(ORF2,        ecdf(unlist(distrib.permute$b.distrib2)),        ecdf(unlist(distrib.synon$b.distrib2)),  "plus_2")
  ,bootResults(com.ORF0,    ecdf(unlist(distrib.permute$com.b.distrib0)),    ecdf(unlist(distrib.synon$com.b.distrib0)),  "com_plus_0")
  ,bootResults(com.ORF1,    ecdf(unlist(distrib.permute$com.b.distrib1)),    ecdf(unlist(distrib.synon$com.b.distrib1)),  "com_plus_1")
  ,bootResults(com.ORF2,    ecdf(unlist(distrib.permute$com.b.distrib2)),    ecdf(unlist(distrib.synon$com.b.distrib2)),  "com_plus_2")
  ,bootResults(rev.ORF0,    ecdf(unlist(distrib.permute$rev.b.distrib0)),    ecdf(unlist(distrib.synon$rev.b.distrib0)),  "rev_plus_0")
  ,bootResults(rev.ORF1,    ecdf(unlist(distrib.permute$rev.b.distrib1)),    ecdf(unlist(distrib.synon$rev.b.distrib1)),  "rev_plus_1")
  ,bootResults(rev.ORF2,    ecdf(unlist(distrib.permute$rev.b.distrib2)),    ecdf(unlist(distrib.synon$rev.b.distrib2)),  "rev_plus_2")
  ,bootResults(revcom.ORF0, ecdf(unlist(distrib.permute$revcom.b.distrib0)), ecdf(unlist(distrib.synon$revcom.b.distrib0)),   "revcom_plus_0")
  ,bootResults(revcom.ORF1, ecdf(unlist(distrib.permute$revcom.b.distrib1)), ecdf(unlist(distrib.synon$revcom.b.distrib1)),   "revcom_plus_1")
  ,bootResults(revcom.ORF2, ecdf(unlist(distrib.permute$revcom.b.distrib2)), ecdf(unlist(distrib.synon$revcom.b.distrib2)),   "revcom_plus_2")
)
