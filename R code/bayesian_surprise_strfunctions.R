# bayesian_surprise_strfunctions.R
# 11/08/2021
# string functions, manipulate and get stats on string (sequence) composition

# 
seq_stats <- function(seq1){
  
  seq1 <- unlist(strsplit(seq1, split = "-"))  # strings deliminated by "-"
  
  unsym <- length(unique(seq1))
  symcount <- length(seq1)
  
  temp <- paste(seq1, collapse = "")
  temp <- unlist(strsplit(temp, split = ""))
  
  total <- sum(table(temp))
  letter_relfreq <- table(temp)/total
  
  relfreqlet <- as.numeric(letter_relfreq)
  tab_names <- names(letter_relfreq)    
  names(relfreqlet) <- tab_names
  
  symfreqtmp <- table(seq1)/symcount
  symfreq <- as.numeric(symfreqtmp)
  tab_names <- names(symfreqtmp)
  names(symfreq) <- tab_names
  
  elist <- list(unique_symbols=unsym,symbol_count=symcount,totalletters=total,
                relfreqLetters=relfreqlet,relfreqSymbols=symfreq)
  return(elist)
  
}


###
reps <- function(s, n) paste(rep(s, n), collapse = "") # repeat s n times

find_string <- function(string, th = 1, len = floor(nchar(string)/th)) {
  for(k in len:1) {
    pat <- paste0("(.{", k, "})", reps("\\1", th-1))
    r <- regexpr(pat, string, perl = TRUE)
    if (attr(r, "capture.length") > 0) break
  }
  if (r > 0) substring(string, r, r + attr(r, "capture.length")-1) else ""
}


find_string2 <- function(string, th = 2, len = floor(nchar(string)/th)) {
  pat <- paste0(c("(.", "{1,", len, "})", rep("\\1", th-1)), collapse = "")
  r <- regexpr(pat, string, perl = TRUE)
  ifelse(r > 0, substring(string, r, r + attr(r, "capture.length")-1), "")
}


# The function that does the matching:
  find_rep_path <- function(vec, reps) {
    regexp <- paste0(c("(.+)", rep("\\1", reps - 1L)), collapse="")
    match <- regmatches(vec, regexpr(regexp, vec, perl=T))
    substr(match, 1, nchar(match) / reps)  
  }

# find_pattern() by Stephen Henderson on stackoverflow
# https://stackoverflow.com/questions/21020032/algorithm-code-in-r-to-find-pattern-from-any-position-in-a-string
find_pattern <- function(pat){
  
  len=nchar(pat)
  thr=3
  reps=floor(len/3)
  
  # all poss strings up to half length of pattern
  pat=str_split(pat, "")[[1]][-1]
  str.vec=vector()
  for(win in 2:reps)
  {
    str.vec= c(str.vec, zoo::rollapply(data=pat,width=win,FUN=paste0, collapse=""))
  }
  
  # the max length string repeated more than 3 times
  tbl=table(str.vec)
  tbl=tbl[tbl >= 3]
  
  return(tbl[which.max(nchar(names(tbl)))])
  
}
  

# find_motifs()
find_motifs <- function(seqstr){
  
  # get freq counts of the alphabet present in the sequence string (seqstr)
  alphabet <- unlist(strsplit(seqstr, split = ""))  # strings deliminated by ""
  unqsym <- length(unique(alphabet))  # save in list
  symcount <- length(alphabet)
  
  relfreq <- paste(seqstr, collapse = "")
  relfreq <- unlist(strsplit(relfreq, split = ""))
 
  total <- sum(table(relfreq))              # save in list
  letter_relfreq <- table(relfreq)/total  # save in list
  
  seq_stats <- list(unqsym,total,letter_relfreq)
  names(seq_stats) <- c("NoUniqSyms", "TotalSyms", "RelFreqs")
  
  motiflist <- list()  # setup motiflist prior to use
  i <- 1
  tempvec <- seqstr
  motif <- "NULL"
  # While loop finds the longest repeating string "motif"
  while(length(motif)>0){
    motif <- find_pattern(tempvec)
    if(length(motif)==0) {  # sorry! my clumsy way of breaking out of while loop!
      break
    }
    tempvec <- str_remove_all(tempvec, names(motif))  # remove motif from tempvec and search again
    motiflist[i] <- names(motif)  # save it in motiflist
    cat("\nmotif is ",as.character(motiflist[i]), " i is ",i)
    i <- i+1
  }
  
  # Can motifs be broken down further? e.g "abcabca" can become "abc" and "abc" with "a" discarded (as noise?)
  cat("\n")
  
  motifvalues <- list(seq_stats,motiflist)  # join stats with found motifs and return
  return(motifvalues)
}


# motif_sequence() will identify the locations of each motif in the original sequence, note each motif
# is likely appear several times. 
  # We need to create:
  # i.   assign M1 identifier to each motif M1, M2, Mn etc
  # ii.  order and position of every motif using the "M" identifiers
  # iii. some motifs overlap, causing collision issues with some motifs not appearing. This issue
  #      can be repaired when only two conflict.
  # iv. finally return the sequence of motifs.

motif_sequence <- function(seqstr,motlist){
  nseq <- length(unlist(strsplit(seqstr, split = "")))
  locations <- rep(NA, nseq)
  nmotifs <- length(motlist)  # how many motifs?
  locationlist <- 0
  collision <- matrix(0, nrow = nmotifs, ncol = nseq)  # collision matrix when motif locations overlap
  
  ncount <- 0
  
  # first for loop is simply to get a count of how many times each motif appears
  for (i in 1:nmotifs) {
     positions <- str_locate_all(seqstr, as.character(motlist[[i]]))
     positions <- positions[[1]]  # we just need the start positions and how many
     tmpcount <- nrow(positions)
     ncount <- ncount +tmpcount
     #cat("\nTotal is ", as.character(motlist[i]))
  }

  cat("\nWe have",ncount,"occurences of",nmotifs,"unique motifs...")
  
  # for each motif find out where it occurs and label it with an M number
  for (i in 1:nmotifs) {
    positions <- str_locate_all(seqstr, as.character(motlist[[i]])) #returns the start & end positions
    positions <- positions[[1]]  # we just need the start positions and how many
    tmplocation <- positions[,1]
    locationlist <- c(locationlist, tmplocation)
    collision[i,positions[,1]] <- 1

  }
  
 
  # Repair collisions between TWO motifs - cannot repair > 2 collisions
  for(i in 1:nseq){
    if(sum(collision[,i]) >1){  # we found a collision!!!
      temp <- collision[,i]
      colind <- which(temp > 0, arr.ind = TRUE)
      cat("\ncollision detected between motifs", colind," at index ",i)
      if(sum(collision[,i-1]) == 0){  # Is there a free space BEFORE collision index?
        cat("\nTrying to repair....")
        collision[colind[1],i-1] <- 1
        collision[colind[1],i] <- 0
      }else{
        if(sum(collision[,i+1]) == 0){  # Is there a free space AFTER collision index?
          cat("\nTrying to repair....")
          collision[colind[1],i+1] <- 1
          collision[colind[1],i] <- 0 }
      }
        
    }
  }
  
  
  # Now convert collision matrix into a sequence or list of motifs as they appear!
  #cat("\ndim",dim(collision))
  newlist <- rep(0, nseq)
  for(i in 1:nseq){
     temp <- collision[,i]
     temp <- which(temp > 0, arr.ind = TRUE)
     if(length(temp)==0){newlist[i] <- 0}
     else{newlist[i] <- temp}
    
  }
  
  newlist <- newlist[newlist != 0 ] # get rid of all unused (zero) entries
  newlist <- paste("M",as.character(newlist),sep = "")  # this version appends "M" but converts vector to strings
  return(newlist)
  
}







  