# bayesian_surprise_sequences.R
# 10/08/2021
# https://stackoverflow.com/questions/27060453/how-to-build-an-alphabetical-tree-from-a-list-of-words-in-r.

seq2 <- "kvuasvclhihijhijhhhihhjhikjhhjjjihjhijklhisvclfvbsvclfvbsvclfvbsvclfvbjjjkhihhijhijhijhijjjsvclfvbsvclfvbkhihkhihkvuasvclkvuasvclkvuasvcl"
motifs <- find_motifs(seq2)
motiforder <- motif_sequence(seq2,motifs[[2]])  # [2] contains the list of motif names

# experiment for detecting differences in repeating pattern M1...M6

motifs.seq <- TraMineR::seqdef(c(
                  "M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M2",
                  "M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M2-M3-M4-M5-M1-M2-M6",
                  "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M4-M2-M2-M3-M4-M5-M1-M6",
                  "M4-M5-M1-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M2",
                  "M1-M2-M3-M4-M5-M1-M1-M2-M3-M4-M5-M1-M2-M3",
                  "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6"))

seqsample <- sample.int(n = nrow(motifs.seq),size=floor(.5*nrow(motifs.seq)), replace = FALSE)
motifs.seq.train <- motifs.seq[seqsample, ]
motifs.seq.test  <- motifs.seq[-seqsample, ]




#highlight(ex4, dict)

