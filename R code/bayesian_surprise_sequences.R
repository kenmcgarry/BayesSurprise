# bayesian_surprise_sequences.R
# 10/08/2021
# https://stackoverflow.com/questions/27060453/how-to-build-an-alphabetical-tree-from-a-list-of-words-in-r.

seq1 <- "dontloveababcabcaaabaacabdcaacccbacabcdeabloveyouloveyouloveyouloveyoucccdabaabcabcabcabcccloveyouloveyoudabadabadontlovedontlovedontlove"
seq2 <- "kvuasvclhihijhijhhhihhjhikjhhjjjihjhijklhisvclfvbsvclfvbsvclfvbsvclfvbjjjkhihhijhijhijhijjjsvclfvbsvclfvbkhihkhihkvuasvclkvuasvclkvuasvcl"
motifs <- find_motifs(seq2)
motiforder <- motif_sequence(seq2,motifs[[2]])  # [2] contains the list of motif names

# experiment for detecting differences in repeating pattern M1...M6
motiforder3 <- "M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4"
seq.motifs.train1 <- TraMineR::seqdef(c("M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M2",
                  "M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M2-M3-M4-M5-M1-M2-M6",
                  "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M5-M6-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6",
                  "M6-M2-M3-M4-M6-M1-M2-M3-M4-M6-M5-M1-M2-M3-M4-M5-M1-M2",
                  "M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M5-M1-M2-M3",
                  "M4-M5-M1-M2-M3-M4-M5-M1-M2-M3-M4-M2-M2-M3-M4-M5-M1-M6"))

seq.motifs.train2 <- TraMineR::seqdef(c("M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1",
                   "M6-M6-M6-M2-M2-M2-M2-M1-M1-M1"))
                   

seq.motifs3 <- TraMineR::seqdef(motiforder3)

seq.motifs.test1 <- TraMineR::seqdef(c("M3-M6-M4-M5-M1-M2-M3-M4-M5-M6-M1-M5-M6-M2-M3-M4-M1"))
seq.motifs.test2 <- TraMineR::seqdef(c("M6-M6-M6-M2-M2-M2-M2-M1-M1-M1"))
seq.motifs.test22 <- TraMineR::seqdef(c("M2-M2-M2-M2-M1-M1-M1-M6-M6-M6"))
seq.motifs.test23 <- TraMineR::seqdef(c("M6-M6-M6-M2-M2-M2-M2-M1-M1-M1"))
seq.motifs.test24 <- TraMineR::seqdef(c("M6-M6-M6-M6-M1-M2-M1-M2-M1-M2-M1-M2-M1-M2"))

seq.motifs.train2 <- TraMineR::seqdef(motifs.train2)


# IDEA #1
# write the sequence to an HTML file and colour code the motifs found within it.
# https://cran.r-project.org/web/packages/tableHTML/vignettes/tableHTML.html
#library(tableHTML)
#create an html table 
#tableHTML(mtcars)
#and to export in a file
#write_tableHTML(tableHTML(mtcars), file = 'myfile.html')
#library("highlightr")
#library("tibble")
#dict <- tibble(
#  feature = c("daba", "abc","loveyou","ccc","aa","dontlove"),
#  bg_colour = c("lightgrey", "cyan","pink","orange","lightblue","lightgreen"),
#  bold=TRUE)
#highlight(ex4, dict)

