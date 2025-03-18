#Installing BioConductor packages ----
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("BSgenome")
BiocManager::install("BiocStyle")

#Loading Packages----
library(GenomicRanges)
library(Biostrings)
library(BSgenome)

#Finding Documentation---- in vignette form 
vignette()
vignette(package="Biostrings")
vignette("BiostringsQuickOverview")

browseVignettes()


#IRanges----
#Integer ranges, 1-based, from start to end inclusive - everything a vector in R so only plural IRanges
myiranges <- IRanges(start = c(5,20,25), end = c(10,30,40))
#all objects in bioconductor = S4 object - class can be obtained w class() - determines functions it can be used for and what slots it has for storing data 
#data stored in slots accessible w  @ - code using accessor function = more generic, less liable to break in future 
myiranges
class(myiranges)
methods(class="IRanges")
start(myiranges) #accessor functions
end(myiranges)
width(myiranges)
resize(myiranges, 3, fix = "start") #changed irange width to 3, start number stayed fixed
resize(myiranges, 3, fix="end") # changed irange width to 3, end number stayed fixed 

show_iranges <- function(obj) {
  for(i in seq_along(obj))
    cat(rep(" ", start(obj)[i]),
        rep("=", width(obj)[i]),
        "\n", sep="")
}
show_iranges( myiranges)
show_iranges( resize(myiranges, 3, fix="start"))
show_iranges( resize(myiranges, 3, fix = "end"))
show_iranges(flank(myiranges, 2, start = TRUE))
show_iranges(flank(myiranges, 2, start = FALSE))
show_iranges(range(myiranges))
show_iranges(reduce(myiranges))
show_iranges(disjoin(myiranges))
show_iranges(setdiff(range(myiranges), myiranges))
show_iranges(intersect(myiranges[2], myiranges[3]))
show_iranges(union(myiranges[2], myiranges[3]))


#GRanges---- refer to location within genome - need sequence name (chr), and whether strand is + or -
mygranges <- GRanges(
  seqnames = c("chrII", "chrI", "chrI"),
  ranges = IRanges(start = c(5, 20, 25), end = c(10,30,40)),
  strand = c("+", "-", "+")
)

mygranges
class(mygranges)
methods(class = "GRanges")
seqnames(mygranges)
strand(mygranges)
ranges(mygranges)
start(mygranges)
as.data.frame(mygranges)
# GRanges are like vectors:
mygranges[2]
c(mygranges, mygranges)
# Like other vectors, elements may be named
names(mygranges) <- c("foo", "bar", "baz")
mygranges
#Granges can have metadata columns  - a bit like a data frame
mygranges$wobble <- c(10,20,30)
mygranges
mcols(mygranges)
mygranges$wobble
# A handy way to create a GRanges
as("chrI:3-5:+", "GRanges")
# All the operations we saw for IRanges are available for GRanges.
# However most GRanges operations will take strand into account.
# Warning: shift() is a notable exception to this.

resize(mygranges, 3, fix="start")
resize(mygranges, 3, fix="end")





#DNAString / DNAStringSet---- 
#DNAString = like character string of bases - possibly IUPAC ambiguity codes such as “N” for any base.
myseq <- DNAString("gcgctgctggatgcgaccgcgcatgcgagcgcgacctatccggaa")
class(myseq)
methods(class = "DNAString")
reverseComplement(myseq)
translate(myseq)
subseq(myseq, 3, 5) # highlights bases 3 to 5
alphabetFrequency(myseq)
letterFrequency(myseq, c("AT", "GC"))
# DNAString objects behave like vectors.
myseq[3:5]
c(myseq, myseq)
# as() converts between types.
as(myseq, "character")
as("ACGT", "DNAString")

#FInding / counting patterns:
matchPattern("GCG", myseq)
countPattern("GCG", myseq)


#DNAStringSet - often want to work w collection of DNA sequences
myset <- DNAStringSet(list(chrI = myseq, chrII = DNAString("ACGTACGTAA")))
myset
class(myset)
elementNROWS(myset)
seqinfo(myset)
# Since a DNAString is like a vector, a DNAStringSet is has to be like a list.
myset$chrII
# Getting sequences with GRanges
getSeq(myset, mygranges)
getSeq(myset, as("chrI:1-3:+", "GRanges")) #on the + strand
getSeq(myset, as("chrI:1-3:-", "GRanges")) # on the - strand - complements above 
