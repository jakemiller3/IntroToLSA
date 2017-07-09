# IC2S2 LSA Workshop Code
# 2017-07-10
# Jacob Miller, Drexel University: jlm479@drexel.edu 
#    with assistance from Jorge Fresneda, David Gefen, James Endicott, Kai Larsen, and others in the R community

###########################
# Section 1: Set up R.
# If you have never loaded these packages in this R environment, you need to do this to download them.
# tm: Text Mining; 
# LSAfun: LSA functions; 
# Matrix: Sparse matrix package
# RSpectra: Fast SVD
# gplots: fancy plotting
install.packages(c("tm", "LSAfun", "Matrix", "RSpectra", "gplots"), dependencies=TRUE)

library(tm) # This loads them into your environment so you can now use that code.
library(LSAfun)
library(Matrix)
library(RSpectra) 
library(gplots)

#
#
###########################
# Section 2: Set up your directory and test the DirSource command
# machFileDirectory = 'C:/Users/username/Desktop/files' # Windows style
# machFileDirectory = '~/Desktop/files' # Mac/Unix style
machFileDirectory = 'C:/Users/Owner/Desktop/Dropbox/LSA Workshop/files'
doc_source <- DirSource(machFileDirectory, recursive = TRUE)
object.size(doc_source)


#
#
###########################
# Section 3: Load the text files
# 
# We shall now create a corpus in memory, using functions from the tm package.

raw_corpus <- VCorpus(doc_source, readerControl=list(language='en'))

omitDocNames = paste(c(452, 454, 456, 458, 460:470, 518), ".txt", sep="")
omitDocs = which(names(raw_corpus) %in% omitDocNames)
# For non-LSA exploration - unweighted tdm, leaving stopwords
machRawtdm <- TermDocumentMatrix(raw_corpus[-omitDocs], 
	control = list(removePunctuation = TRUE, 
	bounds=list(global=c(2,Inf))))
	
machtdm <- TermDocumentMatrix(raw_corpus[-omitDocs], 
	control = list(removePunctuation = TRUE, 
	# stemming=TRUE,
	stopwords = TRUE, 
	weighting = function(x)  weightTfIdf(x, normalize = FALSE),
	bounds=list(global=c(2,Inf))))
machtdm
# The tdm matrix is very sparse - 99% of cells are zero.

inspect(machtdm[1:20,1:8])

# Convert the slam matrix to a Matrix
sparse_tdm <- Matrix::sparseMatrix(i = machtdm$i, j = machtdm$j, x = machtdm$v, dims = c(machtdm$nrow, machtdm$ncol))
dimnames(sparse_tdm) <- dimnames(machtdm)
sparse_tdm[1:20,1:8]

#
#
###########################
# Section 4: Examine the DTM
#
class(machtdm)
class(sparse_tdm)

machtdm[1:20,1:8]
sparse_tdm[1:20,1:8]

dim(sparse_tdm) # What are the dimensions of the matrix? Rows, Columns
object.size(sparse_tdm) # How many bytes?
colnames(sparse_tdm) # What are the column names?


###########################
# Section 5: Run lsa
#

start=Sys.time() # This is a simple timing mechanism. Creates a variable that saves the current time and...
	space <- svds(sparse_tdm, 100) 
	# 	easyspace = lsa(machtdm, dims=100) # This will give a warning - last dimension calculated is very close to zero size
Sys.time()-start # ... this line shows the now-current time minus the "start" time.

class(space) # A list
str(space) #          with five objects inside it
object.size(space)

tk <- space$d * space$u
dk <- space$d * Matrix::t(space$v)
#Assign names
dimnames(tk) <- list(dimnames(sparse_tdm)[[1]], 1:100)   
dimnames(dk) <- list(1:100, dimnames(sparse_tdm)[[2]])  


dim(tk) # Dimensions of term matrix
dim(dk) # Dimensions of doc matrix
length(space$d) # Number of singular values

tk[1:10,1:10]




###########################
# Section A1: Initial investigation

findFreqTerms(machRawtdm, lowfreq=100) # Words which occur 100+ times

# Find closest terms to a focal term

neighbors("trust", 30, tvectors=tk[,1:100])
neighbors("trust", 30, tvectors=tk[,1:30])
#####
# Section A2:
# plot_neighbors find the nearest neighbors and then runs a PCA on the vectors in this set, and plots them
plot_neighbors("trust", 30, dims=2, tvectors=tk, breakdown=F)
# quartz.save("~/Dropbox/Projects/LSA Workshop/graphics/Mach_trust_plotn_72dpi.png", dpi=72)

# You can rotate the 3D one!
plot_neighbors("trust", 30, dims=3, tvectors=tk, breakdown=F)

###########################
# Section B: Cosines

# Need to zero-pad the names to make the next analysis look right, so we can sort it...
charLength = lapply(colnames(dk), nchar) 
colnames(dk)[which(charLength==5)] = paste("00", colnames(dk)[which(charLength==5)], sep="")
colnames(dk)[which(charLength==6)] = paste("0", colnames(dk)[which(charLength==6)], sep="")

# This runs cosines between all paragraphs, returns a matrix of cosines
start=Sys.time()
paraCos = multicos(sort(colnames(dk)), tvectors=t(dk), breakdown=F)
Sys.time()-start

paraCos[1:10,1:10] # Look at first 10x10 of matrix.

# Save your cosine matrix as a CSV for other apps!
write.csv(paraCos, file="C:/Users/Owner/Desktop/Dropbox/LSA Workshop/Mach_paraCos.csv")

###########################
# Section C1: Heatmaps


# Let's see this, and with a hierarchical cluster.
heatmap(paraCos, symm=T, main="Machiavelli Clustered Paras")

# Now, let's see this in paragraph order. Anything becoming clear?
heatmap(paraCos, symm=T, Rowv=NA,  main="Machiavelli Ordered Paras")
# quartz.save("~/Dropbox/Projects/LSA Workshop/graphics/MachMap_72dpi.png", dpi=72) # Save the image...

#####
# Section C2:
termCos = multicos(findFreqTerms(machRawtdm, lowfreq=100, highfreq=300), tvectors= tk[,1:50], breakdown=F)
heatmap(termCos, symm=T, main="Machiavelli Clustered Frequent Terms")

heatmap.2(termCos, symm=T, main="Heatmap 2: Machiavelli Clustered Frequent Terms")

# quartz.save("~/Dropbox/Projects/LSA Workshop/graphics/Mach_TermMap_50D_72dpi.png", dpi=72)

###########################
# Section D: Some extras

# Simple plot of the singular values
plot(space$d)
rTerms = sample(rownames(tk), 30) # random sample of term names
plot(tk[rTerms,1],tk[rTerms,2])
text(tk[rTerms,1],tk[rTerms,2], labels=rTerms)

# Combine strings into a single vector and find cosines.
costring("trust", "useful perfect", tvectors= tk[,1:50], breakdown=TRUE)
costring("trust queen", "useful sin", tvectors= tk[,1:50], breakdown=TRUE)
costring("trust queen", "elizabeth judgment", tvectors= tk[,1:50], breakdown=TRUE)
costring("trust", "prince souldiours", tvectors= tk[,1:50], breakdown=TRUE)
costring("trust", "prince fear", tvectors= tk[,1:50], breakdown=TRUE)
costring("trust", "prince love", tvectors= tk[,1:50], breakdown=TRUE)
costring("trust", "prince pikes", tvectors= tk[,1:50], breakdown=TRUE)
costring("power", "prince pope", tvectors= tk[,1:50], breakdown=TRUE)
costring("power", "prince souldiours", tvectors= tk[,1:50], breakdown=TRUE)
costring("power", "prince", tvectors= tk[,1:50], breakdown=TRUE)
costring("power", "pope", tvectors= tk[,1:50], breakdown=TRUE)
costring("power", "state fear", tvectors= tk[,1:50], breakdown=TRUE)
costring("power", "state love", tvectors= tk[,1:50], breakdown=TRUE)


termlist = c("souldiours", "power", "prince", "princes", "state", "warre", "pikes", "enemies", "trust", "love", "fear", "pope", "alexander") # Create a list of specified terms to use for particular investigation.
specificTermCos = multicos(termlist, tvectors= tk[,1:50], breakdown=F)
specificTermCos
heatmap(specificTermCos, symm=T,  main="Machiavelli Specific Terms")

plot_neighbors("power", 30, dims=3, tvectors=tk, breakdown=F)


Term_count <-apply(tk,1,sum)
TCT <- t(Term_count)
myTerms <- rownames(tk)