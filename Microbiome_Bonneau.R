library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(phyloseq)


################################################################
###  This is to download and filter the original file AG1.biom
###  Did not work-- but shows input filtering FYI
################################################################

url <- "https://github.com/zdk123/SpiecEasi/blob/master/inst/extdata/study_1925_closed_reference_otu_table.biom"

download.file(url, desfile='AG1.biom' method='wget') ## download method is system dependent




amgut1 <- import_biom('AG1.biom')

## filter taxa that appear in fewer than 30% of the samples

amgut1.filt <- prune_taxa(rowSums(sign(amgut1@otu_table@.Data)) > nsamples(amgut1)*.3, amgut1)




ig.mb <- graph.adjacency(se.mb.amgut$refit, mode='undirected')

vsize <- exp(colMeans(se.mb.amgut$data)/2)+2




## filter edges that we are less confidence of (i.e. by StARS 'merge' estimate)

highconfnet <- se.mb.amgut$refit * (se.mb.amgut$merge[[se.mb.amgut$opt.index]] > .9)

ig.hc <- graph.adjacency(highconfnet, mode='undirected')

am.coord <- layout.fruchterman.reingold(ig.mb)




## plot

par(mfrow=c(1,2))

plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")

plot(ig.hc,  layout=am.coord,  vertex.size=vsize, vertex.label=NA, main="High Confidence edges")

###########################################################
###########################################################
###########################################################

##  This part works
##  Simulated Data-- This takes the study data and does simulations based on the marginals
##  Apparently useful to assess the fit of the actual data below-- I think


data(amgut1.filt)   ##  this is apparently derived from the code above
depths <- rowSums(amgut1.filt)
amgut1.filt.n <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, distr='zinegbin', Sigma=Cor, n=n)

se.est <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)

huge::huge.roc(se.est$path, graph, verbose=FALSE)
stars.pr(getOptMerge(se.est), graph, verbose=FALSE)
# stars selected final network under: se.est$refit


##################################################################################
###  This is the run
###  Using actual data

se.mb.amgut <- spiec.easi(amgut1.filt, method='mb', lambda.min.ratio=1e-2, 
                            nlambda=20, icov.select.params=list(rep.num=50))
se.gl.amgut <- spiec.easi(amgut1.filt, method='glasso', lambda.min.ratio=1e-2,
                            nlambda=20, icov.select.params=list(rep.num=50))
sparcc.amgut <- sparcc(amgut1.filt)
## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc.amgut$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
## Create igraph objects
ig.mb <- graph.adjacency(se.mb.amgut$refit, mode='undirected')
ig.gl <- graph.adjacency(se.gl.amgut$refit, mode='undirected')
ig.sparcc <- graph.adjacency(sparcc.graph, mode='undirected')



## set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

elist.gl <- summary(Matrix::triu(se.gl.amgut$opt.cov*se.gl.amgut$refit, k=1))
elist.mb <- summary(symBeta(getOptBeta(se.mb.amgut), mode='maxabs'))
elist.sparcc <- summary(sparcc.graph*sparcc.amgut$Cor)

hist(elist.sparcc[,3], main="", xlab="edge weights")
hist(elist.gl[,3], add=TRUE, col='red', alpha=.5)
hist(elist.mb[,3], add=TRUE, col='forestgreen', alpha=.4)

dd.gl <- degree.distribution(ig.gl)
dd.mb <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)

plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.35), type='b', 
      ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"), 
        col=c("forestgreen", "red", "black"), pch=1, lty=1)