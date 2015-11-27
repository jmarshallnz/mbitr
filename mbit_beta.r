#' code for generating self-attribution for a beta model.
#' 
#' Hopefully later on we can generalise this to other models.

# grab a model matrix and response
mod.matrix <- function(formula, data) {
  if (!('data.frame' %in% class(data))) stop('data parameter must be of class "data.frame"')
  if (!('formula' %in% class(formula))) stop('formula parameter must be of class "formula"')
  mf = model.frame(formula, data = data)
  tf = terms(mf)
  attr(tf, "intercept") <- 0
  X = model.matrix(tf, data=data)
  y = model.response(mf)
  list(X = X, y = y)
}

fit_beta <- function(y, X) {
  # pretty simple: For each level of Y fit a beta distribution to each column of X
  ns = length(levels(y))
  k = table(mm$y)
  ng = ncol(X)
  p_alpha = matrix(NA, ns, ng)
  p_beta  = matrix(NA, ns, ng)
  rownames(p_alpha) = rownames(p_beta) = levels(y)
  colnames(p_alpha) = colnames(p_beta) = colnames(X)
  for (g in 1:ng) {
    n = tapply(X[,g], y, sum)
    p_alpha[,g] = n
    p_beta[,g]  = k - n
  }
  list(alpha=p_alpha, beta=p_beta)
}

predict_beta <- function(mod, x) {
  # TODO: this assumes this is a model object fit by fit_beta
  # TODO: ATM we spit out the mean prediction
  # mean of a beta dist is alpha / (alpha + beta)
  mu = mod$alpha / (mod$alpha + mod$beta)

  # Predictive distribution for Bernoulli/Beta is just the mean of the Beta posterior
  t(exp(log(mu) %*% x + log(1-mu) %*% (1-x)))
}

## do the stupid thing first. This will fit the same model to any repeat types multiple times
library(readxl)
library(dplyr)
library(RColorBrewer)

d = read_excel("MBITs27-11-15.xlsx")[1:22]

# filter the data
d = d %>% filter(!is.na(Group2), Group2 != "Water", Group2 != "Human") # Group2 != "Herbivore")

# fixup!
d$Group2[d$TypeDetails == "Bovine"] <- "Ruminant"
d$Group2 = as.factor(d$Group2)

d <- data.frame(d)
rownames(d) <- d$MBiT


mm = mod.matrix(Group2 ~ . - Key - MBiT - TypeDetails, data=d)
n = length(mm$y)
p = numeric(n)
ps = matrix(NA, n, length(levels(mm$y)))
for (i in 1:n) {
  # fit beta model
  mod_fit = fit_beta(mm$y[-i], mm$X[-i,])
  # predict model
  mod_pred = predict_beta(mod_fit, mm$X[i,])
  # did we get it right? This assumes apriori that all sources are equally likely for the isolate
  p[i] = mod_pred[mm$y[i]] / sum(mod_pred)
  ps[i,] = mod_pred / sum(mod_pred)
}
colnames(ps) <- levels(mm$y)
boxplot(p ~ mm$y)
mean(p)
most_likely = factor(apply(ps, 1, which.max), labels=levels(mm$y))
table(most_likely, mm$y)

l = by(ps, mm$y, FUN=function(x) { x }, simplify=FALSE)
for (i in seq_along(l)) {
  boxplot(l[[i]], main=names(l)[i])
  print(colMeans(l[[i]]))
}

l = by(ps, mm$y, FUN=function(x) { x }, simplify=FALSE)

pal = brewer.pal(n=9, name="GnBu")
pdf("foo.pdf", width=11, height=8)
par(mfrow=c(3,1), mar=c(5,7,2,1))
for (i in seq_along(l)) {
  o = order.dendrogram(as.dendrogram(hclust(dist(as.matrix(l[[i]])))))
  image(as.matrix(l[[i]])[o,], col=pal, xaxt="n", yaxt="n")
  mtext(levels(mm$y), side=2, at=seq(0,1,length.out=length(levels(mm$y))), las=1, line=0.5)
  # TODO: need to runlength encode this...
  labs = d$MBiT[as.numeric(rownames(l[[i]]))[o]]
  pos  = seq(0,1,length.out=nrow(l[[i]]))
  mtext(labs, side=1, at=pos, las=2, line=0.5, cex=0.3)
  title(main=names(l)[i])
}
dev.off()
