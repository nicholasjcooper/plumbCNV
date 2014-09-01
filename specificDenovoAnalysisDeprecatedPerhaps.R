

specific.denovo.analysis <- function(dir,suffix=54) {
  print(load(cat.path(dir$res,"famstats_trios",suffix,ext="RData"))
  print(load(cat.path(dir$res,"TDT_results",suffix,ext="RData"))

  aff.dn <- (tt1[14,])/(tt1[22,])
  unaff.dn <- (tt1[8,]-tt1[14,])/(tt1[6,]-tt1[22,]) 
  unaff.n <- (tt1[8,]-tt1[14,])
  aff.n <- (tt1[14,])
  unaff.t <- (tt1[6,]-tt1[22,])
  aff.t <- (tt1[22,])
  dont <- (unaff.n==0 & aff.n==0) | is.infinite(aff.dn) | is.infinite(unaff.dn) | unaff.t==0 | unaff.n>unaff.t

  aff.n <- aff.n[!dont]; unaff.n <- unaff.n[!dont]
  aff.t <- aff.t[!dont]; unaff.t <- unaff.t[!dont]
  aff.dn <- aff.dn[!dont]; unaff.dn <- unaff.dn[!dont]

  top <- numeric();for (cc in 1:length(aff.n)) { print(cc); top[cc] <- FET(aff.n[cc],unaff.n[cc],case.d=aff.t[cc],cont.d=unaff.t[cc]) }

  or.dn <- aff.dn/unaff.dn
  mm <- cbind(aff.n,aff.t,aff.dn,unaff.n,unaff.t,unaff.dn,or.dn,top)
 
  ss <-sum.del[rownames(mm),c(1:4,8,16:19)]
  ss$genes <- substr(ss$genes, 1,14)
  return(cbind(ss,mm) )
}