# O2PLS integration of RNA-seq and cell-surface proteomics 

####################################################
pkgs <- c("dplyr", "OmicsPLS", "readxl", "writexl")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(dplyr)
library(OmicsPLS)
library(readxl)
library(writexl)



dot = function(a, b) sum(a * b)  
run_cv  <- FALSE               # TRUE to rerun cross-validation (slow)
#load rna and proteome datasets
rna =  as.data.frame(read_excel("rna.xlsx"))   
cp =  as.data.frame(read_excel("cp.xlsx"))
names(rna)[1] = "GeneName"                     
names(cp)[1] = "ID"


#re-shape the datasets into o2pls-compliant dimensions
rnan = rna[,-1]
cpn = cp[,-1]
rna_names = rna[,1]
cp_names = cp[,1]
rnaT = t(rnan)
cpT = t(cpn)
colnames(rnaT) = rna_names                     
colnames(cpT) = cp_names
n = nrow(rnaT)                                 



#Use OmicsPLS to integrate the two datasets.
rnaTS = scale(rnaT, center = TRUE, scale = TRUE)
cpTS = scale(cpT, center = TRUE, scale = TRUE)

if (run_cv) {
  crossval_o2m_adjR2(rnaTS, cpTS, 1:3, 0:3, 0:3, nr_folds = 5)
  crossval_o2m(rnaTS, cpTS, 1:2, 0:2, 0:2, 5)
}

fit = o2m(rnaTS, cpTS, 2,1,1)


#Extract scores and loadings

loadings = loadings(fit)
LX = loadings(fit,"Xjoint")  #joint loading values of RNA-seq 
LY = loadings(fit,"Yjoint")  #joint loading values of proteomics 

LX = data.frame(GeneName = rna_names, V1gene = LX[,1], V2gene = LX[,2])
LY = data.frame(ID = cp_names, `Gene Names` = cp_names, V1cp = LY[,1], V2cp = LY[,2], check.names = FALSE)
write_xlsx(LX, "lv_rnajnt.xlsx")
write_xlsx(LY, "lv_cpjnt.xlsx")
loadings_xorth = loadings(fit,"Xorth")
write_xlsx(data.frame(GeneName = rna_names, orth1 = loadings_xorth[,1]),"lv_rnaorth.xlsx")
loadings_yorth = loadings(fit,"Yorth")
write_xlsx(data.frame(ID = cp_names, orth1 = loadings_yorth[,1]),"lv_cpOrth.xlsx")
scores = scores(fit)
T = scores(fit,"Xjoint") #joint scores of RNA-seq dataset
U = scores(fit,"Yjoint")


samples = rownames(rnaTS)
scores_xorth = scores(fit,"Xorth")
write_xlsx(data.frame(sample = samples, to1 = scores_xorth[,1]),"scores_rnaorth.xlsx")
scores_yorth = scores(fit,"Yorth")
write_xlsx(data.frame(sample = samples, uo1 = scores_yorth[,1]),"scores_cporth.xlsx")

scores_xjnt = scores(fit,"Xjoint")            
write_xlsx(data.frame(sample = samples, t1 = scores_xjnt[,1], t2 = scores_xjnt[,2]),"scores_rnajnt.xlsx")
scores_yjnt = scores(fit,"Yjoint")
write_xlsx(data.frame(sample = samples, u1 = scores_yjnt[,1], u2 = scores_yjnt[,2]),"scores_cpjnt.xlsx")


#Prepare the integrated data for target identification
cov1 = inner_join(LX, LY, by =c("GeneName" = "Gene Names"))

#find correlations
T1 = T[,1]                                     
T2 = T[,2]
U1 = U[,1]
U2 = U[,2]
Cr = apply(cov1, 1, function(row){
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  lvy <- matrix( as.numeric(c(row["V1cp"], row["V2cp"])) , nrow = 2, ncol = 1)
  jointx<-T1*lvx[1]+ T2*lvx[2]
  jointy<-U1*lvy[1]+U2*lvy[2]
  covxy = dot(jointx, jointy) #dot product of jointx and jointy
  varx = dot(jointx, jointx)
  vary = dot(jointy, jointy)
  as.numeric(covxy/sqrt(varx*vary))
})

covcorall = data.frame(gene = cov1$GeneName, ID = cov1$ID, Corr_xyJnt = Cr)


loadings = cov1
lvy1 = loadings[,"V1cp"]                       
lvy2 = loadings[,"V2cp"]
length(lvy1)

#finding joint component of ADAM9
lvy1A9 = as.numeric(cov1[cov1$GeneName== "ADAM9","V1cp"])[1]
lvy2A9 = as.numeric(cov1[cov1$GeneName == "ADAM9","V2cp"])[1]

jntA9 = U1*lvy1A9 + U2*lvy2A9
varA9 = dot(jntA9,jntA9)
jointsy = matrix(nrow = length(lvy1), ncol = n)
cov_y_a9 = matrix(nrow = length(lvy1), ncol = 1)
cor_y_a9 = matrix(nrow = length(lvy1), ncol = 1)

for(i in 1:length(lvy1)){
  #lvy = c(lvy1[i], lvy2[i])
  jnti = U1*lvy1[i] + U2*lvy2[i]
  jointsy[i,] <- jnti
  cova9i = dot(jnti, jntA9)
  cov_y_a9[i,] <- cova9i
  
  vari = dot(jnti,jnti)
  cor_ya9i = cova9i/sqrt(vari*varA9)
  cor_y_a9[i,] <- cor_ya9i
}
cov1 = cbind(cov1,cov_y_a9 )
covcorA9 = cbind(cov1,cor_y_a9 )


#find potential transcriptional targets
bh1.xy = c()
pp1.xy = c()
for (i in 1:nrow(covcorall)){
  r = as.numeric(trimws(covcorall$Corr_xyJnt[i]))
  f =  0.5*log((r+1)/(1-r))
  z = sqrt(n-3)*f
  
  #right-handed p
  p = pnorm(z, lower.tail = FALSE)
  pp1.xy = append(pp1.xy,p)
  
}

bh1.xy = p.adjust(pp1.xy, method = "BH")
cc3 = cbind(covcorall, bh1.xy)
write_xlsx(cc3, "alltargets_bh1xy.xlsx")
cc4 = cc3[cc3$bh1.xy < 0.05,]
write_xlsx(cc4, "tt.xlsx")

sigtt = inner_join(cc4, cp, by ="ID")

if ("FDR" %in% names(sigtt)) sigtt = sigtt[sigtt$FDR < 0.05,] else warning("No FDR column: ttDEP.xlsx is not filtered by differential expression")
write_xlsx(sigtt,"ttDEP.xlsx")

#find potential post-transcriptional targets
ccneg = cc3[cc3$Corr_xyJnt <0,]
ccposns = cc3[cc3$bh1.xy > 0.1,]
ptt = distinct(rbind(ccneg, ccposns))         
dir.create("analysis after integration", showWarnings = FALSE)   
write_xlsx(ptt,"analysis after integration/ptt.xlsx")
sigptt = ptt %>% inner_join(cp,by = "ID")      
if ("FDR" %in% names(sigptt)) sigptt = sigptt %>% filter(FDR < 0.05) else warning("No FDR column: pttDEP.xlsx is not filtered by differential expression")
write_xlsx(sigptt,"pttDEP.xlsx")

#correlation with the sheddase
Yj_lv <- as.data.frame(read_excel("lv_cpjnt.xlsx"))
Yj_sc <- as.data.frame(read_excel("scores_cpjnt.xlsx"))

Yo_lv <- read_excel("lv_cpOrth.xlsx")
Yo_lv = as.data.frame(Yo_lv)
Yo_sc <- as.data.frame(read_excel("scores_cporth.xlsx"))


Yj <- as.matrix(Yj_sc[,c(2,3)]) %*% t(as.matrix(Yj_lv[,c("V1cp","V2cp")]))
Yo <- as.matrix(Yo_sc[,2]) %*% t(as.matrix(Yo_lv[,"orth1"]))

#Yj and Yo are the joint and orthogonal components
#add them to get predicted values

Ypred = Yj + Yo


Yp3 = data.frame(ID = Yj_lv$ID, `Gene Names` = Yj_lv$`Gene Names`, t(Ypred), check.names = FALSE)
pred_cols = paste0("pred_", Yj_sc$sample)
names(Yp3)[-c(1, 2)] = pred_cols
corr.y.A9 = c()
corr.pvalue = c()
a9 = which(Yp3$`Gene Names` == "ADAM9")[1]    
A9 = Yp3[a9,]
A9data = as.numeric(A9[, c(-1, -2)])

for(i in 1:nrow(Yp3)){
  idata = as.numeric(Yp3[i, c(-1, -2)])
  cori = cor.test(A9data, idata, method = "pearson")
  corr.y.A9 = append(corr.y.A9, cori$estimate)
  corr.pvalue = append(corr.pvalue, cori$p.value)
  
}

Yp4 = cbind(Yp3, corr.y.A9,corr.pvalue )


#categorize
#post-transcriptional targets

expr_cols = c(names(cp)[-1], pred_cols)        # control/knockdown and predicted values
ptt <- read_excel("pttDEP.xlsx")
ptt_pred = inner_join(ptt,Yp4, by = "ID" )
ptt_pred2 = ptt_pred[, !(names(ptt_pred) %in% expr_cols)] #remove the control and knockdown expression values
ptt_pred2 = ptt_pred2 %>% rename(corr.pvalue.ya9 = corr.pvalue)

pttneg = ptt_pred2[ptt_pred2$corr.y.A9 < 0,]
pttnegsig = pttneg[pttneg$corr.pvalue.ya9 < 0.05,]
pttneg$corrBH = p.adjust(pttneg$corr.pvalue.ya9, method = "BH")
pttnegsigBH = pttneg[pttneg$corrBH< 0.05,]
write_xlsx(pttnegsigBH,"pttnegsigBH.xlsx")

pttpos = ptt_pred2[ptt_pred2$corr.y.A9 > 0,]
pttpossig = pttpos[pttpos$corr.pvalue.ya9 < 0.05,]
pttpos$corrBH = p.adjust(pttpos$corr.pvalue.ya9, method = "BH")
pttpossigBH = pttpos[pttpos$corrBH< 0.05,]
write_xlsx(pttpossigBH,"pttpossigBH.xlsx")


#transcriptional targets

tt <- as.data.frame(read_excel("ttDEP.xlsx"))
tt_pred = inner_join(tt, Yp4, by = c("ID"))
tt_pred2 = tt_pred[, !(names(tt_pred) %in% expr_cols)]   

ttneg = tt_pred2[tt_pred2$corr.y.A9 < 0,]
ttnegsig = ttneg[ttneg$corr.pvalue < 0.05,]
ttneg$corrBH = p.adjust(ttneg$corr.pvalue, method = "BH")
ttnegsigBH = ttneg[ttneg$corrBH< 0.05,]
write_xlsx(ttnegsigBH, "ttnegsigBH.xlsx")

ttpos = tt_pred2[tt_pred2$corr.y.A9 > 0,]
ttpossig = ttpos[ttpos$corr.pvalue < 0.05,]
ttpos$corrBH = p.adjust(ttpos$corr.pvalue, method = "BH")
ttpossigBH = ttpos[ttpos$corrBH < 0.05,]       
write_xlsx(ttpossigBH,"ttpossigBH.xlsx")