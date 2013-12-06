####################################################################################################
#
#   Functional Data Analysis (FDA) - Case study on Diphthong/Hiatus contrast in European Spanish 
#
#	Michele Gubian (PhD)
#	Centre for Language and Speech Technology
#	Radboud University Nijmegen
#	email: m.gubian@let.ru.nl, michele.gubian@gmail.com
#	website on FDA: http://lands.let.ru.nl/FDA
#
#	Licensed under GPLv3 (See http://www.gnu.org/licenses/gpl-3.0.html)
#
#   This script allows one to reproduce the data analysis based on the output of FDA described in the paper:
#   "Using Functional Data Analysis for investigating multidimensional dynamic phonetic contrasts" by M. Gubian, F. Torreira and L. Boves. 
#   A draft version of the paper is available in the paper/ directory.
#
####################################################################################################


library(fda)
library(lattice)

root_dir = '/home/gubian/Experiments/FDA/lana/' # your root dir
plots_dir = paste(root_dir,'plots/',sep='')
data_dir =  paste(root_dir,'data/',sep='')
scripts_dir =  paste(root_dir,'scripts/',sep='')


DH_data = read.csv(file = paste(data_dir,"DH_data_FDA.csv",sep=''))
speakers = levels(DH_data$spk)

################## Analysis of features one by one ########################

# duration

vdur_class.lm = lm(vdur ~ class,  data = DH_data)
summary(vdur_class.lm) #R^2 = .45

vdur_spk.lm = lm(vdur ~ spk,  data = DH_data)
summary(vdur_spk.lm)  #R^2 = .09

# speaker-specific t-tests
# the alternative hypothesis is set according to the sign of the class coefficient in the lm
# p-values will be adjusted at the end for all features jointly
vdur_class.p = data.frame(spk=speakers,p=rep(NA,length(speakers)),row.names=speakers)
for (spk in vdur_class.p$spk) {
    vdur_class.p[spk,'p'] = t.test( vdur ~ class,  data = DH_data[which(DH_data$spk==spk),] , alternative = 'less'  )$p.value
}


# f0_s1

f0_s1_class.lm = lm(f0_s1 ~ class,  data = DH_data)
summary(f0_s1_class.lm) #R^2 =  .015

f0_s1_spk.lm = lm(f0_s1 ~ spk,  data = DH_data)
summary(f0_s1_spk.lm) #R^2 = .74

f0_s1_class.p =  data.frame(spk=speakers,p=rep(NA,length(speakers)),row.names=speakers)
for (spk in f0_s1_class.p$spk) {
    f0_s1_class.p[spk,'p'] = t.test( f0_s1 ~ class,  data = DH_data[which(DH_data$spk==spk),],alternative='greater')$p.value 
}

# f0_s2

f0_s2_class.lm = lm(f0_s2 ~ class,  data = DH_data)
summary(f0_s2_class.lm) # R^2 = 0.35

f0_s2_spk.lm = lm(f0_s2 ~ spk,  data = DH_data)
summary(f0_s2_spk.lm)$adj.r.squared #R^2 = .09



f0_s2_class.p =  data.frame(spk=speakers,p=rep(NA,length(speakers)),row.names=speakers)
for (spk in f0_s2_class.p$spk) {
    f0_s2_class.p[spk,'p'] = t.test( f0_s2 ~ class,  data = DH_data[which(DH_data$spk==spk),],alternative='less')$p.value 
}


# F12_s1

F12_s1_class.lm = lm(F12_s1 ~ class,  data = DH_data)
summary(F12_s1_class.lm) #R^2 = 0.34

F12_s1_spk.lm = lm(F12_s1 ~ spk,  data = DH_data)
summary(F12_s1_spk.lm)  #R^2 = 0.17


F12_s1_class.p =  data.frame(spk=speakers,p=rep(NA,length(speakers)),row.names=speakers)
for (spk in F12_s1_class.p$spk) {
    F12_s1_class.p[spk,'p'] = t.test( F12_s1 ~ class,  data = DH_data[which(DH_data$spk==spk),],alternative='less')$p.value 
}

# F12_s2

F12_s2_class.lm = lm(F12_s2 ~ class,  data = DH_data)
summary(F12_s2_class.lm)  #R^2 = 0.01696

F12_s2_spk.lm = lm(F12_s2 ~ spk,  data = DH_data)
summary(F12_s2_spk.lm) #R^2 = 0.40



F12_s2_class.p =  data.frame(spk=speakers,p=rep(NA,length(speakers)),row.names=speakers)
for (spk in F12_s2_class.p$spk) {
    F12_s2_class.p[spk,'p'] = t.test( F12_s2 ~ class,  data = DH_data[which(DH_data$spk==spk),],alternative='less')$p.value 
}


# adjust p-values for multiple comparisons

feature_class.p.adjust = data.frame(matrix(p.adjust(c(vdur_class.p$p,f0_s1_class.p$p,f0_s2_class.p$p,F12_s1_class.p$p,F12_s2_class.p$p)),ncol=5),row.names=speakers)
colnames(feature_class.p.adjust) = c("vdur","f0_s1","f0_s2","F12_s1","F12_s2")
feature_class.p.adjust = round(feature_class.p.adjust,3)
# latex-oriented formatting
feature_class.p.adjust[feature_class.p.adjust < 1e-3] = "$< 0.001$"
feature_class.p.adjust = cbind(speakers,feature_class.p.adjust)
 write(t(as.matrix(feature_class.p.adjust)),file="",ncolumns=6,sep = " & ")


################## Combined features analysis ########################

# correlation among (vdur, f0_s2, F12_s1)
library(Hmisc)
png(paste(plots_dir,'correlation_tree.png',sep=''))
plot(varclus(as.matrix(DH_data[,c('vdur','f0_s2','F12_s1')])),cex = 1.3) # ,labels=c('d',expression(s[2]^{f[0]}),expression(s[1]^{F[1-2]})) unfortunately expressions are not supported by plot.varclus
dev.off()
varclus(as.matrix(DH_data[,c('vdur','f0_s2','F12_s1')]))


library(languageR)
collin.fnc(DH_data,c(6,9,10))$cnumber 
# 15.6


# ordinary PCA on  (vdur, f0_s2, F12_s1)
vdur_f0_s2_F12_s1.pca = prcomp( ~  vdur + f0_s2 + F12_s1, data = DH_data, center=TRUE, scale=TRUE)
summary(vdur_f0_s2_F12_s1.pca)
as.data.frame(vdur_f0_s2_F12_s1.pca$rotation)
#           PC1        PC2        PC3
#vdur   -0.5940892  0.3652714  0.7166832
#f0_s2  -0.5451132 -0.8379965 -0.0247665
#F12_s1 -0.5915316  0.4053870 -0.6969590

# change the sign of PC1
DH_data = cbind(DH_data,vdur_f0_s2_F12_s1.pca$x)
DH_data$PC1 = - DH_data$PC1

# ordinary PCA on all five features
all5.pca =  prcomp( ~  vdur + f0_s1 + f0_s2 + F12_s1 + F12_s2, data = DH_data, center=TRUE, scale=TRUE)
summary(all5.pca)
as.data.frame(all5.pca$rotation)
# basically, PC1 is the same as above, PC4 and PC5 are very similar to resp. PC2 and PC3 above, and PC2 and PC3 are similar to each other and involve mainly f0_s1 and F12_s2. With some approximation, there is a subdivision of the 5-dimensional space into two subspaces, (vdur, f0_s2, F12_s1) and (f0_s1, F12_s2.). 

# generalised linear model on PCA scores

# see (Baayen, 2008)
library(Design)
DH_data.dd = datadist(DH_data)
options(datadist = 'DH_data.dd')

class_PCA123_int.lrm = lrm(class ~ (PC1 + PC2 + PC3)^2,  data = DH_data, x=T,y=T)
class_PCA123_int.lrm
# lots of not significant terms. Use bootstrap validation and backward elimination to eliminate variables
validate(class_PCA123_int.lrm,bw=TRUE, B=200)
# or just
fastbw(class_PCA123_int.lrm)
# survived varibles: PC1, PC3, PC1:PC3 (no intercept)
# use glm for convenience
summary(class_0PC13_int.glm <- glm(class ~ 0 +  PC1 * PC3,  data = DH_data, family = 'binomial'))

probs = binomial()$linkinv(fitted(class_0PC13_int.glm))
somers2(probs, as.numeric(DH_data$class) - 1)
# Dxy = 0.86


# center and scale features
DH_data$vdur_CS = scale(DH_data$vdur) 
DH_data$f0_s1_CS = scale(DH_data$f0_s1) 
DH_data$f0_s2_CS = scale(DH_data$f0_s2) 
DH_data$F12_s1_CS = scale(DH_data$F12_s1) 
DH_data$F12_s2_CS = scale(DH_data$F12_s2) 

# scatterplot of centered and scaled vdur and F12_s1, by speaker
png(paste(plots_dir,'scatter_F12_s1_CS_vdur_CS_speaker.png',sep=''))
xyp = xyplot(F12_s1_CS ~ vdur_CS | spk,  data = DH_data,xlab = 'norm. d',ylab=expression(paste('norm.',s[1]^{F[1-2]})), groups= class, 
pch  = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
col =  sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
scales = list(cex=1.5),
#panel = panel.superpose,
#panel.groups = function(...) {
panel = function(...) {
    panel.xyplot(...)
    panel.abline(h=0,lty=2,col='grey')
    panel.abline(v=0,lty=2,col='grey')
    }
)
update	(xyp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)
	    ),
        as.table=TRUE        
    )
dev.off()

# scatterplot of centered and scaled vdur and F12_s1, global
png(paste(plots_dir,'scatter_F12_s1_CS_vdur_CS.png',sep=''))
xyp = xyplot(F12_s1_CS ~ vdur_CS,  data = DH_data,xlab = 'norm. d',ylab=expression(paste('norm.',s[1]^{F[1-2]})), groups= class, 
pch  = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
col =  sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
scales = list(cex=1.5),
panel = function(...) {
    panel.xyplot(...)
    panel.abline(h=0,lty=2,col='grey')
    panel.abline(v=0,lty=2,col='grey')
    }
)
update	(xyp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)
	    ),
        as.table=TRUE        
    )
dev.off()


# plot predicted probability of H  from model class_0PC13_int.glm, translated back into vdur_CS F12_s1_CS coordinates using PCA loadings
cRP = colorRampPalette(c(color[['d']],"white",color[['h']]))
grid = expand.grid( vdur_CS     = do.breaks(c(-3,3),60),
                    F12_s1_CS   = do.breaks(c(-3,3),60))
grid$logit =   with(grid, 
    class_0PC13_int.glm$coefficients['PC1'] * (-1) * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC1'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC1'] * F12_s1_CS) + # (-1)* is because of the manual sign swap operated on PC1
    class_0PC13_int.glm$coefficients['PC3'] * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC3'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC3'] * F12_s1_CS) 
  +class_0PC13_int.glm$coefficients['PC1:PC3']  * (-1) * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC1'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC1'] * F12_s1_CS) * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC3'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC3'] * F12_s1_CS)
    )
grid$prob = with(grid, 1/(1 + exp(-logit)) ) # inverse logit

png(paste(plots_dir,'prob_vdur_CS_F12_s1_CS.png',sep=''))
levelplot(prob~vdur_CS*F12_s1_CS,data=grid,at = do.breaks(c(0,1), 30), col.regions = cRP,
    panel = function(...) {
    panel.levelplot(...)
    panel.points(x=DH_data$vdur_CS,y=DH_data$F12_s1_CS,col='white',cex=1.5,
        pch = sapply(DH_data$class, function(x) symbol[[x]],USE.NAMES = FALSE)
        )
    },
    xlab = list(label='norm. d',cex=1.5), ylab=list(label=expression(paste('norm.',s[1]^{F[1-2]})),cex=1.5),
    scales = list(cex=1.5),
    par.settings=list(
        par.xlab.text = list(cex=1.3),
        par.ylab.text = list(cex=1.3),
        axis.text = list(cex=1.3)
        ),
    colorkey = list(labels=list(cex=1.3)),
) 
dev.off()





# TREES
library(party)

ct  = list()
for (spk in speakers) {
    ct[[spk]] = ctree(class ~  vdur + f0_s2 + F12_s1, data = DH_data[which(DH_data$spk == spk),])
    #png(paste(plots_dir,'ctree_vdur_f0_s2_F12_s1_',spk,'.png',sep=''))
    #plot( ct[[spk]] )
    #dev.off()
    print (paste( spk, sum( DH_data$class[which(DH_data$spk == spk)] !=  Predict(ct[[spk]]) ), length(Predict(ct[[spk]])) ,   
                       sum( DH_data$class[which(DH_data$spk == spk)] !=  Predict(ct[[spk]]) )/ length(Predict(ct[[spk]])) ,
                       ct[[spk]]@tree$psplit$variableName, ct[[spk]]@tree$psplit$splitpoint        ))
}

# use all 5 features
ct5 = list()
for (spk in speakers) {
    ct5[[spk]] = ctree(class ~  vdur + f0_s1 + f0_s2 + F12_s1 + F12_s2, data = DH_data[which(DH_data$spk == spk),])
    png(paste(plots_dir,'ctree_all5_',spk,'.png',sep=''))
    plot( ct5[[spk]] )
    dev.off()
    print (paste( spk, sum( DH_data$class[which(DH_data$spk == spk)] !=  Predict(ct5[[spk]]) ), length(Predict(ct[[spk]])) ,   
                       sum( DH_data$class[which(DH_data$spk == spk)] !=  Predict(ct5[[spk]]) )/ length(Predict(ct[[spk]])) ,
                       ct5[[spk]]@tree$psplit$variableName, ct5[[spk]]@tree$psplit$splitpoint        ))
}
# no difference in prediction, nor in splits that determine a class subdivision.

# global tree for all speakers
CT =  ctree(class ~  vdur + f0_s2 + F12_s1 , data = DH_data)

png(paste(plots_dir,'ctree_vdur_f0_s2_F12_s1_all_speakers.png',sep=''))
plot(CT)
dev.off()
print (paste( 'all speakers' , sum( DH_data$class !=  Predict(CT) ), length(Predict(CT)) ,   
                       sum( DH_data$class !=  Predict(CT) )/ length(Predict(CT)) ))



