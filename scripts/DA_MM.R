# This script contains the data analysis performed unsing linear mixed models.

library(fda)
library(lattice)
library(lme4)
#library(nlme)

root_dir = '/home/gubian/Experiments/FDA/lana/' # your root dir
plots_dir = paste(root_dir,'plots/',sep='')
data_dir =  paste(root_dir,'data/',sep='')
scripts_dir =  paste(root_dir,'scripts/',sep='')


DH_data = read.csv(file = paste(data_dir,"DH_data_FDA.csv",sep=''))
speakers = levels(DH_data$spk)



# graphic parameters
color = list ( 'd' = 'blue','h' = 'orangered')
symbol = list( 'd' = 'D','h' = 'H')
lty = list( 'd' =1,'h' = 3)
lwd  = list( 'd' =1,'h' = 2)
name = list( 'd' ='Diphthong','h' = 'Hiatus')
strip.background.col = trellis.par.get('strip.background')$col[4]
dev.off()






# write.csv(DH_data, file = "DH_data.csv",row.names = FALSE)
DH_data = read.csv(file = "DH_data.csv")

#> colnames(DH_data)
# [1] "filename"  "spk"       "male"      "class"     "ldur"      "vdur"     
# [7] "dur"       "f0_s1"     "f0_s2"     "F12_s1"    "F12_s2"    "pca_s1"   
#[13] "pca_s2"    "pca_s3"    "vdur_CS"   "f0_s1_CS"  "f0_s2_CS"  "F12_s1_CS"
#[19] "F12_s2_CS"


######### Analysis of duration

# vowel sequence duration vdur (d in the text)

# full mixed model:
summary(vdur_class_1classspk.lmer <- lmer(vdur ~ class + (1 + class | spk), data = DH_data))
# null model, only random effects
summary(vdur_1_1spk.lmer <- lmer(vdur ~ 1 + (1 | spk), data = DH_data))
# model with class only as random effect
summary(vdur_1_1classspk.lmer <- lmer(vdur ~ 1 + (1 + class | spk), data = DH_data))
# model with class only as fixed effect
summary(vdur_class_1spk.lmer <- lmer(vdur ~ class + (1 | spk), data = DH_data))

# series of log-likelihood tests to justify the model complexity (model selection)
anova(vdur_class_1classspk.lmer,vdur_class_1spk.lmer) # random slope justified
anova(vdur_class_1classspk.lmer,vdur_1_1classspk.lmer) # fixed effect class justified
# similar tests starting from the null model
anova(vdur_1_1classspk.lmer,vdur_1_1spk.lmer)
anova(vdur_class_1spk.lmer,vdur_1_1spk.lmer)

# Gain in explained variance (a pseudo R^2, see book "Analyzing Linguistic Data: A Practical Introduction to Statistics using R" by Baayen)
# based on the correlation  between fitted and actual values of the dependent variable.
# The comparison is between the full model and the null model. 
vdur_class_1classspk.pseudoR2 = 1 - (cor(fitted(vdur_1_1spk.lmer),DH_data$vdur) / cor(fitted(vdur_class_1classspk.lmer),DH_data$vdur) )^2

# Reduction in std. dev. due to the introduction of class as fixed effect.
# This time the comparison is between the full model and the model without fixed factor (see Baayen book for details).

# Extract std. dev. of the random slope for class 
# (summary( model.lmer ) shows this at the section "Random effects:", row "classh", column "Std.Dev")

vdur_class_1classspk.stddev_reduction = 1 - ( attr(VarCorr(vdur_class_1classspk.lmer)$spk,'stddev')["classh"] / attr(VarCorr(vdur_1_1classspk.lmer)$spk,'stddev')["classh"] )

# Diagnostics for the full model
shapiro.test(residuals(vdur_class_1classspk.lmer)) 
# not passed, but outliers are not many
qqnorm(residuals(vdur_class_1classspk.lmer))
qqline(residuals(vdur_class_1classspk.lmer))
# check whether residuals depend on the fitted values
plot(fitted(vdur_class_1classspk.lmer),residuals(vdur_class_1classspk.lmer))
summary(lm(residuals(vdur_class_1classspk.lmer) ~ fitted(vdur_class_1classspk.lmer)))
# they don't
xyplot(residuals(vdur_class_1classspk.lmer) ~ fitted(vdur_class_1classspk.lmer) |  spk, data = DH_data)
# no sign of heteroscedasticity

detach("package:lme4")
library(nlme)


## alternative with lme from package nlme
# see book "LINEAR MIXED MODELS A Practical Guide Using Statistical Software" by Brady T. West, Kathleen B. Welch and Andrzej T. Galecki
# and companion website http://www-personal.umich.edu/~bwest/almmussp.html

vdur_class_spk.gd = groupedData(vdur ~ class | spk, data = DH_data)
#plot(vdur_class_spk.gd, display = "spk")

summary(vdur_class_1classspk.lme <- lme(vdur ~ class, random = ~ class, method="REML", data = vdur_class_spk.gd))
summary(vdur_1_1classspk.lme <- lme(vdur ~ 1, random = ~ class, method="REML", data = vdur_class_spk.gd))
summary(vdur_class_1spk.lme <- lme(vdur ~ class, random = ~ 1, method="REML", data = vdur_class_spk.gd))
summary(vdur_1_1spk.lme <- lme(vdur ~  1, random = ~  1, method="REML", data = vdur_class_spk.gd))

vdur_class_1classspk.pseudoR22 = 1 - (cor(fitted(vdur_1_1spk.lme),DH_data$vdur) / cor(fitted(vdur_class_1classspk.lme),DH_data$vdur) )^2
# same as vdur_class_1classspk.pseudoR2

vdur_class_1classspk.stddev_reduction2 = 1 - (sqrt(diag(getVarCov(vdur_class_1classspk.lme)))["classh"] / sqrt(diag(getVarCov(vdur_1_1classspk.lme)))["classh"])
# same as vdur_class_1classspk.stddev_reduction


plot(vdur_class_1classspk.lme,resid(., type="p") ~ fitted(.) | factor(class))
plot(vdur_class_1classspk.lme,resid(., type="p") ~ fitted(.) | factor(class), id=.05,idLabels=as.character(DH_data$filename))

plot(vdur_class_1classspk.lme,vdur ~ fitted(.) | factor(class))


# check heteroshedasticity
# see http://www-personal.umich.edu/~bwest/chapter5_R_final.txt , Model 5.3
summary(vdur_class_1classspk_Het.lme <- lme(vdur ~ class, random = ~ class, weights = varIdent(form = ~1 | class), method="REML", data = vdur_class_spk.gd))
anova(vdur_class_1classspk_Het.lme,vdur_class_1classspk.lme) # heteroshedasticity w.r.t. class not justified


detach("package:nlme")
library(lme4)

# plot speaker-specific coefficients
vdur_class_1classspk.coef = coef(vdur_class_1classspk.lmer)
png(paste(plots_dir,'vdur_class_1classspk_coef.png',sep=''))
par(mar=c(5.1,6.6,2.1,2.1),mgp = c(4, 1, 0))
plot(vdur_class_1classspk.coef$spk,type='n',las=1,xlab=expression(beta[0] + u[{list(0,j)}]),ylab=expression(beta[1] + u[{list(1,j)}]),cex.lab=2.5,cex.axis=1.8)
for (spk in speakers) {
    text(vdur_class_1classspk.coef$spk[spk,],label=spk,cex=1.8)
}
points(vdur_class_1classspk.lmer@fixef[1],vdur_class_1classspk.lmer@fixef[2],pch = 17, col='red', cex=2.5)
dev.off()

# boxplot of vdur and by-speaker predicted models

png(paste(plots_dir,'vdur_spk_box.png',sep=''))
bwp = bwplot( vdur ~ class | spk, data = DH_data,  ylab = 'vowel seq duration (ms)',
            panel = function(...,subscripts) { # see (Baayen, 2008)
                panel.bwplot(...)
                speaker = as.character(DH_data[subscripts[1], "spk"]) 
                coefs = as.numeric(unlist(coef(vdur_class_1classspk.lmer)$spk[speaker,]))
                # boxplots are aligned at x=1 and x=2, while coefs refer to values x=0 and x=1, thus the line should be y = coefs[2]*(x-1) + coefs[1]
                #panel.abline(coefs[1]-coefs[2],coefs[2],col='blue',lty=2)
                panel.xyplot(1:2,c(coefs[1],(coefs[1]+coefs[2])),type='b',col='red',lty=2,pch=15)
            }
)
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()


########## The same procedure applied to f0_s1, f0_s2, F12_s1 and F12_s2.

######### Analysis of f0_s1


# full mixed model:
summary(f0_s1_class_1classspk.lmer <- lmer(f0_s1 ~ class + (1 + class | spk), data = DH_data))
# null model, only random effects
summary(f0_s1_1_1spk.lmer <- lmer(f0_s1 ~ 1 + (1 | spk), data = DH_data))
# model with class only as random effect
summary(f0_s1_1_1classspk.lmer <- lmer(f0_s1 ~ 1 + (1 + class | spk), data = DH_data))
# model with class only as fixed effect
summary(f0_s1_class_1spk.lmer <- lmer(f0_s1 ~ class + (1 | spk), data = DH_data))

# series of log-likelihood tests to justify the model complexity (model selection)
anova(f0_s1_class_1classspk.lmer,f0_s1_class_1spk.lmer) # random slope justified
anova(f0_s1_class_1classspk.lmer,f0_s1_1_1classspk.lmer) # fixed effect class NOT justified
# similar tests starting from the null model
anova(f0_s1_1_1classspk.lmer,f0_s1_1_1spk.lmer)
anova(f0_s1_class_1spk.lmer,f0_s1_1_1spk.lmer)

# Gain in explained variance (a pseudo R^2, see book "Analyzing Linguistic Data: A Practical Introduction to Statistics using R" by Baayen)
# based on the correlation  between fitted and actual values of the dependent variable.
# The comparison is between the full model and the null model. 
f0_s1_class_1classspk.pseudoR2 = 1 - (cor(fitted(f0_s1_1_1spk.lmer),DH_data$f0_s1) / cor(fitted(f0_s1_class_1classspk.lmer),DH_data$f0_s1) )^2

# Reduction in std. dev. due to the introduction of class as fixed effect.
# This time the comparison is between the full model and the model without fixed factor (see Baayen book for details).

# Extract std. dev. of the random slope for class 
# (summary( model.lmer ) shows this at the section "Random effects:", row "classh", column "Std.Dev")

f0_s1_class_1classspk.stddev_reduction = 1 - ( attr(VarCorr(f0_s1_class_1classspk.lmer)$spk,'stddev')["classh"] / attr(VarCorr(f0_s1_1_1classspk.lmer)$spk,'stddev')["classh"] )

# Diagnostics for the full model
shapiro.test(residuals(f0_s1_class_1classspk.lmer)) 
# not passed, but outliers are not many
qqnorm(residuals(f0_s1_class_1classspk.lmer))
qqline(residuals(f0_s1_class_1classspk.lmer))
# check whether residuals depend on the fitted values
plot(fitted(f0_s1_class_1classspk.lmer),residuals(f0_s1_class_1classspk.lmer))
summary(lm(residuals(f0_s1_class_1classspk.lmer) ~ fitted(f0_s1_class_1classspk.lmer)))
# they don't
xyplot(residuals(f0_s1_class_1classspk.lmer) ~ fitted(f0_s1_class_1classspk.lmer) |  spk, data = DH_data)
# no sign of heteroscedasticity

detach("package:lme4")
library(nlme)

f0_s1_class_spk.gd = groupedData(f0_s1 ~ class | spk, data = DH_data)
summary(f0_s1_class_1classspk.lme <- lme(f0_s1 ~ class, random = ~ class, method="REML", data = f0_s1_class_spk.gd))

# check heteroshedasticity
# see http://www-personal.umich.edu/~bwest/chapter5_R_final.txt , Model 5.3
summary(f0_s1_class_1classspk_Het.lme <- lme(f0_s1 ~ class, random = ~ class, weights = varIdent(form = ~1 | class), method="REML", data = f0_s1_class_spk.gd))
anova(f0_s1_class_1classspk_Het.lme,f0_s1_class_1classspk.lme) # heteroshedasticity w.r.t. class justified
# effect small: sigma(H) = 0.85 sigma(H)

detach("package:nlme")
library(lme4)


# plot speaker-specific coefficients
f0_s1_class_1classspk.coef = coef(f0_s1_class_1classspk.lmer)
png(paste(plots_dir,'f0_s1_class_1classspk_coef.png',sep=''))
par(mar=c(5.1,6.6,2.1,2.1),mgp = c(4, 1, 0))
plot(f0_s1_class_1classspk.coef$spk,type='n',las=1,xlab=expression(beta[0] + u[{list(0,j)}]),ylab=expression(beta[1] + u[{list(1,j)}]),cex.lab=2.5,cex.axis=1.8)
for (spk in speakers) {
    text(f0_s1_class_1classspk.coef$spk[spk,],label=spk,cex=1.8)
}
points(f0_s1_class_1classspk.lmer@fixef[1],f0_s1_class_1classspk.lmer@fixef[2],pch = 17, col='red', cex=2.5)
dev.off()


# boxplot of f0_s1 and by-speaker predicted models

png(paste(plots_dir,'f0_s1_spk_box.png',sep=''))
bwp = bwplot( f0_s1 ~ class | spk, data = DH_data,  ylab = expression(s[1]^{f[0]}),
            panel = function(...,subscripts) { # see (Baayen, 2008)
                panel.bwplot(...)
                speaker = as.character(DH_data[subscripts[1], "spk"]) 
                coefs = as.numeric(unlist(coef(f0_s1_class_1classspk.lmer)$spk[speaker,]))
                # boxplots are aligned at x=1 and x=2, while coefs refer to values x=0 and x=1, thus the line should be y = coefs[2]*(x-1) + coefs[1]
                #panel.abline(coefs[1]-coefs[2],coefs[2],col='blue',lty=2)
                panel.xyplot(1:2,c(coefs[1],(coefs[1]+coefs[2])),type='b',col='red',lty=2,pch=15)
            }
)
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()

# plot speaker-specific contours as predicted by the mixed model
# note: it requires to load fd objects computed in FDA.R

# see xyplot.ts
table_plot = expand.grid(class = c('d','h'),spk = speakers,stringsAsFactors = FALSE)
curves = matrix(nrow = length(t_f0),ncol = nrow(table_plot))
for (i in 1:nrow(table_plot)) {
    curve = f0_pcafd$meanfd +
        (f0_s1_class_1classspk.coef$spk[table_plot$spk[i],1] + (table_plot$class[i]=='h') * f0_s1_class_1classspk.coef$spk[table_plot$spk[i],2] ) * f0_pcafd$harmonics[1] 
    curves[,i] = eval.fd(t_f0,curve)
}

png(paste(plots_dir,'f0_s1_class_1classspk_curves.png',sep=''))
xyp =	xyplot(
	ts(data=curves,start=t_f0[1],deltat=t_f0[2]-t_f0[1]),
	screens=table_plot$spk,
	col = sapply(table_plot$class, function(x) color[[x]],USE.NAMES = FALSE),
	lty = sapply(table_plot$class, function(x) lty[[x]],USE.NAMES = FALSE),
	lwd = sapply(table_plot$class, function(x) lwd[[x]],USE.NAMES = FALSE),
	layout = c(3,3),
	xlab = 'time (ms)',
	ylab = 'F0 (norm. semitones)',
	default.scales = list(relation='same',cex=1.0),
	panel = function(col,lty,lwd,...) {
	    panel.superpose.plain(col=col,lty=lty,lwd=lwd,...)
	    panel.abline(v=reg$land[2],lty=2,col='black',lwd=1)
	}
	)
update	(xyp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)
	    ),
    	key = list(
    	    space = 'top',
    	    lines = list(
    		    col = as.character(unlist(color)),
    		    lty = as.numeric(unlist(lty)),
    		    lwd = as.numeric(unlist(lwd))
    	        ),
    	    text = list(
    		    lab = as.character(unlist(name)),
    		    cex = 1.2
    		    )
    	),
        as.table=TRUE
    	)
dev.off()




######### Analysis of f0_s2


# full mixed model:
summary(f0_s2_class_1classspk.lmer <- lmer(f0_s2 ~ class + (1 + class | spk), data = DH_data))
# null model, only random effects
summary(f0_s2_1_1spk.lmer <- lmer(f0_s2 ~ 1 + (1 | spk), data = DH_data))
# model with class only as random effect
summary(f0_s2_1_1classspk.lmer <- lmer(f0_s2 ~ 1 + (1 + class | spk), data = DH_data))
# model with class only as fixed effect
summary(f0_s2_class_1spk.lmer <- lmer(f0_s2 ~ class + (1 | spk), data = DH_data))

# series of log-likelihood tests to justify the model complexity (model selection)
anova(f0_s2_class_1classspk.lmer,f0_s2_class_1spk.lmer) # random slope justified
anova(f0_s2_class_1classspk.lmer,f0_s2_1_1classspk.lmer) # fixed effect class justified
# similar tests starting from the null model
anova(f0_s2_1_1classspk.lmer,f0_s2_1_1spk.lmer)
anova(f0_s2_class_1spk.lmer,f0_s2_1_1spk.lmer)

# Gain in explained variance (a pseudo R^2, see book "Analyzing Linguistic Data: A Practical Introduction to Statistics using R" by Baayen)
# based on the correlation  between fitted and actual values of the dependent variable.
# The comparison is between the full model and the null model. 
f0_s2_class_1classspk.pseudoR2 = 1 - (cor(fitted(f0_s2_1_1spk.lmer),DH_data$f0_s2) / cor(fitted(f0_s2_class_1classspk.lmer),DH_data$f0_s2) )^2

# Reduction in std. dev. due to the introduction of class as fixed effect.
# This time the comparison is between the full model and the model without fixed factor (see Baayen book for details).

# Extract std. dev. of the random slope for class 
# (summary( model.lmer ) shows this at the section "Random effects:", row "classh", column "Std.Dev")

f0_s2_class_1classspk.stddev_reduction = 1 - ( attr(VarCorr(f0_s2_class_1classspk.lmer)$spk,'stddev')["classh"] / attr(VarCorr(f0_s2_1_1classspk.lmer)$spk,'stddev')["classh"] )

# Diagnostics for the full model
shapiro.test(residuals(f0_s2_class_1classspk.lmer)) 
# not passed, but outliers are not many
qqnorm(residuals(f0_s2_class_1classspk.lmer))
qqline(residuals(f0_s2_class_1classspk.lmer))
# check whether residuals depend on the fitted values
plot(fitted(f0_s2_class_1classspk.lmer),residuals(f0_s2_class_1classspk.lmer))
summary(lm(residuals(f0_s2_class_1classspk.lmer) ~ fitted(f0_s2_class_1classspk.lmer)))
# they don't
xyplot(residuals(f0_s2_class_1classspk.lmer) ~ fitted(f0_s2_class_1classspk.lmer) |  spk, data = DH_data)
# some sign of speaker-dependent heteroscedasticity, but not consistent, i.e. sometimes D has greater errors, sometimes H does.

detach("package:lme4")
library(nlme)

f0_s2_class_spk.gd = groupedData(f0_s2 ~ class | spk, data = DH_data)
summary(f0_s2_class_1classspk.lme <- lme(f0_s2 ~ class, random = ~ class, method="REML", data = f0_s2_class_spk.gd))

# check heteroshedasticity
# see http://www-personal.umich.edu/~bwest/chapter5_R_final.txt , Model 5.3
summary(f0_s2_class_1classspk_Het.lme <- lme(f0_s2 ~ class, random = ~ class, weights = varIdent(form = ~1 | class), method="REML", data = f0_s2_class_spk.gd))
anova(f0_s2_class_1classspk_Het.lme,f0_s2_class_1classspk.lme) # heteroshedasticity w.r.t. class not justified


detach("package:nlme")
library(lme4)

f0_s2_class_1classspk.coef = coef(f0_s2_class_1classspk.lmer)
png(paste(plots_dir,'f0_s2_class_1classspk_coef.png',sep=''))
par(mar=c(5.1,6.6,2.1,2.1),mgp = c(4, 1, 0))
plot(f0_s2_class_1classspk.coef$spk,type='n',las=1,xlab=expression(beta[0] + u[{list(0,j)}]),ylab=expression(beta[1] + u[{list(1,j)}]),cex.lab=2.5,cex.axis=1.8)
for (spk in speakers) {
    text(f0_s2_class_1classspk.coef$spk[spk,],label=spk,cex=1.8)
}
points(f0_s2_class_1classspk.lmer@fixef[1],f0_s2_class_1classspk.lmer@fixef[2],pch = 17, col='red', cex=2.5)
dev.off()


# boxplot of f0_s2 and by-speaker predicted models

png(paste(plots_dir,'f0_s2_spk_box.png',sep=''))
bwp = bwplot( f0_s2 ~ class | spk, data = DH_data,  ylab = expression(s[2]^{f[0]}),
            panel = function(...,subscripts) { # see (Baayen, 2008)
                panel.bwplot(...)
                speaker = as.character(DH_data[subscripts[1], "spk"]) 
                coefs = as.numeric(unlist(coef(f0_s2_class_1classspk.lmer)$spk[speaker,]))
                # boxplots are aligned at x=1 and x=2, while coefs refer to values x=0 and x=1, thus the line should be y = coefs[2]*(x-1) + coefs[1]
                #panel.abline(coefs[1]-coefs[2],coefs[2],col='blue',lty=2)
                panel.xyplot(1:2,c(coefs[1],(coefs[1]+coefs[2])),type='b',col='red',lty=2,pch=15)
            }
)
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()

# plot speaker-specific contours as predicted by the mixed model
# note: it requires to load fd objects computed in FDA.R

# see xyplot.ts
table_plot = expand.grid(class = c('d','h'),spk = speakers,stringsAsFactors = FALSE)
curves = matrix(nrow = length(t_f0),ncol = nrow(table_plot))
for (i in 1:nrow(table_plot)) {
    curve = f0_pcafd$meanfd +
        (f0_s2_class_1classspk.coef$spk[table_plot$spk[i],1] + (table_plot$class[i]=='h') * f0_s2_class_1classspk.coef$spk[table_plot$spk[i],2] ) * f0_pcafd$harmonics[2] 
    curves[,i] = eval.fd(t_f0,curve)
}

png(paste(plots_dir,'f0_s2_class_1classspk_curves.png',sep=''))
xyp =	xyplot(
	ts(data=curves,start=t_f0[1],deltat=t_f0[2]-t_f0[1]),
	screens=table_plot$spk,
	col = sapply(table_plot$class, function(x) color[[x]],USE.NAMES = FALSE),
	lty = sapply(table_plot$class, function(x) lty[[x]],USE.NAMES = FALSE),
	lwd = sapply(table_plot$class, function(x) lwd[[x]],USE.NAMES = FALSE),
	layout = c(3,3),
	xlab = 'time (ms)',
	ylab = 'F0 (norm. semitones)',
	default.scales = list(relation='same',cex=1.0),
	panel = function(col,lty,lwd,...) {
	    panel.superpose.plain(col=col,lty=lty,lwd=lwd,...)
	    panel.abline(v=reg$land[2],lty=2,col='black',lwd=1)
	}
	)
update	(xyp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)
	    ),
    	key = list(
    	    space = 'top',
    	    lines = list(
    		    col = as.character(unlist(color)),
    		    lty = as.numeric(unlist(lty)),
    		    lwd = as.numeric(unlist(lwd))
    	        ),
    	    text = list(
    		    lab = as.character(unlist(name)),
    		    cex = 1.2
    		    )
    	),
        as.table=TRUE
    	)
dev.off()




######### Analysis of F12_s1



# full mixed model:
summary(F12_s1_class_1classspk.lmer <- lmer(F12_s1 ~ class + (1 + class | spk), data = DH_data))
# null model, only random effects
summary(F12_s1_1_1spk.lmer <- lmer(F12_s1 ~ 1 + (1 | spk), data = DH_data))
# model with class only as random effect
summary(F12_s1_1_1classspk.lmer <- lmer(F12_s1 ~ 1 + (1 + class | spk), data = DH_data))
# model with class only as fixed effect
summary(F12_s1_class_1spk.lmer <- lmer(F12_s1 ~ class + (1 | spk), data = DH_data))

# series of log-likelihood tests to justify the model complexity (model selection)
anova(F12_s1_class_1classspk.lmer,F12_s1_class_1spk.lmer) # random slope justified
anova(F12_s1_class_1classspk.lmer,F12_s1_1_1classspk.lmer) # fixed effect class justified
# similar tests starting from the null model
anova(F12_s1_1_1classspk.lmer,F12_s1_1_1spk.lmer)
anova(F12_s1_class_1spk.lmer,F12_s1_1_1spk.lmer)

# Gain in explained variance (a pseudo R^2, see book "Analyzing Linguistic Data: A Practical Introduction to Statistics using R" by Baayen)
# based on the correlation  between fitted and actual values of the dependent variable.
# The comparison is between the full model and the null model. 
F12_s1_class_1classspk.pseudoR2 = 1 - (cor(fitted(F12_s1_1_1spk.lmer),DH_data$F12_s1) / cor(fitted(F12_s1_class_1classspk.lmer),DH_data$F12_s1) )^2

# Reduction in std. dev. due to the introduction of class as fixed effect.
# This time the comparison is between the full model and the model without fixed factor (see Baayen book for details).

# Extract std. dev. of the random slope for class 
# (summary( model.lmer ) shows this at the section "Random effects:", row "classh", column "Std.Dev")

F12_s1_class_1classspk.stddev_reduction = 1 - ( attr(VarCorr(F12_s1_class_1classspk.lmer)$spk,'stddev')["classh"] / attr(VarCorr(F12_s1_1_1classspk.lmer)$spk,'stddev')["classh"] )

# Diagnostics for the full model
shapiro.test(residuals(F12_s1_class_1classspk.lmer)) 
# not passed, but outliers are not many
qqnorm(residuals(F12_s1_class_1classspk.lmer))
qqline(residuals(F12_s1_class_1classspk.lmer))
# check whether residuals depend on the fitted values
plot(fitted(F12_s1_class_1classspk.lmer),residuals(F12_s1_class_1classspk.lmer))
summary(lm(residuals(F12_s1_class_1classspk.lmer) ~ fitted(F12_s1_class_1classspk.lmer)))
# they don't
xyplot(residuals(F12_s1_class_1classspk.lmer) ~ fitted(F12_s1_class_1classspk.lmer) |  spk, data = DH_data)
# some sign of speaker-dependent heteroscedasticity, but not consistent, i.e. sometimes D has greater errors, sometimes H does.

detach("package:lme4")
library(nlme)

F12_s1_class_spk.gd = groupedData(F12_s1 ~ class | spk, data = DH_data)
summary(F12_s1_class_1classspk.lme <- lme(F12_s1 ~ class, random = ~ class, method="REML", data = F12_s1_class_spk.gd))

# check heteroshedasticity
# see http://www-personal.umich.edu/~bwest/chapter5_R_final.txt , Model 5.3
summary(F12_s1_class_1classspk_Het.lme <- lme(F12_s1 ~ class, random = ~ class, weights = varIdent(form = ~1 | class), method="REML", data = F12_s1_class_spk.gd))
anova(F12_s1_class_1classspk_Het.lme,F12_s1_class_1classspk.lme) # heteroshedasticity w.r.t. class justified, H varies more than D.
# speaker-dep coefficients:
coef(F12_s1_class_1classspk_Het.lme)
coef(F12_s1_class_1classspk.lme)
# are almost identical. The effect of heteroshedasticity is small, sigma(H) = 1.28 sigma(D).

detach("package:nlme")
library(lme4)

# plot speaker-specific coefficients
F12_s1_class_1classspk.coef = coef(F12_s1_class_1classspk.lmer)
png(paste(plots_dir,'F12_s1_class_1classspk_coef.png',sep=''))
par(mar=c(5.1,6.6,2.1,2.1),mgp = c(4, 1, 0))
plot(F12_s1_class_1classspk.coef$spk,type='n',las=1,xlab=expression(beta[0] + u[{list(0,j)}]),ylab=expression(beta[1] + u[{list(1,j)}]),cex.lab=1.5,cex.lab=2.5,cex.axis=1.8)
for (spk in speakers) {
    text(F12_s1_class_1classspk.coef$spk[spk,],label=spk,cex=1.8)
}
points(F12_s1_class_1classspk.lmer@fixef[1],F12_s1_class_1classspk.lmer@fixef[2],pch = 17, col='red', cex=2.5)
dev.off()



# boxplot of F12_s1 and by-speaker predicted models

png(paste(plots_dir,'F12_s1_spk_box.png',sep=''))
bwp = bwplot( F12_s1 ~ class | spk, data = DH_data,  ylab = expression(s[1]^{F[1-2]}),
            panel = function(...,subscripts) { # see (Baayen, 2008)
                panel.bwplot(...)
                speaker = as.character(DH_data[subscripts[1], "spk"]) 
                coefs = as.numeric(unlist(coef(F12_s1_class_1classspk.lmer)$spk[speaker,]))
                # boxplots are aligned at x=1 and x=2, while coefs refer to values x=0 and x=1, thus the line should be y = coefs[2]*(x-1) + coefs[1]
                #panel.abline(coefs[1]-coefs[2],coefs[2],col='blue',lty=2)
                panel.xyplot(1:2,c(coefs[1],(coefs[1]+coefs[2])),type='b',col='red',lty=2,pch=15)
            }
)
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()


# plot speaker-specific contours as predicted by the mixed model
# note: it requires to load fd objects computed in FDA.R

# see xyplot.ts

table_plot_F12 = expand.grid(class = c('d','h'),F12=1:2,spk = speakers)
n_col_plot = 3
n_combined_rows_plot = length(speakers)/n_col_plot # 3 combined row means a row with F1 and f2 one under the other


panel.order = c()
for (c in 1:n_combined_rows_plot) {
    # select rows with spk in the current block (plot combined row)
    table_plot_F12_subset = subset(table_plot_F12,spk %in% levels(spk)[(c-1)*n_col_plot+(1:n_col_plot)])
    panel.order = c(panel.order,as.integer(row.names(table_plot_F12_subset[with(table_plot_F12_subset,order(-F12,spk,class)),])))
}
table_plot_F12 = table_plot_F12[panel.order,]
curves_F12 = matrix(nrow = length(t_F12),ncol = nrow(table_plot_F12))
for (i in 1:nrow(table_plot_F12)) {
    curve = F12_pcafd$meanfd[,table_plot_F12$F12[i]] +
        (F12_s1_class_1classspk.coef$spk[table_plot_F12$spk[i],1] + (table_plot_F12$class[i]=='h') * F12_s1_class_1classspk.coef$spk[table_plot_F12$spk[i],2] ) * F12_pcafd$harmonics[1,table_plot_F12$F12[i]] 
     curves_F12[,i] = eval.fd(t_F12,curve)
}

png(paste(plots_dir,'F12_s1_class_1classspk_curves.png',sep=''))
xyp =	xyplot(
	ts(data=curves_F12,start=t_F12[1],deltat=t_F12[2]-t_F12[1]),
	screens=with( table_plot_F12, paste(spk,', F',F12,sep='')),
    col = sapply(table_plot_F12$class, function(x) color[[x]],USE.NAMES = FALSE),
	lty = sapply(table_plot_F12$class, function(x) lty[[x]],USE.NAMES = FALSE),
	lwd = sapply(table_plot_F12$class, function(x) lwd[[x]],USE.NAMES = FALSE),
    layout = c(3,6),
	xlab = 'norm. time (ms)',
	ylab = 'norm. barks',
	default.scales = list(relation='same',cex=1.0),
)
update	(xyp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)
	    ),
    	key = list(
    	    space = 'top',
    	    lines = list(
    		    col = as.character(unlist(color)),
    		    lty = as.numeric(unlist(lty)),
    		    lwd = as.numeric(unlist(lwd))
    	        ),
    	    text = list(
    		    lab = as.character(unlist(name)),
    		    cex = 1.2
    		    )
    	    ),
        between = list(x = 0.5, y = c(0,0.5))
    	)
dev.off()




######### Analysis of F12_s2



# full mixed model:
summary(F12_s2_class_1classspk.lmer <- lmer(F12_s2 ~ class + (1 + class | spk), data = DH_data))
# null model, only random effects
summary(F12_s2_1_1spk.lmer <- lmer(F12_s2 ~ 1 + (1 | spk), data = DH_data))
# model with class only as random effect
summary(F12_s2_1_1classspk.lmer <- lmer(F12_s2 ~ 1 + (1 + class | spk), data = DH_data))
# model with class only as fixed effect
summary(F12_s2_class_1spk.lmer <- lmer(F12_s2 ~ class + (1 | spk), data = DH_data))

# series of log-likelihood tests to justify the model complexity (model selection)
anova(F12_s2_class_1classspk.lmer,F12_s2_class_1spk.lmer) # random slope justified
anova(F12_s2_class_1classspk.lmer,F12_s2_1_1classspk.lmer) # fixed effect class justified
# similar tests starting from the null model
anova(F12_s2_1_1classspk.lmer,F12_s2_1_1spk.lmer)
anova(F12_s2_class_1spk.lmer,F12_s2_1_1spk.lmer)

# Gain in explained variance (a pseudo R^2, see book "Analyzing Linguistic Data: A Practical Introduction to Statistics using R" by Baayen)
# based on the correlation  between fitted and actual values of the dependent variable.
# The comparison is between the full model and the null model. 
F12_s2_class_1classspk.pseudoR2 = 1 - (cor(fitted(F12_s2_1_1spk.lmer),DH_data$F12_s2) / cor(fitted(F12_s2_class_1classspk.lmer),DH_data$F12_s2) )^2

# Reduction in std. dev. due to the introduction of class as fixed effect.
# This time the comparison is between the full model and the model without fixed factor (see Baayen book for details).

# Extract std. dev. of the random slope for class 
# (summary( model.lmer ) shows this at the section "Random effects:", row "classh", column "Std.Dev")

F12_s2_class_1classspk.stddev_reduction = 1 - ( attr(VarCorr(F12_s2_class_1classspk.lmer)$spk,'stddev')["classh"] / attr(VarCorr(F12_s2_1_1classspk.lmer)$spk,'stddev')["classh"] )

# Diagnostics for the full model
shapiro.test(residuals(F12_s2_class_1classspk.lmer)) 
# not passed, but outliers are not many
qqnorm(residuals(F12_s2_class_1classspk.lmer))
qqline(residuals(F12_s2_class_1classspk.lmer))
# check whether residuals depend on the fitted values
plot(fitted(F12_s2_class_1classspk.lmer),residuals(F12_s2_class_1classspk.lmer))
summary(lm(residuals(F12_s2_class_1classspk.lmer) ~ fitted(F12_s2_class_1classspk.lmer)))
# they don't
xyplot(residuals(F12_s2_class_1classspk.lmer) ~ fitted(F12_s2_class_1classspk.lmer) |  spk, data = DH_data)
# some sign of speaker-dependent heteroscedasticity, but not consistent, i.e. sometimes D has greater errors, sometimes H does.

detach("package:lme4")
library(nlme)

F12_s2_class_spk.gd = groupedData(F12_s2 ~ class | spk, data = DH_data)
summary(F12_s2_class_1classspk.lme <- lme(F12_s2 ~ class, random = ~ class, method="REML", data = F12_s2_class_spk.gd))

# check heteroshedasticity
# see http://www-personal.umich.edu/~bwest/chapter5_R_final.txt , Model 5.3
summary(F12_s2_class_1classspk_Het.lme <- lme(F12_s2 ~ class, random = ~ class, weights = varIdent(form = ~1 | class), method="REML", data = F12_s2_class_spk.gd))
anova(F12_s2_class_1classspk_Het.lme,F12_s2_class_1classspk.lme) # heteroshedasticity w.r.t. class not justified

detach("package:nlme")
library(lme4)

# plot speaker-specific coefficients
F12_s2_class_1classspk.coef = coef(F12_s2_class_1classspk.lmer)
png(paste(plots_dir,'F12_s2_class_1classspk_coef.png',sep=''))
par(mar=c(5.1,6.6,2.1,2.1),mgp = c(4, 1, 0))
plot(F12_s2_class_1classspk.coef$spk,type='n',las=1,xlab=expression(beta[0] + u[{list(0,j)}]),ylab=expression(beta[1] + u[{list(1,j)}]),cex.lab=2.5,cex.axis=1.8)
for (spk in speakers) {
    text(F12_s2_class_1classspk.coef$spk[spk,],label=spk,cex=1.8)
}
points(F12_s2_class_1classspk.lmer@fixef[1],F12_s2_class_1classspk.lmer@fixef[2],pch = 17, col='red', cex=2.5)
dev.off()



# boxplot of F12_s2 and by-speaker predicted models

png(paste(plots_dir,'F12_s2_spk_box.png',sep=''))
bwp = bwplot( F12_s2 ~ class | spk, data = DH_data,  ylab = expression(s[2]^{F[1-2]}),
            panel = function(...,subscripts) { # see (Baayen, 2008)
                panel.bwplot(...)
                speaker = as.character(DH_data[subscripts[1], "spk"]) 
                coefs = as.numeric(unlist(coef(F12_s2_class_1classspk.lmer)$spk[speaker,]))
                # boxplots are aligned at x=1 and x=2, while coefs refer to values x=0 and x=1, thus the line should be y = coefs[2]*(x-1) + coefs[1]
                #panel.abline(coefs[1]-coefs[2],coefs[2],col='blue',lty=2)
                panel.xyplot(1:2,c(coefs[1],(coefs[1]+coefs[2])),type='b',col='red',lty=2,pch=15)
            }
)
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()

# plot speaker-specific contours as predicted by the mixed model
# note: it requires to load fd objects computed in FDA.R

# see xyplot.ts

table_plot_F12 = expand.grid(class = c('d','h'),F12=1:2,spk = speakers)
n_col_plot = 3
n_combined_rows_plot = length(speakers)/n_col_plot # 3 combined row means a row with F1 and f2 one under the other


panel.order = c()
for (c in 1:n_combined_rows_plot) {
    # select rows with spk in the current block (plot combined row)
    table_plot_F12_subset = subset(table_plot_F12,spk %in% levels(spk)[(c-1)*n_col_plot+(1:n_col_plot)])
    panel.order = c(panel.order,as.integer(row.names(table_plot_F12_subset[with(table_plot_F12_subset,order(-F12,spk,class)),])))
}
table_plot_F12 = table_plot_F12[panel.order,]
curves_F12 = matrix(nrow = length(t_F12),ncol = nrow(table_plot_F12))
for (i in 1:nrow(table_plot_F12)) {
    curve = F12_pcafd$meanfd[,table_plot_F12$F12[i]] +
        (F12_s2_class_1classspk.coef$spk[table_plot_F12$spk[i],1] + (table_plot_F12$class[i]=='h') * F12_s2_class_1classspk.coef$spk[table_plot_F12$spk[i],2] ) * F12_pcafd$harmonics[2,table_plot_F12$F12[i]] 
     curves_F12[,i] = eval.fd(t_F12,curve)
}

png(paste(plots_dir,'F12_s2_class_1classspk_curves.png',sep=''))
xyp =	xyplot(
	ts(data=curves_F12,start=t_F12[1],deltat=t_F12[2]-t_F12[1]),
	screens=with( table_plot_F12, paste(spk,', F',F12,sep='')),
    col = sapply(table_plot_F12$class, function(x) color[[x]],USE.NAMES = FALSE),
	lty = sapply(table_plot_F12$class, function(x) lty[[x]],USE.NAMES = FALSE),
	lwd = sapply(table_plot_F12$class, function(x) lwd[[x]],USE.NAMES = FALSE),
    layout = c(3,6),
	xlab = 'norm. time (ms)',
	ylab = 'norm. barks',
	default.scales = list(relation='same',cex=1.0),
)
update	(xyp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)
	    ),
    	key = list(
    	    space = 'top',
    	    lines = list(
    		    col = as.character(unlist(color)),
    		    lty = as.numeric(unlist(lty)),
    		    lwd = as.numeric(unlist(lwd))
    	        ),
    	    text = list(
    		    lab = as.character(unlist(name)),
    		    cex = 1.2
    		    )
    	    ),
        between = list(x = 0.5, y = c(0,0.5))
    	)
dev.off()






save.image('DA_MM.RImage')


# generalised linear mixed model on the final glm in DA.R
class_0PCA13_int_0PCA13.lmer = lmer(class ~ 0 + PC1 * PC3 + (0 + PC1 + PC3 | spk), data = DH_data, family = binomial)
summary(class_0PCA13_int_0PCA13.lmer)
class_0PCA13_int_0PCA13.coef = coef(class_0PCA13_int_0PCA13.lmer)

# plot predicted probability of H  from model class_0PCA13_int_0PCA13.lmer applying per-speaker corrections, translated back into vdur_CS F12_s1_CS coordinates using PCA loadings


cRP = colorRampPalette(c(color[['d']],"white",color[['h']]))
grid = expand.grid( vdur_CS     = do.breaks(c(-3,3),60),
                    F12_s1_CS   = do.breaks(c(-3,3),60))
lp = list() # contains levelplots

png(paste(plots_dir,'prob_vdur_CS_F12_s1_CS_speakers.png',sep=''))
for (spk in speakers) {
    # build the speaker-specific glmm predicting function
    grid$logit =   with(grid, 
    class_0PCA13_int_0PCA13.coef$spk[spk,'PC1'] * (-1) * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC1'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC1'] * F12_s1_CS) + # (-1) is because of the manual sign swap operated on PC1
    class_0PCA13_int_0PCA13.coef$spk[spk,'PC3'] *        (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC3'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC3'] * F12_s1_CS) +
    class_0PCA13_int_0PCA13.coef$spk[spk,'PC1:PC3'] * (-1) * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC1'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC1'] * F12_s1_CS) * (vdur_f0_s2_F12_s1.pca$rotation['vdur','PC3'] * vdur_CS + vdur_f0_s2_F12_s1.pca$rotation['F12_s1','PC3'] * F12_s1_CS)
        )
    grid$prob = with(grid, 1/(1 + exp(-logit)) )  # inverse logit
    # plot
    index_spk = which(DH_data$spk == spk)    
    lp[[spk]] = levelplot(prob~vdur_CS*F12_s1_CS|spk,data=grid,at = do.breaks(c(0,1), 30), col.regions = cRP,
    xlab = list(label='norm. d',cex=1.5), ylab=list(label=expression(paste('norm.',s[1]^{F[1-2]})),cex=1.5),    
    par.settings=list(
        par.xlab.text = list(cex=1.3),
        par.ylab.text = list(cex=1.3),
        axis.text = list(cex=1.3),
        strip.background = list(col=strip.background.col),
        layout.heights=list(strip=1)
        ),
    colorkey = list(labels=list(cex=1.3)),
    strip = strip.custom(par.strip.text = list(cex = 1)),
    panel = function(...) {
        panel.levelplot(...)
        panel.points(x=DH_data$vdur_CS[index_spk],y=DH_data$F12_s1_CS[index_spk],
            pch = sapply(DH_data$class[index_spk], function(x) symbol[[x]],USE.NAMES = FALSE),
            col='white',cex=1
            )
        },
    layout = c(3,3), skip = speakers != spk, as.table=TRUE
    )
    print(lp[[spk]], more = spk != tail(speakers,1))
} 
dev.off()



