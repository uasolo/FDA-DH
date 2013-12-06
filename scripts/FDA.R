#       Licensed under GPLv3


library(fda)
library(lattice)

root_dir = '/home/gubian/Experiments/FDA/lana/' # root dir for this experiment
plots_dir = paste(root_dir,'plots/',sep='')
data_dir =  paste(root_dir,'data/',sep='')
scripts_dir =  paste(root_dir,'scripts/',sep='')

# use pca.fd version from package fda_2.2.5.tar.gz or earlier (you find a copy in the scripts/ dir)
source(paste(scripts_dir,'pca.fd.R',sep=''))
# this is a modified version of the landmarkreg() command 
source(paste(scripts_dir,'landmarkreg.nocurve.R',sep=''))
# this is a slightly modified version of the plot.pca.fd() command,
# but you can also use the standard one.
source(paste(scripts_dir,'plot.pca.fd.corr.R',sep=''))


# filename, speaker, class, durations ('DH' stands for diphthong-hiatus) 
DH_data = read.csv(file = paste(data_dir,"DH_data.csv",sep=''))
n_items = dim(DH_data)[1] # 365
speakers = levels(DH_data$spk) # "AM" "CC" "CP" "CT" "FM" "MC" "MM" "MS" "NM"

# f0 and formants
# all raw contours start 30ms before the beginning of /l/
# samples are every 5ms
# get rid of the first 6 samples (= start at beginning of /l/) for f0
# start at beginning of vowel sequence /ja/ or /i.a/ for formants

time_list = list() # for f0
f0_list = list()
F1bark_list = list()
F2bark_list = list()
len_f0 = c() # number of samples 
len_F1 = c() 
dur_f0 = c() # in ms
dur_F1 = c()
vowel_time = list() # for formants
vowel_time0 = list() # starting from 0

for (i in 1:n_items) {
    data_sample = read.table(paste(data_dir,'pitch_and_formants/',DH_data$filename[i],".pitch",sep=''),h=T)
    time_list[[i]] = data_sample$time[7:length(data_sample$time)] 
    time_list[[i]][1] = 0 # correct possible approx errors
    len_f0 = c(len_f0, length( data_sample$time) - 6) # get rid of the first 6 samples
    dur_f0 = c(dur_f0, time_list[[i]][len_f0[i]])
    sample = data_sample$f0[7:length(data_sample$time)]
    sample = 12 * logb(sample, base = 2)
    f0_list[[i]] = sample - mean(sample) # f0 is in normalized st
    t_vowel = which ( DH_data$ldur[i] <= data_sample$time &
    DH_data$ldur[i] +  DH_data$vdur[i] >= data_sample$time  )
    sample = data_sample$f1bark[t_vowel]
    F1bark_list[[i]] = sample - mean(sample) # formants are in normalized barks
    sample = data_sample$f2bark[t_vowel]
    F2bark_list[[i]] = sample - mean(sample)
    vowel_time[[i]] = time_list[[i]][   which ( DH_data$ldur[i] <= time_list[[i]] & 
    DH_data$ldur[i] +  DH_data$vdur[i] >= time_list[[i]] )  ] 
    vowel_time0[[i]] = vowel_time[[i]] - vowel_time[[i]][1]
    len_F1 = c(len_F1, length(vowel_time0[[i]]))
    dur_F1 = c(dur_F1, vowel_time0[[i]][len_F1[i]])
}

# set some graphic parameters
# this is convenient in order to get the label assignments right when using the groups argument in xyplot, splom etc.
color = list ( 'd' = 'blue','h' = 'orangered')
symbol = list( 'd' = 'D','h' = 'H')
lty = list( 'd' =1,'h' = 3)
lwd  = list( 'd' =1,'h' = 2)
name = list( 'd' ='Diphthong','h' = 'Hiatus')
# change color of trips in lattice plots
strip.background.col = trellis.par.get('strip.background')$col[4]; dev.off() # may open a blank plot window; just close it.

################## Analysis of f0 contours ########################

######## Some exploratory plots

# display raw data

subsamp= runif(n_items) < 0.2 # randomly select 20% of the dataset
# two versions, separated by class or not. Uncomment accordingly.
png(paste(plots_dir,'f0_raw.png',sep=''))
#png(paste(plots_dir,'f0_raw_nocol.png',sep=''))
i=1
plot(time_list[[i]],f0_list[[i]],type = 'n',xlim=c(0,400),ylim=c(-4,6),xlab='time (ms)',ylab='F0 (norm. semitones)',las=1,main = '',cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
#    lines(time_list[[i]],f0_list[[i]],col = 'black', lty=1,lwd=1)
    lines(time_list[[i]],f0_list[[i]],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()

# f0 speaker variability
# select two speakers, show that speaker variability in f0 contours is larger than class variability
    
for (spk in c('CC','MM')) {
    for (class in c('d','h')) {
        index = which(DH_data$spk == spk & DH_data$class == class)
        png(paste(plots_dir,'f0_raw_',spk,'_',class,'.png',sep=''))
        plot(time_list[[index[1]]],f0_list[[index[1]]],type = 'n',xlim=c(0,400),ylim=c(-6,6),xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
        for (i in index) { 
	        lines(time_list[[i]],f0_list[[i]],col = 1)
        }
        dev.off()
    }
}

################## Smoothing ######################################

mean_dur_f0 = mean(dur_f0)

# GCV for smoothing

n_knots_vec <- seq(4,trunc(median(len_f0)/2),4) # explore from 4 knots up to around half the number of samples
loglam_vec <- seq(-4,10,2) # explore lambda from 10^(-4) to 10^10
gcv_err <- array(dim=c(n_items,length(loglam_vec),length(n_knots_vec)),dimnames=c('items','lambda','knots'))
i_sample <- sample(1:n_items,50) # a data subset, to save computation time

# compute GCV error for all (n_knots, lambda) combinations on the i_sample curves, store it in gcv_err (may take some minutes)
for (k in 1:length(n_knots_vec)) {
    for (l in 1:length(loglam_vec)) {
        norm_rng <- c(0,mean_dur_f0)
        knots <- seq(0,mean_dur_f0,length.out = n_knots_vec[k])
        Lfdobj <- 3 # 2 + order of derivative expected to be used. We need velocity to apply Xu's principle, thus order = 1
        norder <- 5 # 2 + Lfdobj 
        nbasis <- length(knots) + norder - 2 # a fixed relation about B-splines
        basis <- create.bspline.basis(norm_rng, nbasis, norder, knots)
        fdPar <- fdPar(basis, Lfdobj, 10^(loglam_vec[l]))
        for (i in i_sample) {
            # apply linear time normalization
            t_norm = (time_list[[i]] / dur_f0[i]) * mean_dur_f0
            gcv_err[i,l,k] = smooth.basis(t_norm,f0_list[[i]],fdPar)$gcv
        }
    }
}    

# compute log of the median of gcv_err for each n_knots and lambda combination
gcv_log_err = log( apply(gcv_err,2:3,median, na.rm=T)) # median to protect from outliers

# plot log GCV errors on a grid
png(paste(plots_dir,'GCV_log_err_f0.png',sep='')) 
col.l <- colorRampPalette(c('blue', 'white'))(30)
levelplot( gcv_log_err, scales=list(y=list(at=1:length(n_knots_vec),labels=n_knots_vec,cex=1.5), x=list(at=1:length(loglam_vec),labels=sapply(loglam_vec,function(x) eval(substitute(expression(10^y) ,list(y=x)) )),cex=1.5) ),xlab = list(label = expression(lambda), cex=2), ylab = list(label= 'k',cex=2),col.regions=col.l,
colorkey=list(label=list(at=-6:-3,label=sapply(-6:-3,function(x) eval(substitute(expression(10^y) ,list(y=x)) )),cex=1.5)),
#aspect = "iso", shrink = c(0.7, 1),
#colorkey=T
)
dev.off()

# min GCV error is in:
argmin = which(gcv_log_err==min(gcv_log_err), arr.ind=TRUE) # rows are lambda indices, cols are n_knots indices
# arg min log lambda is:
10^(loglam_vec[argmin[1]])    
# arg min n_knots is:
n_knots_vec[argmin[2]]

# Inspection of gcv_log_err for f0 shows that:
# min estimated gcv error is obtained basically at the highest number of knots and at low lambda.
# However, the captured detail looks too much (overfitting) for the forthcoming analysis.
# So, lambda and n_knots will be chosen by combining eye inspection of some curves (code below) and the guidance of the GCV figure (above).
# (See the paper for details)


for (loglam in c(-4,2,6,10)) {
    for (n_knots in c(8,28)) {
		lambda = 10^(loglam) 
		norm_rng <- c(0,mean_dur_f0)
		knots <- seq(0,mean_dur_f0,length.out = n_knots)
		Lfdobj <- 3
		norder <- 5
		nbasis <- length(knots) + norder - 2
		basis <- create.bspline.basis(norm_rng, nbasis, norder, knots)
		fdPar <- fdPar(basis, Lfdobj,lambda)
		i=150 # select here a random curve
		t_norm = (time_list[[i]] / dur_f0[i]) * mean_dur_f0
		y_fd = smooth.basis(t_norm,f0_list[[i]],fdPar)$fd
		png(paste(plots_dir,'f0_fit_loglam',loglam,'n_knots',n_knots,'.png',sep=''))
		plot(y_fd,xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='red',lwd=3,ylim=c(-3,5))
		points(t_norm,f0_list[[i]],pch=20)
		#legend('topleft',legend= c(eval(substitute(expression(lambda == 10^x),list(x = loglam))),eval(substitute(expression(k == x),list(x = n_knots))) ),cex=2  )
        dev.off()
    }
}

######## Smoothing using prior knowledge

# From the following paper:
#	author = {Yi Xu and Xuejing Sun},
#	title  = {Maximum speed of pitch change and how it may relate to speech},
#	journal = {J. Acoust. Soc. Am.},
# 	volume = {111},
# 	number = {3},
# 	year = {March 2002},
# 	pages = {1399--1413},

# Xu's empirical equations for rising tones max speed (st and st/s):
# ave. speed = 10.8 + 5.6 * excursion
# max speed  = 12.4 + 10.5 * excursion

# compare solutions in terms of max controllable speed in rising pitch gesture
# so that if we encounter speed values much higher than those found by Xu we consider them
# not realistic, or at least not interesting for us, since they cannot originate from a controlled gesture. 

# looking at the curve fitting plots with the respective velocity curves we can see that in the 
# case of min GCV there are too many 'micro' rising gestures that exhibit peak velocities well above the 
# emirical limits proposed by Xu, whereas the complexity compromise solution looks ok. 

# Plot curves on their original time axis, i.e. not on linearly normalized axis.
# Average and max speed are read on the plots directly by visual inspection.
# Average speed is considered the slope of the f0 contour along a rising gesture.
# Max speed is read directly as the peak of the f0 first derivative curve.


i=150 # select here a random curve
loglam = 6 # change values here 
lambda = 10^(loglam) 
n_knots = 8 # and here
rng_i = range(time_list[[i]])
knots <- seq(rng_i[1],rng_i[2],length.out = n_knots)
Lfdobj <- 3
norder <- 5
nbasis <- length(knots) + norder - 2
basis_i <- create.bspline.basis(rng_i, nbasis, norder, knots)
fdPar_i <- fdPar(basis_i, Lfdobj,lambda)
y_fd = smooth.basis(time_list[[i]],f0_list[[i]],fdPar_i)$fd
png(paste(plots_dir,'f0_orig_time_loglam',loglam,'n_knots',n_knots,'.png',sep=''))
plot(y_fd,xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='red',lwd=3,ylim=c(-3,5))
legend('topleft',legend= c(eval(substitute(expression(lambda == 10^x),list(x = loglam))),eval(substitute(expression(k == x),list(x = n_knots))) ),cex=2  )
dev.off()
# plot first derivative (note: time axis is in ms, should convert to s in order to get st/s on the y axis).
# y(t), t in ms. If T in s, then t = 1000*T. dy(1000*T)/dT = 1000* dy(T)/dT = 1000* dy(t)/dt 
png(paste(plots_dir,'Df0_orig_time_loglam',loglam,'n_knots',n_knots,'.png',sep=''))
plot(1000*y_fd,Lfdobj=1,xlab='time (ms)',ylab='st/s',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='red',lwd=3)
legend('bottomright',legend= c(eval(substitute(expression(lambda == 10^x),list(x = loglam))),eval(substitute(expression(k == x),list(x = n_knots))) ),cex=2  )
dev.off()




# selected values:
lambda = 10^(6) ; n_knots = 8
# build global f0 fd object
norm_rng <- c(0,mean_dur_f0)
knots <- seq(0,mean_dur_f0,length.out = n_knots)
Lfdobj <- 3
norder <- 5
nbasis <- length(knots) + norder - 2
basis <- create.bspline.basis(norm_rng, nbasis, norder, knots)
fdPar <- fdPar(basis, Lfdobj,lambda)
# convenient aliases
basis_f0 = basis
fdPar_f0 = fdPar
# smooth.basis() does not accept different time samples for different curves.
# Thus we create smooth curves one by one on the same basis, store the spline coefficients and compose an fd object at the end.
f0_coefs = matrix(nrow = nbasis, ncol = n_items)
for (i in 1:n_items) {
    t_norm = (time_list[[i]] / dur_f0[i]) * mean_dur_f0
    f0_coefs[,i] = c(smooth.basis(t_norm,f0_list[[i]],fdPar)$fd$coefs)
}
f0_fd = fd(coef=f0_coefs, basisobj=basis)
# curves are linearly time normalized, their duration is mean_dur_f0

# plot the curves
png(paste(plots_dir,'f0_lin.png',sep=''))
#png(paste(plots_dir,'f0_lin_nocol.png',sep=''))
plot(c(0,mean_dur_f0),c(-4,6),type='n',xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
    #lines(f0_fd[i],col = 'black', lty=1,lwd=1)
    lines(f0_fd[i],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()


# this is how the B-spline basis looks
basis_fd = fd(diag(1,nbasis),basis)
png(paste(plots_dir,'B-splines.png',sep=''))
plot(norm_rng,c(0,1),type='n',xlab='time (s)', ylab = '', las=1,cex.axis=1.3,cex.lab=1.3)
for (b in 1:nbasis) {
    lines(basis_fd[b],col='red',lty=2)
}
points(knots,rep(0,n_knots),pch=19,col='blue')
dev.off()


# let us use these B-splines to represent the i-th curve
i=150 # select here a random curve
t_norm = seq(0,mean_dur_f0,length.out = length(f0_list[[i]])) 
y = f0_list[[i]]
y_fd = smooth.basis(t_norm,y,fdPar)$fd
# this is how the splines combine (sum) to approximate the given curve samples
png(paste(plots_dir,'B-splines_smoothing.png',sep=''))
plot(y_fd,lwd=2,col='red',xlab='time (s)', ylab = 'norm. st', las=1,cex.axis=1.3,cex.lab=1.3)
points(t_norm,y,pch=20,col='black')
for (b in 1:nbasis) {
    lines(y_fd$coefs[b] * basis_fd[b],col='red', lty=2)
}
dev.off()

################## Landmark Registration ###############################

# Use landmarkreg.nocurve(), a modified version of the landmarkreg() command. 
# It places knots according to de Boor's theorem, i.e. at landmark positions.
# It operates only on the landmark positions, not on the curves.
# It provides only the time warping curves, which have to be applied to the curves later on.
# It provides also relative rate curves (not used here).

# landmark matrix: 
# one internal landmark: end of /l/
# landmarkreg.nocurve() requires also beginning and end of the token to be included in the landmark matrix.
# because in general the total duration may differ (not in this case, since linear registration already occured).

land = matrix(nrow = n_items, ncol = 3) # one internal landmark + begin and end
for (i in 1:n_items) {
	land[i,] = c(0,DH_data$ldur[i]/dur_f0[i],1) * mean_dur_f0
} 
reg = landmarkreg.nocurve(land, nhknots = n_knots) 
# nhknots are the used for the representation of h(t), not for the actual time warping
# other arguments are left at default values, since in this case registration is easy, having only one landmark (see command code for details).
# Registration may take some minutes.

# fd object for registered f0 contours
f0reg_coefs =  matrix(nrow = nbasis, ncol = n_items)
reg_fdPar = fdPar(basis, Lfdobj,1e-12) # lambda small, since smoothing already occurred
# reg$hfunmat is a matrix whose i-th column contains the time samples h(x) for the i-th curve,
# where x (reg$x) are regularly spaced time samples along the registered time axis and h() is the time warping function 
# that returns the original time axis points.
for (i in 1:n_items) {
	h_i = reg$hfunmat[,i]
	f0reg_coefs[,i] = c(smooth.basis(reg$x, eval.fd(h_i,f0_fd[i]),reg_fdPar)$fd$coefs)
}
f0reg_fd = fd(coef=f0reg_coefs, basisobj=basis)

# Graphical parameters for landmark labels: place a label in the middle of every interval.
landlab = c('/l/', '/ja/ or /i.a/') 
at_land = c() # position of the label along the time axis 
for (i in 1:(length(reg$land)-1)) {
    at_land = c(at_land, mean(reg$land[i:(i+1)]))
}

#png(paste(plots_dir,'f0_reg.png',sep=''))
png(paste(plots_dir,'f0_reg_nocol.png',sep=''))
plot(c(0,mean_dur_f0),c(-4,6),type='n',xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
    lines(f0reg_fd[i],col = 'black', lty=1,lwd=1)
    #lines(f0reg_fd[i],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
abline(v=reg$land[2],lty=2,lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
#legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()



# inverse h(t)
h_inv_list = list()
for (i in 1:n_items) {
	rng <- range(reg$land)
	steps <- seq(rng[1],rng[2],len=50)
	knots <- as.numeric(eval.fd(steps,reg$warpfd[i]))
	# rounding error quick fixes
	knots[1] = 0
	knots[length(steps)] = rng[2]
	norder <- 4
	nbasis <- length(knots) + norder - 2
	basis <- create.bspline.basis(rng, nbasis, norder, knots)
	Lfdobj <- 2
	lambda <- 10^1
	fdPar <- fdPar(basis, Lfdobj, lambda)
	h_inv_list[[i]] <- smooth.basis(knots,steps,fdPar)$fd
}

subsamp_small= (1:n_items)[runif(n_items) < 0.02] # select 2% of the dataset

# plot h(t) for some curves, show alignment 
png(paste(plots_dir,'h_sample.png',sep=''))
plot(reg$warpfd[subsamp_small],lty=1,lwd=1,las=1,xlab='reg. time (ms)',ylab='lin. norm. time (ms)',cex.axis=1.3,cex.lab=1.3,col='black')
for (i in subsamp_small) {
points(eval.fd(land[i,2],h_inv_list[[i]]),land[i,2],col='red',pch=19,cex=1)
}
abline(v=reg$land[2],lty=2,col='black',lwd=1)
dev.off()



# plot some linearly registered curves, show landmark position
png(paste(plots_dir,'registration_lin.png',sep=''))
plot(range(reg$land),c(-4,6),type='n',xlab='lin. norm. time (ms)',ylab='F0 (norm. semitones)',las=1,ylim=c(-4,6),cex.lab=1.3,cex.axis=1.3)
for (i in subsamp_small) {
    lines(f0_fd[i],lty=1,col=1)
    
	points(land[i,2],eval.fd(land[i,2],f0_fd[i]),col='red', pch=19,cex=1.3)
}
dev.off()
# plot the same curves after registration
png(paste(plots_dir,'registration_land.png',sep=''))
plot(range(reg$land),c(-4,6),type='n',xlab='reg. time (ms)',ylab='F0 (norm. semitones)',las=1,ylim=c(-4,6),cex.lab=1.3,cex.axis=1.3) 
for (i in subsamp_small) {
    lines(f0reg_fd[i],lty=1,col=1)
	t_reg = eval.fd(land[i,2],h_inv_list[[i]])
	points(t_reg,eval.fd(t_reg,f0reg_fd[i]),col='red', pch=19,cex=1.3)
}
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
dev.off()


################## Functional PCA on f0 contours ########################


y_fd = f0reg_fd # alias
# usually a good solution is obtained by setting the same lambda and knots (thus basis) used for smoothing
lambda_pca    <- lambda
pcafdPar  <- fdPar(basis_f0, 2, lambda_pca)
f0_pcafd <- pca.fd(y_fd, nharm=3, pcafdPar) # first three PCs
#f0_pcafd = f0_pcafd # alias

# Invert the sign of PC2, in order to facilitate global analysis, i.e. all features will increase from D to H.
# This is not necessary in general, but it can help interpreting results when multiple dimensions are involved.
# (Functional) PCA sings can always be inverted, since they have no intrinsic meaning.
f0_pcafd$harmonics$coefs[,2] = - f0_pcafd$harmonics$coefs[,2]
f0_pcafd$scores[,2] = - f0_pcafd$scores[,2]

# plot PC curves
plot.pca.fd.corr(f0_pcafd,xlab = 'time (ms)',ylab='F0 (norm. semitones)',land = reg$land , nx=40,plots_dir = plots_dir, basename = 'PCA_f0reg_',height=480)

# plot PC scores 
png(paste(plots_dir,'PCsplom_f0reg.png',sep=''))
splom(f0_pcafd$scores ,
groups=DH_data$class,
# in lattice plot functions, the following complex sapply() expression is necessary
# in order to get the order of groups graphical parameters right.
pch  = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
col  = sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
cex=0.8, varnames= c(expression(s[1]),expression(s[2]),expression(s[3])) )
dev.off()

# plot only the first two PC scores
# grouped by class
png(paste(plots_dir,'PCscatter_f0reg.png',sep=''))
xyplot(f0_pcafd$scores[,2] ~  f0_pcafd$scores[,1] , cex=1.5,
xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2), 
 groups= DH_data$class,
pch  = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
col  = sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
,scales = list(cex=1.5)
)
dev.off()

# the original PCA output, i.e. no classes
png(paste(plots_dir,'PCscatter_f0reg_allblack.png',sep=''))
xyplot(f0_pcafd$scores[,2] ~  f0_pcafd$scores[,1] , cex=1,
xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2), 
col  = 'black',pch=20,scales = list(cex=1.5)
)
dev.off()

# PC scores by class and speaker
# see http://tolstoy.newcastle.edu.au/R/e2/help/07/09/24852.html for the use of panel 
png(paste(plots_dir,'PCscatter_f0reg_speaker.png',sep=''))
xyp = xyplot(f0_pcafd$scores[,2] ~  f0_pcafd$scores[,1] | DH_data$spk , groups= DH_data$class,
xlab = list(label=expression(s[1]),cex=1.5), ylab= list(label=expression(s[2]),cex=1.5),cex=1,
	col = sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
	pch = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
panel = panel.superpose,
panel.groups = function(...) {
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

# boxplots for PC scores 
# (manually put s1 or s2, scores[,1 or 2] and s[1] or s[2] accordingly)

png(paste(plots_dir,'s1_f0reg_spk_box.png',sep=''))
bwp = bwplot(  f0_pcafd$scores[,1] ~ DH_data$class | DH_data$spk, ylab = expression(s[1]))
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()


# plot class-specific mean curves

t_f0 = reg$x
f0_pcafd = f0_pcafd 

#png(paste(plots_dir,'f0_mean.png',sep=''))
png(paste(plots_dir,'f0_lm_f0_s2.png',sep=''))
plot(c(0,mean_dur_f0),c(-3,5),type='n',xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
#for (class in names(color)) {
#    lines(f0_pcafd$meanfd + mean(f0_pcafd$scores[which(DH_data$class == class),1]) * f0_pcafd$harmonics[1] + mean(f0_pcafd$scores[which(DH_data$class == class),2]) * f0_pcafd$harmonics[2],col = color[[class]], lwd = lwd[[class]], lty = lty[[class]])
#}

# or use values from linear model f0_s2_class.lm (see script DA.R)
mean_lm_f0_s2 = list(d = -3.8, h = 3.8)
for (class in names(color)) {
    lines(f0_pcafd$meanfd + mean_lm_f0_s2[[class]]  * f0_pcafd$harmonics[2],col = color[[class]], lwd = lwd[[class]], lty = lty[[class]])
}
abline(v=reg$land[2],lty=2)
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()

# plot class- and speaker-specific mean curves
# see xyplot.ts
png(paste(plots_dir,'s2_f0_mean_spk.png',sep=''))
table_plot = expand.grid(class = c('d','h'),spk = speakers,stringsAsFactors = FALSE)

curves = matrix(nrow = length(t_f0),ncol = nrow(table_plot))
for (i in 1:nrow(table_plot)) {
    curve = f0_pcafd$meanfd +
        # choose which PC you want to include
	    #mean(f0_pcafd$scores[which(DH_data$class == table_plot$class[i] & DH_data$spk == table_plot$spk[i]),1]) * f0_pcafd$harmonics[1] #+
	    mean(f0_pcafd$scores[which(DH_data$class == table_plot$class[i] & DH_data$spk == table_plot$spk[i]),2]) * f0_pcafd$harmonics[2] 
    curves[,i] = eval.fd(t_f0,curve)
}
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


# store PC scores in DH_data_FDA
DH_data_FDA = DH_data
DH_data_FDA$f0_s1 = f0_pcafd$scores[,1]
DH_data_FDA$f0_s2 = f0_pcafd$scores[,2]





## FPCA-based reconstruction example (6 plots)

png(paste(plots_dir,'mean','.png',sep=''))
plot(f0_pcafd$meanfd,xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='black',lwd=3,ylim=c(-4,5))
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
dev.off()


png(paste(plots_dir,'PC1','.png',sep=''))
plot(f0_pcafd$harmonics[1],xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='black',lwd=3,)
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
dev.off()


png(paste(plots_dir,'PC2','.png',sep=''))
plot(f0_pcafd$harmonics[2],xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='black',lwd=3,)
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
dev.off()

i = 122
png(paste(plots_dir,'reconstr_mean','.png',sep=''))
plot(f0_pcafd$meanfd,xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='black',lwd=3,ylim=c(-4,5))
lines(f0reg_fd[i],lwd=2, lty=2)
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
legend('topleft',legend=c('original','reconstruction'),lty=c(2,1),lwd=c(2,3))
dev.off()



png(paste(plots_dir,'reconstr_mean_PC1','.png',sep=''))
plot(f0_pcafd$meanfd + f0_pcafd$scores[i,1] * f0_pcafd$harmonics[1] ,xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='black',lwd=3,ylim=c(-4,5))
lines(f0reg_fd[i],lwd=2, lty=2)
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
legend('topleft',legend=c('original','reconstruction'),lty=c(2,1),lwd=c(2,3))
dev.off()


png(paste(plots_dir,'reconstr_mean_PC1_PC2','.png',sep=''))
plot(f0_pcafd$meanfd + f0_pcafd$scores[i,1] * f0_pcafd$harmonics[1] + f0_pcafd$scores[i,2] * f0_pcafd$harmonics[2],xlab='time (ms)',ylab='F0 (norm. semitones)',main = '',las=1,cex.axis=1.5,cex.lab=1.5,col='black',lwd=3,ylim=c(-4,5))
lines(f0reg_fd[i],lwd=2, lty=2)
abline(v=reg$land[2],lty=2,col='black',lwd=1)
axis(3,tick=F,at=at_land, labels=landlab,cex.axis=1.5)
legend('topleft',legend=c('original','reconstruction'),lty=c(2,1),lwd=c(2,3))
dev.off()


################## Analysis of formant contours ########################

######## Some exploratory plots

# display raw data

# F1
png(paste(plots_dir,'F1_raw.png',sep=''))
i=1
plot(vowel_time0[[i]],F1bark_list[[i]],type = 'n',xlim=c(0,200),ylim=c(-4,4),xlab='time (ms)',ylab='F1 (norm. barks)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
    lines(vowel_time0[[i]],F1bark_list[[i]],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()


# F2
png(paste(plots_dir,'F2_raw.png',sep=''))
i=1
plot(vowel_time0[[i]],F2bark_list[[i]],type = 'n',xlim=c(0,200),ylim=c(-4,4),xlab='time (ms)',ylab='F2 (norm. barks)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
    lines(vowel_time0[[i]],F2bark_list[[i]],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()


################## Smoothing ###########################################

mean_dur_F1 = mean(dur_F1)
# build global formants fd object
# use lambda and n_knots from f0 smoothing (or repeat parameter selection).
norm_rng <- c(0,mean_dur_F1)
knots <- seq(0,mean_dur_F1,length.out = n_knots)
Lfdobj <- 3
norder <- 5
nbasis <- length(knots) + norder - 2
basis_F12 <- create.bspline.basis(norm_rng, nbasis, norder, knots)
fdPar_F12 <- fdPar(basis_F12, Lfdobj,lambda)

# store spline coefficients. Note: F12_coefs has a further dimension with respect to f0_coefs, because formants are 2dim trajectories.
F12_coefs = array(dim = c(nbasis,n_items,2))
for (i in 1:n_items) {
    t_norm = (vowel_time0[[i]] / dur_F1[i]) * mean_dur_F1
    F12_coefs[,i,1] = c(smooth.basis(t_norm,F1bark_list[[i]],fdPar_F12)$fd$coefs)
    F12_coefs[,i,2] = c(smooth.basis(t_norm,F2bark_list[[i]],fdPar_F12)$fd$coefs)
}
F12_fd = fd(coef=F12_coefs, basisobj=basis_F12)


# plot the curves
png(paste(plots_dir,'F1_lin.png',sep=''))
plot(c(0,mean_dur_F1),c(-4,4),type='n',xlab='time (ms)',ylab='F1 (norm. barks)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
     lines(F12_fd[i,1],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()

png(paste(plots_dir,'F2_lin.png',sep=''))
plot(c(0,mean_dur_F1),c(-4,4),type='n',xlab='time (ms)',ylab='F2 (norm. barks)',main = '',las=1,cex.axis=1.5,cex.lab=1.5)
for (i in (1:n_items)[subsamp]) {
     lines(F12_fd[i,2],col = color[[DH_data$class[i]]], lty=lty[[DH_data$class[i]]], lwd=lwd[[DH_data$class[i]]])
}
legend('topleft',legend=unlist(name),col=unlist(color),lwd=unlist(lwd),lty=unlist(lty),cex=1.5)
dev.off()


# Note: no landmark registration, since we have no landmarks inside the vowel cluster.


################## Functional PCA on formant contours ########################

y_fd = F12_fd # alias
# usually a good solution is obtained by setting the same lambda and knots (thus basis) used for smoothing
lambda_pca    <- lambda
pcafdPar  <- fdPar(basis_F12, 2, lambda_pca)
F12_pcafd <- pca.fd(y_fd, nharm=3, pcafdPar) # first three PCs

# plot PC curves
plot.pca.fd.corr(F12_pcafd,xlab = 'norm. time',ylab='norm. barks',land = NULL ,nx=40,plots_dir = plots_dir, basename = 'PCA_F12_',height=480)

# plot PC scores grouped by class
png(paste(plots_dir,'PCscatter_F12.png',sep=''))
xyplot(F12_pcafd$scores[,2] ~  F12_pcafd$scores[,1] , cex=1.5,
xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2), 
 groups= DH_data$class,
pch  = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
col  = sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
,scales = list(cex=1.5)
)
dev.off()


# PC scores by class and speaker
# see http://tolstoy.newcastle.edu.au/R/e2/help/07/09/24852.html for the use of panel 
png(paste(plots_dir,'PCscatter_F12_speaker.png',sep=''))
xyp = xyplot(F12_pcafd$scores[,2] ~  F12_pcafd$scores[,1] | DH_data$spk , groups= DH_data$class,
xlab = list(label=expression(s[1]),cex=1.5), ylab= list(label=expression(s[2]),cex=1.5),cex=1,
	col = sapply(levels(DH_data$class), function(x) color[[x]],USE.NAMES = FALSE),
	pch = sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),
panel = panel.superpose,
panel.groups = function(...) {
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


# boxplots for PC scores 
# (change s1 and s2 accordingly)

png(paste(plots_dir,'s2_F12_spk_box.png',sep=''))
bwp = bwplot(  F12_pcafd$scores[,2] ~ DH_data$class | DH_data$spk, ylab = expression(s[2]))
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)),
        as.table=TRUE
	    )
dev.off()



# class-dependent F12 mean curves 
mean_lm_F12_s1 = list(d = -2.3, h = 2.3)


t_F12 = seq(0,mean_dur_F1,length.out = 50)
table_plot_F12 = expand.grid(class = c('d','h'),F12=2:1)
curves_F12 = matrix(nrow = length(t_F12),ncol = nrow(table_plot_F12))
for (i in 1:nrow(table_plot_F12)) {
    coefs = F12_pcafd$meanfd$coefs[,1,table_plot_F12$F12[i]] +
        mean(F12_pcafd$scores[which(DH_data$class == table_plot_F12$class[i]),1]) * F12_pcafd$harmonics$coefs[,1,table_plot_F12$F12[i]]
+       mean(F12_pcafd$scores[which(DH_data$class == table_plot_F12$class[i]),2]) * F12_pcafd$harmonics$coefs[,2,table_plot_F12$F12[i]]
  #      mean_lm_F12_s1[[table_plot_F12$class[i]]]  *  F12_pcafd$harmonics$coefs[,1,table_plot_F12$F12[i]]
     curves_F12[,i] = eval.fd(t_F12,fd(coefs,basis_F12))
}

png(paste(plots_dir,'F12_mean.png',sep=''))
#png(paste(plots_dir,'F12_lm_F12_s1.png',sep=''))
xyp =	xyplot(
	ts(data=curves_F12,start=t_F12[1],deltat=t_F12[2]-t_F12[1]),
	screens=with( table_plot_F12, paste('F',F12,sep='')),
    col = sapply(table_plot_F12$class, function(x) color[[x]],USE.NAMES = FALSE),
	lty = sapply(table_plot_F12$class, function(x) lty[[x]],USE.NAMES = FALSE),
	lwd = sapply(table_plot_F12$class, function(x) lwd[[x]],USE.NAMES = FALSE),
    layout = c(1,2),
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
    	    )
    	)
dev.off()

# plot class- and speaker-specific mean curves
# see xyplot.ts

table_plot_F12 = expand.grid(class = c('d','h'),F12=1:2,spk = speakers)
n_col_plot = 3
n_combined_rows_plot = length(speakers)/n_col_plot # 3 combined row means a row with F1 and F2 one under the other


panel.order = c()
for (c in 1:n_combined_rows_plot) {
    # select rows with spk in the current block (plot combined row)
    table_plot_F12_subset = subset(table_plot_F12,spk %in% levels(spk)[(c-1)*n_col_plot+(1:n_col_plot)])
    panel.order = c(panel.order,as.integer(row.names(table_plot_F12_subset[with(table_plot_F12_subset,order(-F12,spk,class)),])))
}
# select which mean scores you want to add, and change plot name accordingly
table_plot_F12 = table_plot_F12[panel.order,]
curves_F12 = matrix(nrow = length(t_F12),ncol = nrow(table_plot_F12))
for (i in 1:nrow(table_plot_F12)) {
    coefs = F12_pcafd$meanfd$coefs[,1,table_plot_F12$F12[i]] +
   #      mean(F12_pcafd$scores[which(DH_data$class == table_plot_F12$class[i] & DH_data$spk == table_plot_F12$spk[i]),1]) * F12_pcafd$harmonics$coefs[,1,table_plot_F12$F12[i]] 
        mean(F12_pcafd$scores[which(DH_data$class == table_plot_F12$class[i] & DH_data$spk == table_plot_F12$spk[i]),2]) * F12_pcafd$harmonics$coefs[,2,table_plot_F12$F12[i]] 
     curves_F12[,i] = eval.fd(t_F12,fd(coefs,basis_F12))
}

png(paste(plots_dir,'s2_F12_mean_spk.png',sep=''))
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


# store PC scores in DH_data_FDA
DH_data_FDA$F12_s1 = F12_pcafd$scores[,1]
DH_data_FDA$F12_s2 = F12_pcafd$scores[,2]

# save a richer table with metadata and FDA-based features
write.csv(DH_data_FDA,file=paste(data_dir,'DH_data_FDA.csv',sep=''),row.names =FALSE)

################## Analysis of duration ########################

######## Some exploratory plots


png(paste(plots_dir,'vdur_box.png',sep=''))
bwp = bwplot(  DH_data$vdur ~ DH_data$class, ylab = 'vowel seq duration (ms)')
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3)
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) name[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)) 
	    )
dev.off()


png(paste(plots_dir,'vdur_spk_box.png',sep=''))
bwp = bwplot(  DH_data$vdur ~ DH_data$class | DH_data$spk, ylab = 'vowel seq duration (ms)')
update	(bwp, par.settings=list(
	        par.xlab.text = list(cex=1.3),
	        par.ylab.text = list(cex=1.3),
	        strip.background = list(col=strip.background.col)            
            ),
        scales=list(x=list(labels=sapply(levels(DH_data$class), function(x) symbol[[x]],USE.NAMES = FALSE),cex=1.3),y=list(cex=1.3)), 
        as.table=TRUE
	    )

dev.off()

save.image(paste(scripts_dir,'FDA.RImage',sep=''))
