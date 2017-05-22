### R code from vignette source 'McKayProject.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: McKayProjectHold.Rnw:399-547
###################################################
library(BoomSpikeSlab) #for spike and slab
library(glmnet)        #for LASSO
library(tictoc)        #for timing

## Data
# When it comes to image processing, Matlab is king (even for just loading in images).  For this reason,
#  we loaded the images and constructed the entire dictionary D in Matlab and saved the output as
#  a manageable csv file.
D    <-read.csv('pain_crops/dictionary.csv',head=FALSE) #all vectorized images
d    <-read.csv('pain_crops/label.csv',head=FALSE)      #labels (matrix)
d    <-as.vector(as.matrix(d))                          #organize as vector
imdim<-c(217,163) #image size (CHANGE ACCORDING TO NEEDS)

# We need to (two norm) normalize our vectorized images for spike and slab/LASSO
D    <-scale(D,center=FALSE,scale=apply(D,2,function(x) sqrt(sum(x^2))))#two norm

## Testing (Cross-Validation)
# We will randomly select 4 of the 7 options for each subject and see how well the classification 
#  strategies work at classifying the remaining three
#
# As well, we need to construct a classification function to assign an identification based on our
#  models.  We will use one inspired by Wright et al.
classify <-function(y,X,b,labels){
  picks<-unique(labels) #which to pick from
  L    <-length(picks)  #number of picks
  resi <-matrix(1000,L,1)   #residuals for each pick (1000 to see errors)
  for(j in 1:L){
    temB <-matrix(0,length(b),1) #zero except b from class j
    inds <-which(labels==picks[j]) #get indexes
    temB[inds]<-b[inds]            #get class betas
    resi[j]<-sqrt(sum(y-X%*%temB)) #find residuals
  }
  return(picks[which.min(resi)])
}

faces<-unique(d)  #subject identifications
J    <-10         #number of trials
class<-matrix(0,length(faces),J) #hold classifications
lass <-class                    #hold LASSO's
actu <-class                    #hold actual classes
sptm <-matrix(0,length(faces),J) #hold spike, slab time
latm <-matrix(0,length(faces),J) #hold LASSO time
for(j in 1:J){    #loop through trials
  hold <-rep(0,4*length(faces),1)         #hold onto indexes
  for(k in 1:length(faces)){              #loop through subjects
    inds <-sample(7,4)                    #temporary indexes from subjects
    inds <-inds+(k-1)*7                   #adjust for subject location in vector
    hold[((k-1)*4+1):(k*4)]<-inds         #update hold
  }
  
  # Dictionary construction
  X    <-D[,hold] #"trained" dictionary
  temx <-d[hold]  #"trained" labels (just in case)
  Y    <-D[,-hold]# testing  dictionary
  temy <-d[-hold] # testing  labels (just in case)
  #(Note that we should always have 4 of each id for temx, 3 for temy)
  
  # At this point, if we ran this code it would take a while (neither algorithm is fast).
  #  Therefore, we will shorten our testing sizes to make it manageable.  
  samp <-sample(1:3,length(faces),replace=TRUE) #sample 1 of 3 options per class
  sind <-3*(0:(length(faces)-1))+samp
  S    <-length(sind)
  
  # Formatting for X (lm.spike)
  XX   <-data.frame(X)
  # Number of iterations for spike and slab
  niter<-80;
  
  for(k in sind){
    cat('\n')
    ii   <-which(k==sind)
    cat(round(k/(dim(Y)[2]),4))
    y    <-Y[,k]  # test image
    tic.clearlog() #start timer
    tic()
    tems <-lm.spike(y~0+.,data=XX,niter=niter) #no intercept
    toc(log=TRUE) #Spike timer
    temp <-tic.log(format=FALSE)
    tic.clearlog()
    sptm[ii,j]<-mean(unlist(lapply(temp, function(x) x$toc - x$tic)))
    sbet <-as.numeric(tems$beta[niter,])
    tic() 
    teml <-cv.glmnet(X,y,standardize=FALSE,alpha=1,nfold=15,intercept=FALSE) #no intercept
    toc(log=TRUE) #LASSO timer
    temp <-tic.log(format=FALSE)
    tic.clearlog()
    latm[ii,j]<-mean(unlist(lapply(temp, function(x) x$toc - x$tic)))
    lbet <-as.numeric(coef(teml))[-1] #intercept, though 0, still given
    
    #Classifications
    class[ii,j]<-classify(y,X,sbet,temx)
    lass [ii,j]<-classify(y,X,lbet,temx)
    actu [ii,j]<-temy[k]
  }
}

# pdf('spikecoef.pdf')
# temp=which(temx==5)
# hold=sbet
# nold=hold;told=matrix(NA,length(sbet),1);
# nold[temp]=NA
# told[temp]=sbet[temp]
# plot(sbet,col='red',lwd=1.5,main='Spike and Slab Coefficients',
#      xlab='Coefficient For Each Column',ylab='Values',type='l',
#      ylim=c(-.35,.35))
# lines(told,col='darkgreen',lwd=2)
# grid(10,10)
# dev.off()
# pdf('lasscoef.pdf')
# temp=which(temx==5)
# hold=lbet
# nold=hold;told=matrix(NA,length(lbet),1);
# nold[temp]=NA
# told[temp]=lbet[temp]
# plot(lbet,col='red',lwd=1.5,main='LASSO Coefficients',
#      xlab='Coefficient For Each Column',ylab='Values',type='l',
#      ylim=c(-.35,.35))
# lines(told,col='darkgreen',lwd=2)
# grid(10,10)
# dev.off()

# Organize Results
sprite<-colMeans(class==actu)
larite<-colMeans(lass ==actu)
spface<-matrix(0,length(faces),1)
laface<-matrix(0,length(faces),1)
sptt  <-colMeans(sptm) #times
latt  <-colMeans(latm) #times
for(k in faces){
  ii   <-which(k==faces)
  ind  <-actu==k
  spface[ii]<-mean(class[ind]==k)
  laface[ii]<-mean(lass [ind]==k)
}

# Generate Example Image
esti <-as.matrix(XX)%*%sbet; #image from last iteration
dim(esti)<-imdim #reformat
dim(y)<-imdim #reformat actual
toge <-cbind(y,0,0,0,esti)#put togeter 
png('ExampleConstruction.png')
image(t(toge[nrow(toge):1,]),axes=FALSE)
dev.off()

#barplot(rbind(t(spface),t(laface))*100,beside=T,
#names.arg=LETTERS[1:12],ylim=c(0,100),main='Image By Image Result',
#ylab='Classification %',xlab='Image',legend=c('Spike','LASSO'),
#args.legend=list(x='topright',cex=.8,inset=c(0,-.15),horiz=T))


