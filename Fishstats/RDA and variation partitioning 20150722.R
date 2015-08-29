library(vegan)

setwd('g:/FCfish_data')
#loading fish data in spring 2014
fcf<-read.csv('FCfish_biomass_spring2014.csv',row.names=1)# have done
# standardising the community data to decrease the assymmatics 
fcf.biol<-decostand(log1p(fcf),method='hellinger',MARGIN=1)
# loading the envdata in spring
env<-read.csv('Env_data0305spring.csv',row.names=1)
names(env)
# Principal Coordinates of Neighbourhood Matrix of Fangcheng fish community
sites <- read.csv('sites_lonlat20150224.csv',row.names=1)
head(sites)
#pcnm cannot have the sequence or site number

#a PCoA on a neighbour matrix will (typically) produce more eigenvectors
#relative to the same analysis on a standard distance matrix. 
#All resulting eigenvectors with positive eigenvalues may be used 
#as a new set of explanatory, spatial variables in either a multiple
#regression approach (for univariate response data) or a multivariate constrained analysis. 
site.pcnm <- pcnm(dist(sites))

site.pcnm$vectors
# Write CSV in R
write.csv(site.pcnm$vectors, file = "site.pcnm.csv")
site.pcnm <-read.csv('site.pcnm.csv',row.names=1)
View(site.pcnm)
mod<-varpart(fcf.biol,~.,site.pcnm,data=env)
mod
showvarparts(2)
plot(mod)
#warning message of the collinearity detected, so to remove this
library(perturb)
x<-colldiag(env,scale=TRUE,center=FALSE,add.intercept=T)

print(x,dec.places=3,fuzz=NULL,fuzzchar=".")
#PO4,and CHLAa have a large proportion  50% or more
data.matrix(env)# convert a data frame to a numeric matrix
env<-env[,1:8]
env
x2 <- colldiag(site.pcnm,scale=TRUE,center=FALSE,add.intercept=T)
print(x2,dec.places=3,fuzz=NULL,fuzzchar=".")
# no collinearity detected
#reruning the varpart command
mod<-varpart(fcf.biol,~.,site.pcnm,data=env)
mod
#env and site.pcnm may have collinearlity

xb<-read.csv("ENV_data.csv")
View(xb)
#use stepwise VIF selection to solve collinearity
#stepwise VIF function used below
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  for(val in names(in_frame)){
    form_in<-formula(paste(val,' ~ .'))
    vif_init<-rbind(vif_init,c(val,VIF(lm(form_in,data=in_frame,...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(names(in_frame))
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      
      for(val in names(in_dat)){
        form_in<-formula(paste(val,' ~ .'))
        vif_add<-VIF(lm(form_in,data=in_dat,...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}
# check the collinearity ,and improve and select the proper variables

vif_func(in_frame=xb,thresh=10,trace=T)
View(xb)
names(xb)
myvars <-c("DO", "salt","pH","PO4","CHLa","PCNM1",
           "PCNM2","PCNM3","PCNM4","PCNM5",
           "PCNM6" ,"PCNM8","PCNM9","PCNM10" , "PCNM11")
newdata<-xb[,myvars]
newdata
env<- newdata[,c("DO", "salt","pH","PO4","CHLa")]
site.pcnm<-newdata[,c("PCNM1","PCNM2","PCNM3","PCNM4",
                      "PCNM5","PCNM6" ,"PCNM8","PCNM9","PCNM10" , "PCNM11")]
mod<-varpart(fcf.biol,~.,site.pcnm,data=env)
mod
showvarparts(2)
plot(mod,digits=2)
print(x3)
#fraction all[a+b+c]
rda.all<-rda(fcf.biol~.,data=Env_data2)
#fraction[a]
rda.env.geo<-rda(fcf.biol~depth+transp+t+DO+salt+pH+CHLa+Condition(PCNM1+PCNM2),Env_data2)
#fraction[b]
rda.geo.env<-rda(fcf.biol~ Condition(depth+transp+t+DO+salt+pH+CHLa)+PCNM1+PCNM2,Env_data2)
#for completeness, let's define also models iwth simple (marginal)effect of each explanatory variable
# fraction [a+b]
rda.env <- rda(fcf.biol~depth++transp+t+DO+salt+pH+CHLa,data=Env_data2)
# fractions [b+c]:
rda.geo <-rda(fcf.biol~PCNM1+PCNM2,data=Env_data2)
## fraction [a+b+c]
RsquareAdj (rda.all)
##fraction[a]
RsquareAdj (rda.env.geo)
## fraction [c]
RsquareAdj (rda.geo.env)
#Now, let's see the same, using function varpart 
varp <- varpart(fcf.biol,~ sed+depth+transp+t+DO+salt+pH+DIN+PO4+CHLa,~PCNM1+PCNM2+PCNM3+PCNM4+PCNM5+PCNM6+PCNM7,Env_data2)
anova(rda.env.geo)
anova(rda.geo.env)
anova(rda.all)
summary(rda.all)#all explainary variables
summary(rda.env.geo)# environment
summary(rda.geo.env)# geography
vegandocs
