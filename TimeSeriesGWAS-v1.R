#-------------------------------------------------
# This program is developed for the Borevitz lab
# You may use it without any warranty
# Distribution is not permitted
# Copyright - Riyan Cheng 2016
# Edited by Tim Brown and Jared Streich 2017-09
#-------------------------------------------------

#
# This program processed raw phenotype data from the TraitCapture pipeline and performs time-based GWAS
#
# Notes:
# 1) This program runs one data file at a time
# 2) Assumes from input a TraitCapture pipeline output plant area csv file
#     Row 1: plantid; row 2: pot number
#     Column 1: timestamp 
# 3) Plant layout and design file is supplied as a csv named in the format "expid-traitcapture-db-import-Plants.csv" about design
#    e.g. "BVZ0063-traitcapture-db-import-Plants.csv"
# 4) Manually select days 'ud' to reliably smooth if needed
# 5) Set "rootPath" variable below to the parent folder for where your code and data are. 
#    All input files should be in "rootPath/input-data" folder. 
#    Output files will go in the "rootPath/ouput" folder


#############################################################
# basic functions #
#############################################################
strf<- function(strList, pos=1){
   strTmp<- NULL
   for(i in 1:length(strList)){
      strTmp<- c(strTmp,strList[[i]][pos])
   }
   strTmp
}

# smoothing span
nf<- function(n){
   5+(30-4)/(40-5)*(n-5)
}

# remove outliers
drop1f<- function(dat, min.n=5, span=NA){
# dat: data.frame(x=, y=)
# min.n: min no. of useable observations
# span: pass to 'loess'
   drop<- NA
   predicted<- NA
   if(sum(!is.na(dat$y)) >= min.n){
      n<- sum(!is.na(dat$y))
      if(missing(span) || is.na(span))
         spn<- nf(n)/n
      lf<- loess(y ~ x, data=dat, na.action=na.exclude, degree=2, span=spn,
         control=loess.control(surface="direct"))
      rm(n)
      lfp<- predict(lf, dat$x)

      rss<- Inf
      cv<- 0
      for(i in 1:length(dat$y)){
         if(is.na(dat$y[i])) next

         .dtt<- dat; .dtt$y[i]<- NA
         .n<- sum(!is.na(.dtt$y))
         if(missing(span) || is.na(span))
            spn<- nf(.n)/.n
         .lf<- loess(y ~ x, data=.dtt, na.action=na.exclude, degree=2, span=spn,
               control=loess.control(surface="direct"))
         .lfp<- predict(.lf, .dtt$x)
         .sd<- sd(dat$y-.lfp, na.rm=TRUE) #too stringent to use .dtt
         .cv<- abs(dat$y-.lfp)[i]/(.sd + sqrt(.Machine$double.eps))
         if(FALSE){
            # large but minimize RSS
            .sd<- sd(.dtt$y-.lfp, na.rm=TRUE)
            if(.cv > qnorm(1-0.05/.n) && .sd < rss){
               lfp<- .lfp
               drop<- i
               rss<- .sd
            }
         }else{
            if(.cv > max(cv,qnorm(1-0.05/.n))){
               lfp<- .lfp
               drop<- i
               cv<- .cv
            }
         }

         rm(.n)
      }
   }
   if(!is.na(drop)){
      dat$y[drop]<- NA
      predicted<- lfp
   }
   list(data=dat, predicted=predicted, drop=drop)
}

# fit loess per day and then select the max among middle 50%
procDayDataf<- function(pdat, dl=10, min.n=0, span=NA){
# pdat: time points by row, samples by column
#    except the first column timestamp
# dl: extract day string with length of dl
# min.n: min number of observations in a day; or, no data
   itv<- table(substring(pdat$time,15,16))
      itv<- names(itv)[itv > mean(itv)/3]
      itv<- sort(as.numeric(itv),decreasing=FALSE)
      itv<- min(diff(itv))
      itv<- seq(0,60-1e-8,by=itv)
   hrs<- table(substring(pdat$time,12,13))
      hrs<- names(hrs)[hrs > mean(hrs)/3]
      hrs<- sort(as.numeric(hrs),decreasing=FALSE)
   nr<- length(hrs)
   tt<- c(paste(sprintf("%02d",min(hrs)),"_00",sep=""),paste(sprintf("%02d",max(hrs)),"_30",sep=""))
   tm<- paste(sprintf("%02d",rep(hrs,rep(length(itv),nr))), sprintf("%02d",rep(itv,nr)), sep="_")
      tm<- tm[match(tt[1],tm):match(tt[2],tm)]
   day<- sapply(as.character(pdat$time), substr, 1, dl)
   ud<- table(day)
      ud<- names(ud)[ud > mean(ud)/3]
      ud<- sort(ud, decreasing=FALSE)
   dayTime<- paste(rep(ud,rep(length(tm),length(ud))), rep(tm,length(ud)), sep="_")
   dt<- matrix(NA, nrow=length(dayTime), ncol=ncol(pdat)-1)
   rownames(dt)<- dayTime
   colnames(dt)<- colnames(pdat)[-1]
   curTime <- Sys.time()
   timeElapsed = Sys.time()-curTime
   for(j in 2:ncol(pdat)){
     curInd<-j
     curDiff<-(Sys.time()-curTime)
     timeElapsed = timeElapsed + curDiff
     cat(paste("Now on plant", curInd, "of", (ncol(pdat))-1), "\n Time to process last plant:",curDiff, "\n Time passed:", (timeElapsed/60), " | Estimated time left:", ((ncol(pdat) - j)* curDiff)/60,"min. \n\n")
     curTime <- Sys.time()
     if((j-1)%%10 + 1 == 1) cat(j-1) else cat(".")
      for(i in 1:length(ud)){
         #cat(".")
         dat<- pdat[day==ud[i],]
         if(nrow(dat) < min.n){
            #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
            next
         }
         tmTmp<- sapply(as.character(dat$time), substr, max(dl)+2, max(dl)+6)
         dat<- dat[match(tm,tmTmp),]
         datTmp<- data.frame(x=1:nrow(dat), y = dat[, j])
         rng<- c(min(pdat[,-1],na.rm=TRUE), max(pdat[,-1],na.rm=TRUE))

         idx1<- !is.na(datTmp$y)
         #   idx1[idx1]<- datTmp$y[idx1] > 0 # remove 0's?
         if(sum(idx1) < min.n){
            #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
            next
         }
         .dat<- datTmp
            .dat$y[!idx1]<- NA
         while(TRUE){
            if(sum(!is.na(.dat$y)) < min.n){
               #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
               break
            }else{
               dp1<- drop1f(.dat, min.n=min.n, span=span)
               if(is.na(dp1$drop)) break
               .dat<- dp1$data
            }
         }
         # add points that are close
         ytmp<- abs(datTmp$y-dp1$predict)
         idx<- !(ytmp > sd(.dat$y-dp1$predict, na.rm=TRUE)*qnorm(1-0.05/2) + sqrt(.Machine$double.eps))
            idx[is.na(idx)]<- FALSE
            idx<- idx & is.na(.dat$y)
         .dat$y[idx]<- datTmp$y[idx]
         rm(ytmp, idx)

         idx2<- !is.na(.dat$y)
         idx<- idx1 & idx2
         if(sum(idx) >= min.n){
            n<- sum(!is.na(.dat$y))
            # use "interpolate" to avoid unrealistic predicted values
            lf<- loess(y ~ x, data=.dat, na.action=na.exclude, degree=2, span=nf(n)/n,
               control=loess.control(surface="interpolate"))
            lfp<- predict(lf, .dat$x)
               lfp[lfp<rng[1]]<- rng[1]
               lfp[lfp>rng[2]]<- rng[2]
            dt[match(paste(ud[i], tm, sep="_"), rownames(dt)), j-1]<- lfp
            #cat("[",j-1,",",i,"]: ", dt[704:707,134], "\n",sep="")
         }else{
            #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
            next
         }
         #Sys.sleep(1)
      }
   }
   cat("\n")
   invisible(dt)
}

# day data & smoothing
dtf<- function(dayDat, day, ud, dl=10){
   # one data point from each day
   ddt<- matrix(NA,nrow=length(ud),ncol=ncol(dayDat))
   rownames(ddt)<- ud
   colnames(ddt)<- colnames(dayDat)
   for(i in 1:length(ud)){
      idx<- sapply(rownames(dayDat), substr, 1, dl) == ud[i]
      for(j in 1:ncol(ddt)){
         smr<- summary(dayDat[idx,colnames(ddt)[j]])
         # 'median'; better than a time point, e.g. 12:00 pm?
         ddt[i,j]<- smr["Median"]
      }
   }

   # smooth ddt
   adt<- ddt
   for(j in 1:ncol(adt)){
      dtTmp<- data.frame(x=1:nrow(ddt), y=ddt[,j])
      n<- sum(!is.na(dtTmp$y))
      if(n < 5){
         next
      }else{
         lf<- loess(y ~ x, data=dtTmp, na.action=na.exclude, degree=2, span=nf(n)/n,
              control=loess.control(surface="interpolate"))
      }
      rm(n)
      adt[,j]<- predict(lf, dtTmp$x)
   }
   adt[adt<0]<- NA

   # growth rate
   gr<- apply(adt, 2, diff)

   # relative growth rate
   rgr<- gr/adt[-nrow(adt),]

   list(ddt=ddt, adt=adt, gr=gr, rgr=rgr)
}

# quality control
qcf<- function(pdat, missing.pr=0.25, change.pr=0.25){
# missing.pr: max missing proportion
# change.pr: min change ratio
   # remove days with too many missing data
   cn<- colnames(pdat)
   pdat<- pdat[apply(is.na(pdat),1,mean)<missing.pr,]
   # remove plants with too many missing data
   ii<- apply(is.na(pdat),2,sum)
      ii<- (1:length(ii))[ii < nrow(pdat)/2]
   length(ii)
   # remove plants with little change or the control
   dmr<- apply(apply(pdat[,ii],2,range,na.rm=TRUE),2,diff)
   ii<- ii[dmr>IQR(apply(pdat[,ii],2,median,na.rm=TRUE))*change.pr]
   pdt<- pdat[,ii]
   colnames(pdt)<- cn[ii]

   pdt
}

# prepare for myScan
datProcf<- function(pdat){
   eid<- sapply(rownames(pdat),as.character)
   uid<- unique(eid)
   pdt.<- pdat[match(uid,eid),]
   nC<- rep(1,length(uid))
   for(id in uid){
      idx<- is.element(eid,id)
      if(sum(idx) < 2) next
      tmp<- pdat[idx,]
      #if(all(is.na(tmp))) cat(id, ' ')
      pdt.[match(id,uid),]<- apply(tmp,2,mean,na.rm=TRUE)
      nC[match(id,uid)]<- nrow(tmp)
   }
   pdt<- as.data.frame(pdt.)
   colnames(pdt)<- colnames(pdat)

   list(pdat=pdt, n=nC)
}

# genome scan
myScan<- function(pdat,gdat,gmap,day,cdl){
# pdat: phenotype data
# gdat: genotype data
# gmap: genetic/physical map
# day: which column of pdat
# cdl: experimental environment
####-------------------------------------------
   if(missing(cdl)){
      pdatTmp<- datProcf(pdat)
      pdt<- pdatTmp$pdat
      pdt$y<- pdt[,day]
      idx<- !is.na(pdt$y)
      pdt<- pdt[idx,]
      n<- pdatTmp$n[idx]

      rm(pdatTmp,idx)
   }else{
      uc<- unique(cdl)
      pdt<- n<- c<- NULL
      for(u in uc){
         pdatTmp<- datProcf(pdat[cdl==u,])
         pdtTmp<- pdatTmp$pdat
         pdtTmp$y<- pdtTmp[,day]
         idx<- !is.na(pdtTmp$y)
         pdt<- rbind(pdt, pdtTmp[idx,])
         n<- c(n, pdatTmp$n[idx])
         c<- c(c, rep(u,sum(idx)))

         rm(pdatTmp,pdtTmp,idx)
      }
      cdl<- c
      rm(uc,c,u)
   }
   idx<- match(rownames(pdt),rownames(gdat))
   pdt<- pdt[!is.na(idx),]
   n<- n[!is.na(idx)]
   if(!missing(cdl))
      cdl<- cdl[!is.na(idx)]
   idx<- idx[!is.na(idx)]
   gdt<- gdat[idx,]

   cat("Day ", day, " of ",nrow(pdt)," days | ", date(),"\n", sep="")

   vc<- vector("list",5)
   for(j in 1:5){
      idx<- match(rownames(gdt),rownames(gm[[j]]$AA))
      if(any(is.na(idx))) stop("something wrong...")
      v<- list(
         AA = gm[[j]]$AA[idx,idx],
         DD = NULL,
         AD = NULL,
         HH = NULL,
         MH = NULL,
         EE = diag(1/n)
      )
      if(missing(cdl)){
         vc[[j]]<- estVC(y=pdt$y,v=v)
      }else{
         vc[[j]]<- estVC(y=pdt$y,x=cdl,v=v)
      }
   }
   #### genome scan
   lrt<- est<- NULL
   for(j in 1:5){
      idx<- is.element(colnames(gdt),gmap$snp[gmap$chr==j])
      gdtTmp<- gdt[,idx]
         gdtTmp<- as.matrix(gdtTmp)
      if(missing(cdl)){
         sc<- scanOne(y=pdt$y, gdat=gdtTmp, vc=vc[[j]], intcovar=NULL, minorGenoFreq=0.05)
      }else{
         sc<- scanOne(y=pdt$y, x=cdl, gdat=gdtTmp, vc=vc[[j]], intcovar=NULL, minorGenoFreq=0.05)
      }
      lrt<- c(lrt,sc$p)
      estTmp<- matrix(unlist(sc$par),nrow=length(sc$par),byrow=TRUE)
            rownames(estTmp)<- names(sc$par)
            colnames(estTmp)<- names(sc$par[[1]])
      est<- rbind(est, estTmp)
   }

   list( lrt = lrt, est = est)
###--------------------------------
}

plotf<- function(lrt,gmap,main=""){
   m<- as.matrix(lrt[,-c(1:3)])/2/log(10)
   max(m,na.rm=TRUE) #6.70324
   cv<- qchisq(1-0.05/nrow(m),1)/2/log(10)
   br<- c(0, cv/4*1:5, Inf)
   cex<- seq(0,1,length=length(br))[-1]
      cex<- cex[as.numeric(cut(m,breaks=br))]
      cex<- matrix(cex,ncol=ncol(m))
   rbPal <- colorRampPalette(c('yellow','red'))
   col <- rbPal(length(br))[-1]
      col[5:6]<- "blue"
      col<- col[as.numeric(cut(m,breaks=br))]
      col<- matrix(col,ncol=ncol(m))

   par(mar=c(3.75,3.75,2.5,5),mgp=c(2.5,1,0))
   plot(0,0,xlim=c(1,nrow(m)),ylim=c(1,ncol(m)),type="n",xaxt="n",main=main,xlab="Chromosome",ylab="Day")
   ii<- c(FALSE,diff(gmap$chr)!=0)
      ii<- c(0,(1:length(ii))[ii]-1,nrow(gmap))
   h<- ii[-length(ii)]+diff(ii)/2
   axis(1,h,labels=1:5,tick=FALSE)
   axis(1,ii+0.5,labels=FALSE)
   idx<- m<br[2]
   m2<- m; m2[idx]<- NA
   cl2<- col; cl2[idx]<- NA
   cx2<- cex; cx2[idx]<- NA
   for(j in 1:ncol(m)){
      points(1:nrow(m), rep(j,nrow(m)), pch=20, cex=cx2[,j], col=cl2[,j])
   }
   idx<- m<br[5]
   m2<- m; m2[idx]<- NA
   cl2<- col; cl2[idx]<- NA
   cx2<- cex; cx2[idx]<- NA
   for(j in 1:ncol(m)){
      points(1:nrow(m), rep(j,nrow(m)), pch=20, cex=cx2[,j], col=cl2[,j])
   }
   cl<- rbPal(length(br))[-1]
      cl[5:6]<- "blue"
   cx<- seq(0,1,length=length(br))[-1]
   brt<- round(br[-length(br)],2)
   points(rep(nrow(m)*1.05,length(cl)),1:length(cl)*2,col=cl,cex=cx,pch=20,xpd=T)
   text(rep(nrow(m)*1.05,length(cl)),1:length(cl)*2,brt,pos=4,cex=0.75,font=c(1,1,1,1,4,1),xpd=TRUE)
   #text(nrow(m)*1.05,15,paste("0.05 significance threshold is ",round(br[-length(br)],1)[5],sep=""),pos=4,cex=0.75,xpd=T,srt=90)
}



##########################
# plotting data functions
#
library(jpeg)
if(file.exists("../kangaroo.jpeg"))
  img<- readJPEG("../kangaroo.jpeg") else img<- NA

plotDatf<- function(pdat,day,dayDat,ddt,cntNa,img,type=c("pdf","jpg"),outputDirImagesName){
  # outputDirImagesName = Pass the full filepath and root filename for the output images
  # Example input: E:/a_data/temp/GWAS/outputs/images/BVZ0039-GC03L-C01~fullres-area
  # Example: output: E:/a_data/temp/GWAS/outputs/images/BVZ0039-GC03L-C01~fullres-area-cover.jpg
  
  outputDirImage = "E:/a_data/temp/GWAS/outputs/images/BVZ0039-GC03L-C01~fullres-area-cover"
  type<- match.arg(type)
  if(type=="jpg")
    jpeg(paste(outputDirImage,"-Cover.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  # ---------------------------------------------
  plot.new()
  text(0.5, 0.80, "NOTE")
  text(0.5, 0.75,"************************************************************************")
  text(0.5, 0.65, paste("THIS GRAPHICALLY DISPLAYS DATA AND PREDICTION\n\n",
                        "FOR THE PURPOSE OF DIAGNOSIS",sep=""))
  if(!is.na("img")) rasterImage(img,0.35,0.30,0.65,0.50)
  text(0.5, 0.55,"************************************************************************")
  text(0.5, 0.20, format(Sys.time(), "%A, %B %d, %Y"))
  text(1, 0.10, "Warning: you are running the sofware without any warranty", cex=0.95, adj=1)
 # text(1, 0.15, "_________________", adj=1)
  text(1, -0.0, paste("Copyright \uA9 R. Cheng 2015",sep=""), font=3, cex=0.75, adj=1, xpd=TRUE)
  text(1, -0.03, paste("Research funded by an ARC Linkage Grant with the Borevitz Lab ",sep=""), font=3, cex=0.75, adj=1, xpd=TRUE)
  text(1, -0.06, paste("Research School of Biology, ANU Plant Sciences, Canberra Australia 2017.",sep=""), font=3, cex=0.75, adj=1, xpd=TRUE)
  text(1, -0.09, paste("Additional funding from CoE Plant Energy Biology & The Australian Plant Phenomics Facility",sep=""), font=3, cex=0.75, adj=1, xpd=TRUE)
  text(1, -0.12, paste("http://borevitzlab.anu.edu.au | http://traitcapture.org",sep=""), font=3, cex=0.75, adj=1, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
  
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-AllScatter.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  dTmp<- sort(unique(day), decreasing=FALSE)
  udTmp<- dTmp[(1:length(dTmp))%%3==1]
  tTmp<- rownames(dayDat)
  tmIdx<- match(tTmp, sapply(as.character(pdat$time),substr,1,16))
  rng<- range(pdat[tmIdx,-1], na.rm=TRUE)
  matplot(pdat[tmIdx,-1], ylim=rng, type="p", cex=0.25, xaxt="n",
          xlab="Time", ylab="Observed", main=phenoFileName)
  axis(1, at=match(udTmp,day[tmIdx]), tck=-0.010, labels=FALSE)
  text(match(udTmp,day[tmIdx]), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
  
  # raw data + smooth curve per day for each plant
  fmt<- (ncol(pdat)-2)%/%12+1
  fmt<- (fmt-1)%/%10+1
  fmt<- paste("%0",fmt,"d",sep="")
  if(type=="pdf") par(mfrow=c(4,3), mar=c(3.0,3.0,2.0,1), mgp=c(2,1,0))
  for(j in 2:ncol(pdat)){
    cat("Now on image", j, "of", ncol(pdat),"\n")
    if(type=="jpg"){
      jpeg(paste(dirname(outputDirImagesFilename),"/plotsByPlantID/",basename(outputDirImagesFilename),"-",colnames(pdat)[j],".jpg",sep=""),
           width = 720, height = 840, quality=100, res=100)
      par(mar=c(3.0,3.0,2.0,1), mgp=c(2,1,0))
      # ---------------------------------------------
    }
    plot(1:length(tTmp), rep(NA,length(tTmp)), ylim=rng, type="n", xaxt="n", cex=0.25,
         cex.axis=0.65, xlab="Time", ylab="Observed",
         main=paste("Plant ID ", colnames(pdat)[j], sep=""))
    axis(1, at=match(udTmp,day[tmIdx]), tck=-0.015, labels=FALSE)
    text(match(udTmp,day[tmIdx]), rep(min(rng)-diff(rng)*0.125,length(udTmp)),
         labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.5, xpd=TRUE)
    for(d in dTmp){
      dIdx<- grep(d,tTmp)
      tIdx<- match(tTmp[dIdx], sapply(as.character(pdat$time),substr,1,16))
      points(dIdx, pdat[tIdx,][,j], ylim=rng, type="p", cex=0.25)
      #yTmp<- dayDat[dIdx, match(colnames(pdat)[j], colnames(dayDat))]
      # the above has trouble in case duplicate plant IDs
      yTmp<- dayDat[dIdx, j-1]
      idx<- dIdx[!is.na(yTmp)]
      yTmp<- yTmp[!is.na(yTmp)]
      lines(idx, yTmp, col=2)
    }
    if(type=="jpg") dev.off()
    # ---------------------------------------------
  }
  
  udTmp<- ud[(1:length(ud))%%3==1]
  
  #1) predicted
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-Predict1.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  rng<- range(ddt$ddt, na.rm=TRUE)
  matplot(ddt$ddt, type="l", xaxt="n", xlab="Day", ylab="Predicted",
          main="Prediction for All Plants (Per Day)")
  axis(1, at=match(udTmp,rownames(ddt$ddt)), tck=-0.010, labels=FALSE)
  text(match(udTmp,rownames(ddt$ddt)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
  #2) smooth over days (one data point for each data)
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-Predict2.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  rng<- range(ddt$adt, na.rm=TRUE)
  matplot(ddt$adt, type="l", xaxt="n", xlab="Day", ylab="Predicted",
          main="Smoothed Prediction for All Plants")
  axis(1, at=match(udTmp,rownames(ddt$adt)), tck=-0.010, labels=FALSE)
  text(match(udTmp,rownames(ddt$adt)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
  #3) remove data points with too many missing data
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-Predict3.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  rng<- range(ddt$ddt, na.rm=TRUE)
  idx<- apply(is.na(ddt$ddt),2,sum)
  ii<- (1:length(idx))[idx < nrow(ddt$ddt)/4]
  #& remove plants with little change or the control
  dmr<- apply(apply(ddt$ddt[,ii],2,range,na.rm=T),2,diff)
  ii<- ii[dmr>median(dmr)/10]
  matplot(ddt$ddt[,ii],type="l", xaxt="n", xlab="Day", ylab="Predicted",
          main="Prediction for Plants with Missing Proportion < 25%")
  axis(1, at=match(udTmp,rownames(ddt$ddt)), tck=-0.010, labels=FALSE)
  text(match(udTmp,rownames(ddt$ddt)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
  #4) accumulative number of missing data points
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-CumulNA.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  rng<- range(cntNa, na.rm=TRUE)
  if(rng[2]==0) rng<- c(0,1)
  matplot(cntNa, ylim=rng, type="l", xaxt="n", xlab="Day", ylab="Cumulative Number of NAs",
          main="Missing Data Pettern for All Plants (Per Day)")
  axis(1, at=match(udTmp,rownames(cntNa)), tck=-0.010, labels=FALSE)
  text(match(udTmp,rownames(cntNa)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  tbl<- table(cntNa[nrow(cntNa),])
  text(rep(nrow(cntNa)-0.5,length(tbl)), as.numeric(names(tbl)), tbl, cex=0.75, col=2, pos=4, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
  #5) growth rate
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-GrowthRate.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  rng<- range(ddt$gr, na.rm=TRUE)
  matplot(ddt$gr, type="l", xaxt="n", xlab="Day", ylab="Growth Rate", main="")
  axis(1, at=match(udTmp,rownames(ddt$gr)), tck=-0.010, labels=FALSE)
  text(match(udTmp,rownames(ddt$gr)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  if(type=="jpg") dev.off()
  #6) relative growth rate
  if(type=="jpg")
    jpeg(paste(outputDirImagesName,"-GrowthRateRel.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
  par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
  # ---------------------------------------------
  rng<- range(ddt$rgr, na.rm=TRUE)
  matplot(ddt$rgr, type="l", xaxt="n", xlab="Day", ylab="Relative Growth Rate", main="")
  axis(1, at=match(udTmp,rownames(ddt$rgr)), tck=-0.010, labels=FALSE)
  text(match(udTmp,rownames(ddt$rgr)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
       labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
  # ---------------------------------------------
  if(type=="jpg") dev.off()
}

###############################################################################
# process data #
#############################################################


### IMPORTANT: Uncomment these 2 lines if these two packages aren't installed:
#  install.packages("QTLRel")  # Required for GWAS
# install.packages("jpeg")  # Required for output

library(QTLRel)
library(jpeg)

###

# SET YOUR LOCAL FILE PATH HERE
#   * Source data is expected in the "input-data" folder
#   * Output data will go to folder: rootPath/outputs
# If you have two input files for phenotype data, like L and R CSV data, uncomment the "phenoFile1" line and put your second file there
#   - The current code supperts two cameras (e.g. -area.csv for a Left and a Right camera, but if you have more cameras tyou can add additional lines
## IMPORTANT you must use "/" not "\" for filepaths, even on windows

rootPath <- "E:/a_data/projects/GWAS/BVZ0073"

########################

### Folder paths
setwd(rootPath)
dataDir<- paste(rootPath,"/input-data",sep="",collapse="/")
outputDir <-paste(rootPath,"/outputs/",sep="",collapse="/")
outputDirImages <-paste(outputDir,"images/",sep="",collapse="/")
outputDirGraphs <-paste(outputDir,"graphs/",sep="",collapse="/")


# create folders (if unavailable) for ouput
for(f in c(outputDir,outputDirGraphs,outputDirImages, paste(outputDirImages,"plotsByPlantID",sep=""))){
  if(!dir.exists(f)) dir.create(f)
  rm(f)
}

## Load genome data
# NOTE: make sure expInfo$EcotypeID are for the 473x250k genotype data

load("input-data/geno473x250k.RData") #Assumes genome file is in the root folder
gmap<- phyMap; colnames(gmap)<- c("snp","chr","dist")



# Paths to Phenotype -Area file(s)
# Files must have an equal number of rows. Assumption is that the plant-ID data is all in one plantInfo file (see below)
# Make sure dataset with pot# 1 is the phenoFile and right side of chamber is phenoFile1
phenoFile <- paste(dataDir,"/","BVZ0073-GC36L-RGB01~fullres-area.csv",sep="",collapse="/")
phenoFile1 <- paste(dataDir,"/","BVZ0073-GC36R-RGB01~fullres-area.csv",sep="",collapse="/")


# Set filenames/paths and prep data
phenoFileName<- basename(phenoFile) 
phenoFileName <- sub('\\..*$', '', basename(phenoFile)) # this is the filename without the suffix (used for writing output filenames later)
expID<- strsplit(phenoFileName,"-")[[1]][1] 
plantInfo<- paste(dataDir,"/", expID,"-traitcapture-db-import.csv",sep="",collapse="/")
expInfo<- read.csv(plantInfo, check.names=FALSE)
head(expInfo)

#Create root paths with filenames for data output (easier to do here then recreate them 100 times later)
outputDirFilename <- paste(outputDir,phenoFileName,sep="")
outputDirImagesFilename<-paste(outputDirImages,phenoFileName,sep="")
outputDirGraphsFilename<-paste(outputDirGraphs,phenoFileName,sep="")

trayNum<- strf(strsplit(as.character(expInfo$Tray),"[A-Z]"))
trayNum<- as.integer(trayNum)


# Load data
#
pdat1<- read.csv(phenoFile, na.string=c("Na", "NaN"), header=TRUE, check.names=FALSE, skip=0)

# If there are multiple phenotype csv's (e.g. Chamber right and Left, etc) then we import both of them and then concatenate the two files)
if (exists("phenoFile1")){
  pdat2<- read.csv(phenoFile1, na.string=c("Na", "NaN"), header=TRUE, check.names=FALSE, skip=0)
  
  # Make sure the column names are numbeed from 161 onwards rather then restarting at 1
  # Assumes the first data files has 160 columns and timestamps are the same for both -area files
  if (colnames(pdat2)[2] == 1) {
    #remove timestamp (assumption is that this is the same times as in pdat1)
    pdat2<-pdat2[,-1] # delete column 1
    totCol <-ncol(pdat2)+160
    colnames(pdat2) <- seq(161,ncol(pdat2)+160)
    }

  #concatenate the two data files
    pdat<-cbind(pdat1,pdat2)
  } else{
    pdat<-pdat1
}


##############################################
# Run initial analysis and data cleaning

pdat<- pdat[apply(!is.na(pdat),1,sum,na.rm=TRUE)>ncol(pdat)*2/3,] # remove NA rows
pdat<- pdat[apply(pdat[,-1] > 0,1,sum,na.rm=TRUE)>ncol(pdat)*(2/3),] # remove 0 rows 
pid<- sapply(expInfo$PlantID[match(colnames(pdat),expInfo$"Pot")],as.character)
pid[1]<- colnames(pdat)[1]
pid[is.na(pid)]<- paste("Pot",colnames(pdat)[is.na(pid)],sep="")
colnames(pdat)<- pid
rm(pid)


##  *** AFTER HERE THINGS TAKE A LONG TIME TO RUN ***

#  Plan for at least ~ 1min /plant on a big machine for the _procDayDataf_ function

# process data by fitting smoothing curves
#
dl<- 10
min.n<- 7
date()
dayDat<- procDayDataf(pdat, dl=dl, min.n=min.n, span=NA)
date()

#
# one data point / day
#
day<- sapply(as.character(pdat$time), substr, 1, 10)
ud<- sort(unique(day), decreasing=FALSE)
# ddt, adt, gr & rgr
ddt<- dtf(dayDat, day, ud)

cntNa<- apply(is.na(ddt$ddt), 2, cumsum)

#
# save results
#
# smoothed curve over days
write.csv(ddt$adt,file=paste(outputDirFilename,"-smoothed.csv",sep=""))
# growth rate
write.csv(ddt$gr,file=paste(outputDirFilename,"-growthRate.csv",sep=""))

# candidates for manual checking
rownames(ddt$adt)[1:min(nrow(ddt$adt),25)]
pdt<- qcf(ddt$adt[1:min(nrow(ddt$adt),25),], missing.pr=0.25, change.pr=1/3)
ex<- expInfo[is.element(expInfo$PlantID,colnames(ddt$adt)),]
exLst<- ex[!is.element(ex$PlantID,colnames(pdt)),]
rm(pdt,ex)

write.csv(exLst, file=paste(outputDirFilename,"_checkList.csv",sep=""), row.names=FALSE)

save.image(paste(outputDirFilename,"-workspace.RData",sep="")) #Save workspace



#
# Plot Data after initial data cleaning
# (This takes a few minutes)

plotDatf(pdat=pdat, day=day, dayDat=dayDat, ddt=ddt, cntNa=cntNa, img=img, type="jpg", outputDirImagesFilename)
if(FALSE){
  cvt<- paste("convert -compress Zip -quality 100", outputDirImagesFilename, "_Cover.jpg", sep="")
  cvt<- paste(cvt, outputDirImagesFilename, "_AllScatter*.jpg", sep="")
  cvt<- paste(cvt, outputDirImagesFilename, "_Idv*.jpg", sep="")
  cvt<- paste(cvt, outputDirImagesFilename, "_Prediction_*.jpg", sep="")
  cvt<- paste(cvt, outputDirImagesFilename, "_CumulNA.jpg", sep="")
  cvt<- paste(cvt, outputDirImagesFilename, "_GrowthRate*.jpg", sep="")
  cvt<- paste(cvt, " ", phenoFileName, "Tmp.pdf", sep="")
  system(cvt)
}

pdf(paste(outputDirGraphsFilename,".pdf",sep=""), title=phenoFileName, height=9)
plotDatf(pdat=pdat, day=day, dayDat=dayDat, ddt=ddt, cntNa=cntNa, img=img, type="pdf")
dev.off()

# Back up the pdat file so we have a copy since it gets modified in the GWAS and the initial plotting, etc won't work with the modified file
pdatBak = pdat


#q("no")


###################################################################
# GWAS #
###################################################################

# GWAS can be run on shifted or unshifted growth curves. 
# In "shifted", all growth cuvres start when the plant reaches a set pixel size.
# In "unshifted" analysis is run starting with the first timepoint for each plant

###################################
# GWAS without shifting
#

# scanning
pdat<- qcf(ddt$adt, missing.pr=0.25, change.pr=0.25)
   pdat<- t(pdat)
rownames(pdat)<- expInfo$EcotypeID[match(rownames(pdat),expInfo$PlantID)]
days<- colnames(pdat)
lrt<- matrix(NA, nrow=nrow(gmap), ncol=ncol(pdat))
   colnames(lrt)<- days
   lrt<- cbind(gmap, lrt)
   cat("Running genome scan, this will take a while")
   curTime<-Sys.time()
   timeElapsed <- 0
   
for(j in 1:ncol(pdat)){
  #time reporting
  curDiff<-(Sys.time()-curTime)
  timeElapsed = timeElapsed + curDiff
  cat(paste("Now on day", j, "of", ncol(pdat)), "\n Time to process last plant:",curDiff/60, "\n Time passed:", timeElapsed/60, " | Estimated time left:", ((ncol(pdat) - j)* curDiff),"min. \n\n")
  curTime <- Sys.time()
#
     sc<- myScan(pdat,gdat,gmap,j)
   lrt[match(names(sc$lrt),lrt$snp),days[j]]<- sc$lrt
}
   
rm(j,sc)
write.csv(lrt,file=paste(outputDirFilename,"_LRT.csv",sep=""),row.names=FALSE)

# plotting
jpeg(paste(outputDirGraphsFilename,".jpg",sep=""),quality=100,height=600,width=1600,res=100)
   plotf(lrt,gmap,main=phenoFileName)
dev.off()

# peak SNPs
cv<- qchisq(1-0.05/nrow(lrt),1)/2/log(10)
lrtTmp<- lrt[,-c(1:3)]/2/log(10)
   lrtTmp[is.na(lrtTmp)]<- 0
pks<- apply(lrtTmp>cv,1,any)
if(any(pks)){
   idx<- match(rownames(pdat),rownames(gdat))
      idx<- idx[!is.na(idx)]
   gd<- gdat[idx,pks]
   idx<- match(rownames(gdat)[idx],expInfo$EcotypeID)
   gdt<- cbind(expInfo[idx,1:3],gd)
      idx<- is.element(gdt$EcotypeID,rownames(pdat))
      gdt<- cbind(gdt[,1:3],excluded=!idx,gdt[,-c(1:3)])
      colnames(gdt)[-c(1:4)]<- colnames(gdat)[pks]
   ii<- apply(gdt[,-c(1:3)],1,paste,collapse="")
      ii<- order(ii)
   write.csv(gdt[ii,],file=paste(outputDirFilename,"_peakSNPs.csv",sep=""),row.names=FALSE)
}

#
# End GWAS without Time-Shifting
###################################


###################################
# GWAS with time shifting
#

# Reload original version of pdat
pdat = pdatbak

# normalization by shiting time
mmx<- apply(pdat,1,min,na.rm=TRUE)
   mmx<- quantile(mmx,0.99)
ii<- apply(pdat,1,which.min)
cat("Running normalization")
curTime<-Sys.time()
timeElapsed <- 0
for(i in 1:nrow(pdat)){
   tmp<- pdat[i,]
   idx<- (1:length(tmp))[tmp<mmx]
   if(length(idx)==length(tmp)){
      ii[i]<- NA
   }else if(length(idx)>0){
      if(any(!is.na(idx))){
         ii[i]<- max(ii[i],max(idx,na.rm=TRUE))
      }
   }
   rm(i,tmp,idx)
}
table(ii,useNA="ifany")
cv<- 13 #use the above information to make a choice 
pdat.s<- matrix(-Inf,nrow=sum(ii<=cv,na.rm=TRUE),ncol=ncol(pdat)-cv+1)
rownames(pdat.s)<- rownames(pdat)[ii<=cv & !is.na(ii)]
cnt<- 0
for(i in 1:nrow(pdat)){
   if(is.na(ii[i])) next
   if(ii[i]>cv) next
   cnt<- cnt+1
   # larger than mmx the next day
   idx<- 1:ncol(pdat.s)+(ii[i]-1)
   pdat.s[cnt,]<- pdat[i,idx]
   rm(i,idx)
}
rm(mmx,ii,cv,cnt)
sum(is.infinite(pdat.s)) # 0
dim(pdat.s)

# scanning with shifting
lrt.s<- matrix(NA, nrow=nrow(gmap), ncol=ncol(pdat.s))
   lrt.s<- cbind(gmap, lrt.s)
paste("Running Genome Scan")
curTime<-Sys.time()
timeElapsed <- 0
for(j in 1:ncol(pdat.s)){
  curInd<-j
  curDiff<-(Sys.time()-curTime)
  timeElapsed = timeElapsed + curDiff
  cat(paste("Now on day", curInd, "of", ncol(pdat)), "\n Time to process last plant:",curDiff, "min. \n Time elapsed:", timeElapsed/60, " | Estimated time left:", ((ncol(pdat) - j)* curDiff) * 60,"min. \n\n")
  curTime <- Sys.time()
   sc<- myScan(pdat.s,gdat,gmap,j)
   lrt.s[match(names(sc$lrt),lrt.s$snp),j+ncol(gmap)]<- sc$lrt
}
rm(j,sc)
write.csv(lrt.s, file=paste(outputDirFilename,"_shifted_LRT.csv",sep=""), row.names=FALSE)

# plotting
jpeg(paste(outputDirGraphsFilename,"_shifted.jpg",sep=""),quality=100,height=600,width=1600,res=100)
   plotf(lrt.s,gmap,main=phenoFileName)
dev.off()

# peak SNPs
cv<- qchisq(1-0.05/nrow(lrt.s),1)/2/log(10)
lrtTmp<- lrt.s[,-c(1:3)]/2/log(10)
   lrtTmp[is.na(lrtTmp)]<- 0
pks<- apply(lrtTmp>cv,1,any)
if(any(pks)){
   idx<- match(rownames(pdat.s),rownames(gdat))
      idx<- idx[!is.na(idx)]
   gd<- gdat[idx,pks]
   idx<- match(rownames(gdat)[idx],expInfo$EcotypeID)
   gdt<- cbind(expInfo[idx,1:3],gd)
      idx<- is.element(gdt$EcotypeID,rownames(pdat.s))
      gdt<- cbind(gdt[,1:3],excluded=!idx,gdt[,-c(1:3)])
      colnames(gdt)[-c(1:4)]<- colnames(gdat)[pks]
   ii<- apply(gdt[,-c(1:3)],1,paste,collapse="")
      ii<- order(ii)
   write.csv(gdt[ii,],file=paste(outputDirFilename,"_shifted_peakSNPs.csv",sep=""),row.names=FALSE)
}

###################################
# End GWAS with Time-Shifting
###################################

q("no")



#################################################
# the end #
###########