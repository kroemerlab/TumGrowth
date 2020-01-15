################################################################################################
### compModLGpw: compute LG in piecewise model
### compModLG: compute LG: linear model
### formatLGpw: format pairwise comparison table given the Reference - back transformation to be added soon
### prepDiagLG: prep data for plotDiag function
### .makeOnepw: compute contrasts given matrix K
################################################################################################

compModLGpw<-function(lgmat,bfco=0.1,tpcut=NA,adjust=TRUE){
  
  lgdf=lgmat$Df
  lgdf$out=FALSE
  lgdf=groupedData(Resp~Tp|Id/Grp,lgdf)
  ltps=sort(unique(lgdf$Tp))
  
  if(is.na(tpcut) | !(tpcut>min(ltps) & tpcut<max(ltps))){
    tp1<-function(x,bp) ifelse(x <= bp, x, bp)
    tp2<-function(x, bp) ifelse(x > bp, x - bp,0)
    
    optfct<-function(bp,lgdf){
      mod=lmer(Resp~tp1(Tp,bp)+tp2(Tp,bp)+Grp+tp1(Tp,bp):Grp+tp2(Tp,bp):Grp+(tp1(Tp,bp)|Id)+(0+tp2(Tp,bp)|Id), data=lgdf,subset=!out)
      AIC(mod)
    }
    search.range=range(ltps)+c(1,-1)*c(sum(diff(ltps)[1:2]),sum(rev(diff(ltps))[1:2]))/2
    search.grid <- sort(unique(subset(lgdf, Tp > search.range[1] &
                                        Tp<search.range[2], "Tp", drop=TRUE)))
    
    tpcut<-optimize(optfct,interval=range(search.grid),lgdf)$minimum
    
 #   res <- sapply(search.grid, optfct,lgdf)
  #  plot(search.grid, res, type="l")
 #   tpcut <- search.grid[which.min(res)]
  }
  
  
  lgdf$Tp1=lgdf$Tp
  lgdf$Tp1[lgdf$Tp1>tpcut]=max(lgdf$Tp1[lgdf$Tp1<=tpcut])
  lgdf$Tp2=lgdf$Tp-max(lgdf$Tp1)
  lgdf$Tp2[lgdf$Tp2<0]=0
  
  ### check if a mouse with at least  3 measurements in each segments
  if(min(max(tapply(lgdf$Tp2==0,lgdf$Id,sum)),max(tapply(lgdf$Tp2!=0,lgdf$Id,sum)))<3) 
    return(compModLG(lgmat,bfco=bfco,adjust=adjust))
  
  lgdf$out=FALSE
  lfms=list(as.formula(Resp~Tp1+Tp2+Grp+Tp1:Grp+Tp2:Grp+(1|Id)),
            as.formula(Resp~Tp1+Tp2+Grp+Tp1:Grp+Tp2:Grp+(1|Id)+(0+Tp1|Id)+(0+Tp2|Id)),
            as.formula(Resp~Tp1+Tp2+Grp+Tp1:Grp+Tp2:Grp+(Tp1|Id)+(0+Tp2|Id)),
            as.formula(Resp~Tp1+Tp2+Grp+Tp1:Grp+Tp2:Grp+(Tp1+Tp2|Id)))
  lfms2=list(as.formula(Resp~Tp1+Tp2+Grp+(1|Id)),
             as.formula(Resp~Tp1+Tp2+Grp+(1|Id)+(0+Tp1|Id)+(0+Tp2|Id)),
             as.formula(Resp~Tp1+Tp2+Grp+(Tp1|Id)+(0+Tp2|Id)),
             as.formula(Resp~Tp1+Tp2+Grp+(Tp1+Tp2|Id)))

  amods=lapply(lfms,function(fm) (try(lmer(fm, data=lgdf,subset=!out,REML=TRUE),TRUE)))
  lok=which(!sapply(amods,function(x) "try-error"%in%class(x)))
  minAIC=lok[which.min(sapply(amods[lok],AIC))]
  
  ga=lmerTest::lmer(lfms[[minAIC]], data=lgdf,subset=!out,REML=TRUE)
  
  iresid=residuals(ga,scale=TRUE)
  out=outlierTest(lm(iresid~1),cutoff = bfco)
  if(any(out$signif)){
    lgdf$out[which(rownames(lgdf)%in%names(out$bonf.p))]=T
    amods=lapply(lfms,function(fm) (try(lmer(fm, data=lgdf,subset=!out,REML=TRUE),TRUE)))
    lok=which(!sapply(amods,function(x) "try-error"%in%class(x)))
    minAIC=lok[which.min(sapply(amods[lok],AIC))]
    
    ga=lmerTest::lmer(lfms[[minAIC]], data=lgdf,subset=!out,REML=TRUE)
    
  }
  ###################################
  ## Compute global tests
  ga2=lmerTest::lmer(lfms2[[minAIC]], data=lgdf,subset=!out,REML=TRUE) ## avoid refit in anova
  gaml=refitML(ga)
  gaml1=update(gaml,.~.-Tp1:Grp)
  gaml2=update(gaml,.~.-Tp2:Grp)
  ga2ml=refitML(ga2)
  ga21=update(ga2ml,.~.-Grp)
  ga22=update(ga2ml,.~.-Tp1)
  ga23=update(ga2ml,.~.-Tp2)
  
  
  
  mf1=rbind(do.call("cbind",unclass(anova(ga, ddf = "Kenward-Roger")[4:5,c(5,3,4,6)])),
            do.call("cbind",unclass(anova(ga2, ddf = "Kenward-Roger")[1:3,c(5,3,4,6)])))
  mf2=rbind(do.call("cbind",unclass(anova(ga)[4:5,c(5,3,4,6)])),
            do.call("cbind",unclass(anova(ga2)[1:3,c(5,3,4,6)])))
  
  mlrt=rbind("Tp1:Grp"=unlist(unclass(anova(gaml,gaml1)[2,6:8])),
             "Tp2:Grp"=unlist(unclass(anova(gaml,gaml2)[2,6:8])),
             Tp1=unlist(unclass(anova(ga2ml,ga22)[2,6:8])),
             Tp2=unlist(unclass(anova(ga2ml,ga23)[2,6:8])),
             Grp=unlist(unclass(anova(ga2ml,ga21)[2,6:8])))
  
  mwald=rbind(car:::Anova.mer(ga,test="Chisq")[4:5,],car:::Anova.mer(ga2)[1:3,])
  
  antab=cbind(apply(mf1,1,function(x) sprintf("%.2f  (d.f.=%d/%.1f), p<%s",x[1],x[2],x[3],.myf(x[4]))),
              apply(mf2,1,function(x) sprintf("%.2f  (d.f.=%d/%.1f), p<%s",x[1],x[2],x[3],.myf(x[4]))),
              apply(mlrt,1,function(x) sprintf("%.2f  (d.f.=%d), p<%s",x[1],x[2],.myf(x[3]))),
              apply(mwald,1,function(x) sprintf("%.2f  (d.f.=%d), p<%s",x[1],x[2],.myf(x[3]))))
  dimnames(antab)=list(c("Time(1) x Treat","Time(2) x Treat","Time(1)","Time(2)","Treat" ),
                       c("F test (KR)","F test (S)","LR test","Wald test"))
  
  ###################################
  ## Pairwise comparison
  ct=contrMat(table(lgdf$Grp),type="Tukey")[,-1,drop=F]
  jt1=cbind(0,0,0,ct,matrix(0,nrow=nrow(ct),ncol=ncol(ct)))
  jt21=cbind(0,0,0,matrix(0,nrow=nrow(ct),ncol=ncol(ct)),ct,matrix(0,nrow=nrow(ct),ncol=ncol(ct)))
  jt22=cbind(0,0,0,matrix(0,nrow=nrow(ct),ncol=ncol(ct)),matrix(0,nrow=nrow(ct),ncol=ncol(ct)),ct)
  
  ## testing slopes only
  jtab1=.makeOnepw(ga,jt21)
  jtab1=cbind("Segment"=1,jtab1)
  jtab2=.makeOnepw(ga,jt22)
  jtab2=cbind("Segment"=2,jtab2)
  
  jtab=rbind(jtab1,jtab2)
  
  
  ##############################
  # Format output
  tab=summary(ga)$coef
  tab=data.frame(tab)
  tab[,5]=sapply(tab[,5],function(x) sprintf("%.5f",x))
  for(i in 1:4) tab[,i]=sapply(tab[,i],function(x) sprintf("%.3f",x))
  colnames(tab)=c("Coef","SE","DF","Tstat","Pvalue")
  
  ltps=ltp1s=floor(min(lgdf$Tp)-1):ceiling(max(lgdf$Tp)+1)
  ltp1s[ltp1s>max(lgdf$Tp1)]=max(lgdf$Tp1)
  ltp2s=ltps-ltp1s
  
  newdf=lapply(1:length(ltps),function(x) 
    do.call("cbind",LSmeans(ga,effect="Grp",at=list(Tp1=ltp1s[x],Tp2=ltp2s[x]),adjust=adjust)[1:2]))
  newdf=do.call("rbind",newdf)
  newdf=newdf[order(newdf$grid.Grp,newdf$grid.Tp1+newdf$grid.Tp2),]
  names(newdf)=c("fit","se.fit","df","t.stat","p.value","Grp","Tp1","Tp2")
  newdf$Tp=newdf$Tp1+newdf$Tp2
  newdf$Grp=factor(newdf$Grp,levels=levels(lgdf$Grp))
  
  return(list(model=ga,coefTab=tab,pairTab=jtab,AnovaTab=antab,weightCoef=NULL,corrCoef=NULL,
              data=lgdf,pred=newdf,Resp=lgmat$Resp,Trans=lgmat$Trans,Type='piecewise'))
}

################################################################################################
compModLG<-function(lgmat,bfco=0.1,adjust=TRUE){
  
  lgdf=lgmat$Df
  lgdf=groupedData(Resp~Tp|Id/Grp,lgdf)
  lgdf$out=FALSE
  
  # gmod=gls(Resp~Tp*Grp, data=lgdf[!lgdf$out,],correlation=corCAR1(form=~Tp|Id))
  # jt1=cbind(0,0,ct,matrix(0,nrow=nrow(ct),ncol=ncol(ct)))
  # jt2=cbind(0,0,matrix(0,nrow=nrow(ct),ncol=ncol(ct)),ct)
  
  lfms=list(as.formula(Resp~Tp+Grp+Tp:Grp+(1|Id)),
            as.formula(Resp~Tp+Grp+Grp:Tp+(1|Id)+(0+Tp|Id)),
            as.formula(Resp~Tp+Grp+Tp:Grp+(Tp|Id)))
  lfms2=list(as.formula(Resp~Tp+Grp+(1|Id)),
             as.formula(Resp~Tp+Grp+(1|Id)+(0+Tp|Id)),
             as.formula(Resp~Tp+Grp+(Tp|Id)))
  
  amods=sapply(lfms,function(fm) try(lmer(fm, data=lgdf,subset=!out,REML=TRUE),TRUE))
  lok=which(!sapply(amods,function(x) "try-error"%in%class(x)))
  minAIC=lok[which.min(sapply(amods[lok],AIC))]
  
  ga=lmerTest::lmer(lfms[[minAIC]], data=lgdf,subset=!out,REML=TRUE)
#  lme(Resp ~ Tp + Grp + Tp:Grp,random=~Id,correlation=corCAR1(form=~Tp), data=lgdf)
  iresid=residuals(ga,scale=T)
  out=outlierTest(lm(iresid~1),cutoff = bfco)
  if(any(out$signif)){
    lgdf$out[which(rownames(lgdf)%in%names(out$bonf.p))]=T
    amods=sapply(lfms,function(fm) try(lmer(fm, data=lgdf,subset=!out,REML=TRUE),TRUE))
    lok=which(!sapply(amods,function(x) "try-error"%in%class(x)))
    minAIC=lok[which.min(sapply(amods[lok],AIC))]
    
    ga=lmerTest::lmer(lfms[[minAIC]], data=lgdf,subset=!out,REML=TRUE)
    
  }
  
  ###################################
  ## Compute global tests
  ga2=lmerTest::lmer(lfms2[[minAIC]], data=lgdf,subset=!out,REML=TRUE) ## avoid refit in anova
  ga2ml=refitML(ga2)
  ga21=update(ga2ml,.~.-Grp)
  ga22=update(ga2ml,.~.-Tp)
  
  mlrt=rbind("Tp:Grp"=unlist(unclass(anova(refitML(ga),ga2ml)[2,6:8])),
             Tp=unlist(unclass(anova(ga2ml,ga22)[2,6:8])),
             Grp=unlist(unclass(anova(ga2ml,ga21)[2,6:8])))
  
  mf1=rbind(unlist(unclass(anova(ga, ddf = "Kenward-Roger")[3,c(5,3,4,6)])),
            matrix(unlist(unclass(anova(ga2, ddf = "Kenward-Roger")[,c(5,3,4,6)])),nrow=2))
  mf2=rbind(unlist(unclass(anova(ga)[3,c(5,3,4,6)])),
           matrix(unlist(unclass(anova(ga2)[,c(5,3,4,6)])),nrow=2))
  
  mwald=rbind(car:::Anova.mer(ga,test="Chisq")[3,],car:::Anova.mer(ga2)[1:2,])
  
  antab=cbind(apply(mf1,1,function(x) sprintf("%.2f  (d.f.=%d/%.1f), p<%s",x[1],x[2],x[3],.myf(x[4]))),
              apply(mf2,1,function(x) sprintf("%.2f  (d.f.=%d/%.1f), p<%s",x[1],x[2],x[3],.myf(x[4]))),
              apply(mlrt,1,function(x) sprintf("%.2f  (d.f.=%d), p<%s",x[1],x[2],.myf(x[3]))),
              apply(mwald,1,function(x) sprintf("%.2f  (d.f.=%d), p<%s",x[1],x[2],.myf(x[3]))))
  dimnames(antab)=list(c("Time x Treat","Time","Treat" ),c("F test (KR)","F test (S)","LR test","Wald test"))
  
  ###################################
  ## Pairwise comparison
  ct=contrMat(table(lgdf$Grp),type="Tukey")[,-1,drop=F]
  jt1=cbind(0,0,ct,matrix(0,nrow=nrow(ct),ncol=ncol(ct)))
  jt2=cbind(0,0,matrix(0,nrow=nrow(ct),ncol=ncol(ct)),ct)
  ## testing slopes only
  jtab=.makeOnepw(ga,jt2)
  
  
  ##############################
  # Format output
  tab=summary(ga)$coef
  tab=data.frame(tab)
  tab[,5]=sapply(tab[,5],function(x) sprintf("%.5f",x))
  for(i in 1:4) tab[,i]=sapply(tab[,i],function(x) sprintf("%.3f",x))
  colnames(tab)=c("Coef","SE","DF","Tstat","Pvalue")
  
  
  ltps=ltp1s=floor(min(lgdf$Tp)-1):ceiling(max(lgdf$Tp)+1)
  newdf=lapply(ltps,function(x) 
    do.call("cbind",LSmeans(ga,effect="Grp",at=list(Tp=x),adjust=adjust)[1:2]))
  newdf=do.call("rbind",newdf)
  newdf=newdf[order(newdf$grid.Grp,newdf$grid.Tp),]
  names(newdf)=c("fit","se.fit","df","t.stat","p.value","Grp","Tp")
  newdf$Grp=factor(newdf$Grp,levels=levels(lgdf$Grp))
  
  return(list(model=ga,coefTab=tab,pairTab=jtab,AnovaTab=antab,weightCoef=NULL,corrCoef=NULL,
              data=lgdf,pred=newdf,Resp=lgmat$Resp,Trans=lgmat$Trans,Type='linear'))
}


#####################################################################
.makeOnepw<-function(x,K){
  
  if(packageVersion("doBy") >= '4.6.0'){
    jtab2 = doBy::linest(x,L=K,adjust=TRUE)[[1]]
  }else{
    jtab2=doBy:::LSmeans(x,K=K,adjust=TRUE)[[1]]
  }
  selb = grep('se|std.error',colnames(jtab2),value=T)
  tlb =  grep('statistic|t.stat',colnames(jtab2),value=T)
  
  upc=jtab2[,"estimate"]+jtab2[,selb]*qt(0.975,df=jtab2[,"df"])
  loc=jtab2[,"estimate"]-jtab2[,selb]*qt(0.975,df=jtab2[,"df"])
  
  jtab=data.frame(do.call("rbind",strsplit(rownames(jtab2)," - ")),
                  cbind(jtab2[,c(tlb,selb,"df","p.value","estimate"),drop=FALSE],loc,upc),stringsAsFactors=F)
  names(jtab)=c("Largest","Smallest","Tstat","SE","Df","Pvalue","Contrast","Lower","Upper")
  for(i in which(jtab$Contrast<0)){
    jtab[i,c("Largest","Smallest")]=jtab[i,c("Smallest","Largest")]
    jtab[i,c("Lower","Upper")]=-jtab[i,c("Upper","Lower")]
    jtab[i,c("Tstat","Contrast")]=-jtab[i,c("Tstat","Contrast")]
  }
  rownames(jtab)=NULL
  return(jtab)
}

################################################################################################
formatLGpw<-function(lgres,ref='All',padjust="holm",backtrans=F){
  
  if(padjust=='no') padjust="none"
  
  ibtfct<-function(x) x;collab="Difference"
  # if(lgres$Trans=="Log" & backtrans){btfct<-function(x) exp(x);collab='Ratio'}
  #  if(lgres$Trans=="Sqrt" & backtrans) btfct<-function(x) x^2
  if(length(ref)==1) if(ref=='All') ref=levels(lgres$data$Grp)
  ref=ref[ref%in%levels(lgres$data$Grp)]
  if(length(ref)==0) ref=levels(lgres$data$Grp)
  pwt=lgres$pairTab
  if(length(ref)>0)  pwt=pwt[which(pwt[,"Largest"]%in%ref | pwt[,"Smallest"]%in%ref),,drop=F]
  if(length(ref)==1){
    if(any(pwt[,"Largest"]%in%ref)){
      for(i in which(pwt[,"Largest"]%in%ref)){
        pwt[i,c("Largest","Smallest")]=pwt[i,c("Smallest","Largest")]
        pwt[i,c("Lower","Upper")]=-pwt[i,c("Upper","Lower")]
        pwt[i,c("Tstat","Contrast")]=-pwt[i,c("Tstat","Contrast")]
      }
    }
  }
  
  dcts=apply(as.matrix(pwt[,c("Contrast","Lower","Upper")]),1,function(x)
    sprintf("%.3f [%.3f;%.3f]",ibtfct(x[1]),ibtfct(x[2]),ibtfct(x[3])))
  
  top=data.frame(pwt[,c("Largest","Smallest"),drop=F],Contrast=dcts,Df=sprintf('%.2f',pwt[,"Df"]),
                 Pvalue=.myf(pwt$Pval),PvalueAdj=.myf(p.adjust(pwt$Pval,padjust)),stringsAsFactors=FALSE)
  if("Segment"%in%colnames(pwt)){
    lsegs=tapply(lgres$pred$Tp,(lgres$pred$Tp2>0)+1,function(x) paste(range(x),collapse="-"))[c("1","2")]
#    print(str(lsegs))
    tsegs=as.vector(lsegs[as.character(pwt$Segment)])
    top=data.frame(TimeSeg= tsegs,pwt[,c("Largest","Smallest"),drop=F],Contrast=dcts,Df=sprintf('%.2f',pwt[,"Df"]),
                   Pvalue=.myf(pwt$Pval),PvalueAdj=.myf(p.adjust(pwt$Pval,padjust)),stringsAsFactors=FALSE)
  }
  
   if(length(ref)==1){
    top=top[,names(top)!="Smallest"]
    colnames(top)[colnames(top)=="Largest"]=paste("Comp to ",ref)
  }
  rownames(top)=NULL
 return(top)
#  return(list(tab=top,grps=levels(lgres$data$Grp)))
  
}


################################################################################################
prepDiagLG<-function(lgres){
  
  gcols=tapply(lgres$data$color,lgres$data$Grp,unique)
  resp=lgres$Resp
  ibtfct<-btfct(lgres$Trans,TRUE)
  # function(x) x
  # if(lgres$Trans=="Log"){btfct<-function(x) exp(x)}
  # if(lgres$Trans=="SqRt") btfct<-function(x) x^2
  # if(lgres$Trans=="CuRt") btfct<-function(x) x^3
  # 
  f<-function(x) {min(which( x*10^(0:20)==floor(x*10^(0:20)) ))} 
  
  ###############
  ## format data
  tmpdata=lgres$data[!lgres$data$out,c("Resp","Resp_ori","Grp","Tp","Id")]
  #  names(tmpdata)=c("y","Grp","x","Id")
  ndigy=ceiling(median(sapply(tmpdata$Resp_ori,f)))+1
  ndigye=ceiling(median(sapply(tmpdata$Resp,f)))+1
  tmpdata$Id=factor(tmpdata$Id,levels=unique(tmpdata$Id))
  tmpdata$color=gcols[tmpdata$Grp]
  tmpdata$Resp=round(tmpdata$Resp,ndigye)
  tmpdata$Fit=round(fitted(lgres$model,type="response"),ndigye)
  tmpdata$Resp_ori=round(tmpdata$Resp_ori,ndigy)
  tmpdata$Fite=round(ibtfct(fitted(lgres$model,type="response")),ndigy)
  tmpdata$resid=round(unclass(residuals(lgres$model,scale=T)),ndigye+2)
  lso=order(order(tmpdata$resid))
  tmpdata$qt=round(qnorm(ppoints(length(lso)))[lso],3)
  
  ###############
  ## format predicted data
  tmppred=lgres$pred[,c("Tp","Grp","fit","se.fit","df")]
  names(tmppred)[1]="x"
  tmppred$x=as.numeric(tmppred$x)
  tmppred$color=gcols[tmppred$Grp]
  tmppred$y=round(tmppred$fit,ndigye)
  tmppred$ymin=round(tmppred$fit-qt(.975,tmppred$df)*tmppred$se.fit,ndigye)
  tmppred$ymax=round(tmppred$fit+qt(.975,tmppred$df)*tmppred$se.fit,ndigye)
  tmppred$ybt=round(ibtfct(tmppred$fit),ndigy)
  tmppred$ybtmin=round(ibtfct(tmppred$fit-qt(.975,tmppred$df)*tmppred$se.fit),ndigy)
  tmppred$ybtmax=round(ibtfct(tmppred$fit+qt(.975,tmppred$df)*tmppred$se.fit),ndigy)
  tmppred$Grp=factor(tmppred$Grp,levels=levels(tmpdata$Grp))
  
  ###################################
  ## acf: autocorrelation to be added
  # iresid=residuals(lgres$model,scale=T)
  # acfres=do.call("cbind",unclass(acf(iresid,plot=FALSE))[c("lag","acf")])
  
  ###############
  ## format axes
  limtp=pretty(tmpdata$Tp)
  limyfit=pretty(c(tmppred$ybtmin,tmppred$ybtmax,tmpdata$Resp_ori,tmpdata$Fite)) ## original response/backT data
  limyfit0=pretty(c(tmppred$ymin,tmppred$ymax,tmpdata$Resp,tmpdata$Fit)) ## model fit
  limxqq=range(c(tmpdata$resid,tmpdata$qt))
  limxqq=floor(limxqq[1]):ceiling(limxqq[2])
  
  class(tmpdata)=class(tmppred)= "data.frame"
  return(list(data=tmpdata,pred=tmppred,limtp=limtp,limyfit=limyfit,limyfit0=limyfit0,limxqq=limxqq,Resp=lgres$Resp,Trans=lgres$Trans,Type=lgres$Type))
}

# ##################################################################################
# ### Additional function for hockey stick-> yet to be done
# .hockey <- function(x,alpha1,beta1,beta2,brk,eps=diff(range(x))/100)
#   ## alpha1 is the intercept of the left line segment
#   ## beta1 is the slope of the left line segment
#   ## beta2 is the slope of the right line segment
#   ## brk is location of the break point
#   ## 2*eps is the length of the connecting quadratic piece
#   
#   ## reference: Bacon & Watts "Estimating the Transition Between
#   ## Two Intersecting Straight Lines", Biometrika, 1971
#   
#   ## Original function coded by Mary Lindstrom
#   ## <lindstro@biostat.wisc.edu> and taken from
# ## S-NEWS Archive (Mon, 24 Apr 2000) available
# ## from (http://lib.stat.cmu.edu/s-news/Burst/15642). 
# {
#   x1 <- brk-eps
#   x2 <- brk+eps
#   b <- (x2*beta1-x1*beta2)/(2*eps)
#   cc <- (beta2-b)/(2*x2)
#   a <- alpha1+beta1*x1-b*x1-cc*x1^2
#   alpha2 <- (a + b*x2 + cc*x2^2) - beta2*x2
#   
#   lebrk <- (x <= x1)
#   gebrk <- (x >= x2)
#   eqbrk <- (x > x1 & x < x2)
#   
#   result <- rep(0,length(x))
#   result[lebrk] <- alpha1 + beta1*x[lebrk]
#   result[eqbrk] <- a + b*x[eqbrk] + cc*x[eqbrk]^2
#   result[gebrk] <- alpha2 + beta2*x[gebrk]
#   result
# }
# 
# #--------------------------------------------------------------------------------------------required for mergeScans
# getInitPars <- function(x,y)
# {
#   lowLeg <- (x <= quantile(x,probs=0.5))
#   init.lsfit <- lm(y~x, subset=lowLeg)$coef
#   b0 <- init.lsfit[1]; b1 <- init.lsfit[2];
#   init.lsfit2 <- lm(y~x, subset=!lowLeg)$coef
#   b2 <- init.lsfit2[2];
#   b3=-(init.lsfit[1]-init.lsfit2[1])/(init.lsfit[2]-init.lsfit2[2])
#   return(c(b0, b1, b2, b3))
# }
# 
# findThreshold <- function(x,y)
# {
#   n <- length(x)
#   y=y[order(x)]
#   x=x[order(x)]
#   pars <- getInitPars(x=x,y=y)
#   ###########################################           
#   #           if(pars[1]>10)pars[1]<-7
#   #           if(pars[4]>200)pars[4]<-190
#   #           print(pars);
#   ######################################
#   alpha1 <- pars[1]; beta1 <- pars[2]; beta2 <- pars[3];
#   brk <- pars[4];
#   
#   fit1 <- nls(y ~ .hockey(x=x,alpha1,beta1,beta2,brk),
#               start=list(alpha1=alpha1,beta1=beta1,beta2=beta2,brk=brk),
#               control=nls.control(maxiter=100))
#   resid.scale <- sqrt(fit1$m$deviance()/(n-4))
#   yhat <- fit1$m$fitted()
#   Coeff <- fit1$m$getPars()
#   # cat("\nCoefficients are: ",Coeff,"\n")
#   brk <- Coeff[4]; beta2 <- Coeff[3]; beta1 <- Coeff[2]; alpha1 <- Coeff[1];
#   Saturated <- y >= -(x - alpha1*beta1 - brk*(beta1*beta1 + 1))/beta1
#   # Saturated <- Saturated | (x >= brk)
#   high.row <- c(1:length(x))[Saturated]
#   High.fit <- alpha1 + beta1*x[high.row]
#   return(list(x=x,y=y,yhat=yhat,high.row=high.row,resid.scale = resid.scale,
#               High.fit=High.fit))
# }

