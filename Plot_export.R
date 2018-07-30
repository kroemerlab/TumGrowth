transsvg<-function(inam){
  doc=scan(inam,what="raw",sep="\n")
  for(k in grep(">Gst</text",doc)) doc[[k]]="<g>"
  for(k in grep(">Gend</text",doc)) doc[[k]]="</g>"
  write(doc,file=inam,sep="\n")
}

##################################################################################################
exportTC<-function(h1,dogrps=TRUE,type='tc'){
  
  if(length(type)==0) type='tc'
  se=('mese'%in%type)
  
  xtk=h1$xlim
  ytk=h1$ylim
  
  par(xpd=F)
  re=plot(range(xtk),range(ytk),pch=16,cex=0,axes=F,ylim=range(ytk),xlim=range(xtk),xlab="Time",ylab=h1$Title)
  idf=h1$df # Now initialized outside if statement
  
  if(dogrps)  graphics:::text(mean(xtk)*.99,mean(ytk[-1])*.98,"Gst")
  axis(1,at=xtk,pos=ytk[1]-range(ytk)*.02,lwd=2,tick=F)
  axis(1,at=xtk,pos=ytk[1]-range(ytk)*.02,lwd=2,labels = rep("",length(xtk)))
  if(dogrps)    graphics:::text(mean(xtk)*.99,mean(ytk[-1])*.98,"Gend")
  if(dogrps)    graphics:::text(mean(xtk)*.9,mean(ytk[-1])*.8,"Gst")
  axis(2,at=ytk,pos=xtk[1],las=2,lwd=2,tick=F)
  axis(2,at=ytk,pos=xtk[1],lwd=2,labels = rep("",length(ytk)))
  if(dogrps)    graphics:::text(mean(xtk)*.92,mean(ytk[-1])*.92,"Gend")
  
  
  if(any(c('mese','mesd')%in%type)){
    tmppred=h1$sesd
    tmppred$ymin=round(tmppred$y-ifelse(se,1,0)*tmppred$se-ifelse(se,0,1)*tmppred$sd,3)
    tmppred$ymax=round(tmppred$y+ifelse(se,1,0)*tmppred$se+ifelse(se,0,1)*tmppred$sd,3)
    if(any(is.na(tmppred$se))) tmppred$ymin[is.na(tmppred$se)]=tmppred$ymax[is.na(tmppred$se)]=tmppred$y[is.na(tmppred$se)]
    re=tapply(1:nrow(tmppred),tmppred$Grp,function(x){
      x=x[order(tmppred$x[x])]
      xl=tmppred$x[x];yl=tmppred$ymin[x];yh=tmppred$ymax[x]
      if(dogrps)  graphics:::text(mean(xl)*.95,mean(yh)*.78,"Gst")
      lines(xl,tmppred$y[x],col=tmppred$color[x][1],lwd=par("lwd")*1.5,lty=2)
      lines(xl,tmppred$ymin[x],col=tmppred$color[x][1],lwd=par("lwd")*1.5,lty=3)
      lines(xl,tmppred$ymax[x],col=tmppred$color[x][1],lwd=par("lwd")*1.5,lty=3)
      polygon(c(xl,rev(xl),xl[1]),c(yl,rev(yh),yl[1]),border=NA,
              col=adjustcolor(tmppred$color[x][1],alpha.f=0.3))
      if(dogrps)  graphics:::text(mean(xl)*.95,mean(yh)*.58,"Gend")
    })
  }
  
  if(any(c('tc')%in%type)){
  for(igrp in unique(idf$Grp)){
    if(dogrps)  graphics:::text(mean(xtk)*.95,mean(ytk[-1])*.78,"Gst")
    cols=idf$color[idf$Grp==igrp][1]
    for(imid in unique(idf$Id[idf$Grp==igrp])){
      x=idf$x[idf$Id==imid]
      y=idf$y[idf$Id==imid]
      id=idf$Id[idf$Id==imid]
      lines(x,y,col=cols,lwd=2)
    }
    if(dogrps)  graphics:::text(mean(xtk)*.942,mean(ytk[-1])*.912,"Gend")
  }
  }
  legend(mean(xtk[1:2]),max(ytk),unique(idf$Grp),ncol=length(unique(idf$Grp)),
         col=tapply(idf$color,idf$Grp,unique),bty="n",lwd=par()$lwd)
  par(xpd=T)
  
  
}

##################################################################################################
### Export CS boxplot - dogrps=TRUE, will add tags for each group of points, useful when export to SVG
exportCS<-function(p1,dogrps=TRUE,cexpt=1.2,strict=TRUE){
  
  idf=p1$m
  ylim=p1$ylim
  yaxt=p1$yaxt
  ylab=p1$title
  if(!strict) ylim=range(ylim,yaxt)
  rm(list="p1")
  xlim=c(0.4,nlevels(idf$Grp)+.6)
  re=boxplot(idf$Resp~idf$Grp,pch=16,cex=0,axes=F,lwd=par('lwd'),lty=1,ylab=ylab,
             border=tapply(idf$color,idf$Grp,unique)[levels(idf$Grp)],medlwd=par('lwd'),
             ylim=ylim,xlim=xlim)
  
  segments(xlim[1],min(ylim),xlim[2],min(ylim))
  re2=beeswarm(idf$Resp~idf$Grp,pch=16,at=1:nlevels(idf$Grp),pwcex = rep(cexpt,nrow(idf)),
               pwcol=idf$color,spacing=cexpt*1.01,
               corral = 'wrap',corralWidth = .5,add = T,do.plot=F)
  for(i in levels(idf$Grp)){
    ire=re2[!is.na(re2[,1]) & re2$x.ori==i,]
    if(dogrps)     graphics:::text(mean(ire[,1]),mean(ire[,2]),"Gst")
    points(ire[,1],ire[,2],cex=cexpt,pch=ire$pch,col=ire$col)
    if(dogrps)    graphics:::text(mean(ire[,1]),mean(ire[,2]),"Gend")
  }
  if(dogrps)   graphics:::text(nlevels(idf$Grp)*.52,mean(ylim[-1])*.8,"Gst")
  for(i in 1:nlevels(idf$Grp))
    axis(1,at=i,labels = paste(re$names,"\n(",re$n,")",sep="")[i],tick=FALSE,pos=ylim[1])
  if(dogrps)   graphics:::text(nlevels(idf$Grp)*.52,mean(ylim[-1])*.92,"Gend")
  
  if(dogrps)   graphics:::text(nlevels(idf$Grp)*.49,mean(ylim[-1])*.8,"Gst")
  axis(2,at=yaxt,pos=xlim[1],las=2,lwd=2)
  if(dogrps)   graphics:::text(nlevels(idf$Grp)*.42,mean(ylim[-1])*.92,"Gend")
  
}

##################################################################################################
exportKM<-function(p1,lwd=2,cexpt=1.5,dogrps=TRUE){
  
  akms=p1$m
  plot(range(p1$limtp),c(0,100),cex=0,
       xlab="Time",xlim=range(p1$limtp),
       ylab=p1$title,ylim=c(0,100),axes=F)
  if(dogrps)   graphics:::text(mean(p1$limtp),52,"Gst")
  axis(1,at=p1$limtp)
  if(dogrps)   graphics:::text(mean(p1$limtp),45,"Gend")
  
  if(dogrps)   graphics:::text(mean(p1$limtp),50,"Gst")
  axis(2,at=pretty(0:100),las=2)
  if(dogrps)   graphics:::text(mean(p1$limtp),48,"Gend")
  
  ticks=unique(unlist(lapply(akms,function(x) x$x0)))
  ticks=ticks[ticks>0]
  if(dogrps)   graphics:::text(mean(p1$limtp),47,"Gst")
  axis(1,at=ticks,lab=rep("",length(ticks)),pos=0,col.ticks="grey",lwd.ticks=par('lwd')*1.5,lwd=0)
  if(dogrps)   graphics:::text(mean(p1$limtp),55,"Gend")
  for(x in akms) {
    if(dogrps)   graphics:::text(mean(x$x),mean(x$y),"Gst")
    lines(x$x,x$y,col=x$color[1],lwd=lwd)
    if(!is.na(cexpt) & cexpt>0 & any(x$x0>0))
      points(x$x[x$x0>0],x$y[x$x0>0],col=x$color[1],pch=16,cex=cexpt)
    if(dogrps)   graphics:::text(mean(x$x)*.9,mean(x$y)*.8,"Gend")
  }
  
  legend('topright',sapply(akms,function(x) x$Grp[1]),xjust=0,
         ncol=1,bty="n",col=sapply(akms,function(x) x$color[1]),lwd =lwd*1.5,seg.len =1)
  
  
}
