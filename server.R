require(rCharts)
require(shiny)
library(shinyFiles)

source("helpers.R") ## common fcts

source("File_fct.R") ## file I/O
source("Plot_diag.R") ## generic diagnostic plto fct
source("Line_chart.R") ## line chart panel
source("Long_lmer.R")  ## longitudinal panel
source("Surv_km.R") ## survival panel
source("Cross_sect.R") ## cross-sectional panel
source("Plot_export.R") ## reporting functions


options(shiny.maxRequestSize = 200*1024^2)

shinyServer(function(input, output, session) {
  

  mypalette<-function(x) c("#61B270", "#397FB9","#D05043" ,"#8475B5", "#878787", "#EB7B36", "#A5691B" ,"#1B7E76",
                           '#C994C7','#E151A3','#D61967','#980043','#FC9272','#F75A40','#D72322','#A50F15',
                           '#9EBCDA','#8C87BF','#894FA3','#810F7C','#9ECAE1','#5DA4D0','#2C7CBA','#08519C',
                           '#A1D99B','#63BB6D','#2D954D','#006D2C','#BDBDBD','#8A8A8A','#5D5D5D','#252525',
                           '#FEC44F','#F88B22','#D65808','#993404','#C79141','#A5691B','#DFC27D','#7E4808','#543005')
  
  dataParams <- reactiveValues()
  dataParams$lastT <- c()
  
  values <- reactiveValues(
    file1 = NULL
  )
  
  observe({
    input$clearFile1
    input$dataInput
    values$file1 <- NULL
  })
  
  observe({
    values$file1 <- input$browse
  })
  ######## format the side bar
  lcols <- reactiveValues(data = NULL)
  fileData <- reactive({
    if(input$dataInput==1) inFile<-list(Ori=input$sampleData,Where=paste("./",input$sampleData,sep=""))
    if(input$dataInput==2){
  #    print(input$browse)
      shiny:::validate(
        need(try(!is.null(values$file1)), "Please select a data set")
      )
      inFile <- list(Ori=values$file1$name,Where=values$file1$datapath)
    }
    inFile$Ori=gsub(".*/","",gsub("\\.[A-Za-z]+$","",inFile$Ori))
    if(nchar(inFile$Ori)<2) inFile$Ori=paste("TG_",Sys.Date(),format(Sys.time(), "%H:%M"),sep="")
      return(inFile)
  })
  
  dat <- reactive({
    ifile=fileData()
    if (is.null(ifile$Where)) return(NULL)
      cat("Loading",ifile$Where,"\n")
#      print(ifile)
      shiny:::flushReact()
      lcols$data<-NULL
      setday0=suppressWarnings(as.numeric(input$setday0))
      res<-try(suppressWarnings(loadFile(ifile$Where,ndigit=as.numeric(input$ndigits),imputezer=input$imputezer,trim=input$trim, setday0=setday0,
               trimzer=input$trimzer,exclzer=input$exclzer,set2surv=input$set2surv)),TRUE)
      if(class(res)%in%"try-error")
        res=paste("Oups, something went wrong with the format of",ifile$Ori)
      res
  })
  
  mresp <- reactive({cdat=dat();if(!is.list(cdat)) return(NULL);cdat$Resp})
  mresptr <- reactive({cdat=dat();if(!is.list(cdat)) return(NULL);cdat$RespTr})
  mgrps <- reactive({cdat=dat();if(!is.list(cdat)) return(NULL);levels(cdat$dataM$Grp)})
  l2Excl <- reactive({cdat=dat();if(!is.list(cdat)) return(NULL);cdat$Excl})
  censvalt<-reactive({cdat=dat();if(!is.list(cdat)) return(NULL);sort(unique(cdat$data$Tp))})
  
  ## add censoring on time for TFS
  
  censvalrespos<-reactive({
    if(is.null(input$responsekm)) return(NULL)
    cdat=dat();if(!is.list(cdat)) return(NULL)
    if(!input$responsekm %in%cdat$Resp) return(NULL)
    vraw=na.omit(cdat$data[,input$responsekm])
    vraw=round(c(quantile(vraw,.05),max(vraw)),as.numeric(input$ndigits))
    vraw=pretty(vraw,100)
    unique(round(vraw,as.numeric(input$ndigits)))
    })
  
  censvalresptfs<-reactive({
    if(is.null(input$responsekm)) return(NULL)
    cdat=dat();if(!is.list(cdat)) return(NULL)
    if(!input$responsekm %in%cdat$Resp) return(NULL)
    vraw=na.omit(cdat$data[,input$responsekm])
    vraw=round(c(quantile(vraw,.95),min(vraw)),as.numeric(input$ndigits))
    vraw=pretty(vraw,100)
    print(vraw)
    unique(round(vraw,as.numeric(input$ndigits)))
  })
  

  refkm <- eventReactive(input$varskm,input$varskm)
#  reflg <- eventReactive(eventExpr=modelLGpw(),valueExpr=modelLGpw()$grps)
  refcs <- eventReactive(input$varscs,input$varscs)
  
  ###### Update colors
  observeEvent(input$updateCols,{
    num2<-mgrps()
    newcols=sapply(1:length(num2),function(x) eval(parse(text=paste("input$colCustom",x,sep=""))))
    names(newcols)=num2
    lcols$data<-newcols
  })
  
  colGrps<-reactive({
    if(!is.null(lcols$data)) return(lcols$data)
    num2<-mgrps()
    newcols=mypalette()[1:length(num2)]
    names(newcols)=num2
    newcols
  })
  ##################################
  ###### Format ui components
  output$groupscols<-renderUI({
    nums2<-mgrps()
    if(is.null(nums2)) return(NULL)
    listcols0=mypalette()
    if(!is.null(lcols$data)){
      listcols1=lcols$data
      if(length(listcols1)==length(nums2)) if(all(names(listcols1)==nums2)) listcols0=listcols1
    }
    grpcols=list()
    for(i in 1:length(nums2))
      grpcols[[i]]=colourInput(paste("colCustom",i,sep=""), nums2[i], palette = "limited",value=listcols0[i],
                               allowedCols=unique(mypalette()))
    grpcols
  })
  
  ####### data upload
  output$choicesdatexp<-renderUI({
    nums2=mresp()
    if(is.null(nums2)) return(NULL)
    if(length(nums2)==1) return(NULL)
    radioButtons("responsedatexp",NULL, choices=c('All',nums2),selected='All',inline = T)
  })
  output$choicesdat<-renderUI({
    nums2=mresp()
    if(is.null(nums2)) return(NULL)
    radioButtons("responsedat",NULL, choices=nums2,selected=nums2[1],inline = T)
  })
  
  output$groupsdat<-renderUI({
    nums2<-mgrps()
    if(is.null(nums2)) return(NULL)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("groupsdat",NULL, choices=list,selected=list,inline = T)
  })

  ####### line charts
  output$choices<-renderUI({
    nums2=mresp()
    if(is.null(nums2)) return(NULL)
    radioButtons("response", 'Response', choices=nums2,selected=nums2[1],inline = T)
  })
  
  output$choicestrans<-renderUI({
    trtypes=('None')
    if(is.null(input$response)) return(NULL)
    nums2=mresptr()
    if(paste(input$response,".log",sep="") %in% nums2)  trtypes=c(trtypes,'Log')
    if(paste(input$response,".sqrt",sep="") %in% nums2)  trtypes=c(trtypes,'SqRt')
    if(paste(input$response,".curt",sep="") %in% nums2)  trtypes=c(trtypes,'CuRt')
    radioButtons("trans", "Transformation", choices=trtypes,selected='None',inline = TRUE)
  })
  
  output$boxes<-renderUI({
    nums2<-mgrps()
    if(is.null(nums2)) return(NULL)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("vars",'Treatment groups', choices=list,selected=list,inline = T)
  })
  
  ####### longitudinal
  output$choiceslg<-renderUI({
    nums2=mresp()
    if(is.null(nums2)) return(NULL)
    radioButtons("responselg", "Response", choices=nums2,selected=nums2[1],inline = TRUE)
  })
  
  output$choiceslgtrans<-renderUI({
    trtypes=('None')
    if(is.null(input$responselg)) return(NULL)
    nums2=mresptr()
    if(paste(input$responselg,".log",sep="") %in% nums2)  trtypes=c(trtypes,'Log')
    if(paste(input$responselg,".sqrt",sep="") %in% nums2)  trtypes=c(trtypes,'SqRt')
    if(paste(input$responselg,".curt",sep="") %in% nums2)  trtypes=c(trtypes,'CuRt')
    radioButtons("translg", "Transformation", choices=trtypes,selected='None',inline = TRUE)
  })
  
  output$boxeslg<-renderUI({
    nums2<-mgrps()
    if(is.null(nums2)) return(NULL)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("varslg", "Treatment groups", choices=list,selected=list,inline = TRUE)
  })
  
  output$grplgvar<-renderUI({
    pwt=modelLG()
    if(is.null(pwt)) return(NULL)
    nums2=levels(pwt$data$Grp)
    if(length(nums2)>2) nums2=c('All',nums2)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("radiolg", label =NULL,inline = TRUE,choices =unname(nums2), selected = nums2[1])
  })
  
  ####### survival
  output$choiceskm<-renderUI({
    nums2=mresp()
    if(is.null(nums2)) return(NULL)
    radioButtons("responsekm", 'Response', choices=nums2,selected=nums2[1],inline = TRUE)
  })
  
  output$boxeskm<-renderUI({
    nums2<-mgrps()
    if(is.null(nums2)) return(NULL)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("varskm", 'Treatment groups', choices=list,selected=list,inline = TRUE)
  })
  
  output$radiokm<-renderUI({
    nums2=refkm()
    if(is.null(nums2)) return(NULL)
    nums2=c('None',nums2)
    radioButtons("radiokm", NULL, choices=nums2,selected=nums2[1],inline=TRUE)
  })
  
  # TFS
  output$slidekmui2<-renderUI({
    valr=censvalresptfs()
  #  print(valr)
    if(is.null(valr)) return(NULL)
    sliderInput(inputId="slidekm2",label=NULL,
                min=min(valr),max=max(valr),value=min(valr),step =diff(valr)[1])
  })
  
  output$sliderkmtui2<-renderUI({
    valt=censvalt()
    if(is.null(valt)) return(NULL)
    args<- list(inputId="sliderkmt2", label=NULL, ticks=TRUE, value=length(valt)-1,
                min=1,max=length(valt),step=1)
    htmlslider1  <- do.call('sliderInput', args)
    htmlslider1$children[[2]]$attribs[['data-values']] <- paste(valt, collapse=',')
    htmlslider1
  })
  
  # OS
  output$slidekmui<-renderUI({
    valr=censvalrespos()
    if(is.null(valr)) return(NULL)
    sliderInput(inputId="slidekm",label=NULL,
                min=min(valr),max=max(valr),value=max(valr),step =diff(valr)[1])
  })

  output$sliderkmtui<-renderUI({
    valt=censvalt()
    if(is.null(valt)) return(NULL)
    args<- list(inputId="sliderkmt", label=NULL, ticks=TRUE, value=length(valt)-1,
                min=1,max=length(valt),step=1)
    htmlslider  <- do.call('sliderInput', args)
    htmlslider$children[[2]]$attribs[['data-values']] <- paste(valt, collapse=',')
    htmlslider
  })

    ####### cross sectional
  output$choicescs<-renderUI({
    nums2=mresp()
    if(is.null(nums2)) return(NULL)
    radioButtons("responsecs", 'Response', choices=nums2,selected=nums2[1],inline = TRUE)
  })
  
  output$choicescstrans<-renderUI({
    trtypes=('None')
    if(is.null(input$responselg)) return(NULL)
    nums2=mresptr()
    if(paste(input$responsecs,".log",sep="") %in% nums2)  trtypes=c(trtypes,'Log')
    if(paste(input$responsecs,".sqrt",sep="") %in% nums2)  trtypes=c(trtypes,'SqRt')
    if(paste(input$responsecs,".curt",sep="") %in% nums2)  trtypes=c(trtypes,'CuRt')
    radioButtons("transcs", "Transformation", choices=trtypes,selected='None',inline = TRUE)
  })
  
  output$boxescs<-renderUI({
    nums2<-mgrps()
    if(is.null(nums2)) return(NULL)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("varscs", 'Treatment groups', choices=list,selected=list,inline = TRUE)
  })
  
  ## slider
  output$slidecsui<-renderUI({
    cdat=dat()
    if(!is.list(cdat)) return(NULL)
    df=cdat$data
    df=df[df$Id%in%cdat$dataM$Id[cdat$dataM$Use],]
    valt=sort(unique(df$Tp))
    #print(valt)
    if(length(valt)==1){ ## only one time points!
      valt=c(valt,valt+1)
      args       <- list(inputId="slidercs", label=NULL, ticks=TRUE, value=0,min=1,max=1,step=1)
      htmlslider  <- do.call('sliderInput', args)
      htmlslider$children[[2]]$attribs[['data-values']] <- paste(args$valt, collapse=',')
      return(htmlslider)
    }
    
    valt0=range(tapply(df$Tp,df$Id,max,na.rm=T))
    if(length(valt)==2) valt0=valt
    if(valt0[1]==valt0[2]) valt0[1]=max(valt[valt<valt0[1]])
    args       <- list(inputId="slidercs", label=NULL, ticks=TRUE, value=match(valt0,valt)-1,
                       min=1,max=length(valt),step=1)
    htmlslider  <- do.call('sliderInput', args)
    htmlslider$children[[2]]$attribs[['data-values']] <- paste0(valt, collapse=',');
    htmlslider
  })
  
  ## pairwise cmparisons
  output$grpcsvar<-renderUI({
    nums2=c('All',refcs())
    if(is.null(nums2)) return(NULL)
    list<-list()
    for(i in nums2) list[[i]]<-i
    checkboxGroupInput("grpcsvar", label =NULL,inline = T,choices =unname(nums2), selected = 'All')
  })

  #############################################################################################
  #################### Calculations
  
  ###### longitudinal
  modelLG<-eventReactive(input$goButton,{
    
    if(is.null(input$responselg) | is.null(input$varslg)) return(NULL)
    if(any(!input$responselg%in%mresp()) | any(!input$varslg%in%mgrps() )) return(NULL)
    cdat<-dat()
    if(is.null(cdat)) return(NULL)
    objres=getLGmat(cdat,resp=input$responselg,lgrps=input$varslg,
             gcols=colGrps()[levels(cdat$dataM$Grp)],trans=input$translg)
    bfco=ifelse(!input$showPanelLG2,0,as.numeric(gsub('p<','',input$bfcolg)))
    if(!input$showPanelLG1)  res=compModLG(objres,bfco=bfco)
    if(input$showPanelLG1)  res=compModLGpw(objres,bfco=bfco,tpcut=as.numeric(input$tpcut))
    res
  })
  

  modelLGpw<-reactive({
    if(is.null(input$responselg)) return(NULL)
  #  if(any(!input$responselg%in%mresp()) | any(!input$varslg%in%mgrps() )) return(NULL)
    objres=modelLG()
    if(is.null(objres)) return(NULL)
    formatLGpw(objres,input$radiolg,padjust=gsub(" .*","",tolower(input$radiolgadj)))
  })
  
  diagLG<-reactive({
    if(is.null(input$responselg) | is.null(input$varslg)) return(NULL)
    objres=modelLG()
    if(is.null(objres)) return(NULL)
    prepDiagLG(objres)
  })
  
  ###### survival
  ndfKM<-reactive({
    if(is.null(input$responsekm)) return(NULL)
    if(is.null(input$slidekm2) & is.null(input$slidekm)) return(NULL)
    if(any(!input$responsekm%in%mresp()) | any(!input$varskm%in%mgrps() )) return(NULL)
    
    cdat<-dat()
    ltps=sort(unique(cdat$data$Tp))
#    print(input$survTyp)
#    print(input$slidekm)
    if(input$survTyp==1) 
      return(getOSTab(cdat,input$responsekm,input$varskm, gcols=colGrps()[levels(cdat$dataM$Grp)],
                      lastT=ltps[input$sliderkmt+1],lastM=input$slidekm))
    if(input$survTyp==2) 
      return(getTFSTab(cdat,input$responsekm,input$varskm,gcols=colGrps()[levels(cdat$dataM$Grp)],
                       lastT=ltps[input$sliderkmt2+1],firstM=input$slidekm2))
  })
  
  modelKM<-reactive({
    if(is.null(input$responsekm)) return(NULL)
    if(is.null(input$slidekm2) & is.null(input$slidekm)) return(NULL)
    if(input$radiokm=='None') return(NULL)
    objres=ndfKM()
    if(is.null(objres)) return(NULL)
    compKM(objres,ref=input$radiokm,firth=input$radiokmFirth,padjust=gsub(" .*","",tolower(input$radiokmadj)))
  })

  ###### cross sectional
  csectDF<-reactive({
    if(any(!input$responsecs%in%mresp()) | any(!input$varscs%in%mgrps() ) | is.null(input$transcs)) return(NULL)
    cdat<-dat()
    if(is.null(cdat)) return(NULL)
    rangetp=sort(unique(cdat$data$Tp))
    getCSmat(cdat,resp=input$responsecs,lgrps=input$varscs,gcols=colGrps()[levels(cdat$dataM$Grp)],
             trans=input$transcs,rangetp=rangetp[input$slidercs+1],usemax=grepl("Max",input$radiocsMax))
  })
  

  modelCS<-reactive({
    objres<-csectDF()
    if(is.null(objres)) return(NULL)
    bfco=ifelse(!input$showPanelCS1,0,as.numeric(gsub('p<','',input$bfcocs)))
    compCS(objres,bfco=bfco,checkvar=input$radiocsvar)
  })
  
  modelCSpw<-reactive({
    if(any(!input$responsecs%in%mresp()) | any(!input$varscs%in%mgrps()) | is.null(input$grpcsvar)) return(NULL)
    objres=modelCS()
    if(is.null(objres)) return(NULL)
    formatCSpw(objres,input$grpcsvar,padjust=gsub(" .*","",tolower(input$radiocsadj)))
  })
  
  diagCS<-reactive({
    if(any(!input$responsecs%in%mresp()) | any(!input$varscs%in%mgrps()) | is.null(input$grpcsvar)) return(NULL)
    objres=modelCS()
    if(is.null(objres)) return(NULL)
    prepDiagCS(objres)
  })
  
  
  ############################################################################
  ########## output
  ### data upload
  output$design <- DT::renderDataTable({
    cdat=dat()
    if(!is.list(cdat)) return(as.table(matrix(cdat,dimnames=list("Error"," "))))
    
    lmids=cdat$dataM$Id[cdat$dataM$Use]
    lgdf=cdat$data[cdat$data$Id%in%lmids,]
    ltps=sort(unique(cdat$data$Tp))
    if(!is.null(input$responsedatexp))
      if(input$responsedatexp!='All') lgdf=lgdf[!is.na(cdat$data[,input$responsedatexp]),]
    lgdf$Grp=factor(cdat$dataM[lgdf$Id,]$Grp)
    tab=table(lgdf$Grp,factor(lgdf$Tp,levels=ltps))
#    print(as.matrix(tab))
#    return(tab)
    return(as.data.frame.matrix(tab))
  },options = list(paging = FALSE,searching = FALSE,autoWidth = TRUE,info=FALSE),server=TRUE)
  
  
  output$filetableshort <- DT::renderDataTable({
    if(is.null(input$responsedat) | is.null(input$groupsdat)) return(NULL)
    if(!any(input$responsedat%in%mresp()) | !any(input$groupsdat%in%mgrps())) return(NULL)
    cdat=dat()
    if(!is.list(cdat)) return(as.table(matrix(cdat,dimnames=list("Error"," "))))
    
    top=cbind(cdat$dataM[cdat$data$Id,c("Id","Grp","Use","Surv"),],cdat$data[,c("Tp",input$responsedat)])
    top=top[top$Grp%in%input$groupsdat,] 
    
    ltps=sort(unique(top$Tp))
    add=do.call('rbind',tapply(1:nrow(top),top[,"Id"],function(x)
      round(top[x[match(ltps,top[x,'Tp'])],input$responsedat],2)))
    colnames(add)=ltps
    tow=cbind(Use=tapply(as.character(top$Use),top$Id,unique),
              Surv=tapply(as.character(top$Surv),top$Id,unique),
              Grp=tapply(as.character(top$Grp),top$Id,unique),
              add)
    rownames(tow)=tapply(as.character(top$Id),top$Id,unique)
    tow=tow[order(factor(tow[,"Grp"],levels=levels(top$Grp))),]
    tow
  },options = list(paging = FALSE,searching = FALSE,autoWidth = TRUE,info=FALSE),server=TRUE,
  selection = list(mode = 'multiple', selected = l2Excl()))

  ###### line charts
  output$plottcs<-renderChart({
    if(any(!input$response%in%mresp()) | any(!input$vars%in%mgrps() | is.null(input$trans))) return(NULL)
    cdat=dat()
    if(is.null(cdat)) return(NULL)
    lgdata=getLGmat(cdat,resp=input$response,lgrps=input$vars,
                    gcols=colGrps()[levels(cdat$dataM$Grp)],trans=input$trans)
    h1=plotLineC(lgdata,type=input$tcfplottyp,
                 force2zero=input$showPanel2,defzero=as.numeric(input$tcdef),
                 miny=as.numeric(input$tcminy),maxy=as.numeric(input$tcmaxy))$plot
    h1$addParams(dom = 'plottcs')
    return(h1)
  })
  
  ###### longitudinal
  output$plotdiaglg<-renderChart({
    if(is.null(input$responselg)) return(NULL)
    fittedLG<-diagLG()
    if(is.null(fittedLG)) return(NULL)
    #    print(str(fittedLG$data))
    p1=plotDiag(fittedLG,typplot=input$radiodiaglg)
    p1$addParams(dom = 'plotdiaglg')
    return(p1)
  })
  
  output$modelLGsum<-shiny::renderTable(modelLG()$coefTab,align="lcccc",include.rownames = T,include.colnames = T)
  output$modelLGeffect<-shiny::renderTable(data.frame(modelLG()$AnovaTab),align="cccc",include.rownames = F,include.colnames = T)
  output$modelLGpw<-DT::renderDataTable({modelLGpw()},
                    options = list(paging = FALSE,searching = FALSE,autoWidth = TRUE,info=FALSE),rownames = FALSE)
  output$modelLGwei<-renderPrint(modelLG()$weightCoef)
  output$modelLGcor<-renderPrint(modelLG()$corrCoef)
  output$outLG<-renderPrint({
    idat=modelLG()$data
    if(!any((idat$out))) return('No outliers')
    idat=idat[idat$out,]
    louts=tapply(1:nrow(idat),factor(idat$Id),function(x) 
      paste(idat$Id[x[1]]," (",idat$Grp[x[1]],"): ",paste(idat$Tp[x],collapse=","),sep=""))
    c("Outliers: ",unname(louts))
  })
  
  ###### survival
  output$plotkm<-renderChart({
    if(is.null(input$responsekm)) return(NULL)
    ndf=ndfKM()
    if(is.null(ndf)) return(NULL)
    p1=plotKM(ndf,shift=as.numeric(input$kmshift),
              lwd=as.numeric(input$kmlwd),
              title=ifelse(ndf$Typ=="OS","Perc. surviving","Perc. tumour-free"))
    p1=p1$plot
    p1$addParams(dom = 'plotkm')
    return(p1)
  })
  
  
  output$sumKM <- shiny::renderTable({
    re=sumIdKM(ndfKM())
#    print(re)
    re
    },align="llllll",include.rownames = F,include.colnames = T)

  
  output$modKM<-DT::renderDataTable({
    if(is.null(input$responsekm) |is.null(input$radiokm)) return(NULL)
    if(input$radiokm=='None') return(NULL)
    data.frame(modelKM()$modTab)},
    options = list(paging = FALSE,searching = FALSE,autoWidth = TRUE,info=FALSE), rownames = FALSE)
  
  output$hrKM<-DT:::renderDataTable({
    if(is.null(input$responsekm) |is.null(input$slidekm) |is.null(input$radiokm)) return(NULL)
    if(input$radiokm=='None') return(NULL)
    modelKM()$hrTab},
    options = list(paging = FALSE,searching = FALSE,autoWidth = TRUE,info=FALSE), rownames = FALSE)

  ###### cross-sectional
  output$plotcs<-renderChart({
    if(is.null(input$responsecs)) return(NULL)
    objres=csectDF()
    p1=plotCS(objres,miny=as.numeric(input$csminy),maxy=as.numeric(input$csmaxy))$plot
    p1$addParams(dom = 'plotcs')
    return(p1)
  })
  
  output$plotdiagcs<-renderChart({
    if(is.null(input$responsecs)) return(NULL)
    fittedCS<-diagCS()
    if(is.null(fittedCS)) return(NULL)
    p1=plotDiag(fittedCS,typplot=input$radiodiagcs)
    p1$addParams(dom = 'plotdiagcs')
    return(p1)
  })
  
  output$sumCS <- shiny::renderTable({
    objres=csectDF()
    if(is.null(objres)) return(NULL)
    v=objres$Df$Id
    if(grepl("range",objres$Par)) v=paste(objres$Df$Id," (",objres$Df$Tp,")",sep="")
    names(v)=NULL
    lso=order(-objres$Df$Tp)
    v=tapply(v[lso],objres$Df$Grp[lso],sort)
    nnamis=tapply(v[lso],objres$Df$Grp[lso],length)
    ltps=tapply(objres$Df$Tp[lso],objres$Df$Grp[lso],function(x){
      itp=range(x)
      if(itp[1]==itp[2]) return(as.character(itp[1]))
      paste(itp,collapse='-')
    })
    re=data.frame(cbind(Group=names(v),N=nnamis,Tp=ltps,Animals=unlist(sapply(v,paste,collapse=" "))))
    return(re)
  },align="llll",include.rownames = F,include.colnames = T)
  
  
  output$ctCS<-DT::renderDataTable(modelCSpw(),
            options = list(paging = F,searching = FALSE,autoWidth = TRUE,info=FALSE), rownames = FALSE)
 output$modCSeffect<-shiny::renderTable(data.frame(modelCS()$AnovaTab),align="lll")


    output$modCS<-shiny::renderTable(data.frame(modelCS()$coefTab),align="llll",include.rownames = T,include.colnames = T)
  output$modCSwei<-renderPrint(modelCS()$weightCoef)
  output$outCS<-renderPrint({
    idat=modelCS()$data
    if(!any((idat$out))) return('No outliers')
    idat=idat[idat$out,]
    louts=tapply(1:nrow(idat),factor(idat$Id),function(x) 
      paste(idat$Id[x[1]]," (",idat$Grp[x[1]],"): ",paste(idat$Tp[x],collapse=","),sep=""))
    c("Outliers: ",unname(louts))
  })
  

  #########################################################################################################################
  ####### Reportting
  
  output$exporttxtDS<-downloadHandler(
    filename = function() paste(fileData()$Ori,"-data.tsv",sep=""),
    content = function(file) {
      what=sapply(downloadFile(dat(),ndigit=as.numeric(input$ndigits),trans=input$exporttrans),paste,collapse='\t')
      cat(what, file = file, sep = "\n")
    }
  )
  
  #############
  output$exporttxtLG <- downloadHandler(
    filename = function() {
      paste(fileData()$Ori,"-LG.", sep = '', switch(
        input$formatLG, PDF = 'pdf', HTML = 'html', DOCX = 'docx'
      ))
    },
    content = function(file) {
      src <- normalizePath('report/reportLG.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'reportLG.Rmd',overwrite=T)
      
      ifile=fileData()$Ori
      objres=modelLG()
      formatpw=modelLGpw()
      plotdata=prepDiagLG(objres)
      pvadj=gsub(" .*","",tolower(input$radiolgadj))
      out <- render('reportLG.Rmd', switch(
        input$formatLG,PDF = pdf_document(), HTML = html_document(), DOCX = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  #############
  output$exporttxtKM <- downloadHandler(
    filename = function() {
      paste(fileData()$Ori,"-KM.", sep = '', switch(
        input$formatKM, PDF = 'pdf', HTML = 'html', DOCX = 'docx'
      ))
    },
    content = function(file) {
      src <- normalizePath('report/reportKM.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'reportKM.Rmd',overwrite=T)
      
      ifile=fileData()$Ori
      objres=ndfKM()
#      print(objres)
      modKM=modelKM()
      p1=plotKM(objres,shift=as.numeric(input$kmshift),lwd=as.numeric(input$kmlwd),retplot=FALSE)
      pvadj=gsub(" .*","",tolower(input$radiokmadj))
      
      out <- render('reportKM.Rmd', switch(
        input$formatKM,PDF = pdf_document(), HTML = html_document(), DOCX = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  #############
  output$exporttxtCS <- downloadHandler(
    filename = function() {
      paste(fileData()$Ori,"-CS.", sep = '', switch(
        input$formatCS, PDF = 'pdf', HTML = 'html', DOCX = 'docx'
      ))
    },
    content = function(file) {
      src <- normalizePath('report/reportCS.Rmd')
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'reportCS.Rmd',overwrite=T)
      
      #########
      ifile=fileData()$Ori
      objres=csectDF()
      csres=modelCS()
      formatpw=modelCSpw()
      diagdata=prepDiagCS(csres)
      retplotcs=plotCS(objres,retplot=FALSE,
                miny=as.numeric(input$csminy),maxy=as.numeric(input$csmaxy))
      pvadj=gsub(" .*","",tolower(input$radiocsadj))
      #########
      out <- render('reportCS.Rmd', switch(
        input$formatCS,PDF = pdf_document(), HTML = html_document(), DOCX = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  #########################################################################################################################
  ####### Download images
  output$downloadCSsvg <- downloadHandler(
    filename <- function() paste(fileData()$Ori,'-CS.svg',sep=""),
    content <- function(file) {
      require(svglite)
      svglite(file, width =as.numeric(input$cswidth), height =as.numeric(input$csheight))
      
      objres=csectDF()
      if(!is.null(objres)){
        p2=plotCS(objres,
                  miny=as.numeric(input$csminy),maxy=as.numeric(input$csmaxy),retplot=FALSE)
        par(mar=c(3.4,4,1,.1),lwd=2,cex=input$cscex,cex.axis=input$cscex,cex.lab=input$cscex)
        exportCS(p2,dogrps=TRUE,cexpt=as.numeric(input$cscexpt))
      }
      dev.off()
      transsvg(file)
    },
    contentType = 'image/svg'
  )
  
  output$downloadCSpng <- downloadHandler(
    filename <- function() paste(fileData()$Ori,'-CS.png',sep=""),
    content <- function(file) {
      png(file, width =as.numeric(input$cswidth)*96, height =as.numeric(input$csheight)*96)
      objres=csectDF()
      if(!is.null(objres)){
        p2=plotCS(objres,
                  miny=as.numeric(input$csminy),maxy=as.numeric(input$csmaxy),retplot=FALSE)
        par(mar=c(3.4,4,1,.1),lwd=2,cex=input$cscex,cex.axis=input$cscex,cex.lab=input$cscex)
        exportCS(p2,dogrps=FALSE,cexpt=as.numeric(input$cscexpt))
      }
      dev.off()
    },
    contentType ='image/png'
  )
  
  ##################################################################
  output$downloadKMsvg <- downloadHandler(
    filename <- function() paste(fileData()$Ori,'-KM.svg',sep=""),
    content <- function(file) {
      require(svglite)
      svglite(file, width =as.numeric(input$kmwidth), height =as.numeric(input$kmheight))
      ndf=ndfKM()
      if(!is.null(ndf)){
        p1=plotKM(ndf,shift=as.numeric(input$kmshift),
                  lwd=as.numeric(input$kmlwd),retplot=FALSE)
        par(mar=c(5.5,5,.5,.5),lwd=2,cex=input$kmcex,cex.axis=input$kmcex,cex.lab=input$kmcex)
        exportKM(p1,lwd=as.numeric(input$kmlwd),cexpt=as.numeric(input$kmcexpt),dogrps=TRUE)
      }
      dev.off()
      transsvg(file)
    },
    contentType ='image/svg'
  )
  
  output$downloadKMpng <- downloadHandler(
    filename <- function() paste(fileData()$Ori,'-KM.png',sep=""),
    content <- function(file) {
      png(file, width =as.numeric(input$kmwidth)*96, height =as.numeric(input$kmheight)*96)
      
      ndf=ndfKM()
      if(!is.null(ndf)){
        p1=plotKM(ndf,shift=as.numeric(input$kmshift),
                  lwd=as.numeric(input$kmlwd),retplot=FALSE)
        par(mar=c(5.5,5,.5,.5),lwd=2,cex=input$kmcex,cex.axis=input$kmcex,cex.lab=input$kmcex)
        exportKM(p1,lwd=as.numeric(input$kmlwd),cexpt=as.numeric(input$kmcexpt),dogrps=FALSE)
      }
      dev.off()
    },
    contentType = 'image/png'
  )
  
  ##################################################################
  output$downloadTCpng <- downloadHandler(
    filename <- function() paste(fileData()$Ori,'-TC.png',sep=""),
    content <- function(file) {
      png(file, width =as.numeric(input$tcwidth)*96, height =as.numeric(input$tcheight)*96)
      
      par(mar=c(3.4,4,1,.1),lwd=2,lend=0,cex=input$tccex,cex.axis=input$tccex,cex.lab=input$tccex)
      cdat=dat()
#      print(colGrps())
      lgdata=getLGmat(cdat,resp=input$response,lgrps=input$vars,gcols=colGrps()[as.character(input$vars)],trans=input$trans)
      h1=plotLineC(lgdata,type=input$tcfplottyp,
                   force2zero=input$showPanel2,defzero=as.numeric(input$tcdef),
                   miny=as.numeric(input$tcminy),maxy=as.numeric(input$tcmaxy))[-1]
      exportTC(h1,dogrps=F,type=input$tcfplottyp)
      dev.off()
    },  contentType = 'image/png'
  )
  
  output$downloadTCsvg <- downloadHandler(
    filename <- function() paste(fileData()$Ori,'-TC.svg',sep=""),
    content <- function(file) {
      require(svglite)
      svglite(file, width =as.numeric(input$tcwidth), height =as.numeric(input$tcheight),standalone = TRUE)
      par(mar=c(3.4,4,1,.1),lwd=2,lend=0,cex=input$tccex,cex.axis=input$tccex,cex.lab=input$tccex)
      cdat=dat()
      lgdata=getLGmat(cdat,resp=input$response,lgrps=input$vars,gcols=colGrps()[as.character(input$vars)],trans=input$trans)
      h1=plotLineC(lgdata,type=input$tcfplottyp,
                   force2zero=input$showPanel2,defzero=as.numeric(input$tcdef),
                   miny=as.numeric(input$tcminy),maxy=as.numeric(input$tcmaxy))[-1]
      exportTC(h1,dogrps=T,type=input$tcfplottyp)
      dev.off()
      transsvg(file)
    },  contentType = 'image/svg'
  )
  
  #############
  
})
