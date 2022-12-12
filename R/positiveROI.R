positiveROI <- function(x){ # x: cycif obj
  this.ct <- as.character(cellTypes(x))
  this.ct <- sub("(Mac,.+) PDL","\\1, PDL",this.ct)
  # this.ct <- sub(",[^,]+$","",this.ct)
  this.ct <- sub(",.+$","",this.ct)
  this.ct[this.ct == "T"] <- "T, other"
  this.ct[this.ct =="Imm"] <- "Imm, other"
  this.ct <- factor(this.ct,levels=c("Tumor","CD8T","T, other","Mac","Imm, other","unknown"))
  slidePlot(this,type="cell_type",cell_type=this.ct,ttl = smpl,uniq.cols=c(brewer.pal(5,"Set3"),"grey70"),ncells=1e5,pch=".")
  
  if(manual){
  cat("Do you want to modify the 'dna.ths'?\n")
  ans <- "init"
  lo <- dna.ths1[i]
  hi <- dna.ths2[i]
  while(ans!="0"){
    ans <- readline(prompt="0:no, 1:low, 2:high, 3:both [0]")
    if(ans=="1"){
      lo <- locator(1)$x
      abline(v=lo,col=4)
    }else if(ans=="2"){
      hi <- locator(1)$x
      abline(v=hi,col=2)
    }else if(ans=="3"){
      bdr <- locator(2)$x
      bdr <- sort(bdr,decreasing=F)
      lo <- bdr[1]
      hi <- bdr[2]
      abline(v=lo,col=4)
      abline(v=hi,col=2)
    }
  }
}