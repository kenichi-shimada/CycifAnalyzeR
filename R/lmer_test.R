#' Violin plots to show protein expressions
#' @export
setGeneric("lmerCT", function(x,...) standardGeneric("lmerCT"))

#' @export
setMethod("lmerCT", "CycifStack",
          function(x,ld_name,cts_type=c("cell_type","cell_type_full"),
                   summarize=FALSE,summarize.by=c("BOR","TimePoint"),
                   order="cBTgP",margins=c(7,5),used.abs,...){
            cr3 <- cr2 <- cr1 <- as.character(pd$BOR)
            cr2[cr1=="PR" | cr1=="SD"] <- "PR_SD"
            cr3[cr1=="SD" | cr1=="PD"] <- "SD_PD"
            cr1 <- factor(cr1,levels=c("PR","SD","PD"))
            cr2 <- factor(cr2,levels=c("PR_SD","PD"))
            cr3 <- factor(cr3,levels=c("PR","SD_PD"))
            pd$CR1 <- cr1
            pd$CR2 <- cr2
            pd$CR3 <- cr3
            
            pd$Best.Response <- factor(pd$Best.Response,levels=c("PR","SD","PD"))
            
            identical(pd$CR1,pd$Best.Response)
            pData(target_data) <- pd
            
(fm1 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy))
summary(fm1)# (with its own print method; see class?merMod % ./merMod-class.Rd
anova(fm1)