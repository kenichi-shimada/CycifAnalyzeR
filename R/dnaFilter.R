#' @export
setGeneric("dnaFilter", function(x,...) standardGeneric("dnaFilter"))
setMethod("dnaFilter", "Cycif",
          function(x,manual=FALSE,ratio=TRUE,n=1000,n1=1000){
            mat <- x@dna
            smp <- x@name

            dna.list <- names(mat)
            if(ratio){
              mat <- as.data.frame(lapply(mat,function(x)log1p(x/mat[[1]])))
            }else{
              mat <- log1p(mat)
            }

            xmax <- ceiling(max(mat)/n)*n
            xmax <- max(mat)
            brks <- seq(0,xmax,length=(n+1))

            dna.ths1 <- rep(0,length(dna.list))
            dna.ths2 <- rep(Inf,length(dna.list))
            names(dna.ths1) <- names(dna.ths2) <- dna.list

            cat(paste0("Filtering ",smp,"...\n"))

            ## all DNA channels
            if(manual){
              a <- hist(m1[used],breaks=brks,main=channel,freq=FALSE)

              loessMod <- loess(a$density[seq(n1)] ~ brks[seq(n1)], span=0.02)
              smoothed <- predict(loessMod)

              lines(smoothed, x=brks[seq(n1)], col=2)

              # points(a$mids[ind.p],a$density[ind.p],pch=20,col=2)
              # points(a$mids[ind.v],a$density[ind.v],pch=20,col=4)

              th <- locator(1)$x
              if(th < 0) th <- 0
              abline(v=th,col=4)
              dna.ths1[i] <- th

              th2 <- locator(1)$x
              if(th2 < 0) th2 <- 0
              abline(v=th2,col=2)
              dna.ths2[i] <- th2
            }else{
              inds.v <- inds.p <- list()
              par(mfcol=c(4,2))
              par(mar=c(3,3,1,1))
              par(bg="white",fg="black")

              for(i in 2:length(dna.list)){
                channel <- dna.list[i]
                m <- mat[[channel]]
                xmax <- max(mat)
                brks <- seq(0,xmax,length=(n+1))
                a <- hist(m,breaks=brks,main=channel,freq=FALSE)

                loessMod <- loess(a$density[seq(n1)] ~ brks[seq(n1)], span=0.02)
                smoothed <- predict(loessMod)
                lines(smoothed, x=brks[seq(n1)], col=2)

                th.up <- quantile(m,.9)

                ind.v <- which(diff(sign(diff(smoothed)))>0  & a$mids[c(-1,-n1)] < th.up)
                ind.p <- which(diff(sign(diff(smoothed)))<0  & a$mids[c(-1,-n1)] < th.up)
                idx <- sort(c(ind.v,ind.p))

                idx.adj <- 13
                idx.hi <- 40

                while(any(diff(idx)<= idx.adj)){
                  k <- which.min(diff(idx))
                  idx <- idx[-c(k,k+1)]
                }

                if(length(ind.v)>0){
                  ind.v <- ind.v
                  tmp <- a$mids[max(ind.v)]
                  if(tmp < 0.25){
                    x.drop.th <- tmp
                  }else{
                    x.drop.th <- 0
                  }
                }else{
                  x.drop.th <- 0
                }
                abline(v=x.drop.th,col=4)
                dna.ths1[i] <- x.drop.th

                if(length(ind.p)>0){
                  if(max(ind.p) > idx.hi){
                    p.max <- max(ind.p)
                    x.max <- a$mids[p.max]
                    y.max <- a$density[p.max]
                    x.half.max <- which.min(abs(a$density - y.max/2)[-seq(p.max)])
                    x.bunch.th <- a$mids[p.max + x.half.max*2]
                    dna.ths2[i] <- x.bunch.th
                    abline(v=x.bunch.th,col=2)
                  }else{
                    dna.ths1[i] <- Inf
                  }
                }

                inds.v[[i]] <- ind.v <- ind.v[ind.v %in% idx]
                inds.p[[i]] <- ind.p <- ind.p[ind.p %in% idx]

              }
            }

            used.cells <- sapply(names(mat),function(channel){
              ind <- (mat[[channel]] > dna.ths1[channel]) + (mat[[channel]] > dna.ths2[[channel]])
              return(ind)
            })

            used.cells[,1] <- 1		# used.cells[,1] <- TRUE

            x@dna_thres <- data.frame(low=dna.ths1,high=dna.ths2)
            x@used_cells <- used.cells
            # validObject(x)
            return(x)
          }
)
