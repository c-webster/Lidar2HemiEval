
#####
# This is the R script that is associated with the matlab script 'Lidar2HemiEval_Prep.m'

# Input is given through running Lidar2HemiEval.m from Matlab.

# This script is hard coded and should not be changed.


################################################################################
### BEGIN



# Check packages
check.packages <- function(pkg){
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  
  if (length(new.pkg)) 
    
    install.packages(new.pkg, dependencies = TRUE)
  
  sapply(pkg, require, character.only = TRUE)
  
}

packages<-c("sp", "rLiDAR", "raster", "rlas", "itcSegment", "pracma")

check.packages(packages)


FindTreesCHM <- function(chm, fws, minht) {
  
  if (class(chm)[1] != "RasterLayer") {
    chm <- raster(chm)
  }
  if (class(fws) != "numeric") {
    stop("The fws parameter is invalid. It is not a numeric input")
  }
  if (class(minht) != "numeric") {
    stop("The minht parameter is invalid. It is not a numeric input")
  }
  w <- matrix(c(rep(1, fws * fws)), nrow = fws, ncol = fws)
  chm[chm < minht] <- NA
  f <- function(chm) max(chm)
  rlocalmax <- focal(chm, fun = f, w = w, pad = TRUE, padValue = NA)
  setNull <- chm == rlocalmax
  XYmax <- SpatialPoints(xyFromCell(setNull, Which(setNull == 
                                                     1, cells = TRUE)), proj4string = crs(chm))                # Edited
  htExtract <- over(XYmax, as(chm, "SpatialGridDataFrame"))
  treeList <- cbind(coordinates(XYmax), htExtract)              # Edited
  colnames(treeList) <- c("x", "y", "height")
  return(treeList)
  
} # end of FindTreesCHM
# FindTreesCHM function taken from blog post by Andrew Sanchez Meador:
# (http://quantitativeecology.org/using-rlidar-and-fusion-to-delineate-individual-trees-through-canopy-height-model-segmentation/)

SegmentCan <- function(chmf,imout){
  
  ### import CHM and associated normalised lidar data
  chm    <- raster(chmf)

  # Setting the fixed window size (fws)
  fws <- 3
  # Set the specified height above ground for the detection break
  minht <- 2
  # Create the individual tree detection list and summarize it
  loc <- FindTreesCHM(chm, fws, minht)
  
  # Set the maxcrown parameter - maximum individual tree crown radius expected
  maxcrown = 5.0
  # Set the exclusion parameter - A single value from 0 to 1 that represents
  # the % of pixel exclusion. E.g. a value of 0.5 will exclude all of the
  # pixels for a single tree that has a height value of less than 50% of
  # the maximum height from the same tree. Default value is 0.3.
  exclusion = 0.3
  
  cat("\n Segmenting canopy... \n \n")
  
  # Compute canopy areas for the individual tree detections
  canopy <- ForestCAS(chm, loc, maxcrown, exclusion)
  
  cat("\n \n Finished segmenting canopy. \n")
  
  # Retrieve the boundary for individual tree detection and canopy area calculation
  boundaryTrees <- canopy[[1]]
  
  # remove small tree crowns
  canopy[[2]] <- canopy[[2]][!(canopy[[2]]$ca < 4),] # get rid of the small shitty ones
  
  # save a figure of the segmented canopy for reference
  png(file.path(imout,paste(fname,'Segmented_Canopy.png',sep="")))
  plot(chm)
  lines(boundaryTrees)
  points(loc,pch=20,cex=0.25)
  dev.off()

  return(canopy)

} # end of SegmentCan


ClassTrees <- function(canopy,lasf,species,biome,branches){
  
  lasdat <- read.las(lasf)
  
  cat("\n Calculating tree crown statistics. \n")
  
  # get dbh data using function in itcSegment
  # calculate crown diameter
  canopy[[2]]$cd_m<-2*sqrt(canopy[[2]]$ca/pi)
  treedim <- dbh(H = canopy[[2]]$z, CA = canopy[[2]]$cd_m, biome = biome)
  
  boundaryTrees <- canopy[[1]]
  
  cat("\n Getting tree and branch info. \n")
  
  
  # empty matrix for classifying branch points (ltc)
  lastc         <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  treex         <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  treey         <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  horzdist      <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  ba            <- matrix(data=-NaN,nrow=length(lasdat$X),ncol=1)
  ltcdf    <- data.frame(lasdat$X,lasdat$Y,lasdat$Z,treex,treey,horzdist,lastc,ba)
  colnames(ltcdf)[1] <- "lasdat_X"
  colnames(ltcdf)[2] <- "lasdat_Y"
  colnames(ltcdf)[3] <- "lasdat_Z"
  
  ltcdf$lasdat_Z[ltcdf$lasdat_Z < 1.0] = NaN # don't bother with branches within 1 metre of the ground
  ltcdf <- ltcdf[complete.cases(ltcdf$lasdat_Z), ]

  # empty matrix for tree trunk and crown statistics (dbh)
  treeptsX      <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsY      <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsH_m    <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsCA_m2  <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  treeptsdbh_cm <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  # treeptscd_m   <- matrix(data=NaN,nrow=length(boundaryTrees$Trees),ncol=1)
  dbhdf    <- data.frame(treeptsX,treeptsY,treeptsH_m,treeptsdbh_cm,treeptsCA_m2)


  # classify points in lidar cloud based on what tree boundary it is in
  for (i in seq_along(1:length(boundaryTrees$Trees))){
    
    pix <- boundaryTrees@polygons[[i]]@plotOrder[1]
    outline <- boundaryTrees@polygons[[i]]@Polygons[[pix]]@coords
    
    # find which chm maxima is in that polygon because loc and boundaryTrees don't match
    tpts_in <- point.in.polygon(canopy[[2]]$x,canopy[[2]]$y,outline[,1],outline[,2], mode.checked = FALSE)
    
    if (length(which(tpts_in == 1) > 0)) {
      dbhdf$treeptsX[i]      <- canopy[[2]]$x[which(tpts_in==1)]
      dbhdf$treeptsY[i]      <- canopy[[2]]$y[which(tpts_in==1)]
      dbhdf$treeptsH_m[i]    <- canopy[[2]]$z[which(tpts_in==1)]
      dbhdf$treeptsCA_m2[i]  <- canopy[[2]]$ca[which(tpts_in==1)]
      dbhdf$treeptsdbh_cm[i] <- treedim[which(tpts_in==1)]
    }
    
    if (branches == 1){ # only bother looping through and classifying branches if the user wants them
      # now classify the lidar by which polygon it falls in
      if (length(which(tpts_in > 0) > 0)) {
        pts_in <- point.in.polygon(ltcdf$lasdat_X,ltcdf$lasdat_Y,outline[,1],outline[,2], 
                                   mode.checked = FALSE)
        
        ltcdf$lastc[which(pts_in>0)] <- i
        
        ltcdf$treex[which(pts_in>0)]    <- boundaryTrees@polygons[[i]]@Polygons[[pix]]@labpt[1]
        ltcdf$treey[which(pts_in>0)]    <- boundaryTrees@polygons[[i]]@Polygons[[pix]]@labpt[2]
        ltcdf$horzdist[which(pts_in>0)] <- sqrt(((ltcdf$lasdat_X[which(pts_in>0)]-
                                                 ltcdf$treex[which(pts_in>0)])**2) 
                                           + ((ltcdf$lasdat_Y[which(pts_in>0)]-
                                                 ltcdf$treey[which(pts_in>0)])**2))
        
        # assign a branch angle
        if (species == 1){
          ang <- runif(1, 60, 100)
        } else if (species == 2){
          ang <- runif(1, 50, 100)
        } else if (species == 3){
          ang <- runif(1, 45, 75)
        } else if (species == 4){
          ang <- runif(1, 32, 50)
        } else {
          ang <- 90
        }
        
        ltcdf$ba[which(pts_in>0)] <- ang
        rm(ang)
      
      } # end classifying lidar points
      
    } # end if if(branches==1) 
    
  } # end of i in seq_along (looping through polygons)

  ltcdf$horzdist[ltcdf$horzdist < 0.2] = NaN
  
  cat("\n Finished getting tree crown info. \n")
  
  return(list(ltcdf=ltcdf,dbhdf=dbhdf,canopy=canopy,treedim=treedim))
  
} # end of FindClassBranches


GenGrid <- function(dbhdf,canopy,treedim,xlim,ylim,epsgstr){
  
  cat("\n Generating analysis points on grid... \n")
  
  Xmat <- seq(from = xlim[1], to = xlim[2], by = spacing)
  Ymat <- seq(from = ylim[1], to = ylim[2], by = spacing)
  cnew <- meshgrid(Xmat,Ymat)
  xnew <- c(cnew$X)
  ynew <- c(cnew$Y)
  
  rad <- ((treedim/100) / 2)
  
  df <- as.data.frame(canopy[[2]])
  xy <- df[,c(1,2)]
  
  if (epsgstr == '0000'){ # doesn't remove points in trunks
    
    cat("\n \n Not processing grid points for those in trunks. \n")
    
    
  } else { # remove points within trunks
  
    projstr = strcat('+init=EPSG:',epsgstr)

    spdf <- SpatialPointsDataFrame(coords = xy,data=df,proj4string = CRS(projstr))

    sptA <- gBuffer(spdf,byid=TRUE,id=NULL,width=rad,quadsegs = 10,
                    capStyle = "ROUND",joinStyle = "ROUND")
    
    lastcA <- matrix(data=-NaN,nrow=length(xnew),ncol=1)

    for (i in seq_along(1:length(sptA$X))){
      
      outline <- sptA@polygons[[i]]@Polygons[[1]]@coords
      
      pts_in <- point.in.polygon(xnew,ynew,outline[,1],outline[,2], mode.checked = FALSE)
      
      lastcA[which(pts_in==1)] <- 1
      
    }
    
    ptsdf = data.frame(xnew,ynew)
    ptsdf$xnew[which(lastcA>0)] <- NaN
    ptsdf$ynew[which(lastcA>0)] <- NaN
    
    ptsdf <- ptsdf[complete.cases(ptsdf),]

    return(ptsdf)
    
    cat("\n \n Done. \n")
    
  }

} # end of gengrid


# Define input arguments
args <-  commandArgs(TRUE)

# parse input arguments to variables
basefolder <- args[1]
chmf       <- args[2]
lasf       <- args[3]
fname      <- args[4]
trunks     <- as.numeric(args[5])
branches   <- as.numeric(args[6])
gengrid    <- as.numeric(args[7])
spacing    <- as.numeric(args[8])
epsgstr    <- args[9]
biome      <- as.numeric(args[10])
species    <- as.numeric(args[11])

# basefolder <- unlist(strsplit(basefolder,split=files))

cat("\n \n Done. \n")


# define paths for output

  dbh     <- file.path(basefolder,"Data_Surface","DBH",paste(fname,"_treeinfo.txt",sep=""))
  ltc     <- file.path(basefolder,"Data_Surface","LTC",paste(fname,"_lastreeclass.txt",sep=""))
  pts     <- file.path(basefolder,"Data_Points",paste(fname,"_gridpts.txt",sep=""))
  imout   <- file.path(basefolder,"Output_Prep",fname)

  if (trunks == 1 | branches == 1){

    canopy    <- SegmentCan(chmf,imout)
    output    <- ClassTrees(canopy,lasf,species,biome,branches)
    canopynew <- output$canopy
    treedim   <- output$treedim

    dbhdf  <- output$dbhdf[complete.cases(output$dbhdf), ]
    ltcdf  <- output$ltcdf[complete.cases(output$ltcdf), ]

    if (trunks == 1){
      dbhdf <- round(dbhdf,2)
      write.table(dbhdf,dbh,sep = "\t", row.names = FALSE,
                  col.names = TRUE, quote = FALSE)
    }
    if (branches == 1){
      ltcdf <- round(ltcdf,2)
      write.table(ltcdf,ltc,sep = "\t", row.names = FALSE,
                  col.names = TRUE, quote = FALSE)
    }

    if (gengrid == 1) {

      # get dimensions
      lims    <- read.table(file.path(basefolder,"Output_Prep",fname,paste(fname,"_analysisarea.txt",sep="")))
      xlim <- c(as.numeric(lims[1]),as.numeric(lims[2]))
      ylim <- c(as.numeric(lims[3]),as.numeric(lims[4]))

      ptsdf  <- GenGrid(dbhdf,canopynew,treedim,xlim,ylim,epsgstr)

      write.table(format(ptsdf,digits=0,scientific=FALSE),pts,
                  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    }

  }


  
tmpout   <- file.path(basefolder,"temp",'R_finish_confirmation.txt')
write.table(NaN,tmpout,row.names = FALSE, col.names = FALSE)









