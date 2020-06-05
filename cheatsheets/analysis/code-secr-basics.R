library(secr)

traps_csv <- read.table("data/Tost_capthist.txt")


traps_csv <- read.csv("data/Tost_cams3.csv")
write.csv(traps_csv, file = "data/Tost_cams3.csv", row.names = FALSE, col.names = FALSE)
capthist_csv <- read.csv("data/Tost_capthist.csv")

ch <- read.capthist(captfile = "data/Tost_capthist.txt", 
                    binary.usage = FALSE, trapfile = "data/Tost_cams.csv", 
                    detector="count", fmt = "trapID", 
                    trapcovnames = c("Rgd","Topo", "Water", "Winter"))

ch <- read.capthist(captfile = "data/Tost_capthist.txt", 
                    binary.usage = FALSE, trapfile = "data/Tost_cams4.txt", 
                    detector="count", fmt = "trapID")

mask <- make.mask(traps = traps(ch), buffer = 12000, type = "trapbuffer")

plot(mask)
plot(traps(ch), add = TRUE)

m0 <- secr.fit(ch, detectfn="HHN", mask = mask,
               model=list(D~1, lambda0~1, sigma~1))

trapfile <- "data/Tost_cams.csv"
chfile <- "data/Tost_capthist.csv"
traps <- read.traps(file = trapfile, trapID = "X.ID", 
                    detector="count",
                    binary.usage = FALSE, 
                    covnames = c("Rgd","Topo", "Water", "Winter"))
ch <- read.capthist(captfile = chfile,
                    trapfile = trapfile,  
                    binary.usage = FALSE, 
                    detector="count", fmt = "trapID", 
                    trapcovnames = c("Rgd","Topo", "Water", "Winter"))

saveRDS(ch, file = "data/Tost_ch.Rds")


TNN.trapfiles = c(
  "Analysis4paper/Data/Tost_cams.csv",
  "Analysis4paper/Data/Noyon_cams.csv",
  "Analysis4paper/Data/Nemegt_cams.csv"
)
TNN_ch<-read.capthist(captfile = "Analysis4paper/Data/TNN_capthist.csv", 
                      binary.usage = FALSE, trapfile = TNN.trapfiles, 
                      detector="count", fmt = "trapID", 
                      trapcovnames = c("Rgd","Topo", "Water", "Winter"))


my_mask_df <- data.frame(X = c(0,1,0,1), Y = c(0,0,1,1), elevation = c(0,110,80,30))
my_mask <- read.mask(data = my_mask_df, spacing = 1)
covariates(my_mask) <- data.frame(elevation = c(0,110,80,30), temp = c(25,26,36,37))

xx<-make.mask(traps(TNN_ch))
load("Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
load("Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_chTNN_ch

m0 <- secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
               model=list(D~1, lambda0~1, sigma~1))

m1 <- secr.fit(Tost_ch, detectfn="HR", mask=TostMask,
               model=list(D~1, g0~1, sigma~1))

par(mai = c(0.4,0.4,0.4,0.4))
plot(m1, xval = 0:15000, col = "blue", lwd = 3,
     cex.lab=2,ylab="",xlab="",
     yaxt='n',xaxt='n')
title(ylab="Detection prob.", xlab="Distance (m)", line=0, cex.lab=1.4)
plot(m0, xval = 0:15000, add = TRUE, col = "red", lwd = 3)

covariates(TostMask)$elev <- attr(TostMask, "covariates")$stdGC
m2 <- secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
               model=list(D~elev, lambda0~Water, sigma~1))

print(m0)

AIC(m0,m1,m2)

xx <- cbind(TostMask[,1:2], covariates(TostMask)$stdGC)
write.table(xx, file = "data/TostMask.txt", sep = " ", row.names = FALSE)

