# IGT plot
rm(list=ls(all=TRUE))
par(mfrow=c(1,2))

IGTBlockDeck=read.table("IGTBlockXDeckBehavData.txt")

colors = c("orangered","darkgoldenrod1", "blue3", "green2")
pchs = c(21,24,19,17)
lwds = c(3, 3, 5, 5)
xlabs= c("Block", "","","")
ylabs= c("Proportion of Choices from Each Deck", "", "", "")
par(oma=c(0,2,0,0))
par(mar=c(5,4.1,2,2))
for (i in 1:4)  {

  if (i ==1) {
    plot(1:6,IGTBlockDeck[,i],ylim=c(0,0.5), type="b"
       ,col=colors[i],pch=pchs[i], cex = 2, lwd=lwds[i],xlab="Block",ylab="Proportion of Choices from Each Deck"
       ,cex.lab=1.5, cex.axis=1.3, cex.main=2
       ,main="IGT")
  } else {
    lines(1:6,IGTBlockDeck[,i],ylim=c(0,0.5), type="b"
       ,col=colors[i],pch=pchs[i], cex = 2, lwd=lwds[i],xlab="",ylab=""
       ,cex.lab=1.5, cex.axis=1.3, cex.main=2, xaxt="n", yaxt="n"
       ,main="")
  }
}

leg.txt <- c("Deck A","Deck B","Deck C","Deck D")
legend("topleft", leg.txt,col=c("orangered","darkgoldenrod1","blue3","green2"),cex=1.2
	       ,bty="n",lty=1,pch=pchs,pt.bg=rep("white",4),lwd=lwds/2)		#legend

#par(oma=c(0,0,0,0))
par(mar=c(5,2,2,4.1))
#par(mgp=c(0,1,0))

# SGT plot
SGTBlockDeck=read.table("SGTBlockXDeckBehavData.txt")

colors = c("orangered","darkgoldenrod1", "blue3", "green2")
pchs = c(21,24,19,17)
lwds = c(3, 3, 5, 5)

for (i in 1:4)  {
  if (i ==1) {
    plot(1:6,SGTBlockDeck[,i],ylim=c(0,0.5), type="b"
       ,col=colors[i],pch=pchs[i], cex = 2, lwd=lwds[i],xlab="Block",ylab="Proportion of Choices from Each Deck"
       ,cex.lab=1.5, cex.axis=1.3, cex.main=2
       ,main="SGT")
  } else {
    lines(1:6,SGTBlockDeck[,i],ylim=c(0,0.5), type="b"
       ,col=colors[i],pch=pchs[i], cex = 2, lwd=lwds[i],xlab="",ylab=""
       ,cex.lab=1.5, cex.axis=1.3, cex.main=2, xaxt="n", yaxt="n"
       ,main="")
  }
}

leg.txt <- c("Deck A","Deck B","Deck C","Deck D")
legend("topleft", leg.txt,col=c("orangered","darkgoldenrod1","blue3","green2"),cex=1.2
	       ,bty="n",lty=1,pch=pchs,pt.bg=rep("white",4),lwd=lwds/2)		#legend

dev.copy2eps(file="BlockByDeckIGTSGT.eps" )