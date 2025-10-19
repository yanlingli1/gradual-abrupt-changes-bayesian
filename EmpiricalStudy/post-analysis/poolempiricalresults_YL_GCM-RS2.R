library(ggplot2)
library(ggpubr)
library(jtools)

load("../data/datblock_Meaning.Rdata")

load("../Updated_resulttable_rstvar_gcm_empirical_IIV_final_0703_GCM-RS2.Rdata")

plotres=function(pid, id, IDnum){
  Y = datblock[datblock$PID==pid,"Meaning"]
  O = length(Y)
  p_2s = c()
  for(t in 1:O){
    p_2s[t] = as.numeric(resulttable[paste0("midx[",id,",",t,"]"), ][1])-1
  }
  y = c(Y, p_2s*100)
  df = as.data.frame(cbind(1:O, y))
  colnames(df) = c("time","y")
  df$variable = c(rep("Meaning",O),
                  rep("P(S=2)",O))
  
  gg = ggplot(df,aes(x = time, y= y)) + 
    geom_area(df[df$variable=="P(S=2)",], mapping=aes(x=time, y=y), fill="snow2")+
    geom_line(df[df$variable %in% c("Meaning"),], mapping=aes(x = time, y= y, color=variable, linetype=variable))+ 
    #scale_size_manual(values=c(0.2,1,0.2))+
    xlab("time")  +
    scale_y_continuous(name = "Meaning", 
                       sec.axis = sec_axis(~./100, name = "P(S=2)"),
                       limits = c(0, 100))+
    ggtitle(paste("Participant", IDnum))+
    theme_apa()+
    theme(legend.position = "none")
  return(gg)
}

#97->68; 108->74; 98->69; 4->2
p1 = plotres(98,69,1)
p2 = plotres(4,2,2)
p3 = plotres(97,68,3)
p4 = plotres(108,74,4)

pdf("pred.pdf", width=10, height=8)
pp = ggarrange(p1,p2,p3,p4, 
               ncol=2, nrow=2)

annotate_figure(pp, top = text_grob("Meaning Dynamics with Model-implied Regime Probabilities",face = "bold",size=20))
dev.off()

# check entropy values for each participant
pr2 = matrix(NA,108,224)
pr1 = matrix(NA,108,224)
Entropy = c()
for(i in 1:108){
  for(t in 1:224){
    pr2[i,t] = as.numeric(resulttable[paste0("midx[",i,",",t,"]"), ][1])-1
    pr1[i,t] = 1-pr2[i,t]
  }
  Entropy[i] = 1 + (sum(pr2[i,]*log(pr2[i,]),na.rm=T)+sum(pr1[i,]*log(pr1[i,]),na.rm=T))/224/log(2)
}
#which(Entropy>0.8)
highentropyids = c(25, 31, 74, 79, 116, 125, 139)
hist(Entropy, main="(A) Histogram of Participants' Entropy Values from the GCM-RS2 Model")
