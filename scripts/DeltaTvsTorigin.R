

DeltaT=c(10^(-4), 10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3), 10^(4))
# T_origin_model_time=c(0.782507, 
#                       0.782594,
#                       0.707956,
#                       0.750964,
#                       0.578051,
#                       0.195762,
#                       0.040818,
#                       0.006377,
#                       0.000864
#                       )
T=0.5
T_origin_model_time_n_2=c(0.272667, 
                      0.268560,
                      0.268887,
                      0.269333,
                      0.271363,
                      0.195587,
                      0.040803,
                      0.006368,
                      0.000864
)


T_origin_model_time_n_20=c(0.378114, 
                           0.378584,
                           0.376758,
                           0.376904,
                           0.374960 ,
                           0.276482,
                           0.051258,
                           0.007405,
                           0.000971 
)
library(Hmisc)
library(dplyr)
library(grid)
library(gtable)
library(grDevices)
library(gridtext)
dat=data.frame(DeltaT=DeltaT,TOriginModelTime= T_origin_model_time_n_20 )

grob1 <- grobTree(textGrob("Time of origin", x=0.1,  y=0.9, hjust=0,
                           gp=gpar(col="black", fontsize=13, fontface="bold")))

grob2 <- grobTree(textGrob(expression('Effect. pop size 10^5'), x=0.1,  y=0.3, hjust=0,
                           gp=gpar(col="black", fontsize=13, fontface="bold")))


grob3 <- grobTree(textGrob(expression('Sample size 2'), x=0.1,  y=0.2, hjust=0,
                           gp=gpar(col="black", fontsize=13, fontface="bold")))

grob4 <- grobTree(textGrob(expression('Sample size 20'), x=0.1,  y=0.2, hjust=0,
                           gp=gpar(col="black", fontsize=13, fontface="bold")))

ggplot(dat) +
  geom_point(aes(x=log10(DeltaT), y=TOriginModelTime))+
  geom_smooth(aes(x=log10(DeltaT), y=TOriginModelTime))+
  ggtitle("Mean time of MRCA in model time vs log DeltaT")+
  geom_hline(yintercept=T, linetype="dashed", color = "red")+
  labs(x="log base 10 of DeltaT", y= "time of MRCA ")+
  theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=15),
        axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=15),
        plot.title = element_text(size=15, face="bold", vjust=2),
        axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
        axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
        legend.title = element_blank(),
        legend.position = c(.30,.70),
        legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        annotation_custom(grob1)+
  annotation_custom(grob2)+
  annotation_custom(grob4)
  
path_to_save="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/"
ggsave(paste(path_to_save, "/", "MRCAvsLogDeltaTn=20", ".png", sep=""))