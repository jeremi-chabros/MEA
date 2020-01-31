




############
###########
###########
#https://www.psychologie.uni-heidelberg.de/ae/meth/team/mertens/blog/anova_in_r_made_easy.nb.html



#clustering coeff
# gene effect ctx c spike 

#########################################
##########################################
#install.packages(c("afex", "car", "devtools", 
#                  "ggplot2", "lme4", "lmerTest", "lsmeans"))
#############################################
library(car)
library(lme4)
library(lmerTest)
library(lsmeans)
library(devtools)
library(ggplot2)
#3.1.2 ctx FR over development using threshold vs template method
data<-read.csv(file.choose(),header=T) #load data3.1.1.CTXvsHPCmSpikes.csv
head(data)

#mean, median, interquartile range
aggregate(br~group+age, data = data,summary)

my.lmer<-lmer(br ~ age + group  + age*group + (1|id) , data = data)
#my.lmer<-lmer(br ~ age + genotype + (1|id) , data = data)
anova(my.lmer,type = 3)

#check assumptions; normally distributed residuals; equal variance
par(mfrow=c(1,1))
plot(my.lmer)

resid.model<-residuals(my.lmer)
shapiro.test(resid.model) 
par(mfrow=c(1,1))
hist(resid.model) #looks nearly normal
qqnorm(resid.model)
qqline(resid.model)

leveneTest(br~group,data = data)
#need to do mauchlys test for sphericity to test variance between ages
#if significant do greenhouse geiiser coreection
bartlett.test(br~group,data = data)
#if normal; levenes test; if not normal, bartlett test for equal variance between groups (across all ages)


############################
#using afex package
#https://www.psychologie.uni-heidelberg.de/ae/meth/team/mertens/blog/anova_in_r_made_easy.nb.html
#https://github.com/singmann/afex
afexdata<-data[c(4,3,2,1)] #order: ID,group,paired factor,response
library(afex)
head(afexdata)
(fit_nice <- aov_ez("id","br",afexdata,between = "group", within= "age" ,return="nice"))
(fit_all <- aov_ez("id","br",afexdata, within="age", between = "group"))
names(fit_all)
summary(fit_all) #mauchlys test not sig; we can assume sphericity

#check normality
hist(subset(data,group=="HE")$br) #all ages of threshold
hist(subset(data,group=="WT")$br)
hist(subset(data,age=="DIV07"&group=="HE")$br)
hist(subset(data,age=="DIV14"&group=="HE")$br)
hist(subset(data,age=="DIV21"&group=="HE")$br)
hist(subset(data,age=="DIV28"&group=="HE")$br)
hist(subset(data,age=="DIV07"&group=="WT")$br)
hist(subset(data,age=="DIV14"&group=="WT")$br)
hist(subset(data,age=="DIV21"&group=="WT")$br)
hist(subset(data,age=="DIV28"&group=="WT")$br)

qqnorm(subset(data,group=="HE")$br) #all ages of threshold
qqline(subset(data,group=="HE")$br) #all ages of threshold
qqnorm(subset(data,group=="WT")$br) #all ages of threshold
qqline(subset(data,group=="WT")$br) #all ages of threshold
qqnorm(data$br) #all ages of threshold
qqline(data$br) #all ages of threshold
qqnorm(resid.model)
qqline(resid.model)
my.lmer2<-lmer(br ~ age + (1|id) , data = data)
resid.2<-residuals(my.lmer2)
hist(resid.2)
qqnorm(resid.2)
qqline(resid.2)
shapiro.test(resid.2)

my.lm<-lm(br~group, data=data)
anova(my.lm)
resid.g<-residuals(my.lm)
hist(resid.g)
qqnorm(resid.g)
qqline(resid.g)
shapiro.test(resid.g)
par(mfrow=c(2,2))
plot(my.lm)
par(mfrow=c(1,1))

### post hoc non para
install.packages("dunn.test")
library(dunn.test)
pdata<-data[c(1,3)]
uns.data<-unstack(pdata)
dunn.test(uns.data)


#manual calculation of mean and sem for each group
data<-read.csv(file.choose(),header=T) #load data3.1.1.CTXvsHPCmSpikes.csv

toplot<-as.data.frame(aggregate(br~group+age, data = data,summary))
toplot<-toplot[c(1,2)]
meanFRs<-rbind(
  mean(subset(data,age=="DIV07"&group=="HE")$br),
  mean(subset(data,age=="DIV07"&group=="WT")$br),
  mean(subset(data,age=="DIV14"&group=="HE")$br),
  mean(subset(data,age=="DIV14"&group=="WT")$br),
  mean(subset(data,age=="DIV21"&group=="HE")$br),
  mean(subset(data,age=="DIV21"&group=="WT")$br),
  mean(subset(data,age=="DIV28"&group=="HE")$br),
  mean(subset(data,age=="DIV28"&group=="WT")$br))

SEMFRs<-rbind(
  sd(subset(data,age=="DIV07"&group=="HE")$br)/sqrt(length(subset(data,age=="DIV07"&group=="HE")$br)),
  sd(subset(data,age=="DIV07"&group=="WT")$br)/sqrt(length(subset(data,age=="DIV07"&group=="WT")$br)),
  sd(subset(data,age=="DIV14"&group=="HE")$br)/sqrt(length(subset(data,age=="DIV14"&group=="HE")$br)),
  sd(subset(data,age=="DIV14"&group=="WT")$br)/sqrt(length(subset(data,age=="DIV14"&group=="WT")$br)),
  sd(subset(data,age=="DIV21"&group=="HE")$br)/sqrt(length(subset(data,age=="DIV21"&group=="HE")$br)),
  sd(subset(data,age=="DIV21"&group=="WT")$br)/sqrt(length(subset(data,age=="DIV21"&group=="WT")$br)),
  sd(subset(data,age=="DIV28"&group=="HE")$br)/sqrt(length(subset(data,age=="DIV28"&group=="HE")$br)),
  sd(subset(data,age=="DIV28"&group=="WT")$br)/sqrt(length(subset(data,age=="DIV28"&group=="WT")$br)))

meanFRs<-rbind(
  mean(subset(data,age=="DIV07"&group=="ctx")$br),
  mean(subset(data,age=="DIV07"&group=="hpc")$br),
  mean(subset(data,age=="DIV14"&group=="ctx")$br),
  mean(subset(data,age=="DIV14"&group=="hpc")$br),
  mean(subset(data,age=="DIV21"&group=="ctx")$br),
  mean(subset(data,age=="DIV21"&group=="hpc")$br),
  mean(subset(data,age=="DIV28"&group=="ctx")$br),
  mean(subset(data,age=="DIV28"&group=="hpc")$br))

SEMFRs<-rbind(
  sd(subset(data,age=="DIV07"&group=="ctx")$br)/sqrt(length(subset(data,age=="DIV07"&group=="ctx")$br)),
  sd(subset(data,age=="DIV07"&group=="hpc")$br)/sqrt(length(subset(data,age=="DIV07"&group=="hpc")$br)),
  sd(subset(data,age=="DIV14"&group=="ctx")$br)/sqrt(length(subset(data,age=="DIV14"&group=="ctx")$br)),
  sd(subset(data,age=="DIV14"&group=="hpc")$br)/sqrt(length(subset(data,age=="DIV14"&group=="hpc")$br)),
  sd(subset(data,age=="DIV21"&group=="ctx")$br)/sqrt(length(subset(data,age=="DIV21"&group=="ctx")$br)),
  sd(subset(data,age=="DIV21"&group=="hpc")$br)/sqrt(length(subset(data,age=="DIV21"&group=="hpc")$br)),
  sd(subset(data,age=="DIV28"&group=="ctx")$br)/sqrt(length(subset(data,age=="DIV28"&group=="ctx")$br)),
  sd(subset(data,age=="DIV28"&group=="hpc")$br)/sqrt(length(subset(data,age=="DIV28"&group=="hpc")$br)))


meanFRs<-rbind(
  mean(subset(data,age=="DIV07"&group=="template")$br),
  mean(subset(data,age=="DIV07"&group=="threshold")$br),
  mean(subset(data,age=="DIV14"&group=="template")$br),
  mean(subset(data,age=="DIV14"&group=="threshold")$br),
  mean(subset(data,age=="DIV21"&group=="template")$br),
  mean(subset(data,age=="DIV21"&group=="threshold")$br),
  mean(subset(data,age=="DIV28"&group=="template")$br),
  mean(subset(data,age=="DIV28"&group=="threshold")$br))

SEMFRs<-rbind(
  sd(subset(data,age=="DIV07"&group=="template")$br)/sqrt(length(subset(data,age=="DIV07"&group=="template")$br)),
  sd(subset(data,age=="DIV07"&group=="threshold")$br)/sqrt(length(subset(data,age=="DIV07"&group=="threshold")$br)),
  sd(subset(data,age=="DIV14"&group=="template")$br)/sqrt(length(subset(data,age=="DIV14"&group=="template")$br)),
  sd(subset(data,age=="DIV14"&group=="threshold")$br)/sqrt(length(subset(data,age=="DIV14"&group=="threshold")$br)),
  sd(subset(data,age=="DIV21"&group=="template")$br)/sqrt(length(subset(data,age=="DIV21"&group=="template")$br)),
  sd(subset(data,age=="DIV21"&group=="threshold")$br)/sqrt(length(subset(data,age=="DIV21"&group=="threshold")$br)),
  sd(subset(data,age=="DIV28"&group=="template")$br)/sqrt(length(subset(data,age=="DIV28"&group=="template")$br)),
  sd(subset(data,age=="DIV28"&group=="threshold")$br)/sqrt(length(subset(data,age=="DIV28"&group=="threshold")$br)))


msems<-as.data.frame(cbind(meanFRs,SEMFRs))
toplot1<-cbind(toplot,msems)
names(toplot1)[c(3, 4)] <- c("mean","SE")
toplot1



library(ggplot2)
#ref_df <- as.data.frame(summary(ref))
theme_update(text = element_text(size=30))
pd <- position_dodge(0.1)
g4 <- ggplot(toplot1, aes(x=age, y=mean,group=group,colour=group))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.15,position=pd,size=1) +
  geom_line(position=pd,size=1)+
  geom_point(position=pd,size = 3)+
  theme_classic(base_size = 22,base_line_size = 1)+
  ylab('Firing rate (Hz)') +
  xlab('Age') 

  #ylim(c(-0.01, 0.4))

print(g4)
g4 + scale_color_manual(values=c("#FF0000", "#0000FF"))+
     scale_linetype_manual(values=c("dashed", "solid"))

#change colour of lines and type of line # RGB in hexidecimal ###### 00 to FF




#note this is plotting the least squares mean which is the adjustment to the means made in the anova
(ref2 <- lsmeans(fit_all,~age|group))
comps<-contrast(ref2,method="pairwise",adjust="bonferroni") #pairwise comparisons, correcting using bonferroni method
summary(comps) #i.e. post hoc tests using Tukey method with bonferoni correction

(ref3 <- lsmeans(fit_all,~group|age))
comps<-contrast(ref3,method="pairwise",adjust="bonferroni") #pairwise comparisons, correcting using bonferroni method
summary(comps) #i.e. post hoc tests using Tukey method with bonferoni correction
toplot1
aggregate(br~group+age, data = data,summary)











############main effect of age
(fit_all <- aov_ez("id","br",afexdata, within="age"))
names(fit_all)
summary(fit_all) #mauchlys test not sig; we can assume sphericity

(ref4 <- lsmeans(fit_all,~age))
comps<-contrast(ref4,method="pairwise",adjust="bonferroni") #pairwise comparisons, correcting using bonferroni method
summary(comps) #i.e. post hoc tests with bonferoni correction


toplot2<-as.data.frame(aggregate(br~age, data = data,summary))
toplot2<-toplot2[c(1)]
meanFRs<-rbind(
  mean(subset(data,age=="DIV07")$br),
  mean(subset(data,age=="DIV14")$br),
  mean(subset(data,age=="DIV21")$br),
  mean(subset(data,age=="DIV28")$br))

SEMFRs<-rbind(
  sd(subset(data,age=="DIV07")$br)/sqrt(length(subset(data,age=="DIV07")$br)),
  sd(subset(data,age=="DIV14")$br)/sqrt(length(subset(data,age=="DIV14")$br)),
  sd(subset(data,age=="DIV21")$br)/sqrt(length(subset(data,age=="DIV21")$br)),
  sd(subset(data,age=="DIV28")$br)/sqrt(length(subset(data,age=="DIV28")$br)))

msems<-as.data.frame(cbind(meanFRs,SEMFRs))
toplot3<-cbind(toplot2,msems)
names(toplot3)[c(2, 3)] <- c("mean","SE")
toplot3

theme_update(text = element_text(size=30))

# Basic line plot with points
#ggplot(data=toplot3, aes(x=age, y=mean, group=1)) +
#  geom_line()+
#  geom_point()+
#  geom_errorbar() +
  
# Change the line type
#ggplot(data=df, aes(x=dose, y=len, group=1)) +
#  geom_line(linetype = "dashed")+
#  geom_point()
# Change the color
#ggplot(data=df, aes(x=dose, y=len, group=1)) +
#  geom_line(color="red")+
#  geom_point()


pd <- position_dodge(0.1)
g4 <- ggplot(toplot3, aes(x=age, y=mean,group=1))+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=0.15,position=pd,size=1,color="orchid4") +
  geom_line(position=pd,size=1,color="orchid4")+
  geom_point(position=pd,size = 3,color="orchid4")+
  theme_classic(base_size = 22,base_line_size = 1)+
  ylab('??_norm (k)') +
  xlab('Age')
  
#ylim(c(0.0, 2.22))

print(g4)






####################
###################
#################
#non parametric tests
#WT, comparison between ages
#question is; does CC increase over time (is DIV14>DIV07)
wtdata=subset(data, data$genotype=='WT')
head(wtdata[c('br','age')])
uns.wtdata<-unstack(wtdata[c('br','age')])
t1<-wilcox.test(uns.wtdata$DIV07, uns.wtdata$DIV14,
            alternative=  "less", paired=T)
t2<-wilcox.test(uns.wtdata$DIV14, uns.wtdata$DIV21,
                alternative=  "less", paired=T)
t3<-wilcox.test(uns.wtdata$DIV21, uns.wtdata$DIV28,
                alternative=  "less", paired=T)

#DIV28 relative to baseline
p.adjust(c(t1$p.value,t2$p.value,t3$p.value),method = 'fdr')

t4<-wilcox.test(uns.wtdata$DIV14, uns.wtdata$DIV28,
                alternative=  "less", paired=T)
t5<-wilcox.test(uns.wtdata$DIV07, uns.wtdata$DIV28,
                alternative=  "less", paired=T)

p.adjust(c(t4$p.value,t5$p.value),method = 'fdr')

#comp between agees in Het culture
hedata=subset(data, data$group=='HE')
head(hedata[c('br','age')])
uns.hedata<-unstack(hedata[c('br','age')])
t11<-wilcox.test(uns.hedata$DIV07, uns.hedata$DIV14,
                alternative=  "less", paired=T)
t22<-wilcox.test(uns.hedata$DIV14, uns.hedata$DIV21,
                alternative=  "less", paired=T)
t33<-wilcox.test(uns.hedata$DIV21, uns.hedata$DIV28,
                alternative=  "less", paired=T)
t.test(uns.hedata$DIV21, uns.hedata$DIV28,
       alternative=  "greater", paired=T)

#DIV28 relative to baseline
c(t11$p.value,t22$p.value,t33$p.value)
p.adjust(c(t11$p.value,t22$p.value,t33$p.value),method = 'fdr')

t44<-wilcox.test(uns.hedata$DIV14, uns.hedata$DIV28,
                alternative=  "less", paired=T)
t55<-wilcox.test(uns.hedata$DIV07, uns.hedata$DIV28,
                alternative=  "less", paired=T)

p.adjust(c(t44$p.value,t55$p.value),method = 'fdr')

#across all genotypes

head(data[c('br','age')])
uns.data<-unstack(data[c('br','age')])
t111<-wilcox.test(uns.data$DIV07, uns.data$DIV14,
                alternative=  "less", paired=T)
t222<-wilcox.test(uns.data$DIV14, uns.data$DIV21,
                alternative=  "less", paired=T)
t333<-wilcox.test(uns.data$DIV21, uns.data$DIV28,
                alternative=  "less", paired=T)

#CC increased significantly from DIV21 to 28; result is basically the same as the ANOVA
p.adjust(c(t111$p.value,t222$p.value,t333$p.value),method = 'fdr')



data<-read.csv(file.choose(),header = T)
head(data)
###group comparisons
d07data=subset(data, data$age=='DIV07')
d14data=subset(data, data$age=='DIV14')
d21data=subset(data, data$age=='DIV21')
d28data=subset(data, data$age=='DIV28')

par(mfrow=c(2,3))
par(cex=0.7,cex.lab=1.2,cex.axis=1.1)  #font size
plot(1,1,type='n',yaxt='n',xaxt='n',ylab=' ',xlab=' ')
plot(1,1,type='n',yaxt='n',xaxt='n',ylab=' ',xlab=' ')
plot(1,1,type='n',yaxt='n',xaxt='n',ylab=' ',xlab=' ')
boxplot(br~group, data=d14data,main='DIV14',ylab='Burst rate (/min)')
boxplot(br~group, data=d21data,main='DIV21')
boxplot(br~group, data=d28data,main='DIV28')

aggregate(br~group,data = d07data,summary)
aggregate(br~group,data = d14data,summary)
aggregate(br~group,data = d21data,summary)
aggregate(br~group,data = d28data,summary)

head(d14data[c('br','group')])
#t1<-wilcox.test(br~group, data = d14data,alternative=  "two.sided")
t1<-wilcox.test(br~group, data = d14data,alternative=  "greater")
t2<-wilcox.test(br~group, data = d21data,alternative=  "greater")
t3<-wilcox.test(br~group, data = d28data,alternative=  "greater")
c(t1$p.value,t2$p.value,t3$p.value)
p.adjust(c(t1$p.value,t2$p.value,t3$p.value),method = 'fdr')
t1 
t2
t3
t4<-wilcox.test(br~group, data = data,alternative=  "greater")
t4
aggregate(br~group+age,data = data,  summary)

bartlett.test(br~group,data = d14data)
bartlett.test(br~group,data = d21data)
bartlett.test(br~group,data = d28data)

leveneTest(br~group,data = d14data)
leveneTest(br~group,data = d21data)
leveneTest(br~group,data = d28data)

shapiro.test(subset(d14data,d14data$group=='WT')$br) 
hist(subset(d14data,d14data$group=='WT')$br) 
shapiro.test(subset(d21data,d21data$group=='WT')$br) 
hist(subset(d21data,d21data$group=='WT')$br) 
shapiro.test(subset(d28data,d28data$group=='WT')$br)
hist(subset(d28data,d28data$group=='WT')$br) 

shapiro.test(subset(d14data,d14data$group=='HE')$br) 
hist(subset(d14data,d14data$group=='HE')$br) 
shapiro.test(subset(d21data,d21data$group=='HE')$br)
hist(subset(d21data,d21data$group=='HE')$br) 
shapiro.test(subset(d28data,d28data$group=='HE')$br) 
hist(subset(d28data,d28data$group=='HE')$br)


#t tests

t1<-t.test(br~group, data=d14data, alternative='greater',var.equal=TRUE) #False = welch's
t2<-t.test(br~group, data=d21data, alternative='greater',var.equal=TRUE) #False = welch's
t3<-t.test(br~group, data=d28data, alternative='greater',var.equal=TRUE) #False = welch's
c(t1$p.value,t2$p.value,t3$p.value)
p.adjust(c(t1$p.value,t2$p.value,t3$p.value),method = 'fdr')
t.test(br~group, data=data, alternative='greater',var.equal=TRUE) #False = welch's
# non para results confirmed anova - main effect of age, no effects of genotype
#interpret t test differences between genotypes with caution 

data<-read.csv(file.choose(),header=T) #load data3.1.1.CTXvsHPCmSpikes.csv
  head(data)

p1<-subset(data,data$id=='MPT190403_2B')
p2<-subset(data,data$id=='MPT190403_6B')
p3<-subset(data,data$id=='MPT190403_6C')
p4<-subset(data,data$id=='MPT190515_2B')
p5<-subset(data,data$id=='MPT190515_3B')
p6<-subset(data,data$id=='MPT190515_4B')
p7<-subset(data,data$id=='MPT190515_4C')
p8<-subset(data,data$id=='MPT190528_1A')
p9<-subset(data,data$id=='MPT190528_1B')
p10<-subset(data,data$id=='MPT190528_1C')
p11<-subset(data,data$id=='MPT190528_2B')
p12<-subset(data,data$id=='MPT190528_4A')

par(cex=1,cex.lab=1.5,cex.axis=1.3)  #font size
#interaction.plot(data$age,data$group,data$br,col = c("#0000FF","#FF0000"),lty = 1,lwd=2,
#                 ylab = 'Firing rate (Hz)', xlab='Age',bty = 'l',legend = 'F',ylim = c(0,max(data$br)))
interaction.plot(data$age,data$group,data$br,col = c("#0000FF","#FF0000"),lty = 1,lwd=2,
                 ylab = 'Firing rate (Hz)', xlab='Age',bty = 'l',legend = 'F',
                 ylim = c(0,max(data$br)))


g1<-subset(data,data$group=='WT')
g2<-subset(data,data$group=='HE')

#g1<-subset(data,data$group=='ctx')
#g2<-subset(data,data$group=='hpc')

g1<-subset(data,data$group=='template')
g2<-subset(data,data$group=='threshold')

m1<-c(mean(subset(g1,g1$age=='DIV07')$br),mean(subset(g1,g1$age=='DIV14')$br),mean(subset(g1,g1$age=='DIV21')$br),mean(subset(g1,g1$age=='DIV28')$br))
m2<-c(mean(subset(g2,g2$age=='DIV07')$br),mean(subset(g2,g2$age=='DIV14')$br),mean(subset(g2,g2$age=='DIV21')$br),mean(subset(g2,g2$age=='DIV28')$br))

s1<-   c(sd(subset(g1,g1$age=='DIV07')$br/sqrt(length(subset(g1,g1$age=='DIV07')$br))),
         sd(subset(g1,g1$age=='DIV14')$br/sqrt(length(subset(g1,g1$age=='DIV14')$br))),
         sd(subset(g1,g1$age=='DIV21')$br/sqrt(length(subset(g1,g1$age=='DIV21')$br))),
         sd(subset(g1,g1$age=='DIV28')$br/sqrt(length(subset(g1,g1$age=='DIV28')$br))))

s2<-   c(sd(subset(g2,g2$age=='DIV07')$br/sqrt(length(subset(g2,g2$age=='DIV07')$br))),
         sd(subset(g2,g2$age=='DIV14')$br/sqrt(length(subset(g2,g2$age=='DIV14')$br))),
         sd(subset(g2,g2$age=='DIV21')$br/sqrt(length(subset(g2,g2$age=='DIV21')$br))),
         sd(subset(g2,g2$age=='DIV28')$br/sqrt(length(subset(g2,g2$age=='DIV28')$br))))

#points(jitter(seq(1.9,2.1,length.out = length(p1))),c2bt0$Hurtz,pch=4,col='red',type = 1)
interaction.plot(data$age,data$group,data$br,col = c("white","white"),lty = 1,lwd=2,
                 ylab = 'Firing rate (Hz)', xlab='Age',bty = 'l',legend = 'F',
                 ylim = c(0,max(data$br)),cex.axis=1.3,cex.lab=1.3,cex=1,type = 'b',pch=16)
#try parameter mgp=c(1.2,-1,0) to move axis labels closer

points(p1$age,p1$br,pch=2,col='#AAAAFF',type = 'l',lwd=1.5,lty=2,cex=1.5) #190403;het
points(p1$age,p1$br,pch=2,col='#AAAAFF',type = 'p',lwd=2,lty=2,cex=1.5) #190403;het
jitvec<-c(1,2,3,4)
points(jitvec+0.03,p2$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec+0.03,p2$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het
points(jitvec+0.06,p3$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec+0.06,p3$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het

points(jitvec-0.06,p4$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec-0.06,p4$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het
points(jitvec-0.03,p5$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec-0.03,p5$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het
points(jitvec-0.09,p6$br,pch=2,col='#AAAAFF',type = 'l',lwd=1.5,lty=2,cex=1.5) #190403;het
points(jitvec-0.09,p6$br,pch=2,col='#AAAAFF',type = 'p',lwd=2,lty=2,cex=1.5) #190403;het

points(jitvec+0.09,p7$br,pch=2,col='#AAAAFF',type = 'l',lwd=1.5,lty=2,cex=1.5) #190403;het
points(jitvec+0.09,p7$br,pch=2,col='#AAAAFF',type = 'p',lwd=2,lty=2,cex=1.5) #190403;het
points(jitvec+0.12,p8$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec+0.12,p8$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het
points(jitvec-0.12,p9$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec-0.12,p9$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het

points(jitvec+0.15,p10$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec+0.15,p10$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het
points(jitvec-0.15,p11$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec-0.15,p11$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het
points(jitvec+0.18,p12$br,pch=1,col='#FFAAAA',type = 'l',lwd=1.5,lty=2,cex=1) #190403;het
points(jitvec+0.18,p12$br,pch=1,col='#FFAAAA',type = 'p',lwd=2,lty=2,cex=1) #190403;het

#means
points(c(1,2,3,4),m1,pch=16,col='#FF0000',type='b', lty = 1,lwd=2)
points(c(1.1,2.1,3.1,4.1),m2,pch=17,col='#0000FF',type='b', lty = 1,lwd=2) 

arrows(c(1,2,3,4), m1-s1, c(1,2,3,4), m1+s1, length=0.05, angle=90, code=3,col = '#FF0000',lwd=2)
arrows(c(1.1,2.1,3.1,4.1), m2-s2, c(1.1,2.1,3.1,4.1), m2+s2, length=0.05, angle=90, code=3,col = '#0000FF',lwd=2)


####density  just sems
interaction.plot(data$age,data$group,data$br,col = c("white","white"),lty = 1,lwd=2,
                 ylab = 'N. active channels', xlab='Age',bty = 'l',legend = 'F', 
                 ylim = c(0,35),cex.axis=1.2,cex.lab=1.3,cex=1,type = 'o',pch=16)

#
points(c(1,2,3,4),m1,pch=16,col='#FF0000',type='o', lty = 2,lwd=2)
points(c(1.1,2.1,3.1,4.1),m2,pch=16,col='#0000FF',type='o', lty = 1,lwd=2) 
arrows(c(1,2,3,4), m1-s1, c(1,2,3,4), m1+s1, length=0.05, angle=90, code=3,col = '#FF0000',lwd=2)
arrows(c(1.1,2.1,3.1,4.1), m2-s2, c(1.1,2.1,3.1,4.1), m2+s2, length=0.05, angle=90, code=3,col = '#0000FF',lwd=2)

###########
#all groups
data<-read.csv(file.choose(),header=T) #load data3.1.1.CTXvsHPCmSpikes.csv
head(data)

m3<-c(mean(subset(data,data$age=='DIV07')$br),mean(subset(data,data$age=='DIV14')$br),
      mean(subset(data,data$age=='DIV21')$br),mean(subset(data,data$age=='DIV28')$br))
s3<-   c(sd(subset(data,data$age=='DIV07')$br/sqrt(length(subset(data,data$age=='DIV07')$br))),
         sd(subset(data,data$age=='DIV14')$br/sqrt(length(subset(data,data$age=='DIV14')$br))),
         sd(subset(data,data$age=='DIV21')$br/sqrt(length(subset(data,data$age=='DIV21')$br))),
         sd(subset(data,data$age=='DIV28')$br/sqrt(length(subset(data,data$age=='DIV28')$br))))
interaction.plot(data$age,data$group,data$br,col = c("white","white"),lty = 1,lwd=2,
                 ylab = 'Firing rate (Hz)', xlab='Age',bty = 'l',legend = 'F',
                 ylim = c(0,max(data$br)),cex.axis=1.3,cex.lab=1.3,cex=1,type = 'b',pch=16)

points(c(1,2,3,4),m3,pch=16,col='orchid4',type='b', lty = 1,lwd=2)

arrows(c(1,2,3,4), m3-s3, c(1,2,3,4), m3+s3, length=0.05, angle=90, code=3,col = 'orchid4',lwd=2)


##############
toplot1
############
##############
#mean and sd and sem for genotype main effect


toplot2<-as.data.frame(aggregate(br~genotype, data = data,summary))
toplot2<-toplot2[c(1)]

meanFRs<-rbind(
  mean(subset(data,group=="HE")$br),
  mean(subset(data,group=="WT")$br))
 
SEMFRs<-rbind(
  sd(subset(data,group=="HE")$br)/sqrt(length(subset(data,group=="HE")$br)),
  sd(subset(data,group=="WT")$br)/sqrt(length(subset(data,group=="WT")$br)))

  SDFRs<-rbind(
    sd(subset(data,group=="HE")$br),
    sd(subset(data,group=="WT")$br))
  
  msems<-as.data.frame(cbind(meanFRs,SEMFRs,SDFRs))
  toplot3<-cbind(toplot2,msems)
  names(toplot3)[c(2, 3, 4)] <- c("mean","SE","SD")
  
  
  
##### sample size calculation based on DIV28
  muA=meanFRs[1]
  muB=meanFRs[2]
  kappa=3/9
  sdA=SDFRs[1]
  sdB=SDFRs[2]
  alpha=0.05
  beta=0.20
  (nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
  ceiling(nA) # 85
  z=(muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)
  (Power=pnorm(z-qnorm(1-alpha)))
  ## Note: Rosner example on p.303 is for 2-sided test.
  ## These formulas give the numbers in that example
  ## after dividing alpha by 2, to get 2-sided alpha.
  
 # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-1-Sided
  
  #manual calculation of mean and sem for each group
  toplot<-as.data.frame(aggregate(br~group+age, data = data,summary))
  toplot<-toplot[c(1,2)]
  meanFRs<-rbind(
    mean(subset(data,age=="DIV28"&group=="HE")$br),
    mean(subset(data,age=="DIV28"&group=="WT")$br))
  
  SEMFRs<-rbind(
    sd(subset(data,age=="DIV28"&group=="HE")$br)/length(subset(data,age=="DIV28"&group=="HE")$br),
    sd(subset(data,age=="DIV28"&group=="WT")$br)/length(subset(data,age=="DIV28"&group=="WT")$br))
  
  SDFRs<-rbind(
    sd(subset(data,age=="DIV28"&group=="HE")$br),
    sd(subset(data,age=="DIV28"&group=="WT")$br))
  
  msems<-as.data.frame(cbind(meanFRs,SEMFRs,SDFRs))
  toplot3<-cbind(toplot2,msems)
  names(toplot3)[c(2, 3, 4)] <- c("mean","SE","SD")
  
  muA=meanFRs[2]
  muB=meanFRs[1]
  kappa=3/9
  sdA=SDFRs[2]
  sdB=SDFRs[1]
  alpha=0.05
  beta=0.20
  (nA=(sdA^2+sdB^2/kappa)*((qnorm(1-alpha)+qnorm(1-beta))/(muA-muB))^2)
  ceiling(nA) # 85
  nB <- kappa*nA
  ceiling(nB)
  z=(muA-muB)/sqrt(sdA^2+sdB^2/kappa)*sqrt(nA)
  (Power=pnorm(z-qnorm(1-alpha)))
  ## Note: Rosner example on p.303 is for 2-sided test.
  ## These formulas give the numbers in that example
  ## after dividing alpha by 2, to get 2-sided alpha.
  
  
  
  
  
  
  
  
  
  
  ##################
  #correlating variables to see if variability in basic firing can explain them
  data<-read.csv(file.choose(),header=T)
  #cordata<-data[c(1,2,3)]
  #cordata<-subset(data,data$age=="DIV28")
  cordata<-subset(data,data$age=="DIV28"|data$age=="DIV21")
  #cordata<-cordata[c(1,2,3)]
  cordata<-cordata[c(1,2)]
  head(cordata)
  pairs(cordata,lower.panel = NULL,cex=2,pch=1,lwd=2,cex.labels = 3,cex.axis=2)

  #use spearman due to heteoscedasiticity
  #all recordings, n = 48
  spear<-cor(cordata,method = "spearman")
  spear
 spear['clustering','firing']
 spt1<-cor.test(cordata$clustering,cordata$firing,method='spearman')
 spt2<-cor.test(cordata$clustering,cordata$channels,method='spearman')
 spt3<-cor.test(cordata$channels,cordata$firing,method='spearman')
 
  c(spt1$p.value,spt2$p.value,spt3$p.value)
  
  cor(cordata,method = "pearson")
  
  lm1<-lm(clustering~firing, data = cordata)
anova(lm1)  
plot(clustering~firing,data = cordata,cex=1.5,pch=1,lwd=2,cex.lab=1.5,cex.axis=1.5,bty='l')
abline(lm1,col='red')

lm2<-lm(clustering~channels, data = cordata)
anova(lm2)
plot(clustering~channels,data = cordata,cex=1.5,pch=1,lwd=2,cex.lab=1.5,cex.axis=1.5,bty='l')
abline(lm2,col='red')

lm3<-lm(channels~firing, data = cordata)
anova(lm3)
plot(channels~firing,data = cordata,cex=1.5,pch=1,lwd=2,cex.lab=1.5,cex.axis=1.5,bty='l')
abline(lm3,col='red')

##clust and density and sttc
cordata<-subset(data,data$age=="DIV28"|data$age=="DIV21"|data$age=="DIV14")
cordata<-subset(data,data$age=="DIV28"|data$age=="DIV21")
cordata<-subset(data,data$age=="DIV28")
cordata<-subset(data,data$age=="DIV21")
cordata<-subset(data,data$age=="DIV14")
#cordata<-data[c(1,2)]
cordata<-cordata[c(1,2,3)]
pairs(cordata,lower.panel = NULL,cex=2,pch=1,lwd=2,cex.labels = 3,cex.axis=2)

#use spearman due to heteoscedasiticity
#all recordings, n = 48
spear<-cor(cordata,method = "spearman")
spear
pear<-cor(cordata,method = "pearson")
pear
spear['Clustering','Density']
spt1<-cor.test(cordata$Density,cordata$Clustering,method='spearman')
spt2<-cor.test(cordata$Density,cordata$Clustering,method='pearson')
c(spt1$p.value,spt2$p.value)

spt1<-cor.test(cordata$Density,cordata$Clustering,method='spearman')
spt2<-cor.test(cordata$Density,cordata$STTC,method='spearman')
spt3<-cor.test(cordata$Clustering,cordata$STTC,method='spearman')
c(spt1$p.value,spt2$p.value,spt3$p.value)


cor(cordata,method = "pearson")

lm1<-lm(clustering~firing, data = cordata)
anova(lm1)  
plot(clustering~firing,data = cordata,cex=1.5,pch=1,lwd=2,cex.lab=1.5,cex.axis=1.5,bty='l')
abline(lm1,col='red')

lm2<-lm(clustering~channels, data = cordata)
anova(lm2)
plot(clustering~channels,data = cordata,cex=1.5,pch=1,lwd=2,cex.lab=1.5,cex.axis=1.5,bty='l')
abline(lm2,col='red')

lm3<-lm(channels~firing, data = cordata)
anova(lm3)
plot(channels~firing,data = cordata,cex=1.5,pch=1,lwd=2,cex.lab=1.5,cex.axis=1.5,bty='l')
abline(lm3,col='red')
