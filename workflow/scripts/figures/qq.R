library(dplyr)
library(data.table)
library(ggplot2)




########### Loading data needed ###########
dat = fread(snakemake@input[[1]])




########### Formating data ###########
dat = dat[!duplicated(dat$ID), ] 

dat= dat %>% arrange(desc(LOG10P)) # sample 

dat= dat %>% ungroup() %>% arrange(desc(LOG10P)) %>% mutate(exp1= -log10(1:length(LOG10P)/length(LOG10P))) # expected, under the null hypothesis are uniformly distributed between 0 and 1




########### Calculate lambda ###########
chisq= qchisq(1-10**-dat$LOG10P, 1)
lambda_gc= median(chisq)/qchisq(0.5,1)
cat("Labda: ", lambda_gc)




########### QQ plot ###########
p1= ggplot(dat, aes(exp1, LOG10P)) +
	  geom_point(size= 0.4, color= "#E69F00") +
	  geom_abline(intercept = 0, slope = 1, alpha = .5) +
	  labs(colour="") +
	  theme_bw( ) +
	  xlab('Expected (-log10(p-value))') +
	  ylab('Observed (-log10(p-value))') +
	  theme(legend.position= 'bottom') +
	  geom_text(aes(6, 0), label= paste("lambda", "==", round(lambda_gc, 2)), size= 11/.pt, parse= T)




########### Save figure ###########
ggsave(snakemake@output[[1]],p1)
