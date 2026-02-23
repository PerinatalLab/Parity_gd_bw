library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)

parsed_opts <- list(
  make_option("--PRS",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--phe",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--cov",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--quantile_plot",
              action = "store",
              default = NA,
              type = 'character') 
)


if (interactive()) {
  opt <- list(PRS='results/PGS/gd_Mother_parity0/PGS.txt',
              phe='results/replication/gd_Mother_parity0/phe.tsv',
              cov='results/replication/gd_Mother_parity0/cov.tsv',
              quantile_plot="results/PGS/gd_Mother_parity0/qplot.png"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}
#incremental R² and β per SD

phe <- fread(opt$phe)
cov <- fread(opt$cov) %>%
  inner_join(phe,by=c("FID","IID"))

df <- fread(opt$PRS) %>%
  inner_join(cov,by=c("FID","IID"))

df$PRS_z <- scale(df$PRS)



# baseline model with covariates
mod0 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df)
mod1 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PRS_z, data=df)

cat("Effect size per SD PRS:",summary(mod1)$coef["PRS_z",])
r2_0 <- summary(mod0)$r.squared
r2_1 <- summary(mod1)$r.squared
delta_r2 <- r2_1 - r2_0
cat("Incremental R2:", delta_r2, "\n")


##=====================PRS quantile
get_quantile <- function(x, num.quant, quant.ref){
  quant <- as.numeric(cut(x,
                          breaks = unique(quantile(
                            x, probs = seq(0, 1, 1 / num.quant)
                          )),
                          include.lowest = T))
  if(is.null(quant.ref)){
    quant.ref <- ceiling(num.quant / 2)       #define reference factor
  }
  quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
  return(quant)
}




lr <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df)
lr_res <- resid(lr)
lr_res_zscore <- (lr_res-mean(lr_res))/sd(lr_res)
df$lr_res_zscore <- lr_res_zscore

df$quantile_group <- get_quantile(df$PRS_z,num.quant=10, quant.ref=1)
q975 <- qnorm(0.975)

lr <- summary(lm(lr_res_zscore ~ -1 + quantile_group, data=df))
lr_m <- lr$coefficients[,1]
lr_sd <- lr$coefficients[,2]
lr_counts <- table(df[,"quantile_group"])
clr <- cbind(c(1:10),lr_m,lr_m+q975*lr_sd,lr_m-q975*lr_sd,lr_counts)
colnames(clr) <- c("Quantile","Coef","CI.U","CI.L","N")

myplot <- ggplot(clr,aes(x=Quantile,y=Coef,ymin=CI.L,ymax=CI.U))+
  geom_point(position=position_dodge(0.5))+
  geom_errorbar(width=0.1,position=position_dodge(0.5)) +
  xlab("PGS Quantiles")+
  ylab("Effect (95% CI)")+
  scale_x_continuous(breaks=seq(1,10,1)) +
  theme_bw() +
  theme(legend.title = element_blank())

ggsave(opt$quantile_plot,myplot,height=4,width=4)

