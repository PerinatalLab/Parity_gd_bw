## interaction effect plot

## model: GD ~ INT + PRS + parity + prs*parity
## if parity=0, GD ~ INT + PRS, so plot a line with intercept=INT, slope=effect of PRS
## if parity=1, GD ~ INT + PRS + parity + prs*parity, so  plot a line with intercept = INT + effect of parity, slope = effect of PRS + effect of interaction

library(data.table)
library(dplyr)
library(ggplot2)


phe='results/replication/gd_Mother_parityall/phe.tsv'
cov='results/replication/gd_Mother_parityall/cov.tsv'
prs="results/PGS/gd_Mother_effectSNPnoInteraction_parityall/PGS.txt"
genome="Maternal"
traits="GD"
outfile="results/PGS/gd_Mother_effectSNPnoInteraction_parityall/interaction_plot.png"
outfile2="results/PGS/gd_Mother_effectSNPnoInteraction_parityall/interaction_plot2.png"

  phe <- fread(phe)
  cov <- fread(cov)
  prs <- fread(prs)
  df <- prs %>%
    mutate(PRS_z=scale(as.numeric(PRS), scale = F)) %>%
    inner_join(phe,by=c("FID","IID"))%>%
    inner_join(cov,by=c("FID","IID")) %>%
    filter(sex_assigned_at_birth==1 | sex_assigned_at_birth==2)
 

  s <- coef(summary(mod))
  int0 <- s["(Intercept)","Estimate"]
  slope0 <- s["PRS_z","Estimate"]
  
  int1 <- s["(Intercept)","Estimate"] + s["PARITET1","Estimate"]
  slope1 <- s["PRS_z","Estimate"] + s["PARITET1:PRS_z","Estimate"]
  
  interactionpvalue = formatC(s["PARITET1:PRS_z","Pr(>|t|)"], format = "e", digits = 2)
  
  plotdata = c()
  plotdata = as.data.frame(rbind(plotdata, c(int0,int1,slope0,slope1,interactionpvalue,traits,genome)))
  names(plotdata) <- c("int0","int1","slope0","slope1","interactionpvalue","traits","genome")

  p <- plotdata %>% ggplot() +
    coord_cartesian(
      xlim = c(round(min(df$PRS_z),0), round(max(df$PRS_z),0)),
      ylim = c(265, 290)
    )+
    geom_hline(yintercept = 0, alpha = .5, color = "black") 
    xlab(paste("Standardized polygenic genetic score", paste(genome, "genome",sep =" "), sep="\n")) +
    ylab(paste(traits)) +
    geom_abline(
      aes(intercept = as.numeric(plotdata$int0),
          slope    = as.numeric(plotdata$slope0),
          color    = "= 0"),
      linewidth = 1.5
    ) +
    geom_abline(
      aes(intercept = as.numeric(plotdata$int1),
          slope    = as.numeric(plotdata$slope1),
          color    = "> 0"),
      linewidth = 1.5
    )+
    annotate("text",x = max(df$PRS_z, na.rm = TRUE),y = 290,hjust = 1,label = paste0("PGS × Parity p = ", interactionpvalue), size = 3) +
    theme_bw()+
    theme(
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=11, color="black"),
      axis.title = element_text(size=11, color="black"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(size = 11, color="black"),
      panel.border = element_rect(color ="black", fill=NA, linewidth =0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    scale_color_manual(
      name   = "Parity",
      values = c("= 0" = "#008837", "> 0" = "#7b3294")
    )
  
  
  ggsave(filename=outfile,plot=p)

 

### add confidence interval
  
  df$PARITET <- as.factor(df$PARITET)
  df$PRS_z <- as.numeric(df$PRS_z)
  
  fit <- lm(phe ~ PRS_z*PARITET,data=df)
  
  
  newdata <- expand.grid(
    PRS_z = seq(min(fit$model$PRS_z), max(fit$model$PRS_z), length.out = 100),
    PARITET = levels(fit$model$PARITET)
  )
  
  
  pred <- predict(
    fit,
    newdata = newdata,
    interval = "confidence"
  )
  
  pred_df <- cbind(newdata, as.data.frame(pred))
  
  #pred_df <- pred_df[order(pred_df$PARITET, pred_df$PRS_z), ]
  pred_df$PARITET <- factor(
    pred_df$PARITET,
    levels = c("0", "1"),
    labels = c("= 0", "> 0")
  )
  pp <- ggplot(pred_df, aes(x = PRS_z, y = fit, color = PARITET, fill = PARITET)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
    scale_colour_manual(name = "PARITET", values=c("= 0"="#008837","> 0"="#7b3294"))+
    scale_fill_manual(name = "PARITET",values = c("= 0" = "#008837", "> 0" = "#7b3294"))+
    xlab(paste("Standardized polygenic genetic score", paste(genome, "genome",sep =" "), sep="\n")) +
    ylab("Gestational duration") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=11, color="black"),
      axis.title = element_text(size=11, color="black"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(size = 11, color="black"),
      panel.border = element_rect(color ="black", fill=NA, size=0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA))+
    annotate("text",x = max(df$PRS_z, na.rm = TRUE),y = 290,hjust = 1,label = paste0("PGS × Parity p = ", interactionpvalue), size = 3) 
  ggsave(outfile2,plot=pp)
  
