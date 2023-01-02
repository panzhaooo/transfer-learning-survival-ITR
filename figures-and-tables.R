library(ggplot2)
setwd("")

# (semi)parametric models 
load("survITR1.RData")
load("survITR2.RData")
load("survITR3.RData")
load("survITR4.RData")

n.MonteCarlo <- 350

# double robustness
data1 <- data.frame(value = c(result1$Vd$Vd_naive, # Vd
                              result1$Vd$Vd_ipsw,
                              result1$Vd$Vd_cwipw,
                              result1$Vd$Vd_cwor,
                              result1$Vd$Vd_ort,
                              result1$Vd$Vd_acw,
                              result1$pcd$pcd_naive, # PCD
                              result1$pcd$pcd_ipsw,
                              result1$pcd$pcd_cwipw,
                              result1$pcd$pcd_cwor,
                              result1$pcd$pcd_ort,
                              result1$pcd$pcd_acw,
                              result1$Vhat$Vhat_naive, # Vhat
                              result1$Vhat$Vhat_ipsw,
                              result1$Vhat$Vhat_cwipw,
                              result1$Vhat$Vhat_cwor,
                              result1$Vhat$Vhat_ort,
                              result1$Vhat$Vhat_acw),
                    type = rep(c("True Value", "PCD", "Estimated Value"), each = n.MonteCarlo * 6),
                    Estimators = rep(rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo), times = 3),
                    sce = "O:T / S:T, A:T, C:T")
data2 <- data.frame(value = c(result2$Vd$Vd_naive, # Vd
                              result2$Vd$Vd_ipsw,
                              result2$Vd$Vd_cwipw,
                              result2$Vd$Vd_cwor,
                              result2$Vd$Vd_ort,
                              result2$Vd$Vd_acw,
                              result2$pcd$pcd_naive, # PCD
                              result2$pcd$pcd_ipsw,
                              result2$pcd$pcd_cwipw,
                              result2$pcd$pcd_cwor,
                              result2$pcd$pcd_ort,
                              result2$pcd$pcd_acw,
                              result2$Vhat$Vhat_naive, # Vhat
                              result2$Vhat$Vhat_ipsw,
                              result2$Vhat$Vhat_cwipw,
                              result2$Vhat$Vhat_cwor,
                              result2$Vhat$Vhat_ort,
                              result2$Vhat$Vhat_acw),
                    type = rep(c("True Value", "PCD", "Estimated Value"), each = n.MonteCarlo * 6),
                    Estimators = rep(rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo), times = 3),
                    sce = "O:W / S:T, A:T, C:T")
data3 <- data.frame(value = c(result3$Vd$Vd_naive, # Vd
                              result3$Vd$Vd_ipsw,
                              result3$Vd$Vd_cwipw,
                              result3$Vd$Vd_cwor,
                              result3$Vd$Vd_ort,
                              result3$Vd$Vd_acw,
                              result3$pcd$pcd_naive, # PCD
                              result3$pcd$pcd_ipsw,
                              result3$pcd$pcd_cwipw,
                              result3$pcd$pcd_cwor,
                              result3$pcd$pcd_ort,
                              result3$pcd$pcd_acw,
                              result3$Vhat$Vhat_naive, # Vhat
                              result3$Vhat$Vhat_ipsw,
                              result3$Vhat$Vhat_cwipw,
                              result3$Vhat$Vhat_cwor,
                              result3$Vhat$Vhat_ort,
                              result3$Vhat$Vhat_acw),
                    type = rep(c("True Value", "PCD", "Estimated Value"), each = n.MonteCarlo * 6),
                    Estimators = rep(rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo), times = 3),
                    sce = "O:T / S:W, A:W, C:W")
data4 <- data.frame(value = c(result4$Vd$Vd_naive, # Vd
                              result4$Vd$Vd_ipsw,
                              result4$Vd$Vd_cwipw,
                              result4$Vd$Vd_cwor,
                              result4$Vd$Vd_ort,
                              result4$Vd$Vd_acw,
                              result4$pcd$pcd_naive, # PCD
                              result4$pcd$pcd_ipsw,
                              result4$pcd$pcd_cwipw,
                              result4$pcd$pcd_cwor,
                              result4$pcd$pcd_ort,
                              result4$pcd$pcd_acw,
                              result4$Vhat$Vhat_naive, # Vhat
                              result4$Vhat$Vhat_ipsw,
                              result4$Vhat$Vhat_cwipw,
                              result4$Vhat$Vhat_cwor,
                              result4$Vhat$Vhat_ort,
                              result4$Vhat$Vhat_acw),
                    type = rep(c("True Value", "PCD", "Estimated Value"), each = n.MonteCarlo * 6),
                    Estimators = rep(rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo), times = 3),
                    sce = "O:W / S:W, A:W, C:W")

data_final <- rbind(data1, data2, data3, data4)
data_final$value[data_final$value > 5] <- NA # remove extreme estimated values for better plots
data_final$Estimators <- factor(data_final$Estimators, levels = c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"))
hline_value <- data.frame(type = c("True Value", "PCD", "Estimated Value"), y.value = c(Vmax, 1, Vmax))
p_final <- ggplot(data_final, aes(Estimators, value, fill = Estimators)) +
             geom_boxplot() + 
             facet_grid(type ~ sce, scales = "free_y") +
             geom_hline(data = hline_value, aes(yintercept = y.value), linetype = "twodash", col = "blue") +
             theme(legend.position = "bottom",
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank()) +
             labs(y=NULL, x=NULL) + guides(fill = guide_legend(nrow = 1))
p_final

ggsave(plot = p_final,
       filename = "figures/double_robust.pdf",
       device = "pdf",
       width = 6.5,
       height = 8)


# summary table 
apply(result1$res$res.naive, 2, mean)
apply(result1$res$res.ipsw, 2, mean)
apply(result1$res$res.cw.ipw, 2, mean)
apply(result1$res$res.cw.or, 2, mean)
apply(result1$res$res.ort, 2, mean)
apply(result1$res$res.acw, 2, mean)

apply(cbind(result1$Vhat$Vhat_naive,
            result1$Vhat$Vhat_ipsw,
            result1$Vhat$Vhat_cwipw,
            result1$Vhat$Vhat_cwor,
            result1$Vhat$Vhat_ort,
            result1$Vhat$Vhat_acw), 2, sd)


apply(result2$res$res.naive, 2, mean)
apply(result2$res$res.ipsw, 2, mean)
apply(result2$res$res.cw.ipw, 2, mean)
apply(result2$res$res.cw.or, 2, mean)
apply(result2$res$res.ort, 2, mean)
apply(result2$res$res.acw, 2, mean)

apply(cbind(result2$Vhat$Vhat_naive,
            result2$Vhat$Vhat_ipsw,
            result2$Vhat$Vhat_cwipw,
            result2$Vhat$Vhat_cwor,
            result2$Vhat$Vhat_ort,
            result2$Vhat$Vhat_acw), 2, sd)


apply(result3$res$res.naive, 2, mean)
apply(result3$res$res.ipsw, 2, mean)
apply(result3$res$res.cw.ipw, 2, mean)
apply(result3$res$res.cw.or, 2, mean)
apply(result3$res$res.ort, 2, mean)
apply(result3$res$res.acw, 2, mean)

apply(cbind(result3$Vhat$Vhat_naive,
            result3$Vhat$Vhat_ipsw,
            result3$Vhat$Vhat_cwipw,
            result3$Vhat$Vhat_cwor,
            result3$Vhat$Vhat_ort,
            result3$Vhat$Vhat_acw), 2, sd)


apply(result4$res$res.naive, 2, mean)
apply(result4$res$res.ipsw, 2, mean)
apply(result4$res$res.cw.ipw, 2, mean)
apply(result4$res$res.cw.or, 2, mean)
apply(result4$res$res.ort, 2, mean)
apply(result4$res$res.acw, 2, mean)

apply(cbind(result4$Vhat$Vhat_naive,
            result4$Vhat$Vhat_ipsw,
            result4$Vhat$Vhat_cwipw,
            result4$Vhat$Vhat_cwor,
            result4$Vhat$Vhat_ort,
            result4$Vhat$Vhat_acw), 2, sd)


# ML models
load("survITR_NP1.RData")
n.MonteCarlo <- 200

data.np1.1 <- data.frame(value = c(Vhat_naive, Vhat_ipsw, Vhat_cwipw, Vhat_cwor, Vhat_ort, Vhat_acw),
                         Estimators = rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo))
data.np1.1$value[data.np1.1$value > 4] <- NA
data.np1.1$Estimators <- factor(data.np1.1$Estimators, levels = c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"))
p1.1 <- ggplot(data.np1.1, aes(Estimators, value, fill = Estimators)) + 
          geom_boxplot() + 
          geom_hline(yintercept = Vmax, linetype = "twodash", col = "blue") +
          theme(legend.position = "bottom",
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank()) +
          labs(x=NULL, y="Estimated Value") + guides(fill = guide_legend(nrow = 1))
p1.1

data.np1.2 <- data.frame(value = c(pcd_naive, pcd_ipsw, pcd_cwipw, pcd_cwor, pcd_ort, pcd_acw),
                         Estimators = rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo))
data.np1.2$Estimators <- factor(data.np1.2$Estimators, levels = c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"))
p1.2 <- ggplot(data.np1.2, aes(Estimators, value, fill = Estimators)) + 
          geom_boxplot() + 
          geom_hline(yintercept = 1, linetype = "twodash", col = "blue") +
          theme(legend.position = "bottom",
                axis.ticks.x=element_blank(),
                axis.text.x = element_blank()) +
          labs(x=NULL, y="PCD") + guides(fill = guide_legend(nrow = 1))
p1.2

data.np1.3 <- data.frame(value = c(Vd_naive, Vd_ipsw, Vd_cwipw, Vd_cwor, Vd_ort, Vd_acw),
                         Estimators = rep(c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"), each = n.MonteCarlo))
data.np1.3$Estimators <- factor(data.np1.3$Estimators, levels = c("Naive", "IPSW", "CW-IPW", "CW-OR", "ORt", "ACW"))
p1.3 <- ggplot(data.np1.3, aes(Estimators, value, fill = Estimators)) + 
          geom_boxplot() + 
          geom_hline(yintercept = Vmax, linetype = "twodash", col = "blue") +
          theme(legend.position = "bottom",
                axis.ticks.x=element_blank(),
                axis.text.x = element_blank()) +
          labs(x=NULL, y="True Value") + guides(fill = guide_legend(nrow = 1))
p1.3

library(ggpubr)
p1 <- ggarrange(p1.1, p1.2, p1.3,
                nrow = 1,
                common.legend = TRUE, legend = "bottom")
p1

ggsave(plot = p1,
       filename = "figures/ML_performance.pdf",
       device = "pdf",
       width = 6.5,
       height = 4.5)


# sample size
load("survITR_NP2.RData")
n.MonteCarlo <- 200
data.np2 <- data.frame(value = c(Vacw1, Vacw2, Vacw3, Vacw4, Vacw5, Vacw6),
                       N = rep(factor(c(0.5, 1, 2, 4, 6, 8)), each = n.MonteCarlo))
p2 <- ggplot(data.np2, aes(x = N, y = value)) +
        geom_boxplot() +
        geom_hline(yintercept = Vmax, linetype = "twodash", col = "blue") +
        theme(axis.ticks.x=element_blank()) +
        labs(x = "Target super population size", y = "Estimated Value (ACW)")
p2

ggsave(plot = p2,
       filename = "figures/ML_sample_size.pdf",
       device = "pdf",
       width = 5,
       height = 3.5)

# summary table
apply(cbind(Vacw1, Vacw2, Vacw3, Vacw4, Vacw5, Vacw6),
      2, mean) - Vmax

apply(cbind(CPacw1, CPacw2, CPacw3, CPacw4, CPacw5, CPacw6),
      2, mean)

apply(cbind(Vacw1, Vacw2, Vacw3, Vacw4, Vacw5, Vacw6),
      2, sd)
apply(cbind(SEacw1, SEacw2, SEacw3, SEacw4, SEacw5, SEacw6),
      2, mean)

