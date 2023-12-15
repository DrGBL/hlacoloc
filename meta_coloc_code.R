library(tidyverse)
library(vroom)
library(susieR)
library(rstanarm)
library(bayestestR)
library(ggpubr)

pathMS<-"ms.tsv"
pathMS_R<-"r2_ms.tsv"
pathEBNA<-"ebna.tsv"
pathEBNA_R<-"r2_ebna.tsv"

pheno1=read.table(pathEBNA,header = TRUE)
pheno2=read.table(pathMS,header = TRUE)
pheno1R=read.table(pathEBNA_R,header = TRUE)
pheno2R=read.table(pathMS_R,header = TRUE)
is_cohort_ld_pheno1=TRUE
is_cohort_ld_pheno2=TRUE
max_iter_susieR=10000

meta_coloc<-function(pheno1,pheno1R,is_cohort_ld_pheno1,
                     pheno2,pheno2R,is_cohort_ld_pheno2,
                     max_iter_susieR,plot_susie=TRUE){
  
  if( length(which(!pheno1$Name %in% pheno2$Name)) > 0 |
      length(which(!pheno2$Name %in% pheno1$Name)) > 0){
    print("Not same alleles in both cohorts.")
  } else {
    #susie pheno1
    susie_pheno1<-susieR::susie_rss(R=as.matrix(pheno1R),
                                    bhat=pheno1$beta,
                                    shat=pheno1$se,
                                    n=max(pheno1$N),
                                    estimate_residual_variance=is_cohort_ld_pheno1,
                                    max_iter=max_iter_susieR)
    
    #susie pheno2
    susie_pheno2<-susieR::susie_rss(R=as.matrix(pheno2R),
                                    bhat=pheno2$beta,
                                    shat=pheno2$se,
                                    n=max(pheno2$N),
                                    estimate_residual_variance=is_cohort_ld_pheno2,
                                    max_iter=max_iter_susieR)
    
    full_final<-inner_join(pheno1 %>% 
                             rename(beta1=beta,
                                    se1=se,
                                    N1=N),
                           pheno2 %>% 
                             rename(beta2=beta,
                                    se2=se,
                                    N2=N)) %>%
      mutate(pip_pheno1=susie_pheno1$pip) %>%
      mutate(pip_pheno2=susie_pheno2$pip) %>%
      mutate(coloc=pip_pheno1*pip_pheno2) %>% 
      mutate(gene=str_extract(Name,"[A-Z0-9]*\\*")) %>%
      mutate(gene=str_replace(gene,"\\*",""))
    
    full_final_summary<-full_final %>%
      group_by(gene) %>% 
      summarize(susie_coloc_prob = 1-prod(1-coloc))
    
    
    
    
    #get the bayesian correlation probability
    #specifically this is the probability that the correlation term is not zero.
    bayes_corr_df<-c()
    posteriors_fit<-c()
    posteriors_ci<-c()
    for(g in (full_final %>% pull(gene) %>% unique())){
      bayes_lm<-stan_glm(pip_pheno2~pip_pheno1, data=full_final %>% filter(gene==g),family = gaussian(link = "identity"))
      
      posteriors_ci<-as.data.frame(as.matrix(bayes_lm)) %>%
        mutate(gene=g) %>%
        rename(alpha=`(Intercept)`) %>%
        rename(beta=pip_pheno1) %>%
        head(n=100) %>%
        bind_rows(posteriors_ci,.)
      
      posteriors_fit<-data.frame(alpha=bayes_lm$coefficients[1],
                                 beta=bayes_lm$coefficients[2],
                                 gene=g) %>%
        bind_rows(posteriors_fit,.)
      
      bayes_corr_df<-data.frame(gene=g,
                                bayes_pd=p_direction(bayes_lm)$pd[2],
                                direction_of_correlation=ifelse(bayes_lm$coefficients[2] >= 0,
                                                                "Correct",
                                                                "Incorrect")) %>%
        bind_rows(bayes_corr_df,.)
    }
    
    #build a summary table
    
    full_final_summary<-full_final_summary %>% 
      left_join(.,bayes_corr_df) %>%
      mutate(meta_colocalization_probability=ifelse(is.na(susie_coloc_prob) | is.na(bayes_pd),
                                                    NA,
                                                    susie_coloc_prob*bayes_pd)) %>%
      arrange(desc(meta_colocalization_probability),desc(susie_coloc_prob))
    
    
    
    if(plot_susie==TRUE){
      #ebna
      reg_coloc<-full_final %>%
        ggplot(aes(x=beta1,y=beta2))+
        geom_point()+
        geom_smooth(method="lm",formula=y~x-1)+
        facet_wrap(~gene, scales="free",ncol=1)+
        theme_bw()+
        ylab("Pheno1 Betas")+
        xlab("Pheno2 Betas")
      
      annotate_df<-c()
      for(g in full_final %>% pull(gene) %>% unique()){
        
        annotate_df<-data.frame(gene=g,
                                x=0.25,
                                y=0.75,
                                text=paste0("P_coloc: ", format(full_final_summary %>% 
                                                                  filter(gene==g) %>% 
                                                                  pull(meta_colocalization_probability), digits=2))) %>%
          bind_rows(annotate_df,)
        
        
      }
      
      pip_coloc<-full_final %>%
        ggplot(aes(x=pip_pheno1,y=pip_pheno2))+
        geom_abline(data=posteriors_ci, aes(intercept=alpha, slope=beta),color = "skyblue", linewidth = 0.2, alpha = 0.25)+
        geom_abline(data=posteriors_fit, aes(intercept=alpha, slope=beta),color = "black", linewidth = 1)+
        geom_point()+
        geom_text(data=annotate_df,aes(x=x,y=y,label=text),inherit.aes = FALSE)+
        facet_wrap(~gene,ncol=1)+
        theme_bw()+
        coord_cartesian(ylim = c(0, 1), xlim = c(0,1))+
        xlab("Pheno1 PIPs")+
        ylab("Pheno2 PIPs")
      
      plot_susie_pip<-ggarrange(reg_coloc,pip_coloc,ncol=2,labels=c("a","b"))
      
      
      
    } else {
      plot_susie_pip<-NA
    }
    
    return(list(meta_colocalization=full_final_summary,
                plot=plot_susie_pip))
    
  }
}
