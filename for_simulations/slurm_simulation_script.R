#!/usr/bin/env Rscript

library(optparse)

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--anc"), type="character", default=NULL,
              help="name of ancestry or cohort", metavar="character"),
  make_option(c("--ids"), type="character", default=NULL,
              help="path to file containing IDs of cohort", metavar="character"),
  make_option(c("--call_folder"), type="character", default=NULL,
              help="path to folder with QCed calls", metavar="character"),
  make_option(c("--af_folder"), type="character", default=NULL,
              help="path to folder with af and counts", metavar="character"),
  make_option(c("--bim_folder"), type="character", default=NULL,
              help="path to bim file", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt",
              help="output file name with path", metavar="character"),
  make_option(c("--n_min_alleles"), type="numeric", default=10,
              help="Minimum number of alleles in a gene for inclusion", metavar="integer"),
  make_option(c("--all_causal_alleles"), type="logical", default=FALSE,
              help="For the causal genes, are all alleles causal? If FALSE (default), then up to 50% of alleles are randomly assigned an effect of 0.", metavar="logical"),
  make_option(c("--all_causal_genes"), type="logical", default=FALSE,
              help="For the causal genes, are all alleles causal? If TRUE (default), then up to 1/3 of gene are randomly assigned an effect of 0.", metavar="logical"),
  make_option(c("--pheno2_mean_shared_factor"), type="numeric", default=0,
              help="Mean of multiplicative factor between genes with shared signal", metavar="numeric"),
  make_option(c("--pheno2_variance_shared_factor"), type="numeric", default=1,
              help="Variance of multiplicative factor between genes with shared signal", metavar="numeric"),
  make_option(c("--gen_var"), type="numeric", default=0.4,
              help="Percentage of Variance explained by genetic markers (on average)", metavar="numeric"),
  make_option(c("--freq_threshold"), type="numeric", default=0.005,
              help="MAF threshold to include an HLA allele in analysis", metavar="numeric"),
  make_option(c("--plot_susie"), type="logical", default=FALSE,
              help="Whether to plot the susie PIPs or not. Useful to visualize the pearson correlation. Default=FALSE.", metavar="logical"),
  make_option(c("--iteration"), type="numeric", default=1,
              help="Iteration number. This has no impact on the simulation, it's just a way to keep track of things in the output.", metavar="numeric"),
  make_option(c("--continuous_trait"), type="logical", default=TRUE,
              help="Whether the traits are continuous (TRUE, default) or binary (FALSE).", metavar="logical"),
  make_option(c("--n_clone"), type="numeric", default=1,
              help="Number of times to clone the cohort, helpful for case control studies simulations where power is limited", metavar="numeric"),
  make_option(c("--susie_L"), type="numeric", default=10,
              help="Maximum number of non-zero effects in the susie regression model. Default = 10", metavar="numeric"),
  make_option(c("--fix_alpha"), type="logical", default=TRUE,
              help="Whether to fix the intercept of PIP regression in STAN to 0 (not recommended) or also estimate it (the default).", metavar="logical")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(tidyverse)
library(vroom)
library(susieR)
library(rstanarm)
library(bayestestR)
library(ggpubr)

#for testing
#can comment out when ready
# opt<-data.frame(anc="sas",
#                 bim_folder="~/Microbiologie/RichardsLab/HLA/UKB_WES_HLA/",
#                 ids="~/Microbiologie/RichardsLab/HLA/ebv_coloc/ukb.sasIDsPCA.txt",
#                 call_folder="~/Microbiologie/RichardsLab/HLA/UKB_WES_HLA/qced_calls/",
#                 af_folder="~/Microbiologie/RichardsLab/HLA/UKB_WES_HLA/",
#                 ld_folder="~/Microbiologie/RichardsLab/HLA/UKB_WES_HLA/")

anc<-opt$anc
ids<-scan(opt$ids,
          what = character())

hla_alleles<-read.table(paste0(opt$bim_folder, "ld_six_digit_", anc, ".bim"), header=FALSE)
genes<-hla_alleles %>%
  dplyr::select(V2) %>%
  rename(allele=V2) %>%
  mutate(gene=str_extract(allele, "^[A-Z0-9]*"))

#ld
ld<-read_tsv(paste0(opt$ld_folder,"ld_six_digit_",anc,".ld.gz"), col_names=FALSE)

#haplotypes for that population
hap<-c()
for(f in 10:60){
  hap<-read_tsv(paste0(opt$call_folder,
                       "hla_df_batch_qced_",f,".tsv.gz"),
                col_types = cols(.default = "c")) %>%
    bind_rows(hap,.)
}
hap<-hap %>%
  filter(ID %in% ids)

#now arrange the haplotypes in a 0-1-2 fashion for linear regression
#first extract all alleles
all_alleles<-hap %>%
  dplyr::select(-c(ID)) %>%
  pivot_longer(cols=everything()) %>%
  filter(value!="0") %>%
  pull(value) %>%
  unique()

#now build the data frame with NA values
hap_lm_ready<-matrix(0,nrow=nrow(hap),ncol=1+length(all_alleles))
colnames(hap_lm_ready)<-c("ID",all_alleles)
hap_lm_ready<-as.data.frame(hap_lm_ready)
hap_lm_ready$ID<-hap$ID

#a special tmp column
hap_concat<-hap %>% 
  bind_cols(unite(data=hap[,-1],col=alleles_concat),.) %>%
  dplyr::select(c(ID, alleles_concat)) %>% 
  mutate(alleles_concat=paste0("_",alleles_concat,"_")) %>%
  mutate(alleles_concat=str_replace_all(alleles_concat,"_","__"))
hap_lm_ready<-left_join(hap_lm_ready,hap_concat,by="ID")

#now loop through the alleles
for(i in colnames(hap_lm_ready %>% dplyr::select(-c("ID","alleles_concat")))){
  print(i)
  hap_lm_ready[,i]<-str_count(hap_lm_ready$alleles_concat,gsub("\\*","\\\\*",paste0("_",i,"_")))
}

#clean the data
hap_lm_ready<-hap_lm_ready %>% dplyr::select(-c(alleles_concat))

#allele and their counts
ac<-read.table(paste0(opt$af_folder, "freq_six_digit_",anc,".frq.counts.gz"), header=TRUE) %>%
  filter(SNP %in% genes$allele) %>%
  mutate(flip=ifelse(A1=="C",TRUE,FALSE)) %>%
  mutate(tmp_C1=C1) %>%
  mutate(A1=ifelse(flip,"A",A1)) %>%
  mutate(A2=ifelse(flip,"C",A2)) %>%
  mutate(C1=ifelse(flip,C2,C1)) %>%
  mutate(C2=ifelse(flip,tmp_C1,C2)) %>%
  dplyr::select(-tmp_C1) %>%
  rename(allele=SNP) %>%
  rename(count=C1) %>%
  dplyr::select(c(allele,count)) %>%
  mutate(freq=count/(2*nrow(hap)))

genes<-genes %>%
  left_join(.,
            ac) %>%
  group_by(gene) %>%
  mutate(freq_adj=count/sum(count))

#for testing purposes, can comment out later
# n_min_alleles=10
# pheno2_mean_shared_factor=0
# pheno2_variance_shared_factor=1
# plot_predicted=TRUE
# gen_var=0.4
# haplotypes=hap
# haplotypes_lm=hap_lm_ready
# freq_threshold=0.005
# genes=genes
# all_causal_alleles=FALSE
# plot_susie=FALSE
# all_causal_genes=FALSE
# continuous_trait=TRUE
# negative_threshold=0.1
# n_clone=10
# susie_L=10
# fix_alpha=FALSE

simulate_hla<-function(n_min_alleles=10,
                       pheno2_mean_shared_factor=0,
                       pheno2_variance_shared_factor=1,
                       plot_predicted=TRUE,
                       gen_var=0.4,
                       haplotypes=hap,
                       haplotypes_lm=hap_lm_ready,
                       freq_threshold=0.005,
                       genes=genes,
                       all_causal_alleles=FALSE,
                       all_causal_genes=FALSE,
                       plot_susie=FALSE,
                       continuous_trait=TRUE,
                       negative_threshold=0.1,
                       n_clone=1,
                       susie_L=10,
                       fix_alpha=FALSE){
  
  #apply frequency and allele number thresholds to obtain list of genes to simulate
  genes_filtered<-genes %>%
    filter(freq>=freq_threshold) %>%
    group_by(gene) %>%
    filter(n()>=n_min_alleles) %>%
    ungroup()
  
  genes_to_keep<-genes_filtered %>%
    dplyr::select(gene) %>%
    distinct()
  
  #obtain the ld matrix corresponding the alleles remaining after the filters above
  alleles_to_keep<-which(genes$allele %in% genes_filtered$allele)
  ld_filtered<-as.matrix(ld[alleles_to_keep,alleles_to_keep])
  
  #simulate gene causal effects
  #here is also where you decide if you want to make some genes non-causal (at most one third for each phenotypes).
  n_shared_causal_pheno2<-sample(c(1:nrow(genes_to_keep)),1)
  
  if(all_causal_genes==TRUE){
    genes_causal<-data.frame(gene=genes_to_keep) %>%
      mutate(var_pheno1=rchisq(nrow(.), 1)) %>%
      mutate(var_pheno2=rchisq(nrow(.), 1)) %>%
      mutate(causal_gene_pheno1=TRUE) %>%
      mutate(causal_gene_pheno2=TRUE) %>%
      mutate(shared_genetics=ifelse(row_number() %in% sample(1:nrow(genes_to_keep),n_shared_causal_pheno2),
                                    TRUE,FALSE)) %>%
      mutate(tmp_factor=rnorm(nrow(.), mean=pheno2_mean_shared_factor, sd=sqrt(pheno2_variance_shared_factor))) %>%
      mutate(factor_shared_causal_pheno2=ifelse(shared_genetics==TRUE,
                                                tmp_factor,
                                                0)) %>%
      mutate(var_pheno2=ifelse(shared_genetics==TRUE,
                               0,
                               var_pheno2)) %>%
      dplyr::select(-tmp_factor)
  } else {
    vec_causal<-c(rep(FALSE,floor(nrow(genes_to_keep)/3)), rep(TRUE, nrow(genes_to_keep)-floor(nrow(genes_to_keep)/3)))
    
    genes_causal<-data.frame(gene=genes_to_keep) %>%
      mutate(var_pheno1=rchisq(nrow(.), 1)) %>%
      mutate(var_pheno2=rchisq(nrow(.), 1)) %>%
      mutate(causal_gene_pheno1=sample(vec_causal,
                                       length(vec_causal),
                                       replace = FALSE)) %>%
      mutate(causal_gene_pheno2=sample(vec_causal,
                                       length(vec_causal),
                                       replace = FALSE)) %>%
      mutate(var_pheno1=ifelse(causal_gene_pheno1==TRUE,
                               var_pheno1,
                               0)) %>%
      mutate(var_pheno2=ifelse(causal_gene_pheno2==TRUE,
                               var_pheno2,
                               0)) %>%
      mutate(shared_genetics=ifelse(row_number() %in% sample(1:nrow(genes_to_keep),n_shared_causal_pheno2) &
                                      causal_gene_pheno1==TRUE &
                                      causal_gene_pheno2==TRUE,
                                    TRUE,FALSE)) %>%
      mutate(tmp_factor=rnorm(nrow(.), mean=pheno2_mean_shared_factor, sd=sqrt(pheno2_variance_shared_factor))) %>%
      mutate(factor_shared_causal_pheno2=ifelse(shared_genetics==TRUE,
                                                tmp_factor,
                                                0)) %>%
      mutate(var_pheno2=ifelse(shared_genetics==TRUE,
                               0,
                               var_pheno2)) %>%
      dplyr::select(-tmp_factor)
  }
  
  #combine with all effect in a dataframe
  effect<-genes_filtered %>%
    left_join(.,genes_causal)
  
  #here is where allele causal effect are sampled from the gene level causal effect variances drawn above
  #the variance used to draw allele level causal effect is further adjusted by the allele frequency
  for(g in genes_to_keep$gene){
    print(g)
    
    eff_tmp<-effect %>% 
      filter(gene==g)
    
    
    eff_tmp$causal_effect_pheno1_tmp<-NA
    eff_tmp$causal_effect_pheno2_tmp<-NA
    
    for(i in 1:nrow(eff_tmp)){
      eff_tmp$causal_effect_pheno1_tmp[i]<-ifelse(eff_tmp$var_pheno1[i]==0,
                                                  0,
                                                  rnorm(1,0,sqrt(eff_tmp$var_pheno1[i]/(eff_tmp$freq_adj[i]*(1-eff_tmp$freq_adj[i])))))
      eff_tmp$causal_effect_pheno2_tmp[i]<-ifelse(eff_tmp$var_pheno2[i]==0,
                                                  eff_tmp$causal_effect_pheno1_tmp[i]*eff_tmp$factor_shared_causal_pheno2[i],
                                                  rnorm(1,0,sqrt(eff_tmp$var_pheno2[i]/(eff_tmp$freq_adj[i]*(1-eff_tmp$freq_adj[i])))))
    }
    
    effect<-effect %>%
      filter(!(allele %in% eff_tmp$allele)) %>%
      bind_rows(.,
                eff_tmp)
    
  }
  
  #if decided that not all alleles are causal, put some allele level effects to 0
  effect<-effect %>%
    mutate(drop_to_zero_pheno1=FALSE,
           drop_to_zero_pheno2=FALSE)
  if(all_causal_alleles==FALSE){
    for(g in genes_to_keep$gene){
      non_causal_pheno1<-sample(0:ceiling((effect %>%
                                           filter(gene==g) %>%
                                           nrow())/5),
                              1)
      row_non_causal_pheno1<-sample(1:ceiling((effect %>%
                                             filter(gene==g) %>%
                                             nrow())),
                                non_causal_pheno1)
      effect<-effect %>%
        group_by(gene) %>%
        mutate(drop_to_zero_pheno1=ifelse(gene==g & (row_number() %in% row_non_causal_pheno1),
                                          TRUE,
                                          drop_to_zero_pheno1)) %>%
        ungroup()
      
      non_causal_pheno2<-sample(0:ceiling((effect %>%
                                           filter(gene==g) %>%
                                           nrow())/5),
                              1)
      row_non_causal_pheno2<-sample(1:(effect %>%
                                     filter(gene==g) %>%
                                     nrow()),
                                    non_causal_pheno2)
      effect<-effect %>%
        group_by(gene) %>%
        mutate(drop_to_zero_pheno2=ifelse(gene==g & (row_number() %in% row_non_causal_pheno2),
                                          TRUE,
                                          drop_to_zero_pheno2)) %>%
        ungroup()
    }
  }
  
  #update the effect data frame with the above data
  #then adjust the overall allele level effects so that each gene has an weighted average of 0
  #weights are proportional to frequency
  #this represents the fact that allele effects at the HLA are always relative to other alleles, hence the average population effect is zero.
  effect<-effect %>%
    mutate(causal_effect_pheno1_tmp=ifelse(drop_to_zero_pheno1==TRUE,
                                           0,
                                           causal_effect_pheno1_tmp)) %>%
    mutate(causal_effect_pheno2_tmp=ifelse(drop_to_zero_pheno1==TRUE &
                                             shared_genetics==TRUE,
                                           0,
                                           causal_effect_pheno2_tmp)) %>%
    mutate(causal_effect_pheno2_tmp=ifelse(drop_to_zero_pheno2==TRUE &
                                             shared_genetics==FALSE,
                                           0,
                                           causal_effect_pheno2_tmp)) %>%
    group_by(gene) %>%
    mutate(weighted_mean_effect1=sum(freq_adj*causal_effect_pheno1_tmp)) %>%
    mutate(weighted_mean_effect2=sum(freq_adj*causal_effect_pheno2_tmp)) %>%
    mutate(causal_effect_pheno1=causal_effect_pheno1_tmp-weighted_mean_effect1) %>%
    mutate(causal_effect_pheno2=causal_effect_pheno2_tmp-weighted_mean_effect2) %>%
    ungroup()
  
  #check that the effects aren't correlated by chance only
  mod_df<-data.frame(gene=genes_to_keep$gene,
                     pval_pheno_effects=rep(NA,nrow(genes_to_keep)),
                     z_pheno_effects=rep(NA,nrow(genes_to_keep)))
  for(i in 1:nrow(mod_df)){
    tmp<-effect %>% filter(gene==genes_to_keep$gene[i])
    mod<-lm(causal_effect_pheno1~causal_effect_pheno2,data=tmp)
    if(nrow(summary(mod)$coefficients)==1){
      mod_df$pval_pheno_effects[i]<-NA
      mod_df$z_pheno_effects[i]<-NA
    } else {
      mod_df$pval_pheno_effects[i]<-summary(mod)$coefficients[2,4]
      mod_df$z_pheno_effects[i]<-summary(mod)$coefficients[2,1]/summary(mod)$coefficients[2,2]
    }
  }
  
  #now simulate individual level phenotypes
  pheno_sims<-data.frame(ID=haplotypes$ID,
                         pheno1=NA,
                         pheno2=NA)
  for(i in 1:nrow(haplotypes)){
    if(i %% 100==0){
      print(paste0(i, " of ", nrow(haplotypes), "."))
    }
    
    alleles_tmp<-data.frame(allele=t(haplotypes[i,-1])) %>%
      left_join(.,effect,
                by = "allele")
    
    pheno_sims$pheno1[i]<-sum(alleles_tmp$causal_effect_pheno1, na.rm=TRUE)
    pheno_sims$pheno2[i]<-sum(alleles_tmp$causal_effect_pheno2, na.rm=TRUE)
  }
  
  #clone the samples to increase power
  if(continuous_trait==FALSE){
    pheno_sims_orig<-pheno_sims
    haplotypes_lm_orig<-haplotypes_lm
    for(i in 1:n_clone){
      pheno_sims<-pheno_sims_orig %>%
        mutate(ID=paste(ID,"_",i)) %>%
        bind_rows(pheno_sims,.)
      haplotypes_lm<-haplotypes_lm_orig %>%
        mutate(ID=paste(ID,"_",i)) %>%
        bind_rows(haplotypes_lm,.)
    }
  }
  
  #add error term so that the causal genes explain roughly 40% of variance (value can be changed above, but 40% is default)
  
  var1<-var(pheno_sims$pheno1)
  var2<-var(pheno_sims$pheno2)
  
  pheno_sims$pheno1_rand<-pheno_sims$pheno1+rnorm(nrow(pheno_sims),0,sqrt(var1*(1-gen_var)/gen_var))
  pheno_sims$pheno2_rand<-pheno_sims$pheno2+rnorm(nrow(pheno_sims),0,sqrt(var2*(1-gen_var)/gen_var))
  
  #if binary trait, turn these into 0 or 1
  if(continuous_trait==FALSE){
    pheno_sims$pheno_invlogit1<-invlogit(pheno_sims$pheno1_rand)
    pheno_sims$pheno_invlogit2<-invlogit(pheno_sims$pheno2_rand)
    pheno_sims$pheno1_binary<-rbernoulli(nrow(pheno_sims),pheno_sims$pheno_invlogit1)
    pheno_sims$pheno2_binary<-rbernoulli(nrow(pheno_sims),pheno_sims$pheno_invlogit2)
    
    
    
    
    #pheno_sims<-pheno_sims %>% arrange(pheno1_rand)
    #pheno_sims$pheno1_binary<-c(rep(0,floor(nrow(pheno_sims)*0.6)), rep(1,nrow(pheno_sims)-floor(nrow(pheno_sims)*0.6)))
    #pheno_sims<-pheno_sims %>% arrange(pheno2_rand)
    #pheno_sims$pheno2_binary<-c(rep(0,floor(nrow(pheno_sims)*0.6)), rep(1,nrow(pheno_sims)-floor(nrow(pheno_sims)*0.6)))
  }
  
  #now run gwas of HLA alleles on both phenotypes
  gwas_res<-matrix(NA,nrow=ncol(haplotypes_lm)-1,ncol=7)
  colnames(gwas_res)<-c("allele", "beta1", "se1", "pval1", "beta2", "se2", "pval2")
  gwas_res<-as.data.frame(gwas_res)
  gwas_res$allele<-colnames(haplotypes_lm[,-1])
  gwas_res<-gwas_res %>%
    filter(allele %in% effect$allele)
  if(continuous_trait==TRUE){
    for(i in 1:nrow(gwas_res)){
      if(i %% 100==0){
        print(paste0(i, " of ", nrow(gwas_res), "."))
      }
      lm_1<-lm(scale(pheno_sims$pheno1_rand)~haplotypes_lm[,gwas_res$allele[i]])
      lm_2<-lm(scale(pheno_sims$pheno2_rand)~haplotypes_lm[,gwas_res$allele[i]])
      
      gwas_res$beta1[i]<-coefficients(lm_1)[2]
      gwas_res$se1[i]<-sqrt(diag(vcov(lm_1)))[2]
      gwas_res$pval1[i]<-summary(lm_1)$coefficients[2,4]
      gwas_res$beta2[i]<-coefficients(lm_2)[2]
      gwas_res$se2[i]<-sqrt(diag(vcov(lm_2)))[2]
      gwas_res$pval2[i]<-summary(lm_2)$coefficients[2,4]
    }
    
    gwas_res<-gwas_res %>%
      mutate(gene=str_extract(allele,"^[A-Z0-9]*")) %>%
      mutate(weight=1/(se1^2+se2^2)) %>%
      left_join(.,genes)
    
    min_p_pheno1<-gwas_res %>%
      arrange(pval1) %>%
      group_by(gene) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      rename(min_pval_pheno1=pval1) %>%
      dplyr::select(c(gene,min_pval_pheno1))
    
    min_p_pheno2<-gwas_res %>%
      arrange(pval2) %>%
      group_by(gene) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      rename(min_pval_pheno2=pval2) %>%
      dplyr::select(c(gene,min_pval_pheno2))
    
    #variances explained by each genes on each phenotype
    varexp<-data.frame(gene=genes_to_keep,
                       r2_1=NA,
                       adj_r2_1=NA,
                       r2_2=NA,
                       adj_r2_2=NA,
                       role="beta")
    for(i in 1:nrow(varexp)){
      alleles_tmp<-effect %>%
        filter(gene==varexp$gene[i]) %>%
        mutate(allele=paste0("haplotypes_lm[,'",
                             allele,
                             "']")) %>%
        pull(allele)
      
      
      form_tmp1<-paste0("scale(pheno_sims$pheno1_rand)~",
                        paste(alleles_tmp,collapse="+"))
      mod_tmp1<-lm(as.formula(form_tmp1))
      varexp$r2_1[i]<-summary(mod_tmp1)$r.squared
      varexp$adj_r2_1[i]<-summary(mod_tmp1)$adj.r.squared
      
      form_tmp2<-paste0("scale(pheno_sims$pheno2_rand)~",
                        paste(alleles_tmp,collapse="+"))
      mod_tmp2<-lm(as.formula(form_tmp2))
      varexp$r2_2[i]<-summary(mod_tmp2)$r.squared
      varexp$adj_r2_2[i]<-summary(mod_tmp2)$adj.r.squared
      
    } 
  } else {
    for(i in 1:nrow(gwas_res)){
      if(i %% 100==0){
        print(paste0(i, " of ", nrow(gwas_res), "."))
      }
      glm_1<-glm(pheno_sims$pheno1_binary~haplotypes_lm[,gwas_res$allele[i]],family=binomial(link="logit"))
      glm_2<-glm(pheno_sims$pheno2_binary~haplotypes_lm[,gwas_res$allele[i]],family=binomial(link="logit"))
      
      gwas_res$beta1[i]<-coefficients(glm_1)[2]
      gwas_res$se1[i]<-sqrt(diag(vcov(glm_1)))[2]
      gwas_res$pval1[i]<-summary(glm_1)$coefficients[2,4]
      gwas_res$beta2[i]<-coefficients(glm_2)[2]
      gwas_res$se2[i]<-sqrt(diag(vcov(glm_2)))[2]
      gwas_res$pval2[i]<-summary(glm_2)$coefficients[2,4]
    }
    
    gwas_res<-gwas_res %>%
      mutate(gene=str_extract(allele,"^[A-Z0-9]*")) %>%
      mutate(weight=1/(se1^2+se2^2)) %>%
      left_join(.,genes)
    
    min_p_pheno1<-gwas_res %>%
      arrange(pval1) %>%
      group_by(gene) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      rename(min_pval_pheno1=pval1) %>%
      dplyr::select(c(gene,min_pval_pheno1))
    
    min_p_pheno2<-gwas_res %>%
      arrange(pval1) %>%
      group_by(gene) %>%
      filter(row_number()==1) %>%
      ungroup() %>%
      rename(min_pval_pheno2=pval2) %>%
      dplyr::select(c(gene,min_pval_pheno2))
    
    #variances explained by each genes on each phenotype
    varexp<-data.frame(gene=genes_to_keep,
                       r2_1=NA,
                       adj_r2_1=NA,
                       r2_2=NA,
                       adj_r2_2=NA,
                       role="beta")
    for(i in 1:nrow(varexp)){
      alleles_tmp<-effect %>%
        filter(gene==varexp$gene[i]) %>%
        mutate(allele=paste0("haplotypes_lm[,'",
                             allele,
                             "']")) %>%
        pull(allele)
      
      
      form_tmp1<-paste0("pheno_sims$pheno1_binary~",
                        paste(alleles_tmp,collapse="+"))
      mod_tmp1<-glm(as.formula(form_tmp1), family=binomial(link="logit"))
      varexp$r2_1[i]<-1 - summary(mod_tmp1)$deviance/summary(mod_tmp1)$null.deviance
      
      form_tmp2<-paste0("pheno_sims$pheno2_binary~",
                        paste(alleles_tmp,collapse="+"))
      mod_tmp2<-glm(as.formula(form_tmp2), family=binomial(link="logit"))
      varexp$r2_2[i]<-1 - summary(mod_tmp2)$deviance/summary(mod_tmp2)$null.deviance
      
    }
  }
  
  
  #susie pheno1
  if(is.na(susie_L)){
    susie_L<-nrow(ld_filtered)
  }
  susie1<-susieR::susie_rss(R=ld_filtered,
                            bhat=left_join(effect %>% dplyr::select(allele),gwas_res)$beta1,
                            shat=left_join(effect %>% dplyr::select(allele),gwas_res)$se1,
                            n=nrow(pheno_sims),
                            estimate_residual_variance=TRUE,
                            L=susie_L)
  
  #susie pheno2
  susie2<-susieR::susie_rss(R=ld_filtered,
                            bhat=left_join(effect %>% dplyr::select(allele),gwas_res)$beta2,
                            shat=left_join(effect %>% dplyr::select(allele),gwas_res)$se2,
                            n=nrow(pheno_sims),
                            estimate_residual_variance=TRUE,
                            L=susie_L)
  
  #combine the pips in the summary stats dataframe
  final_summ_stats<-left_join(effect,gwas_res) %>%
    mutate(pip_pheno1=susie1$pip) %>%
    mutate(pip_pheno2=susie2$pip) %>%
    mutate(coloc=pip_pheno1*pip_pheno2)
  
  
  
  #get the bayesian correlation probability
  #specifically this is the probability that the correlation term is not zero.
  posteriors_fit<-c()
  posteriors_ci<-c()
  if(fix_alpha==FALSE){
    bayes_corr_df<-data.frame("gene"=NA, "bayes_pd"=NA, "direction_of_correlation"=NA, "bayes_pd2"=NA,               
                              "odds_p_map"=NA, "mean_alpha"=NA, "median_alpha"=NA, "low95_alpha"=NA,            
                              "high95_alpha"=NA, "posterior_prob_map"=NA)
    
    for(g in genes_to_keep$gene){
      if(length(which(final_summ_stats %>% filter(gene==g) %>% pull(pip_pheno1) > negative_threshold)) == 0 |
         length(which(final_summ_stats %>% filter(gene==g) %>% pull(pip_pheno2) > negative_threshold)) == 0){
        bayes_corr_df<-data.frame(gene=g, bayes_pd=NA, direction_of_correlation=NA) %>% bind_rows(bayes_corr_df,.)
      } else {
        bayes_lm<-stan_glm(pip_pheno2~pip_pheno1, 
                           data=final_summ_stats %>% filter(gene==g),
                           family = gaussian(link = "identity"))
        
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
                                  odds_p_map=1/p_map(bayes_lm)$p_MAP[2],
                                  mean_alpha=mean(posteriors_fit$alpha),
                                  median_alpha=median(posteriors_fit$alpha),
                                  low95_alpha=quantile(posteriors_fit$alpha,0.025),
                                  high95_alpha=quantile(posteriors_fit$alpha,0.975),
                                  #nullInterval = c(0,1))$pd,
                                  direction_of_correlation=ifelse(bayes_lm$coefficients[2] >= 0,
                                                                  "Correct",
                                                                  "Incorrect")) %>%
          mutate(posterior_prob_map=ifelse(p_map(bayes_lm)$p_MAP[2]==0,1,odds_p_map/(1+odds_p_map))) %>%
          bind_rows(bayes_corr_df,.)
      }
    }
  } else {
    bayes_corr_df<-data.frame("gene"=NA, "bayes_pd"=NA, "direction_of_correlation"=NA, "bayes_pd2"=NA,               
                              "odds_p_map"=NA, "posterior_prob_map"=NA)
    
    for(g in genes_to_keep$gene){
      if(length(which(final_summ_stats %>% filter(gene==g) %>% pull(pip_pheno1) > negative_threshold)) == 0 |
         length(which(final_summ_stats %>% filter(gene==g) %>% pull(pip_pheno2) > negative_threshold)) == 0){
        bayes_corr_df<-data.frame(gene=g, bayes_pd=NA, direction_of_correlation=NA) %>% bind_rows(bayes_corr_df,.)
      } else {
        bayes_lm<-stan_glm(pip_pheno2~pip_pheno1-1, 
                           data=final_summ_stats %>% filter(gene==g),
                           family = gaussian(link = "identity"))
        
        posteriors_ci<-as.data.frame(as.matrix(bayes_lm)) %>%
          mutate(gene=g) %>%
          rename(beta=pip_pheno1) %>%
          head(n=100) %>%
          bind_rows(posteriors_ci,.)
        
        posteriors_fit<-data.frame(beta=bayes_lm$coefficients[1],
                                   gene=g) %>%
          bind_rows(posteriors_fit,.)
        
        bayes_corr_df<-data.frame(gene=g,
                                  bayes_pd=p_direction(bayes_lm)$pd,
                                  odds_p_map=1/p_map(bayes_lm)$p_MAP,
                                  #nullInterval = c(0,1))$pd,
                                  direction_of_correlation=ifelse(bayes_lm$coefficients[1] >= 0,
                                                                  "Correct",
                                                                  "Incorrect")) %>%
          mutate(posterior_prob_map=ifelse(p_map(bayes_lm)$p_MAP==0,1,odds_p_map/(1+odds_p_map))) %>%
          bind_rows(bayes_corr_df,.)
      }
    }
  }
  
  #build a summary table
  summary_res<-final_summ_stats %>% 
    group_by(gene) %>% 
    summarize(susie_coloc_prob = 1-prod(1-coloc),
              var_causal_effect_pheno1=var(causal_effect_pheno1),
              var_causal_effect_pheno2=var(causal_effect_pheno2),
              correlation_p=cor.test(pip_pheno1,pip_pheno2,method="pearson",use = "complete.obs")$p.value) %>%
    left_join(.,genes_causal) %>%
    left_join(.,varexp) %>%
    left_join(.,bayes_corr_df) %>%
    mutate(colocalizes_frequentist=ifelse(correlation_p<0.05 & susie_coloc_prob>0.8,
                                          TRUE,
                                          FALSE)) %>%
    mutate(colocalizes_frequentist=ifelse(is.na(colocalizes_frequentist),
                                          FALSE,
                                          colocalizes_frequentist)) %>%
    mutate(bayesian_prob=ifelse(is.na(susie_coloc_prob) | is.na(posterior_prob_map),
                                  NA,
                                  susie_coloc_prob*posterior_prob_map)) %>%
    arrange(desc(shared_genetics),desc(bayesian_prob),desc(susie_coloc_prob),correlation_p,desc(var_pheno1),desc(var_pheno2)) %>%
    left_join(.,min_p_pheno1) %>%
    left_join(.,min_p_pheno2) %>%
    left_join(.,mod_df)
  
  #effect %>% ggplot(aes(x=causal_effect_pheno1,y=causal_effect_pheno2)) +geom_point()+facet_wrap(~gene)
  
  #plot pips if want to
  if(plot_susie==TRUE){
    
    reg_coloc<-final_summ_stats %>% 
      filter(gene %in% genes$gene) %>%
      ggplot(aes(x=beta1,y=beta2))+
      geom_point()+
      geom_smooth(method="lm",formula=y~x-1)+
      facet_wrap(~gene, scales="free",ncol=1)+
      theme_bw()+
      xlab("Pheno1 Betas")+
      ylab("Pheno2 Betas")
    
    annotate_df<-c()
    for(g in unique(final_summ_stats$gene)){
      
      annotate_df<-data.frame(gene=g,
                              x=0.25,
                              y=0.75,
                              text=paste0("P_coloc: ", format(summary_res %>% 
                                                                filter(gene==g) %>% 
                                                                pull(bayesian_prob), digits=2))) %>%
        bind_rows(annotate_df,)
      
    }
    
    if(fix_alpha==FALSE){
      pip_coloc<-final_summ_stats %>%
        filter(gene %in% genes$gene) %>%
        ggplot(aes(x=pip_pheno1,y=pip_pheno2))+
        geom_abline(data=posteriors_ci %>% filter(gene %in% genes$gene), aes(intercept=alpha, slope=beta),color = "skyblue", linewidth = 0.2, alpha = 0.25)+
        geom_abline(data=posteriors_fit %>% filter(gene %in% genes$gene), aes(intercept=alpha, slope=beta),color = "black", linewidth = 1)+
        geom_point()+
        geom_text(data=annotate_df,aes(x=x,y=y,label=text),inherit.aes = FALSE)+
        facet_wrap(~gene,ncol=1)+
        theme_bw()+
        coord_cartesian(ylim = c(0, 1), xlim = c(0,1))+
        xlab("Pheno1 PIPs")+
        ylab("Pheno2 PIPs")
    } else {
      pip_coloc<-final_summ_stats %>%
        filter(gene %in% genes$gene) %>%
        ggplot(aes(x=pip_pheno1,y=pip_pheno2))+
        geom_abline(data=posteriors_ci %>% filter(gene %in% genes$gene), aes(intercept=0, slope=beta),color = "skyblue", linewidth = 0.2, alpha = 0.25)+
        geom_abline(data=posteriors_fit %>% filter(gene %in% genes$gene), aes(intercept=0, slope=beta),color = "black", linewidth = 1)+
        geom_point()+
        geom_text(data=annotate_df,aes(x=x,y=y,label=text),inherit.aes = FALSE)+
        facet_wrap(~gene,ncol=1)+
        theme_bw()+
        coord_cartesian(ylim = c(0, 1), xlim = c(0,1))+
        xlab("Pheno1 PIPs")+
        ylab("Pheno2 PIPs")
    }
    
    plot_susie_pip<-ggarrange(reg_coloc,pip_coloc,ncol=2,labels=c("a","b"))
    
  } else {
    plot_susie_pip<-NA
  }
  
  #plot_susie_pip
  
  if(plot_susie==TRUE){
    return(list(summary_res=summary_res,
                plot_susie_pip=plot_susie_pip))
  }
  else{
    return(list(summary_res=summary_res))
  }
}

munged_sim<-c()
for(i in 1:50){
  sim_entry<-simulate_hla(n_min_alleles=opt$n_min_alleles,
                          pheno2_mean_shared_factor=opt$pheno2_mean_shared_factor,
                          pheno2_variance_shared_factor=opt$pheno2_variance_shared_factor,
                          gen_var=opt$gen_var,
                          haplotypes=hap,
                          haplotypes_lm=hap_lm_ready,
                          freq_threshold=opt$freq_threshold,
                          genes=genes,
                          all_causal_alleles=opt$all_causal_alleles,
                          all_causal_genes=opt$all_causal_genes,
                          plot_susie=opt$plot_susie,
                          continuous_trait=opt$continuous_trait,
                          n_clone=opt$n_clone,
                          susie_L=opt$susie_L,
                          fix_alpha=opt$fix_alpha)
  
  # sim_entry<-simulate_hla(n_min_alleles=10,
  #                         pheno2_mean_shared_factor=0,
  #                         pheno2_variance_shared_factor=1,
  #                         plot_predicted=TRUE,
  #                         gen_var=0.4,
  #                         haplotypes=hap,
  #                         haplotypes_lm=hap_lm_ready,
  #                         freq_threshold=0.005,
  #                         genes=genes,
  #                         all_causal_alleles=FALSE,
  #                         plot_susie=FALSE,
  #                         all_causal_genes=FALSE,
  #                         continuous_trait=TRUE,
  #                         negative_threshold=0.1,
  #                         n_clone=10,
  #                         susie_L=10,
  #                         fix_alpha=FALSE)
  
  munged_sim<-sim_entry$summary_res %>%
    mutate(iteration=paste0(opt$iteration,".",i)) %>%
    bind_rows(munged_sim,.)
  
  if(opt$plot_susie==TRUE){
   ggsave(paste0(opt$out, "susie_pip_plots_",anc,"_",opt$iteration,".",i,".pdf"),sim_entry$plot_susie)
  }
}

write_tsv(munged_sim, paste0(opt$out, "empirical_df_",anc,"_",opt$iteration, ".tsv.gz"))
