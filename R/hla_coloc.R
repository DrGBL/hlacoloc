#' Perform HLA colocalization
#'
#' @param pheno1 Dataframe of HLA allele associations for phenotype 1. Needs to contain the following columns:
#' * `Name`: the name of the HLA allele in the IMGT-HLA format.
#' * `z`: the z-score of the HLA allele association.
#' * `beta`: the beta (effect sizes) of the HLA allele association (required if z-score not provided).
#' * `se`: the standard error (required if z-score not provided).
#' * `N`: the sample size.
#' @param pheno1R A dataframe or matrix of correlation coefficient (R) or the alleles in pheno1. The rows and columns need to be in the same order as the alleles in pheno1.
#' @param is_cohort_ld_pheno1 Whether the LD matrix for the pheno1 cohort is from the same cohort (TRUE) or from an external reference (FALSE).
#' @param pheno2 Dataframe of HLA allele associations for phenotype 2. Needs to be in the same format as pheno1, including rows in same order.
#' @param pheno2R A dataframe or matrix of correlation coefficient (R) or the alleles in pheno2. The rows and columns need to be in the same order as the alleles in pheno2.
#' @param is_cohort_ld_pheno2 Whether the LD matrix for the pheno2 cohort is from the same cohort (TRUE) or from an external reference (FALSE).
#' @param max_iter_susieR Maximum number of iterations for susieR (default=100).
#' @param plot_susie Whether to plot the results (TRUE/FALSE). This is automatically set to FALSE if the beta of the HLA allele association tests are not available.
#' @param plot_assoc Whether to plot the HLA allele association results. This is automatically set to FALSE if the beta of the HLA allele association tests are not available, or if plot_susie is set to false.
#' @param negative_threshold Minimum susieR posterior inclusion probability needed for both phenotypoes in order to check for colocalization using stanR (default=0.001).
#' @param susie_L Maximum number of alleles with non-zero effect in the susieR model (default=10).
#' @param n_min_alleles Minimum number of alleles required at a gene in order to attempt HLA colocalization at that gene. Genes with less than this threshold will be excluded from the analyses (default=10).
#' @param pheno1_name Name of pheno1, used for plotting only (default="Pheno1").
#' @param pheno2_name Name of pheno2, used for plotting only (default="Pheno2").
#'
#' @importFrom graphics "text"
#' @importFrom stats "gaussian"
#' @importFrom utils "head"
#' @importFrom rlang ".data"
#' @return List with a dataframe containing the colocalization result, and a plot of the results (optional).
#' @export
hla_coloc<-function(pheno1,pheno1R,is_cohort_ld_pheno1=FALSE,
                    pheno2,pheno2R,is_cohort_ld_pheno2=FALSE,
                    max_iter_susieR=100,plot_susie=TRUE,plot_assoc=TRUE,
                    negative_threshold=0.001,
                    susie_L=10,n_min_alleles=10,
                    pheno1_name="Pheno1",
                    pheno2_name="Pheno2"){

  if( length(which(!pheno1$Name %in% pheno2$Name)) > 0 |
      length(which(!pheno2$Name %in% pheno1$Name)) > 0){
    print("Not same alleles in both cohorts.")
  } else if(all(pheno1$Name==pheno2$Name)==FALSE) {
    print("Alleles not in same order in both cohorts.")
  } else {
    if("z" %in% colnames(pheno1)){
      susie_pheno1<-susieR::susie_rss(R=as.matrix(pheno1R),
                                      z=pheno1$z,
                                      n=max(pheno1$N),
                                      estimate_residual_variance=is_cohort_ld_pheno1,
                                      max_iter=max_iter_susieR,
                                      L=susie_L)
    } else {
      susie_pheno1<-susieR::susie_rss(R=as.matrix(pheno1R),
                                      bhat=pheno1$beta,
                                      shat=pheno1$se,
                                      n=max(pheno1$N),
                                      estimate_residual_variance=is_cohort_ld_pheno1,
                                      max_iter=max_iter_susieR,
                                      L=susie_L)
    }

    #susie pheno2
    if("z" %in% colnames(pheno2)){
      susie_pheno2<-susieR::susie_rss(R=as.matrix(pheno2R),
                                      z=pheno2$z,
                                      n=max(pheno2$N),
                                      estimate_residual_variance=is_cohort_ld_pheno2,
                                      max_iter=max_iter_susieR,
                                      L=susie_L)
    } else {
      susie_pheno2<-susieR::susie_rss(R=as.matrix(pheno2R),
                                      bhat=pheno2$beta,
                                      shat=pheno2$se,
                                      n=max(pheno2$N),
                                      estimate_residual_variance=is_cohort_ld_pheno2,
                                      max_iter=max_iter_susieR,
                                      L=susie_L)
    }

    if("beta" %in% colnames(pheno1) & "beta" %in% colnames(pheno2)){
      full_final<-dplyr::inner_join(pheno1 %>%
                                      dplyr::rename(beta1=.data$beta,
                                                    se1=.data$se,
                                                    N1=.data$N),
                                    pheno2 %>%
                                      dplyr::rename(beta2=.data$beta,
                                                    se2=.data$se,
                                                    N2=.data$N)) %>%
        dplyr::mutate(pip_pheno1=susie_pheno1$pip) %>%
        dplyr::mutate(pip_pheno2=susie_pheno2$pip) %>%
        dplyr::mutate(coloc=.data$pip_pheno1*.data$pip_pheno2) %>%
        dplyr::mutate(gene=stringr::str_extract(.data$Name,"[A-Z0-9-]*\\*")) %>%
        dplyr::mutate(gene=stringr::str_replace(.data$gene,"\\*","")) %>%
        dplyr::group_by(.data$gene) %>%
        dplyr::filter(dplyr::n()>=n_min_alleles) %>%
        dplyr::ungroup()
    } else {
      full_final<-data.frame(Name=pheno1$Name) %>%
        dplyr::mutate(pip_pheno1=susie_pheno1$pip) %>%
        dplyr::mutate(pip_pheno2=susie_pheno2$pip) %>%
        dplyr::mutate(coloc=.data$pip_pheno1*.data$pip_pheno2) %>%
        dplyr::mutate(gene=stringr::str_extract(.data$Name,"[A-Z0-9-]*\\*")) %>%
        dplyr::mutate(gene=stringr::str_replace(.data$gene,"\\*","")) %>%
        dplyr::group_by(.data$gene) %>%
        dplyr::filter(dplyr::n()>=n_min_alleles) %>%
        dplyr::ungroup()
    }



    full_final_summary<-full_final %>%
      dplyr::group_by(.data$gene) %>%
      dplyr::summarize(susie_coloc_prob = 1-prod(1-.data$coloc))




    #get the bayesian correlation probability
    #specifically this is the probability that the correlation term is not zero.
    bayes_corr_df<-c()
    posteriors_fit<-c()
    posteriors_ci<-c()
    for(g in (full_final %>% dplyr::pull(.data$gene) %>% unique())){
      dat_tmp<-full_final %>% dplyr::filter(.data$gene==g)
      if(length(which(dat_tmp$pip_pheno1 > negative_threshold)) > 0 &
         length(which(dat_tmp$pip_pheno2 > negative_threshold)) > 0){
        bayes_lm<-rstanarm::stan_glm(pip_pheno2~pip_pheno1,data=dat_tmp,family = gaussian(link = "identity"))

        posteriors_ci<-as.data.frame(as.matrix(bayes_lm)) %>%
          dplyr::mutate(gene=g) %>%
          dplyr::rename(alpha=.data$`(Intercept)`) %>%
          dplyr::rename(beta=.data$pip_pheno1) %>%
          head(n=100) %>%
          dplyr::bind_rows(posteriors_ci,.)

        posteriors_fit<-data.frame(alpha=bayes_lm$coefficients[1],
                                   beta=bayes_lm$coefficients[2],
                                   gene=g) %>%
          dplyr::bind_rows(posteriors_fit,.)

        bayes_corr_df<-data.frame(gene=g,
                                  bayes_pd=bayestestR::p_direction(bayes_lm)$pd[2],
                                  direction_of_correlation=ifelse(bayes_lm$coefficients[2] >= 0,
                                                                  "Correct",
                                                                  "Incorrect")) %>%
          dplyr::bind_rows(bayes_corr_df,.)
      } else {
        bayes_corr_df<-data.frame(gene=g,
                                  bayes_pd=0.5,
                                  direction_of_correlation="Not applicable") %>%
          dplyr::bind_rows(bayes_corr_df,.)
      }
    }

    #build a summary table

    full_final_summary<-full_final_summary %>%
      dplyr::left_join(.,bayes_corr_df) %>%
      dplyr::mutate(hla_colocalization_probability=ifelse(is.na(.data$susie_coloc_prob) | is.na(.data$bayes_pd),
                                                    NA,
                                                    .data$susie_coloc_prob*(.data$bayes_pd-0.5)*2)) %>%
      dplyr::arrange(dplyr::desc(.data$hla_colocalization_probability),dplyr::desc(.data$susie_coloc_prob))



    if(plot_susie==TRUE){
      if(plot_assoc==TRUE & "beta" %in% colnames(pheno1) & "beta" %in% colnames(pheno2)){
        df_reg<-data.frame(gene=full_final %>% dplyr::pull(.data$gene) %>% unique(),
                           eq=NA,
                           min_x=NA,
                           max_x=NA,
                           min_y=NA,
                           max_y=NA)
        for(i in 1:nrow(df_reg)){
          gene_tmp<-(df_reg %>% dplyr::pull(.data$gene))[i]
          mod<-lm(beta2~beta1-1, data=full_final %>% dplyr::filter(.data$gene == gene_tmp))
          df_reg$eq[i]<-as.character(as.expression(substitute(italic(y) == b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                                              list(b = format(unname(coef(mod)[1]), digits = 2),
                                                                   r2 = format(summary(mod)$r.squared, digits = 2)))))
          df_reg$min_x[i]<-min(full_final %>%
                                 dplyr::filter(.data$gene == gene_tmp) %>%
                                 dplyr::pull(.data$beta1))
          df_reg$max_x[i]<-max(full_final %>%
                                 dplyr::filter(.data$gene == gene_tmp) %>%
                                 dplyr::pull(.data$beta1))
          df_reg$min_y[i]<-min(df_reg$min_x[i]*coef(mod)[1],min(full_final %>%
                                                                  dplyr::filter(.data$gene == gene_tmp) %>%
                                                                  dplyr::pull(.data$beta2)))
          df_reg$max_y[i]<-max(df_reg$max_x[i]*coef(mod)[1],max(full_final %>%
                                                                  dplyr::filter(.data$gene == gene_tmp) %>%
                                                                  dplyr::pull(.data$beta2)))
        }

        reg_coloc<-full_final %>%
          ggplot2::ggplot(ggplot2::aes(x=.data$beta1,y=.data$beta2))+
          ggplot2::geom_point()+
          ggplot2::geom_smooth(method="lm",formula=y~x-1)+
          ggplot2::geom_text(ggplot2::aes(label=eq, x=0.8*min_x+0.2*max_x, y=0.1*min_y+0.9*max_y),data=df_reg,parse = TRUE, size=3)+
          ggplot2::facet_wrap(~gene, scales="free",ncol=1)+
          ggplot2::theme_bw()+
          ggplot2::xlab(paste0(pheno1_name," Betas"))+
          ggplot2::ylab(paste0(pheno2_name," Betas"))
      }

      annotate_df<-c()
      for(g in full_final %>% dplyr::pull(.data$gene) %>% unique()){

        annotate_df<-data.frame(gene=g,
                                x=0.25,
                                y=0.75,
                                text_tmp=paste0("P_coloc: ", format(full_final_summary %>%
                                                                  dplyr::filter(.data$gene==g) %>%
                                                                  dplyr::pull(.data$hla_colocalization_probability), digits=2))) %>%
          dplyr::bind_rows(annotate_df,)


      }

      annotate_df<-annotate_df %>%
        dplyr::mutate(text_tmp=ifelse(is.na(.data$text_tmp),0,.data$text_tmp)) %>%
        dplyr::mutate(text_tmp=scales::label_percent()(.data$text_tmp)) %>%
        dplyr::mutate(text=paste0("P_coloc: ", .data$text_tmp)) %>%
        dplyr::select(-.data$text_tmp) %>%
        dplyr::mutate(text=paste0(stringr::str_extract(.data$text,"P_coloc: [0-9]*\\.[0-9]"),"%"))

      pip_coloc<-full_final %>%
        ggplot2::ggplot(ggplot2::aes(x=.data$pip_pheno1,y=.data$pip_pheno2))+
        ggplot2::geom_abline(data=posteriors_ci, ggplot2::aes(intercept=alpha, slope=beta),color = "skyblue", linewidth = 0.2, alpha = 0.25)+
        ggplot2::geom_abline(data=posteriors_fit, ggplot2::aes(intercept=alpha, slope=beta),color = "black", linewidth = 1)+
        ggplot2::geom_point()+
        ggplot2::geom_text(data=annotate_df,ggplot2::aes(x=x,y=y,label=text),inherit.aes = FALSE)+
        ggplot2::facet_wrap(~gene,ncol=1)+
        ggplot2::theme_bw()+
        ggplot2::coord_cartesian(ylim = c(0, 1), xlim = c(0,1))+
        ggplot2::xlab(paste0(pheno1_name," PIPs"))+
        ggplot2::ylab(paste0(pheno2_name," PIPs"))

      if(plot_assoc==TRUE & "beta" %in% colnames(pheno1) & "beta" %in% colnames(pheno2)){
        plot_susie_pip<-ggpubr::ggarrange(reg_coloc,pip_coloc,ncol=2,labels=c("a","b"))
      } else {
        plot_susie_pip<-pip_coloc
      }

    } else {
      plot_susie_pip<-NA
    }

    return(list(hla_colocalization=full_final_summary,
                plot=plot_susie_pip))

  }
}
