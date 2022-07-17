# Note-2-for-DESeq2-time-series-data-analysis

More notes on using LRT to test time-series data. Thanks for the discussion with Jie. 

    swapping the levels of time factor won’t change the LRT results, as if the time variable is a factor, LRT won’t see it as a trajectory analysis but rather a factor analysis (e.g. condition-specific difference at ANY time point). 
    
    subsetting only two time points {t0, ti} in LRT will get different numbers of DE genes at different time point i.  See example below; when only testing the {0, 60} minutes, 7 DE genes found. If including all time points in LRT, it will only find 4 DE genes. This could be the case that the likelihood ratio gets smaller when including time points with no or smaller condition-specific difference. Another possibility this may be happening is that there is a large dependence on t2 independent of the condition. When you add t2, you are adding it in the variable “time” to both the numerator and denominator of the LRT, and that may result in a smaller ratio.
    
    converting the time covariate from a factor / categorical variable to a continuous variable can get different results. The categorical variable in LRT does not consider the slope or trajectory nature, but it can detect genes that contribute a big condition-specific difference at a specific time point while the overall slope may not change. A continuous variable can consider the slope change or trajectory analysis. It’s more like a time-series analysis.  Sometimes we may need to do both. 

        library(“fission”)

        library(“DESeq2”)

        data(“fission”)

        ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)

        dim(ddsTC)

        colData(ddsTC)

        ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)

        resTC <- results(ddsTC)

        head(resTC[order(resTC$padj),], 4)

        # log2 fold change (MLE): strainmut.minute180 

        # LRT p-value: ‘~ strain + minute + strain:minute’ vs ‘~ strain + minute’ 

        # DataFrame with 4 rows and 6 columns

        # baseMean log2FoldChange     lfcSE      stat      pvalue        padj

        #          

        #   SPBC2F12.09c   174.671     -2.6567195  0.752261   97.2834 1.97415e-19 1.33453e-15

        # SPAC1002.18    444.505     -0.0509321  0.204299   56.9536 5.16955e-11 1.74731e-07

        # SPAC1002.19    336.373     -0.3927490  0.573494   43.5339 2.87980e-08 6.48916e-05

        # SPAC1002.17c   261.773     -1.1387648  0.606129   39.3158 2.05137e-07 3.46682e-04

        dim(subset(resTC, padj<=0.05))

        # 4


        # swap the levels of time factor

        colData(ddsTC)$minute = factor(colData(ddsTC)$minute, levels = c(0, 60, 120, 180, 15, 30))

        ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)

        resTC <- results(ddsTC)

        head(resTC[order(resTC$padj),], 4)

        # log2 fold change (MLE): strainmut.minute30 

        # LRT p-value: ‘~ strain + minute + strain:minute’ vs ‘~ strain + minute’ 

        # DataFrame with 4 rows and 6 columns

        # baseMean log2FoldChange     lfcSE      stat      pvalue        padj

        #          

        #   SPBC2F12.09c   174.671      -2.600469  0.634343   97.2834 1.97415e-19 1.33453e-15

        # SPAC1002.18    444.505       0.556944  0.195812   56.9536 5.16955e-11 1.74731e-07

        # SPAC1002.19    336.373       1.246269  0.525847   43.5339 2.87980e-08 6.48916e-05

        # SPAC1002.17c   261.773       1.479279  0.542824   39.3158 2.05137e-07 3.46682e-04

        dim(subset(resTC, padj<=0.05))

        # 4

        # Conclusion: no change when swapping the levels of time factor


        # change time to continous variable to test the slope difference

        colData(ddsTC)$minute = as.numeric(as.character(colData(ddsTC)$minute))

        colData(ddsTC)

        ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ strain + minute)

        resTC <- results(ddsTC)

        head(resTC[order(resTC$padj),], 4)

        dim(subset(resTC, padj<=0.05))

        # 0


        ## subset with only two time points: 0 and one of the following time points 

        for(i in c(15, 30, 60, 120, 180)){

          dds_subset=ddsTC[,ddsTC$minute %in% c(0,i)]; 

          colData(dds_subset)$minute = factor(colData(dds_subset)$minute, levels = c(0,i))

          design(dds_subset) <- as.formula(" ~ strain + minute + strain:minute")

          dds_subset <- DESeq(dds_subset, test="LRT", reduced = ~ strain + minute, quiet=T)

          res <- results(dds_subset)

          message(paste(“DE genes in time”,i,”:”,nrow(subset(res, padj<=0.05))))  

        }


        # DE genes in time 15 : 0

        # DE genes in time 30 : 1

        # DE genes in time 60 : 7

        # DE genes in time 120 : 0

        # DE genes in time 180 : 0
