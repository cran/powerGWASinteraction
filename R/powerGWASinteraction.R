powerGWASinteraction <- function (env, b, maf, cc, nsnps, alpha1=c(1,0.1,0.01,0.001,0.0001), crit=0.05, caseonly=FALSE,
             designinfo=TRUE)
{
    # code for the approximate power calculations for identifying
    # interactions in a two-stage analysis of a genome-wide association
    # study, as described in (the appendix of)
    # C. Kooperberg, M. LeBlanc (2008) Increasing the power of identifying
    # gene x gene interactions in genome-wide association studies.
    # Genetic Epidemiology, In Press.
    # Copyright (C) 2007, 2008 Charles Kooperberg
    #
    #  This program is free software; you can redistribute it and/or modify
    #  it under the terms of the GNU General Public License as published by
    #  the Free Software Foundation; either version 2 of the License, or
    #  (at your option) any later version.
    #
    #  This program is distributed in the hope that it will be useful,
    #  but WITHOUT ANY WARRANTY; without even the implied warranty of
    #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #  GNU General Public License for more details.
    #
    #  The text of the GNU General Public License, version 2, is available
    #  as http://www.gnu.org/copyleft or by writing to the Free Software
    #  Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
    #
    # Charles Kooperberg, 12/3/2007, 10/20/2008, clk@fhcrc.org
    #------------------------------------------------------------------------
    # This revised version also works for SNP x Environment interactions.
    # All SNPs and ENvironmental factors are assumed to be binary
    # The "maf" is therefore not really the traditional minor allele frequency,
    # but it is P(SNP>0) or P(Env>0). For SNPs in HWE, we can relate the maf
    # in this code (which I refer to as "maf") to the traditional minor
    # allele frequency (which I refer to as "MAF") for dominant and recessive
    # models.
    # In particular, for a recessive model P(SNP>0)=MAF^2,
    #                       thus maf=MAF^2 and MAF=sqrt(maf)
    # For a dominant model P(SNP>0)=1-(1-MAF)^2
    #                       thus maf=1-(1-MAF)^2 and MAF=1-sqrt(1-maf)
    #
    # For gene x environment interactions only the SNP is being filtered for.
    #--------------------------------------------------------------------------
    #
    # INPUT parameters
    #
    # ==> b: parameters in a logistic model we assume
    #     logit(P(Y=1|X))=b[1]+b[2]*(SNP1>0)+b[3]*(SNP2>0)+b[4]*(SNP1>0)*(SNP2>0)
    #     or 
    #     logit(P(Y=1|X))=b[1]+b[2]*(SNP1>0)+b[3]*(ENV>0)+b[4]*(SNP1>0)*(ENV>0)
    #

    if(length(b)!=4) stop("b needs to have length 4")

    # ==> maf: c(P(SNP1>0),P(SNP2>0)) or c(P(SNP1>0),P(ENV>0))
    #      if only one minor allele frequency is specified, we assume they're equal

    if(length(maf)==1) maf <- c(maf,maf)

    if(length(maf)!=2) stop("maf needs to have length 1 or 2")

    if(min(maf)<=0 || max(maf)>=1) stop("maf needs to be strictly between 0 and 1")

    # ==> cc:  c(#cases, #controls)
    #      if only one sample size is specified we assume #cases = #controls

    if(length(cc)==1)  cc <- c(cc,cc)

    if(length(cc)!=2) stop("cc needs to have length 1 or 2")

    if(min(cc)<5) stop("cc needs to be at least 5")

    #
    # ==> nsnps: number of SNPs in the GWA

    nsnps <- round(nsnps)

    if(nsnps<5) stop("nsnps needs to be at least 5")

    # ==> alpha1: marginal significance level at which testing at the first stage
    #         is being done. 

    if(min(alpha1) <=0 )stop("alpha1 needs to be strictly larger than 0")
    if(max(alpha1) >1 )stop("alpha1 needs to be smaller or equal to 1")
    if(min(alpha1)<1e-5)warning("for some alpha1 you test very few SNPs at the second stage, the independence\nassumption gets questionable, and the results unreliable")

    #
    # ==> crit: overall Type 1 error - note, if crit>1 this is the expected number
    #       of false positives - thus it is OK two use, for example, crit=3
 
    if(crit<=0) stop("crit needs to be strictly larger than 0")

    #
    # ==> env: are we looking for SNP1 x SNP2 (gene-gene) interactions (FALSE) 
    #      SNP1 x  ENV (gene-environment) interactions (TRUE)
    #
    # ==> caseonly: also provide power if the analyses are carried out using
    #     a case-only analysis. This typically assumes that the two factors
    #     are independent in the population. This is almost certainly NOT true
    #     for gene x gene interactions, it may be true for gene x environment
    #     interactions in some special cases, e.g. a randomized treatment assignment

    if(caseonly && env==FALSE)warning("a case-only analysis for gene x gene interactions is almost certainly wrong")
    #
    # ==> designinfo: should some info about the design be printed?
    #
    # OUTPUT values
    #
    # if caseonly = FALSE
    #     two values - the first one gives the power of identifying the correct
    #     (SNP1 x SNP2) or (SNP1x ENV) interaction (equation (9) in the paper),
    #     the second gives
    #     the expected number of false positives
    # if caseonly = TRUE
    #     four values - the first two are the same as when caseonly = FALSE, the
    #     second two are the same using a case-only test
    #
    zpower <- function (n1,n2,p1,p2, alpha)
    {
        # this routine computes the power of finding a difference
        # between two binomial samples - sample i (i=1,2) with sample
        # size ni and success pi. The Type 1 error (2 sided) is alpha.
        # it's a slight generalization of the R function
        # power.prop.test() for different sample sizes.
        zalpha <- qnorm(alpha/2, lower.tail=FALSE)
        q1 <- 1-p1
        q2 <- 1-p2
        pbar <- (p1+p2)/2
        qbar <- (q1+q2)/2
        r <- n2/n1
        top <- sqrt(r*(p1-p2)*(p1-p2)*n1)-zalpha*sqrt((r+1)*pbar*qbar)
        bottom <- sqrt(r*p1*q1+p2*q2)
        pnorm(top/bottom)
    }
    # the design matrix for the two different SNPs and their interaction
    desmat <- matrix(c(1,1,1,1,0,1,0,1,0,0,1,1,0,0,0,1),ncol=4)

    # the fitted values, exponentiated
    p00 <- exp(desmat%*%as.matrix(b))

    # P(Y=0|X,b)
    q00 <- 1/(1+p00)

    # P(Y=1|X,b)
    p00 <- 1 - q00

    # now convert to conditional probabilities for case-control sampling
    q00 <- q00 * (1-maf[1])^(1-desmat[,2])*maf[1]^desmat[,2]
    p00 <- p00 * (1-maf[1])^(1-desmat[,2])*maf[1]^desmat[,2]
    q00 <- q00 * (1-maf[2])^(1-desmat[,3])*maf[2]^desmat[,3]
    p00 <- p00 * (1-maf[2])^(1-desmat[,3])*maf[2]^desmat[,3]
    q00 <- q00/sum(q00)
    p00 <- p00/sum(p00)

    lalpha1 <- length(alpha1)

    # the power that SNP1 will be selected at stage 1
    z1 <- zpower(cc[1], cc[2], p00[2] + p00[4], q00[2] + q00[4], alpha1)

    # the power that SNP2 will be selected at stage 1
    if(env==FALSE){
       z2 <- zpower(cc[1], cc[2], p00[3] + p00[4], q00[3] + q00[4], alpha1)
    }
    else {
       z2 <- 1
    }

    # the power that both SNPs will be selected at stage 1
    zz <- z1*z2

    # expected standard error for the coefficient of interaction
    v <- sqrt((sum(1/q00)/cc[2]+sum(1/p00)/cc[1]))

    # v-caseonly
    if(caseonly)vcaseonly <- sqrt(sum(1/p00)/cc[1])

    # normalized coefficient of the interaction parameter
    v <- abs(b[4]/v)
    if(caseonly)vcaseonly <- abs(b[4]/vcaseonly)

    # ztable
    if(designinfo){
       kr <- p00/(p00+q00)
       kk <- matrix(0,ncol=3,nrow=3)
       kk[1:2,1] <- kr[1:2]
       kk[1:2,2] <- kr[3:4]
       n00 <- c((1-maf[1])*(1-maf[2]),maf[1]*(1-maf[2]),maf[2]*(1-maf[1]),maf[1]*maf[2])
       kk[3,1] <- sum(kr[1:2]*n00[1:2])/sum(n00[1:2])
       kk[3,2] <- sum(kr[3:4]*n00[3:4])/sum(n00[3:4])
       kk[1,3] <- sum(kr[c(1,3)]*n00[c(1,3)])/sum(n00[c(1,3)])
       kk[2,3] <- sum(kr[c(2,4)]*n00[c(2,4)])/sum(n00[c(2,4)])
       kk[3,3] <- sum(kr*n00)/sum(n00)
       orh <- kk[,2]*(1-kk[,1])/(kk[,1]*(1-kk[,2]))
       orv <- kk[2,]*(1-kk[1,])/(kk[1,]*(1-kk[2,]))
       kkx <- cbind(kk,orh)
       kkx <- rbind(kkx,c(orv,NA))
       kkx <- round(kkx,2)
       kk <- data.frame(kkx[,1:2],"|",kkx[,3],"|",kkx[,4])
       if(env==FALSE){
         names(kk) <- c("SNP2=0","SNP2=1","|","SNP2=any","|","odds ratio")
         row.names(kk) <- c("SNP1=0","SNP1=1","SNP1=any","odds ratio")
       }
       else{
         names(kk) <- c("ENV=0","ENV=1","|","ENV=any","|","odds ratio")
         row.names(kk) <- c("SNP=0","SNP=1","SNP=any","odds ratio")
       }
       
       cat("\nFraction of subjects that are a case in each cell of your design\n\n")
       print(kk)
       kk <- kk[1:3,1:4]
       kk[,1] <- c(n00[1:2],n00[1]+n00[2])
       kk[,2] <- c(n00[3:4],n00[3]+n00[4])
       kk[,4] <- kk[,1]+kk[,2]
       kk[,-3] <- round(kk[,-3]*(cc[1]+cc[2]))
       cat("\nNumber of subjects that are in each cell of your design\n\n")
       print(kk)
       cat("\n")
    }
    

    # the number of other (not interesting) SNPs
    nother <- nsnps-2+env

    # the plausible range of the number of other SNPs being selected
    u1 <- u2 <- u3 <- u4 <- rep(0,lalpha1)
    for(i in 1:lalpha1){
       qb <- (qbinom(.0000001,nother,alpha1[i])-10):(qbinom(1-.0000001,nother,alpha1[i])+10)
       qb <- qb[qb>=0]

        # the following lines yield vectors with different values for different qb
        # the probabilities of seeing that many other SNPs
        db <- dbinom(qb,nother,alpha1)

        # tails - just in case the probabilities don
        ub <- 1-sum(db)
        db[c(1,length(db))] <- db[c(1,length(db))]+ub/2
 
        # number of interactions involved in the selected SNPs
        if(env==FALSE){
           nu <- (qb+2)*(qb+1)/2
        }
        else{
           nu <- qb+1
        }

        # Bonferonni alpha for each interaction
        alpha2 <- -qnorm(0.5*crit/nu)
    
        # power of identifying the "correct" interaction, if SNP1 and SNP2 survived stage 1
        p1 <- pnorm(-alpha2 - v) + 1 - pnorm(alpha2 - v)
        if(caseonly)p1caseonly <- pnorm(-alpha2 - vcaseonly) + 1 - pnorm(alpha2 - vcaseonly)
  
    # integrating out: zz=P(snps survived stage 1), 
    #                  p1=P(identified|snps survived stage 1 and qb identified)
    #                  db=p(qb identified)
       u1[i] <- sum(zz[i]*p1*db)
       if(caseonly)u3[i] <- sum(zz[i]*p1caseonly*db)

    # second term - when SNP1 and SNP2 are not significant, power of any interaction
    # thirs term  - when SNP1 and SNP2 are  significant, power of another interaction
    #
       u2[i] <- (1-zz[i])*crit+zz[i]*crit*sum(db*(nu-1)/nu*(1-p1))
       if(caseonly) u4[i] <- (1-zz[i])*crit+zz[i]*crit*sum(db*(nu-1)/nu*(1-p1caseonly))
    }
    u5 <- round(alpha1*nother+z1+z2,2)
    if(env) u5 <- u5-1
    if(caseonly==FALSE){
       output <- data.frame(alpha1=alpha1,exp.stage1=u5,power=u1,false.pos=u2)
    }
    else{
       output <- data.frame(alpha1=alpha1,exp.stage1=u5,power.cc=u1,false.pos.cc=u2,
                           power.co=u3,false.pos.co=u4)
    }
    output
}
