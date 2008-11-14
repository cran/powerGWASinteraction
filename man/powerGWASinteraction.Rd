\name{powerGWASinteraction}
\alias{powerGWASinteraction}
\title{Power calculations for identifying interactions in GWAS studies}
\description{This function carries out approximate
power calculations for identifying SNP x SNP and SNP x environment interactions in 
genome-wide association (GWAS) studies. It assumes a two-stage analysis, where only
SNPs that are significant at a marginal significance level alpha1 are investgigated
for interactions, and a binary environmental covariate (for SNP x environment interactions.}
\usage{powerGWASinteraction(env, b, maf, cc, nsnps, alpha1, crit, caseonly, designinfo)}
\arguments{
\item{env}{Are you interested in SNP x SNP interactions (FALSE) or SNP x Environment interactions (TRUE).}
\item{b}{Vector of length four: parameters in a logistic regression model. We assume:\cr
logit(P(Y=1|SNPs))=b[1]+b[2]*(SNP1>0)+b[3]*(SNP2>0)+b[4]*(SNP1>0)*(SNP2>0)\cr
or\cr
logit(P(Y=1|SNPs,ENV))=b[1]+b[2]*(SNP1>0)+b[3]*(ENV>0)+b[4]*(SNP1>0)*(ENV>0)\cr
as the model.}
\item{maf}{Scalar or vector of length two, with the probabilities
that the relevant SNP(s) and environmental factors are 1; i.e.: c(P(SNP1>0),P(SNP2>0)) or c(P(SNP1>0),P(ENV>0)).
If only one maf (minor allele frequency) is specified, we assume they're equal.
All SNPs and ENvironmental factors are assumed to be binary.
The "maf" is therefore not really the traditional minor allele frequency,
but it is P(SNP>0). For SNPs in HWE, we can relate the maf
in this code (which we refer to as "maf") to the traditional minor
allele frequency (which we refer to as "MAF") for dominant and recessive models.
In particular, for a recessive model P(SNP>0)=MAF*MAF, thus maf=MAF*MAF and
MAF=sqrt(maf).  For a dominant model P(SNP>0)=1-(1-MAF)*(1-MAF)
    thus maf=1-(1-MAF)*(1-MAF) and MAF=1-sqrt(1-maf).}
\item{cc}{Vector of length 2, number of cases and controls, respectively.
If only one is specified the number of cases and number of controls is assumed
to be the same.}
\item{nsnps}{Number of SNPs in the GWAS.}
\item{alpha1}{Marginal significance level that the SNPs (but not the environmental factor) have to pass to be tested for interactions. This can be a vector, in which case the calculations will be carried out for each alpha1. Reasonable numbers are 1 (all SNPs are tested, thus this is a one-phase analysis)
and numbers in
the order 0.01 to 0.0001, i.e. between 1 in 100 and 1 in 10000 of the SNPs are tested for interactions.}
\item{crit}{Multiple comparisons corrected overall significance level (Family-wise
error) at which interactions are tested. Traditionally this would be 0.05 (which is the default), but it can also be a number larger than 1, in which case
it becomes the number of expected false positives.}
\item{caseonly}{Also provide power if the analyses are carried out using
      a case-only analysis. This typically assumes that the two factors
      are independent in the population. This is almost certainly NOT true
      for gene x gene interactions, it may be true for gene x environment
      interactions in some special cases, e.g. a randomized treatment assignment. Default is FALSE.}
\item{designinfo}{Should some informative information about the selected design be printed. Examining this info protects you against choosing real weird designs. Default is TRUE}
}
\value{A data frame. One row for each alpha1. Column 1: alpha1; column 2: expected number of SNPs that are significant at level alpha1 and will make it to phase 2; column 3:
the power of identifying the correct
        (SNP1xSNP2) or (SNP1xENV) interaction (equation (9) in Kooperberg and LeBlanc (2008);
column 4:
the expected number of false positives;
column 5 and 6: as column 3 and 4, for a case-only analysis, if caseonly=TRUE.}
\references{Kooperberg C, LeBlanc M (2008). Increasing the power of identifying gene x gene interactions in genome-wide association studies.
\emph{Genetic Epidemiology}, \bold{32}, 255-263.}
\author{Charles Kooperberg, \email{clk@fhcrc.org}}
\examples{
powerGWASinteraction(env=FALSE,b=c(-2,0,0,.6),maf=0.3,cc=4000,nsnps=500000,crit=3)
#
# power for SNP x SNP interactions, for a pure epistatic effect
# with OR exp(.6) where both SNPs have P(SNP>0)=0.3, there are
# 4000 cases and 4000 controls, and 500000 SNPs). Provide results
# for an expected number of 3 false psoitives.
#
powerGWASinteraction(env=TRUE,b=c(-2,0,.5,.5),maf=c(0.4,0.5),cc=c(2000,3000),nsnps=500000,crit=0.05,caseonly=TRUE)
#
# power for a SNP x Environment interaction, where the Environmental
# factor has an effect that is enhanced by a SNP. P(SNP>0)=0.4, and
# P(Env>0)=0.5, 2000 cases and 3000 contols, 500000 SNPs, testing
# at a FWER of 0.05, and providing details on a case-only analysis
# as well.}
\keyword{methods}

