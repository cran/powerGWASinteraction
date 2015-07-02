\name{powerGE}
\alias{powerGE}
\title{Power for GxE interactions in genetic association studies}
\description{This routine carries out (analytical, approximate) power calculations for identifying Gene-Environment interactions in Genome Wide Association Studies}
\usage{powerGE(n, power, model, caco, alpha, alpha1, maintain.alpha)}
\arguments{
\item{n}{Sample size: combined number of cases and controls. Note: exactly one of \code{n} and \code{power} should be specified.}
\item{power}{Power: targeted power. Note: exactly one of \code{n} and \code{power} should be specified.}
\item{model}{List specifying the genetic model. This list contains the following objects:
\itemize{
\item \code{prev} Prevalence of the outcome in the \emph{population}. Note that for case-only and empirical Bayes estimators to be valid,
the prevalence needs to be low.
\item \code{pGene} Probability that a \emph{binary} SNP is 1 (i.e. not the minor allele frequency for a three level SNP).
\item \code{pEnv} Frequency of the binary environmental variable.
\item \code{orGE} Odds ratio between the binary SNP and binary environmental variable.
\item \code{beta.LOR} Vector of length three with the odds ratios of the genetic, environmental, and GxE interaction effect, respectively.
\item \code{nSNP} Number of SNPs (genes) being tested.}}
\item{caco}{Fraction of the sample that are cases (default = 0.5).}
\item{alpha}{Overall (family-wise) Type 1 error (default = 0.05).}
\item{alpha1}{Significance level at which testing during the first stage (screening) takes place. If alpha1 = 1, there is no screening.}
\item{maintain.alpha}{Some combinations of screening and GxE testing methods do not maintain the proper Type 1 error. Default is \code{True}: combinations that do not maintain the Type 1 error are not computed. If \code{maintain.alpha} is \code{False} all combinations are computed.}
}

\details{The routine computes power for a variety of two-stage procedures. Five different screening procedures are used:
\itemize{
\item \bold{No screening} All SNPs are tested for interaction
\item \bold{Marginal screening} Only SNPs that are marginally significant at level alpha1 are screened for interaction. See Kooperberg and LeBlanc (2010).
\item \bold{Correlation screening} Only SNPs that are, combined over all cases and controls, associated with the environmental variable at level alpha 1
are screened for interaction. See Murcray et al. (2012).
\item \bold{Cocktail screening} SNPs are screened on the most significant of marginal and correlation screening. See Hsu et al. (2012).
\item \bold{Chi-square screening} SNPs are screened using a chi-square combination of correlation and marginal screening. See Gauderman et al. (2013).}

After screening, the SNPs that pass the screen can be tested using
\itemize{
\item \bold{Case-control} The standard case-control estimator.
\item \bold{Case-only} The case-only estimator.
\item \bold{Empirical Bayes} The empirical Bayes estimator of Mukherjee and Chatterjee (2010).}

If screening took place using the correlation or chi-square screening, the Type 1 error won't be
maintained if the final GxE testing is carried out using either the case-only or empirical Bayes estimator.
See Dai et al. (2012). The cocktail screening maintains the Type 1 family wise error rate, since only those
SNPs that pass on to the second stage using marginal screening will use the case-only or empirical Bayes estimator, the
SNPs that pass on to the second stage using correlation screening will always use the case-control estimator.

When SNP and environment are correlated in the population (i.e. \code{model$orGE} does not equal 1) the case-only estimator does not maintain the Type 1 error.
The empirical Bayes estimator may also have a moderately inflated Type 1 error. When the disease is common either the case-only
estimator or the empirical Bayes estimator also may not estimate the GxE interaction.

Power calculations are described in Kooperberg, Dai, and Hsu (2014). Briefly, for a given genetic model we compute the expected p-values for all
screening statistics. We then use a normal approximation to compute the probability that this SNP passes the screening (e.g., if \code{alpha1} equaled this expected p-value this probability would be exactly 0.5), and combine this with power calculations for the second stage of GxE testing.}
\value{A list with three components.
\item{power}{A 5x3 matrix with estimated power for all testing approaches, only if \code{n} was specified.}
\item{samplesize}{A 5x3 matrix with required sample sizes for all testing approaches, only if \code{power} was specified.}
\item{expected.p}{A 5x3 matrix with the expected p value for the SNP to pass screening. This p-value depends on the sample size, but not on the second stage testing.}
\item{prob.select}{A 5x3 matrix with the probability that the interacting SNP would pass the screening stage. This probability depends on the sample size, but not on the second stage testing.}}
\references{
Dai JY, Kooperberg C, LeBlanc M, Prentice RL (2012).
Two-stage testing procedures with  independent filtering for genome-wide
gene-environment interaction. \emph{Biometrika}, \bold{99}, 929-944.

Gauderman WJ, Zhang P, Morrison JL, Lewinger JP (2013).
Finding novel genes by testing GxE interactions in a genome-wide association study.
\emph{Genetic Epidemiology}, \bold{37}, 603-613.

Hsu L, Jiao S, Dai JY, Hutter C, Peters U, Kooperberg C (2012).
Powerful cocktail methods for detecting genome-wide gene-environment interaction.
\emph{Genetic Epidemiology}, \bold{36}, 183-194. 

Kooperberg C, Dai, JY, Hsu L (2014).
Two-stage procedures for the identification of gene x
environment and gene x gene interactions in
genome-wide association studies. \emph{To appear}.

Kooperberg C, LeBlanc ML (2008). Increasing the power of
identifying gene x gene interactions in genome-wide
association studies. \emph{Genetic Epidemiology}, \bold{32}, 255-263.

Mukherjee B, Chatterjee N (2008). Exploiting gene-environment inde-
pendence for analysis of case-control studies: an empirical
Bayes-type shrinkage estimator to trade-off between bias and efficiency
\emph{Biometrics}, \bold{64}, 685-694.

Murcray CE, Lewinger JP, Gauderman WJ (2009).
Gene-environment interaction in genome-wide association studies.
\emph{American Journalk of Epidemiology}, \bold{169}, 219-226.
}


\author{
Li Hsu \email{lih@fredhutch.org} and
Charles Kooperberg \email{clk@fredhutch.org}.
}


\examples{
mod1 <- list(prev=0.01,pGene=0.2,pEnv=0.2,beta.LOR=log(c(1.0,1.2,1.4)),orGE=1.2,nSNP=10^6)
results <- powerGE(n=20000, model=mod1,alpha1=.01)
print(results)

mod2 <- list(prev=0.01,pGene=0.2,pEnv=0.2,beta.LOR=log(c(1.0,1.0,1.4)),orGE=1,nSNP=10^6)
results <- powerGE(power=0.8, model=mod2,alpha1=.01)
print(results)
}
\seealso{powerGG}
