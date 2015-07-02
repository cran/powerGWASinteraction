\name{powerGG}
\alias{powerGG}
\title{Power for GxG interactions in genetic association studies}
\description{This routine carries out (analytical, approximate) power calculations for identifying Gene-Gene interactions in Genome Wide Association Studies}
\usage{powerGG(n, power, model, caco, alpha, alpha1)}
\arguments{
\item{n}{Sample size: combined number of cases and controls. Note: exactly one of \code{n} and \code{power} should be specified.}
\item{power}{Power: targeted power. Note: exactly one of \code{n} and \code{power} should be specified.}
\item{model}{List specifying the genetic model. This list contains the following objects:
\itemize{
\item \code{prev} Prevalence of the outcome in the \emph{population}. Note that for case-only and empirical Bayes estimators to be valid,
the prevalence needs to be low.
\item \code{pGene1} Probability that the first \emph{binary} SNP is 1 (i.e. not the minor allele frequency for a three level SNP).
\item \code{pGene2} Probability that the first \emph{binary} SNP is 1 (i.e. not the minor allele frequency for a three level SNP).
\item \code{beta.LOR} Vector of length three with the odds ratios of the first genetic, second genetic, and GxG interaction effect, respectively.
\item \code{nSNP} Number of SNPs (genes) being tested.}}
\item{caco}{Fraction of the sample that are cases (default = 0.5).}
\item{alpha}{Overall (family-wise) Type 1 error (default = 0.05).}
\item{alpha1}{Significance level at which testing during the first stage (screening) takes place. If alpha1 = 1, there is no screening.}
}
\details{The routine computes power calculations for a two-stage procedure with marginal screening followed by either case-control or case-only testing.}
\value{A data frame consisting of two numbers: the power for the case-control and case-only approaches if \code{n} is specified or the required combined sample size 
for the case-control and case-only approaches if \code{power} is specified.}
\references{Kooperberg C, LeBlanc M (2008). Increasing the power of identifying gene x gene interactions in genome-wide association studies.
\emph{Genetic Epidemiology}, \bold{32}, 255-263.}
\author{Charles Kooperberg, \email{clk@fredhutch.org}}
\examples{
mod1 <- list(prev=0.05, pGene1=0.3, pGene2=0.3, beta.LOR=c(0,0,.6),nSNP=500000)
powerGG(n=10000,mod=mod1,caco=0.5,alpha=.05,alpha1=.001)
powerGG(power=0.8,mod=mod1,caco=0.5,alpha=.05,alpha1=.001)
}
\seealso{powerGG}
\keyword{methods}

