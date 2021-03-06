% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ScriptHomemadeFunctions.R
\name{comb_mod}
\alias{comb_mod}
\title{Creation of 2-by-2 comparison vectors from factor}
\usage{
comb_mod(FAC)
}
\arguments{
\item{FAC}{A factor containing the modalities that need to be combined.}
}
\value{
A list of two vectors.
}
\description{
Creates a list of two vectors that contains all possible combinations of modalities of the input factor. Useful to perform 2-by-2 comparisons
(e.g. post-hoc tests). Therefore, it gives 2 vectors of length L = (n*n-1)/2, where n is the number of different modalities (i.e. factor levels)
contained in FAC. For example, if the input factor contains 3 modalities A, B and C, comb_mod will give AB, AC and BC (in that order).
}
\note{
The comb_mod function always returns a list of 2 vectors.
}
\examples{
# aa1 <- comb_mod(FAC = mydata$treatment)$mod1
# aa2 <- comb_mod(FAC = mydata$treatment)$mod2
# Le $mod1 ou 2 permet de dire lequel des 2 vecteurs d'output tu souhaites dans l'objet crÃ©Ã©
# (puisque la fonction comb_mod crÃ©Ã© toujours une liste de 2 vecteurs).
# Sans cette prÃ©cision, la fonction crÃ©erait 2 objets identiques (aa1 et aa2):
# 2 liste contenant les 2 vecteurs
}
