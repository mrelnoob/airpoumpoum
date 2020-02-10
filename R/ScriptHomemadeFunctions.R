####################
### MY FUNCTIONS ###
####################


##### Timeline and history #####
#
# 17/01/2020: Creation of the airpoumpoum package
# 17/01/2020: Creation of the functions: se.m() ; IC95 ; rar.rm_v2() ; comb_mod() ;
# 21/01/2020: Creation of the function: sensiNMDS() ;
# 23/01/2020: Creation of the function: posthoc.NMDS() ;
# 27/01/2020: Creation of the function: super_rich.group() ;
# 28/01/2020: Uploading of the package on GitHub --> To install package use: devtools::install_github("mrelnoob/airpoumpoum")
# 29/01/2020: Creation of the function: super_distriplot();
# 31/01/2020: Updates on the super_distriplot() function.
# 04/02/2020: Updates on the super_rich.group() function.





##### ----------------------------------- se.m ----------------------------------- #####

#' Standard Error of a mean
#'
#' @description Computes the standard error of a sample mean. Useful to add error bars on plots or to compute intervals of confidence.
#' @param x a vector
#'
#' @return The standard error of the input vector's mean.
#' @import stats
#' @export
#'
#' @examples myvector <- rnorm(n = 100, mean = 48, sd = 16)
#' se.m(x = myvector)
se.m<-function(x)
{
  sqrt(var(x,na.rm=T))/sqrt(length(x))
}





##### ----------------------------------- IC95 ----------------------------------- #####

#' Ninety-five percent confidence interval
#'
#' @description Computes the 95 percent confidence interval of a mean. Only trustworthy for normally distributed data (or data with n > 30).
#' @param x A vector.
#'
#' @return The range of a confidence interval.
#' @import stats
#' @export
#'
#' @examples myvector <- rnorm(n = 100, mean = 48, sd = 16)
#' IC95(x = myvector)
#' # Gives a result of 2.992398. It means that the mean of myvector
#' # has 95% chance to be comprised in an interval of +/- 2.992398.
IC95 <- function(x){
  qt(p = 0.975, df = (length(x)-1))*sqrt(var(x, na.rm=T))/sqrt(length(x))
}





##### ----------------------------------- rar.rm_v2 ----------------------------------- #####

#' Uncommon species remover
#'
#' @description Enables removing uncommon species from a contingency table (matrix or dataframe) as a function of their abundance
#' (number of obervations per site) and frequency (number of site where the species has been observed). Useful to only keep common and structuring
#' species in community data.
#' @param tableau.AD A contingency table (typically with species as columns and sites/releves as lines).
#' @param NB_REL The minimum number of sites in which the species needs to be oberved to be kept.
#' @param ABUN The minimum number of observations within site required for the species to be kept.
#'
#' @return A contingency table with only the species that are found at least ABUN times in NB_REL different sites.
#' @export
#'
#' @examples #Nope
rar.rm_V2 <- function (tableau.AD, NB_REL = 1, ABUN = 0.00000001) # Fonction de suppression des espèces peu fréquentes
{
  tableau.PA <- data.frame(apply(tableau.AD, c(1, 2), function(x) if (x >= ABUN) 1  else 0)) # Veut dire: créer un data.frame dans lequel tu appliques (apply)
  # la fonction x (si x >= à l'abondance spécifiée, tu inscris 1, sinon tu inscris 0) sur les lignes et colonnes (i.e. c(1, 2)) de tableau.AD
  occur <- apply(tableau.PA, 2, sum) # Veut dire: dans l'objet "occur", tu fais la somme des valeurs par colonne (i.e. 2) du tableau.AD
  tableau.AD.wr <- tableau.AD[, occur >= NB_REL] # Veut dire: créé une table (tableau.AD.wr) ne gardant de tableau.AD que les colonnes (espèces) qui
  # apparaissent dans NB_REL sites (lignes) ou plus
} # NOTE IMPORTANTE: il utiliser en entrée un tableau de contingence, donc avec comme seules variables les espèces !!!





##### ----------------------------------- comb_mod ----------------------------------- #####

#' Creation of 2-by-2 comparison vectors from factor
#'
#' @description Creates a list of two vectors that contains all possible combinations of modalities of the input factor. Useful to perform 2-by-2 comparisons
#' (e.g. post-hoc tests). Therefore, it gives 2 vectors of length L = (n*n-1)/2, where n is the number of different modalities (i.e. factor levels)
#' contained in FAC. For example, if the input factor contains 3 modalities A, B and C, comb_mod will give AB, AC and BC (in that order).
#' @param FAC A factor containing the modalities that need to be combined.
#'
#' @note  The comb_mod function always returns a list of 2 vectors.
#'
#' @return A list of two vectors.
#' @export
#'
#' @examples # aa1 <- comb_mod(FAC = mydata$treatment)$mod1
#' # aa2 <- comb_mod(FAC = mydata$treatment)$mod2
#' # Le $mod1 ou 2 permet de dire lequel des 2 vecteurs d'output tu souhaites dans l'objet créé
#' # (puisque la fonction comb_mod créé toujours une liste de 2 vecteurs).
#' # Sans cette précision, la fonction créerait 2 objets identiques (aa1 et aa2):
#' # 2 liste contenant les 2 vecteurs
comb_mod <- function(FAC){ # FAC étant un facteur (e.g. variable/colonne catégorielle d'un dataframe)
  mod1 <- NULL
  mod2 <- NULL
  n_mod1 <- length(levels(FAC))-1
  for (i in 1:(n_mod1)) # Pour chaque i (i.e. chaque modalité de FAC) moins une...
  {
    mod1 <- c(mod1,rep(levels(FAC)[i],n_mod1)) # Créé un vecteur mod1 dans lequel tu répètes les modalités i, n_mod1 fois !
    n_mod1 <- n_mod1-1
    mod2 <- c(mod2,levels(FAC)[(i+1):length(levels(FAC))]) # Créé un vecteur mod2 en listant toutes les modalités de FAC mais en décalé (i.e. +1), sauf
    # celle de i.
  }
  output <- list(mod1 = mod1,mod2 = mod2)
  return(output)
}




##### ----------------------------------- sensiNMDS ----------------------------------- #####

#' Sensitivity analysis for the selection of species for NMDSs
#'
#' @description This function performs a sensitivity analysis to choose the combination of species to keep in an ecological dataset to optimize the
#' stress in NMDS analyses. This function works with the \code{\link[airpoumpoum:rar.rm_V2]{rar.rm_V2}} function.As a consequence, sensiNMDS() performs a
#' NMDS for each possible combination of species removal in terms of abundance (of a species in a given site) and number of colonized sites (by that species),
#' as in \code{\link[airpoumpoum:rar.rm_V2]{rar.rm_V2}}.
#'
#' @details For time-efficiency reasons, the functions currently only tests for a predefined limited number of combinations of possible abundances (i.e.
#' accounted for by the \strong{ABUN} argument in the \code{\link[airpoumpoum:rar.rm_V2]{rar.rm_V2}} function) with a maximal value set to 40 (i.e. a species
#' has to be observed 40 times to be kept). For additional values, please contact me. \cr
#' The NMDSs are performed using the \code{\link[vegan:metaMDS]{metaMDS}} function in the \strong{vegan} package.
#'
#' @param MYDATA A contingency table (a purely numeric species-site matrix).
#' @param TREATMENT A factor with a length matching the number of samples (i.e. the number of lines of MYDATA). Typically, a factor representing a given
#' treatment associated to the different samples/sites. In fact, this argument is not used in the current version of the function, but is required because
#' it is called by \emph{hidden} and \emph{inactive} parts of the source code (useful for simulation studies). If you do not have such a factor, you can
#' easily generate one (see Example section below).
#' @param NB_DIM The desired number of dimensions kept in the NMDS analyses (>1). Usually, 2 or 3. The number of dimensions kept depends on the quality of
#' the data.
#'
#' @return A graph showing the behaviour of NMDSs' stress for each combination of species removal.
#' @import stats vegan graphics
#' @export
#'
#' @examples  ## For a simulated matrix with 20 lines (sites) and
#' ## 30 columns (species), split in 2 treatments:
#' library(vegan)
#' data(dune)
#' mytreatment <- gl(2, 10)
#' sensiNMDS(MYDATA = dune, TREATMENT = mytreatment, NB_DIM = 2)
sensiNMDS <- function(MYDATA, TREATMENT, NB_DIM){
  val_occ <- c(1:nrow(MYDATA)) # Tout d'abord, on cree des vecteurs possedant la longueur de toutes les combinaisons de NMDS que l'on souhaite tester :
  val_ab <- c(1,2,5,10,15,20,30,40) # Ces 4 lignes de code servent a cela.
  val_occ <- rep(val_occ,length(val_ab)) # Cela revient donc a creer des vecteurs qui comportent toutes les combinaisons possibles d'occurrences et
  # d'abondance
  val_ab <- rep(val_ab,each = nrow(MYDATA))
  n_sp <- NULL
  n_empty_row <- NULL
  R2_test <- NULL
  p_test <- NULL
  stress_test <- NULL

  Facteur <- as.factor(TREATMENT) # Le traitement/la variable dont l'effet sur les compositions m'interesse

  for(i in 1:length(val_occ))
  {
    rel_test <- rar.rm_V2(MYDATA, NB_REL = val_occ[i], ABUN = val_ab[i]) #le releve avec retrait des especes
    if(is.null(dim(rel_test))==TRUE){ #s'il n'y a qu'une espece on considere comme nul et pas utilise
      n_sp[i] <- 0
      n_empty_row[i] <- 0
      R2_test[i] <- NA
      p_test[i] <- NA
      stress_test[i] <- NA

    }else{
      n_sp[i] <- ncol(rel_test) # pour avoir le nombre d'especes
      sumab <- apply(rel_test,1,sum) # pour avoir le nombre de releves vides :
      n_empty_row[i] <- length(sumab)-length(sumab[sumab>0]) # Ici, on doit toujours etre à 0 (sinon, c'est qu'il y a un releve vide)

      if(n_sp[i]>2 & n_empty_row[i]==0){ # pour faire la NMDS, on attends d'avoir au moins 2 especes et 0 releves vides
        NMDS <- metaMDS(rel_test, distance = "bray", k = NB_DIM, trace = FALSE, try = 10, trymax = 50) # Parametres a modifier si besoin (distance, nombre de
        # dimensions etc.)
        stress_test[i] <- NMDS$stress
        env_an <- adonis(rel_test~Facteur)
        R2_test[i] <- env_an$aov.tab$R2[1]
        p_test[i] <- env_an$aov.tab$`Pr(>F)`[1]

      }else{
        R2_test[i] <- NA
        p_test[i] <- NA
        stress_test[i] <- NA
      }
    }
  }

  # Plots de visualisation
  super_stressplot <- plot(x = val_ab,val_occ, cex = 20*stress_test, main = "Stress values for various combinations of species removal",
                           col = ifelse(test = stress_test < 0.05, yes = 2, no = ifelse(test = stress_test<0.2, yes = "orange", no = 1)),
                           pch = ifelse(test = stress_test<0.2, yes = 19, no = 1), # indique que si stress > 20, alors les cercles sont vides (au lieu
                           # d'etre pleins)
                           ylab = "Required number of colonized sites", xlab = "Required number of species observations")
  val_stress <- c(min(stress_test,na.rm=T),mean(stress_test,na.rm=T),max(stress_test,na.rm=T))
  legend(x = "topright", legend = c(paste0("Min stress = ", round(val_stress[1], 2)), paste0("Mean stress = ", round(val_stress[2], 2)),
                                    paste0("Max stress = ", round(val_stress[3], 2)), "stress < 0.05 ~ excellent", "stress < 0.2 ~ bon"), # Cree 5 elements
         # de legende (les 3 de val_stress, et les 2 ecris ensuite)
         pt.cex =  20*c(val_stress,0.05,0.2),  pch=c(1,1,1,19,19), col=c(1,1,1,2,"orange"), bty = "n")

  print(super_stressplot) # Tout ?a pour sortir cette petite figure !

  #plot(x = val_ab,val_occ, cex = 4*R2_test, main = "taille des points relative au R2 de envfit",
  #     col = ifelse(R2_test>0.5, 2, 1), pch = ifelse(R2_test>0.1, 19, 1)) # Les cercles sont pleins si le R? > 10%
  #val_R2 <- c(min(R2_test,na.rm=T),mean(R2_test,na.rm=T),max(R2_test,na.rm=T))
  #legend("topright",legend = c(paste0("R2 = ",round(val_R2,2)), "R2 > 0.5"),
  #       pt.cex = 4*c(val_R2,0.5), pch=c(1,1,1,19), col=c(1,1,1,2), bty = "n")

  #plot(x = val_ab,val_occ, cex = 1, main = "Taille en fonction des p-values",
  #     col = ifelse(test = p_test < 0.05, yes = 2, no = ifelse(test = p_test < 0.08, yes = "orange", no = 1)),
  #     pch = ifelse(test = p_test < 0.08, yes = 16, no = 1))
  #val_pval <- c(min(p_test,na.rm=T),mean(p_test,na.rm=T),max(p_test,na.rm=T))
  #legend("topright",legend = c(paste0("p = ",round(val_pval,2)), "p < 0.08", "p < 0.05"),
  #       pt.cex = 1, pch=c(1,1,1,16,16), col=c(1,1,1,"orange",2), bty = "n")

  #text(val_ab,val_occ, labels = p_test,cex = 0.8) # Si l'on veut pouvoir lire les p-values sur le graphe !
}
# NOTE IMPORTANTE: dans cette ?tude, envfit() donne des p-values au comportement erratique, donc on utilise une ADONIS a la place !





##### ----------------------------------- posthoc.NMDS() ----------------------------------- #####

#' Two-by-two post-hoc tests for NMDSs using PERMANOVA
#'
#' @description This function performs as many NMDSs and PERMANOVAs as there are possible combinations between the input factor levels. In other words, the
#' function will compute (using the \code{\link[airpoumpoum:comb_mod]{comb_mod}} function) all possible combinations of the modalities (levels) contained in
#' \emph{FAC}, and run a separate NMDS and associated PERMANOVA between each pair of those modalities (i.e. treatments) and the community data in \emph{REL}. \cr
#' The post-hoc tests are corrected using \strong{Holm-Bonferroni} method to control the family-wise error rate. By default, the functions only performs
#' NMDSs and PERMANOVAs using \strong{Bray-Curtis} distances, and keep only 2 dimensions in each ordination.
#'
#' @note This function should only be used when a preliminary PERMANOVA detected a significant effect of \emph{FAC} on the centroid position of \emph{REL},
#' and after ensuring that there was no \emph{multivariate spread} in the data groups (see \code{\link[vegan:adonis]{adonis}}).
#'
#' @param FAC A factor containing the different modalities/treatments we want to combine and test.
#' @param REL A contingency table. Typically, a species-site matrix.
#'
#' @return A table containing the p-values and corrected p-values for each PERMANOVA test.
#' @import vegan stats
#' @export
#'
#' @examples ## For a simulated matrix with 20 lines (sites) and
#' ## 30 columns (species), split in 2 treatments:
#' library(vegan)
#' data(dune)
#' mytreatment <- gl(2, 10)
#' aa <- posthoc.NMDS(REL = dune, FAC = mytreatment)
#' aa
#' ## Here, since there is only 2 levels in the factor, there
#' ## will only be one combination and thus, one test.
#' ## It's a silly example allright.
posthoc.NMDS <- function(REL, FAC){
mod1 <- comb_mod(FAC)$mod1
mod2 <- comb_mod(FAC)$mod2

t <- REL
t$treatment <- FAC

pposthoc <- NULL
Fposthoc <- NULL
for (i in 1:length(mod1)) # Pour chaque i (allant de 1 ? la longueur de mod1, ici = 3)
{
    ttt <- t[t$treatment == mod1[i] | t$treatment == mod2[i], ] # Cr?? une table ne gardant que les lignes associ?es aux modalit?s i de
  # mod1 ET de mod2 (e.g. pour i = 1, ?a gardera donc les lignes des 2 premi?res modalit?s de traitement par ordre alphab?tique; pour i = 2, la 1?re
  # modalit? et la 3?me etc.)
  NMDS_2by2 <- metaMDS(comm = ttt[,1:ncol(ttt)-1], distance = "bray", k = 2, try = 10, trymax = 100, trace = FALSE) # Fait une NMDS avec ttt (sans la colonne
  # "treatment", qui n'a rien ? faire dans un tableau de contingence)

  tperm_NMDS2by2 <- adonis(formula = ttt[,1:ncol(ttt)-1]~ttt$treatment, permutations = 999, method = "bray") # Fait le test de permutation sur la NMDS
  pposthoc[i] <- tperm_NMDS2by2$aov.tab$`Pr(>F)`[1] # Remplie la ligne i de pposthoc avec la p-value du test de permutation
  Fposthoc[i] <- tperm_NMDS2by2$aov.tab$F.Model[1] # Idem pour la statistique F.
}

pposthoc_adj <- p.adjust(pposthoc, method = "holm")
data.frame(comp=paste(mod1,mod2), pposthoc, pposthoc_adj, Fposthoc)
}





##### ----------------------------------- super_rich.group() ----------------------------------- #####

#' Computation of species richness per groups
#'
#' @description This functions gives the specific richness (mean and standard deviation) for each \emph{group} of data contained in the argument
#' \code{MYGROUP} (e.g. treatments, types of site, sites etc.). Additionally, \code{super_rich.group()} also returns a dataframe containing all the species
#' data for each \emph{group}, sorted out by decreasing order of abundance. \cr
#' For instance, if \code{MYDATA} is a dataframe containing the abundance (or presence/absence) of 100 species in 30 sites
#' belonging to three groups of equal size named A, B and C (detailed in \code{MYGROUP}); then \code{super_rich.group()} will return the mean
#' \strong{specific richness} of the A, B and C sites, and create 3 different dataframes:
#'
#' * One containing all the species of the 10 sites of group A (ordered by decreasing abundance);
#' * One containing all the species of the 10 sites of group B (ordered by decreasing abundance);
#' * One containing all the species of the 10 sites of group C (ordered by decreasing abundance).
#'
#' Note that the last line of each of these dataframes is \emph{the sum of species abundances}: i.e. the line used to order the dataframes' columns!
#'
#'
#' @param MYGROUP A factor with a length matching the number of rows of \code{MYDATA} (typically, a factor highlighting to which group each line belongs).
#' @param MYDATA A contingency table (typically, a matrix with species as columns and sites/releves as lines).
#'
#' @return Two list objects: one named \code{richness_group} and another named \code{topabun_group}. The former contains the specific richness for each group,
#' the latter contains the subset of data of each group (the aforementioned dataframes).
#' @export
#'
#' @examples library(vegan)
#' data("dune")
#' data("dune.env")
#' # Creates a list containing the species richness per group.
#' aa <- super_rich.group(MYDATA = dune, MYGROUP = dune.env$Use)[[1]]
#' # Creates a list containing the subset dataframes for each group.
#' bb <- super_rich.group(MYDATA = dune, MYGROUP = dune.env$Use)[[2]]
#' # You can also create a list of lists by excluding the
#' # double-brackets: i.e. [[]]:
#' cc <- super_rich.group(MYDATA = dune, MYGROUP = dune.env$Use)
#' aa <- cc[[1]]
#' bb <- cc[[2]]
super_rich.group <- function(MYDATA, MYGROUP){

  MYDATA$treatment <- MYGROUP
  richness_group <- NULL
  topabun_group <- NULL

  for(i in 1:length(unique(MYGROUP))){


    fff <- MYDATA[MYDATA$treatment == unique(MYGROUP)[i],]
    fff1 <- fff[, 1:ncol(fff)-1]
    fff2 <- rar.rm_V2(fff1, NB_REL = 1, ABUN = 1)
    nbl <- nrow(fff2)+1
    fff3 <- fff2
    fff3[nbl,] <- apply(fff3, 2, sum)
    fff4 <- fff3[, order(fff3[nbl,], decreasing = T)]

    topabun_group[[i]] <- data.frame(fff4) # NOTE IMPORTANTE: sur ces 2 lignes de code, en fonction du nombre de [], on ne fait pas du tout la m?me chose !

    tt <- data.frame(apply(X = fff2, c(1,2), function(x) if (x > 0) 1 else 0))
    tt$richness <- apply(X = tt, MARGIN = 1, FUN = sum)


    richness_group$treatment[i] <- as.character(unique(MYGROUP))[i]
    richness_group$mean_richness[i] <- mean(tt$richness)
    richness_group$sd_richness[i] <- sd(tt$richness)

  }

  richness_group <- as.data.frame(richness_group)
  output <- list(richness_group, topabun_group)

  return(output)
}





##### ----------------------------------- super_distriplot() ----------------------------------- #####

#' Creation of a multi-plot of variables distribution
#'
#' @description The \code{super_distriplot()} function plots, in a single viewing window, the distribution of \code{MYVARIABLES} for all the \code{GROUPS} to
#' which they belong. It is very useful to quickly visualize if the \emph{among groups} \strong{normality assumption} is respected or not. The function
#' additionally plots the curve of the \emph{probability density function} of each sampled population, and that of a \strong{Normal} distribution.
#'
#' @note For graphical convenience, the function cannot plot too many distributions in a single window. For that reason, only a maximum of 4 different
#' groups and 5 different variables is allowed as inputs in the function (giving a window with 20 plots). \cr
#' If ou wish to plot more variables and/or groups, please divide your data and run the function several times on each subset. \cr
#' Also, if you wish to plot a single variable divided in different groups, you need to specify \code{MYVARIABLES} that your \emph{single vector} (i.e.
#' your variable) is a dataframe. Otherwise, \code{super_distriplot()} will not be able to write a proper title to the histograms. If you want a proper
#' title, then use a synthax as follows: \code{... MYVARIABLES = as.data.frame(myvector) ...} or, for a \emph{single column}, \code{... MYVARIABLES =
#' as.data.frame(mydata$myvariable) ...} (see also 'examples').
#'
#' @param MYVARIABLES A numeric vector or an array (up to 5) of numeric vectors (typically, variables as columns of a dataframe or matrix).
#' @param GROUPS A factor whose length matches the length of the variables of interest. The maximal number of groups/treatments to divide the observations
#' is limited to 4 (see 'note').
#' @param breaks A value, a vector of values, an algorithm etc. to set the width of the bars in the histograms (for more details on the possible input
#' values, please refer to the help page of \code{\link[graphics:hist]{hist}}).
#'
#' @return A vinwing window with as many histograms as there are combinations between the input variables and groups (up to 20).
#' @import stats graphics
#' @export
#'
#' @examples library(vegan)
#' data("dune")
#' data("dune.env")
#' super_distriplot(MYVARIABLES = dune[,3:6], GROUPS = dune.env[,4], breaks = 5)
#' # For a single variable (e.g. Bellis perennis):
#' super_distriplot(MYVARIABLES = as.data.frame(dune$Bellpere), GROUPS = dune.env[,4], breaks = 5)
super_distriplot <- function(MYVARIABLES, GROUPS, breaks){

  mydataset <- as.data.frame(MYVARIABLES)
  mydataset$treatment <- GROUPS
  mydataset$treatment <- droplevels(mydataset$treatment)

  nb_gr <- length(unique(mydataset$treatment))
  ifelse((ncol(mydataset)-1) == 1, nb_var <- 1, nb_var <- ncol(mydataset)-1)

  par(mfrow=c(nb_gr,nb_var))
  for (j in 1:nb_gr){
    for (i in 1:nb_var){

      tt <- mydataset[mydataset$treatment == levels(mydataset$treatment)[j],i]
      tt <- as.data.frame(tt)

      hist(tt[,1], breaks = breaks, probability = TRUE, xlab = NULL,
           border = "pink", col = "coral3", main = paste("Distribution of", colnames(mydataset)[i], "\nfor", levels(mydataset$treatment)[j]))
      lines(density(tt[,1]), col = "gray 12", lwd = 2)
      f <- function(x){
        dnorm(x = x, mean = mean(tt[,1]), sd = sd(tt[,1]))}
      curve(f, add = TRUE, col = "red", lwd = 2, lty = 2)

    }
  }
}
