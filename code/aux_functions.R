# auxiliar functions


#####################################################################
#VIF FUNCTION. From Zuur et al. 2009 http://www.highstat.com/Books/Book2/HighstatLibV10.R
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  #cat("Correlations of the variables\n\n")
  #tmp_cor <- cor(dataz,use="complete.obs")
  #print(tmp_cor)

  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)

  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}
#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}


#####################################################################
# functions to create posthocs tables
posthoc.commons <-
  function(modelo){emmeans(modelo, list(pairwise ~ commons), adjust = "tukey")}
posthoc.transh<- function(modelo){
  emmeans(modelo, list(pairwise ~ transh), adjust = "tukey")}
posthoc.int <- function(modelo){
  emmeans(modelo, list(pairwise ~ commons | transh), adjust = "tukey")}


tabla_postHoc_commons <- function(modelo){
  posthoc <- emmeans(modelo, list(pairwise ~ commons), adjust = "tukey")
  as.data.frame(posthoc$`pairwise differences of commons`) %>%
    rename(commons = 1) %>%
    kbl(digits = 4) %>%
    kable_paper("hover", full_width = F)
}

tabla_postHoc_transh <- function(modelo){
  posthoc <- emmeans(modelo, list(pairwise ~ transh), adjust = "tukey")
  as.data.frame(posthoc$`pairwise differences of transh`) %>%
    rename(transh = 1) %>%
    kbl(digits = 4) %>%
    kable_paper("hover", full_width = F)
}

tabla_postHoc_int <- function(modelo){
  posthoc <- emmeans(modelo, list(pairwise ~ commons | transh), adjust = "tukey")
  as.data.frame(posthoc$`pairwise differences of commons | transh`) %>%
    rename(commons = 1) %>%
    kbl(digits = 4) %>%
    kable_paper("hover", full_width = F)
}

#####################################################################
# functions to formatAnova
formatAnova <- function(df, yvar){
  df$variable <- yvar
  df$factor <- rownames(df)
  df %>%
    dplyr::select(fvalue = `F value`, p=`Pr(>F)`, variable, factor) %>%
    assign(paste0("anova_", yvar), ., inherits = TRUE)
}
