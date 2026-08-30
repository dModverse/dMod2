## sanitize.R -- input validation / normalization helpers

sanitizeConditions <- function(conditions) {
  
  new <- str_replace_all(conditions, "[[:punct:]]", "_")
  new <- str_replace_all(new, "\\s+", "_")
  return(new)
  
}



sanitizePars <- function(pars = NULL, fixed = NULL) {
  
  # Convert fixed to named numeric
  if (!is.null(fixed)) fixed <- structure(as.numeric(fixed), names = names(fixed))
  
  # Convert pars to named numeric
  if (!is.null(pars)) {
    pars <- structure(as.numeric(pars), names = names(pars))
    # remove fixed from pars
    pars <- pars[setdiff(names(pars), names(fixed))]
  }
    
  
  return(list(pars = pars, fixed = fixed))
  
}




sanitizeData <- function(x, required = c("name", "time", "value"), imputed = c(sigma = NA, lloq = -Inf)) {
  
  all.names <- names(x)
  
  missing.required <- setdiff(required, all.names)
  missing.imputed <- setdiff(names(imputed), all.names)
  
  if (length(missing.required) > 0)
      stop("These mandatory columns are missing: ", paste(missing.required, collapse = ", "))
  
  if (length(missing.imputed) > 0) {
    
    for (n in missing.imputed) x[[n]] <- imputed[n]
    
  }
  
  list(data = x, columns = c(required, names(imputed)))
  
}
