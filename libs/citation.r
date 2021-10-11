citation.text <- function(...) {
    cit = citation(...)
    class(cit) = 'list'
    return(attr(cit[[1]], 'textVersion'))
}

