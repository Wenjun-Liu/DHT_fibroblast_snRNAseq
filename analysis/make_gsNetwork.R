make_gsNetwork <- function(normalisedScores, gsTopology,  colorBy = c("robustZ", "pvalue"),
                           foldGSname = TRUE, foldafter = 2, plotIsolated ){

    # create dummy variable to pass R CMD CHECK
    from <- to <- E <- robustZ <- NULL
    GS2Gene <- get_GSgenelist(gsTopology)
    GS2Gene <- left_join(normalisedScores, GS2Gene, by = "gs_name")
    
    GSlist <- split(GS2Gene[,c("gs_name", "gene_id")], f = GS2Gene$gs_name)
    nGS <- length(GSlist)
    GSname <- names(GSlist)
    
    w <- lapply(seq_len(nGS-1), function(x){
        lapply((x+1):nGS, function(y){
            data.frame(from = GSname[x], to = GSname[y], weight = jacIdex_func(GSlist[[x]]$gene_id, GSlist[[y]]$gene_id))
        })
    })
    
    w <- bind_rows(lapply(w, bind_rows))
    w <- dplyr::filter(w, from != to)
    
    g <- graph.data.frame(dplyr::select(w, from, to), directed = FALSE)
    g <- set_edge_attr(g, "weight", value = w$weight)
    
    GSsize <- reshape2::melt(lapply(GSlist, nrow))
    colnames(GSsize) <- c("size", "from")
    g <- set_vertex_attr(g, "size", index = GSsize$from, value = GSsize$size)
    
    if (colorBy == "robustZ"){
        GScolor <- mutate(GS2Gene, color =  ifelse(robustZ < 0, "Inhibited", "Activated"))
        g <- set_vertex_attr(g, "color", index = GScolor$gs_name, value = GScolor$color)
    }
    
    if (colorBy == "pvalue"){
        GSpvalue <- unique(GS2Gene[,c("gs_name", "pvalue")])
        g <- set_vertex_attr(g, "color", index = GSpvalue$gs_name, value = GSpvalue$pvalue)
    }
    
    if(!plotIsolated){
        removeEdge <- which(E(g)$weight == 0)
        g <-  delete_edges(g, removeEdge)
        IsolatedNode <- which(degree(g) == 0)
        g <- delete_vertices(g, IsolatedNode)
    }
    
    if(foldGSname){
        g <- set_vertex_attr(g, "name", value = vapply(V(g)$name, function(x){ifelse(length(strsplit(x, " ")[[1]]) > foldafter,
                                                                                     str_replace_nth(x, " ", "\n", foldafter),
                                                                                     x)}, character(1))) }
    g
}

get_GSgenelist <- function(gsTopology){
    GStoGene <- lapply(gsTopology, rownames)
    GStoGene <- reshape2::melt(GStoGene)
    colnames(GStoGene) <- c("gene_id", "gs_name")
    GStoGene
}


jacIdex_func <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

# This function was a modified version of the str_replace_nth function from martinctc/textworks
str_replace_nth <- function(x, pattern, replacement, n) {
    g <- gregexpr(pattern, x)[[1]][n]
    s <- scan(text = gsub("[()]", "", pattern),
              sep = "|",
              what = "",
              quiet = TRUE)
    substr(x, g, g) <- replacement[match(substr(x, g, g), s)]
    x
}