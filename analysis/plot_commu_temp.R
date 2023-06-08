plot_community_temp <- function(
        normalisedScores, gsTopology, gsAnnotation = NULL, colorBy = "community", 
        communityMethod = c(
            "louvain", "walktrap", "spinglass", "leading_eigen", 
            "edge_betweenness", "fast_greedy", "label_prop", "leiden"
        ),
        foldGSname = TRUE, foldafter = 2, labelFun = .rm_prefix,
        layout = c(
            "fr", "dh", "gem", "graphopt", "kk", "lgl", "mds", "sugiyama"
        ),
        markCommunity = "ellipse", markAlpha = 0.2, color_lg_title = NULL,
        edgeAlpha = 0.8, scale_edgeWidth = c(0.5, 3), edgeLegend = FALSE, 
        scale_nodeSize = c(3,6), nodeShape = 16,  lb_size = 3,
        lb_color = "black", plotIsolated = FALSE, ...
){
    
    name <- data <- community <- category <- category_n <- weight <- color <-
        size <- Community <- x <- y <- NULL
    communityMethod <- match.arg(communityMethod)
    communityMethod <- paste0("cluster_", communityMethod)
    layout <- match.arg(layout)
    
    ## check if input has required columns
    cols <- colnames(normalisedScores)
    colorBy <- match.arg(colorBy, c(cols, "community"))
    stopifnot("gs_name" %in% cols)
    
    # Make sure the gs topologies are a named list with at least two elements
    stopifnot(length(names(gsTopology)) == length(gsTopology))
    if (length(unique(normalisedScores$gs_name)) < 2) 
        stop("At least 2 gene-sets are required for a network plot")
    gsTopology <- gsTopology[names(gsTopology) %in% normalisedScores$gs_name]
    
    if (is.null(gsAnnotation)) {
        ## if gene-set annotation info is not provided, use built-in KEGG 
        ## pathway annotations
        gsAnnotation_df <- c()
        data("gsAnnotation_df", envir = environment())
        gsAnnotation <- gsAnnotation_df
    } else {
        ## if user provided gene-set annotation, gs_name and category column 
        ## must be in the df
        if (!all(c("gs_name", "category") %in% colnames(gsAnnotation)))
            stop("gsAnnotation must include gs_name and category columns") 
    }
    # create igraph object
    if (colorBy != "community") {
        g <- .make_gsNetwork(
            normalisedScores, gsTopology, colorBy = colorBy,  
            plotIsolated = plotIsolated, labelFun = NULL
        )
    } else {
        g <- .make_gsNetwork(
            normalisedScores, gsTopology, colorBy = NULL,  
            plotIsolated = plotIsolated, NULL
        )
    }
    
    # perform community detection
    comm_method <- get(communityMethod, envir = rlang::ns_env("igraph"))
    comm_result <- comm_method(g, resolution = 0.3)
    comm_df <- data.frame(
        gs_name = comm_result$names, community = comm_result$membership
    )
    
    if (length(intersect(names(gsTopology), gsAnnotation$gs_name)) == 0) {
        warning(
            "Gene-set annotation does not match with topology provided. ",
            "Communities won't be annotated"
        )
        g <- set_vertex_attr(
            g, "Community", index = comm_df$gs_name, 
            value = as.character(comm_df$community)
        )
        silent_commu <- TRUE
    } else {
        ## left_join community detection results to pathway annotation
        comm_df <- left_join(comm_df, gsAnnotation, by = "gs_name")
        
        ## Find the category with highest occurrence for each community
        comm_summary <- group_by(comm_df, community, category)
        
        comm_summary <- summarise(
            comm_summary, category_n = dplyr::n(), .groups = "keep"
        )
        comm_summary <- group_by(comm_summary, community)
        comm_summary <- dplyr::filter(
            comm_summary, category_n == max(category_n)
        )
        
        comm_summary <- drop_na(comm_summary)
        
        ## if there's a tie between two categories for a given community, 
        ## category names are pasted together
        comm_summary <- mutate(
            comm_summary, category = paste(category, collapse  = " &\n")
        )
        comm_summary <- ungroup(comm_summary)
        comm_summary <- distinct(comm_summary, community, category)
        comm_summary <- left_join(
            comm_df[, c("gs_name", "community")], comm_summary, by = "community"
        )
        
        
        
        g <- set_vertex_attr(
            g, "Community", index = comm_summary$gs_name, 
            value = comm_summary$category
        )
        silent_commu <- FALSE
    }
    
    if (foldGSname) {
        nm <-  .str_replace_nth(
            V(g)$name, pattern = " ", replacement = "\n", n = foldafter
        )
        g <- set_vertex_attr(g, "name", value = as.character(nm) )
    }
    
    # Find xy coordinators of nodes based on the layout chosen
    layout_method <- get(
        paste("layout_with_", layout, sep = ""), 
        envir = rlang::ns_env("igraph")
    )
    xy <- layout_method(g)
    
    ## Tidy up node labels if a function has been passed to this argument
    if (!is.null(labelFun)) {
        ## This allows for stadard label replacement, but also for users to
        ## provide more complex methods of string manipulation
        stopifnot(is(labelFun, "function"))
        nm <- vertex_attr(g, "name")
        new_nm <- labelFun(nm)
        stopifnot(length(nm) == length(new_nm))
        g <- set_vertex_attr(g, "name", seq_along(nm), new_nm)
    }
    
    # plot network edges
    pl <- ggraph(g, layout = "manual", x = xy[,1], y = xy[,2]) +
        geom_edge_link(
            alpha = edgeAlpha, aes(width = weight), colour = 'darkgrey'
        ) +
        scale_edge_width_continuous(range = scale_edgeWidth, guide = "none")
    
    # plot node points
    if (!is.null(colorBy)) {
        color <- ifelse(colorBy == "community", "Community", "color")
        pl <- pl + geom_node_point(
            aes(color = !!sym(color), size = size), shape = nodeShape, 
            stroke = 0.5
        ) +
            scale_size(range =  scale_nodeSize, guide = "none")  +
            labs(colour = color_lg_title)
    } else {
        pl <- pl +
            geom_node_point(aes(size = size), shape = nodeShape,stroke = 0.5) +
            scale_size(range =  scale_nodeSize, guide = "none")
    }
    
    if (!is.null(markCommunity)) {
        mark_method <- get(
            paste("geom_mark_", markCommunity, sep = ""), 
            envir = rlang::ns_env("ggforce")
        )
        if (silent_commu) {
            pl <- pl + mark_method(
                aes(x = x, y = y, fill = Community), alpha = markAlpha, ...
            )
        } else {
            pl <- pl + mark_method(
                aes(x = x, y = y, fill = Community, label = Community), 
                color = NA, alpha = markAlpha, ...
            )
        }
    }
    
    pl + geom_node_text(
        aes(label = name), size = lb_size, repel = TRUE, colour = lb_color
    ) 
    
}

.make_gsNetwork <- function(
        normalisedScores, gsTopology, colorBy = NULL, plotIsolated, labelFun
){
    
    # create dummy variable to pass R CMD CHECK
    from <- to <- NULL
    GS2Gene <- .get_GSgenelist(gsTopology)
    GS2Gene <- left_join(
        normalisedScores, GS2Gene, by = "gs_name", multiple = "all"
        ## Setting multiple = "all" requires dplyr > 1.1.0
    )
    
    GSlist <- split(GS2Gene[,c("gs_name", "entrezid")], f = GS2Gene$gs_name)
    nGS <- length(GSlist)
    GSname <- names(GSlist)
    
    w <- lapply(
        seq_len(nGS - 1),
        function(x){
            lapply(
                seq(x + 1, nGS, by = 1),
                function(y) {
                    w <- .jaccard(GSlist[[x]]$entrezid, GSlist[[y]]$entrezid)
                    data.frame(from = GSname[x], to = GSname[y], weight = w)
                }
            )
        }
    )
    
    w <- bind_rows(lapply(w, bind_rows))
    w <- dplyr::filter(w, from != to)
    GSsize <- melt(lapply(GSlist, nrow))
    colnames(GSsize) <- c("size", "from")
    
    g <- graph.data.frame(dplyr::select(w, from, to), directed = FALSE)
    g <- set_edge_attr(g, "weight", value = w$weight)
    g <- set_vertex_attr(g, "size", index = GSsize$from, value = GSsize$size)
    
    if (!is.null(colorBy)) {
        GSvalue <- unique(GS2Gene[,c("gs_name", colorBy)])
        g <- set_vertex_attr(
            g, "color", index = GSvalue$gs_name, value = GSvalue[[colorBy]]
        )
    }
    
    if (!plotIsolated) {
        removeEdge <- which(E(g)$weight == 0)
        g <-  delete_edges(g, removeEdge)
        IsolatedNode <- which(degree(g) == 0)
        g <- delete_vertices(g, IsolatedNode)
    }
    
    if (!is.null(labelFun)) {
        ## This allows for stadard label replacement, but also for users to
        ## provide more complex methods of string manipulation
        stopifnot(is(labelFun, "function"))
        nm <- vertex_attr(g, "name")
        new_nm <- labelFun(nm)
        stopifnot(length(nm) == length(new_nm))
        g <- set_vertex_attr(g, "name", seq_along(nm), new_nm)
    }
    
    g
}

#' @importFrom reshape2 melt
#' @keywords internal
.get_GSgenelist <- function(gsTopology, mapEntrezID = NULL){
    GStoGene <- lapply(gsTopology, rownames)
    GStoGene <- reshape2::melt(GStoGene)
    colnames(GStoGene) <- c("entrezid", "gs_name")
    
    if (all(c("entrezid","mapTo") %in% colnames(mapEntrezID))){
        if (any(GStoGene$entrezid %in% mapEntrezID$entrezid)){
            left_join(
                GStoGene, mapEntrezID[,c("entrezid","mapTo")], by = "entrezid",
                multiple = "all"
            )
        } else {
            warning("None of the Entrez IDs in mapEntrezID mapped to gsTopology.")
            return(GStoGene)
        }
    } else {
        return(GStoGene)
    }
}

#' @keywords internal
.jaccard <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    cmn <- intersect(x, y)
    universe <- union(x, y)
    length(cmn) / length(universe)
}

#' @importFrom stringr str_split
#' @keywords internal
.str_replace_nth <- function(x, pattern, replacement, n) {
    x_list <- str_split(x, pattern = pattern)
    fold <- vapply(x_list, length, integer(1)) > n
    x_list[fold] <- lapply(
        x_list[fold],
        function(x){
            index <- seq(n, length(x), by = n)
            x[index] <- paste(x[index], replacement, sep = "")
            not_index <- setdiff(seq_along(x), index)
            x[not_index] <- paste(x[not_index], " ", sep = "")
            x <- paste(x, collapse = "")
        }
    )
    x_list[!fold] <- lapply(x_list[!fold], paste, collapse = " ")
    x_list <- unlist(x_list)
    trails <- paste0(replacement, "$")
    gsub(trails, "", x_list) # Remove any trailing line breaks
}

#' @keywords internal
.rm_prefix <- function(x) {
    gsub("^[a-z]+\\.", "", x)
}
