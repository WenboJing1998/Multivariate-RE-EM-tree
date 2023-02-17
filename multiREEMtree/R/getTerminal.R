#' Extract the terminal nodes (leaves) of a multiREEMtree, GmultiREEMtree, or multitree object
#'
#' Extract the terminal nodes (leaves) of a multiREEMtree, GmultiREEMtree, or multitree object, and provide summary statistics for
#' each response variable in each node.
#'
#' @param object A multiREEMtree, GmultiREEMtree, or multitree object.
#'
#'
#' @return A list containing
#' \itemize{
#'  \item terminal_nodes - a character vector showing all of the terminal nodes.
#'  \item terminal_nodes_table - a table showing the number of observations in each terminal node.
#'  \item terminal_nodes_where - a vector showing the terminal node to which each observation belongs.
#'  \item terminal_nodes_summary - a matrix summarizing the mean and standard deviation of each
#' response variable in each terminal node.
#'
#' }
#'
#' @export



getTerminal <- function(object){

    if(inherits(object, c("multiREEMtree", "multitree", "GmultiREEMtree"))){

      ResponseVariables <- object$ResponseVariables

      VariableNumber <- length(ResponseVariables)

      tree <- object$Tree

      terminal_nodes <- with(tree, names(table(where))[table(where)  %in%  with(frame, n[var=="<leaf>"])])

      summary_terminal <- matrix(NA, nrow=length(terminal_nodes), ncol=2*VariableNumber+1)
      for(i in 1:length(terminal_nodes)){
        summary_terminal[i, 1] <- table(tree$where)[terminal_nodes[i]]
        summary_terminal[i, 2:(VariableNumber+1)] <-
                 apply(object$Y_original[tree$where == terminal_nodes[i],  ], 2, mean)
        summary_terminal[i, (VariableNumber+2):(2*VariableNumber+1)] <-
                 apply(object$Y_original[tree$where == terminal_nodes[i],  ], 2, sd)
      }
      rownames(summary_terminal) <- terminal_nodes
      colnames(summary_terminal) <- c("n", paste("mean.", ResponseVariables, sep=""), paste("sd.", ResponseVariables, sep=""))


      tree$frame[tree$frame$var=="<leaf>", c("n", "yval2")]

      return(list(terminal_nodes = terminal_nodes,
                  terminal_nodes_table = table(tree$where),
                  terminal_nodes_where = tree$where,
                  terminal_nodes_summary = summary_terminal))


    }else{
      stop("The object should be multiREEMtree, GmultiREEMtree, or multitree.")
    }
}

