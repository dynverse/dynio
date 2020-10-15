#' Copy trajectory from h5ad
#'
#' @param anndata AnnData object
#' @param name_prefix A prefix to add to each created object in the anndata.
#'
#' @importFrom dplyr bind_cols select mutate rename_all
#' @importFrom purrr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom reticulate py_set_item
#' @import dynwrap
#' @importFrom Matrix Matrix
#'
#' @export
from_h5ad <- function(anndata, name_prefix = "") {
  # anndata <- anndata::read_h5ad("/home/rcannood/workspace/jnj/scpipelineviash/src/ti/slingshot/output.h5ad")
  # traj <- readr::read_rds("/home/rcannood/workspace/jnj/scpipelineviash/src/ti/slingshot/output.rds")
  # gimp <- readr::read_rds("/home/rcannood/workspace/jnj/scpipelineviash/src/ti/slingshot/gimp.rds")
  # name_prefix <- "slingshot_"

  obs_names <- python_builtins$list(anndata$obs_names)
  var_names <- python_builtins$list(anndata$var_names)

  # read from X
  expression <- as(anndata$X, "CsparseMatrix")
  rownames(expression) <- obs_names
  colnames(expression) <- var_names
  counts <- expression
  counts@x <- 2^counts@x - 1

  # read from uns slots
  # as.vector removes extra attributes
  milestone_ids <- anndata$uns[paste0(name_prefix, "milestone_ids")] %>% as.vector
  milestone_network <- anndata$uns[paste0(name_prefix, "milestone_network")]
  dimred_ids <- anndata$uns[paste0(name_prefix, "dimred_ids")]

  dimred_segments <- anndata$uns[paste0(name_prefix, "dimred_segments")]
  dimred_segment_progressions <- dimred_segments %>% select(-!!dimred_ids)
  dimred_segment_points <- dimred_segments %>% select(!!dimred_ids) %>% as.matrix

  dimred_milestones <- anndata$uns[paste0(name_prefix, "dimred_milestones")]
  rownames(dimred_milestones) <- milestone_ids
  colnames(dimred_milestones) <- dimred_ids

  # read from obs slots
  progr_names <- paste0(name_prefix, c("from", "to", "percentage"))
  progressions <- anndata$obs %>% select(!!progr_names) %>%
    rename_all(function(x) gsub(name_prefix, "", x, fixed = TRUE)) %>%
    rownames_to_column("cell_id") %>%
    mutate(from = milestone_ids[from], to = milestone_ids[to])

  # read from var slots
  gimp_name <- paste0(name_prefix, "importance")
  if (gimp_name %in% colnames(anndata$var)) {
    gimp <- anndata$var[gimp_name] %>% rownames_to_column("feature_id")
    colnames(gimp) <- gsub(name_prefix, "", colnames(gimp), fixed = TRUE)
  } else {
    gimp <- NULL
  }

  # read from obsm slots
  dimred <- anndata$obsm[paste0(name_prefix, "dimred")]
  rownames(dimred) <- obs_names
  colnames(dimred) <- paste0("comp_", seq_len(ncol(dimred)))

  # recreate trajectory
  traj <-
    wrap_data(
      cell_ids = obs_names
    ) %>%
    add_expression(
      expression = expression,
      counts = counts
    ) %>%
    add_trajectory(
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      progressions = progressions
    ) %>%
    add_dimred(
      dimred = dimred,
      dimred_segment_points = dimred_segment_points,
      dimred_segment_progressions = dimred_segment_progressions,
      dimred_milestones = dimred_milestones
    )

  traj$feature_importances <- gimp

  traj
}
