#' Copy trajectory from h5ad
#'
#' @param anndata AnnData object
#' @param name_prefix A prefix to add to each created object in the anndata.
#'
#' @importFrom dplyr bind_cols select mutate rename_all
#' @importFrom purrr %>%
#' @importFrom tibble rownames_to_column
#' @import dynwrap
#' @importFrom Matrix Matrix
#'
#' @export
from_h5ad <- function(
  anndata,
  trajectory_prefix = "traj",
  dimred_prefix = trajectory_prefix,
  grouping_prefix = trajectory_prefix
) {
  obs_names <- anndata$obs_names
  var_names <- anndata$var_names

  # read from X
  expression <- as(anndata[], "CsparseMatrix")
  counts <- expression
  counts@x <- 2^counts@x - 1

  # read from uns slots
  # as.vector removes extra attributes
  milestone_ids <- anndata$uns[[paste0(trajectory_prefix, "milestone_ids")]] %>% as.vector
  milestone_network <- anndata$uns[[paste0(trajectory_prefix, "milestone_network")]]
  dimred_ids <- anndata$uns[[paste0(dimred_prefix, "dimred_ids")]]

  dimred_segments <- anndata$uns[[paste0(trajectory_prefix, "dimred_segments")]]
  dimred_segment_progressions <- dimred_segments %>% select(-!!dimred_ids)
  dimred_segment_points <- dimred_segments %>% select(!!dimred_ids) %>% as.matrix

  dimred_milestones <- anndata$uns[[paste0(trajectory_prefix, "dimred_milestones")]]
  rownames(dimred_milestones) <- milestone_ids
  colnames(dimred_milestones) <- dimred_ids

  # read from obs slots
  progr_names <- paste0(trajectory_prefix, c("from", "to", "percentage"))
  progressions <- anndata$obs %>% select(!!progr_names) %>%
    rename_all(function(x) gsub(trajectory_prefix, "", x, fixed = TRUE)) %>%
    rownames_to_column("cell_id") %>%
    mutate(from = milestone_ids[from], to = milestone_ids[to])

  # read from var slots
  gimp_name <- paste0(fimp_prefix, "importance")
  if (gimp_name %in% colnames(anndata$var)) {
    gimp <- anndata$var[gimp_name] %>% rownames_to_column("feature_id")
    colnames(gimp) <- gsub(fimp_prefix, "", colnames(gimp), fixed = TRUE)
  } else {
    gimp <- NULL
  }

  # read from obsm slots
  dimred <- anndata$obsm[paste0(dimred_prefix, "dimred")]
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
    ) %>%
    add_feature_importance(
      feature_importance = gimp
    )

  traj
}
