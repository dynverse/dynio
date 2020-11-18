#' Copy trajectory from h5ad
#'
#' @param anndata AnnData object
#' @param layer_name A layer name to read. If `NULL`, will use `ad$X`.
#' @param trajectory_prefix The prefix of a trajectory to read.
#' @param dimred_name A dimred to read, a matrix in `ad$obsm`.
#' @param grouping_name A grouping to read, a column in `ad$obs`.
#' @param fimp_name A feature importance to read, a column in `ad$var`.
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
  layer_name = NULL,
  trajectory_prefix = "traj_",
  dimred_name = "dimred",
  grouping_name = "cluster",
  fimp_name = "importance"
) {
  obs_names <- anndata$obs_names
  var_names <- anndata$var_names

  # create trajectory
  traj <-
    wrap_data(
      cell_ids = obs_names,
      feature_ids = var_names
    )

  # read expression
  expression <-
    if (is.null(layer_name)) {
      as(anndata[], "CsparseMatrix")
    } else {
      anndata$layers[[layer_name]]
    }
  if (!is.null(expression)) {
    counts <- expression
    counts@x <- 2^counts@x - 1

    traj <- traj %>% add_expression(
      counts = counts,
      expression = expression
    )
  }

  # read trajectory
  if (!is.null(trajectory_prefix) && paste0(trajectory_prefix, "milestone_network") %in% anndata$uns_keys()) {

    # read milestone ids
    # as.vector removes extra attributes
    milestone_ids <- anndata$uns[[paste0(trajectory_prefix, "milestone_ids")]] %>% as.vector

    # read milestone network
    milestone_network <- anndata$uns[[paste0(trajectory_prefix, "milestone_network")]]

    # read progressions
    progr_names <- paste0(trajectory_prefix, c("from", "to", "percentage"))
    progressions <- anndata$obs %>% select(!!progr_names) %>%
      rename_all(function(x) gsub(trajectory_prefix, "", x, fixed = TRUE)) %>%
      rownames_to_column("cell_id") %>%
      mutate(from = milestone_ids[from], to = milestone_ids[to])

    # read dimred ids
    dimred_ids <- anndata$uns[[paste0(trajectory_prefix, "dimred_ids")]]

    # read dimred
    dimred <- anndata$obsm[[paste0(trajectory_prefix, "dimred")]]
    if (!is.null(dimred)) {
      rownames(dimred) <- obs_names
      colnames(dimred) <- dimred_ids
    }

    # read dimred segments
    dimred_segments <- anndata$uns[[paste0(trajectory_prefix, "dimred_segments")]]
    if (!is.null(dimred_segments)) {
      dimred_segment_progressions <- dimred_segments %>% select(-!!dimred_ids)
      dimred_segment_points <- dimred_segments %>% select(!!dimred_ids) %>% as.matrix
    } else {
      dimred_segment_progressions <- NULL
      dimred_segment_points <- NULL
    }

    # read dimred milestones
    dimred_milestones <- anndata$uns[[paste0(trajectory_prefix, "dimred_milestones")]]
    if (!is.null(dimred_milestones)) {
      rownames(dimred_milestones) <- milestone_ids
      colnames(dimred_milestones) <- dimred_ids
    }

    # create traj
    traj <- traj %>%
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
  }

  # read grouping
  if (!is.null(grouping_name) && grouping_name %in% anndata$obs_keys()) {
    clus <- anndata$obs[[grouping_name]]
    names(clus) <- obs_names

    traj <- traj %>% add_grouping(
      grouping = clus
    )
  }

  # read feature importances
  if (!is.null(fimp_name) && fimp_name %in% anndata$var_keys()) {
    gimp <- data.frame(
      feature_id = var_names,
      importance = anndata$var[[fimp_name]],
      stringsAsFactors = FALSE
    )

    traj <- traj %>% add_feature_importance(
      feature_importance = gimp
    )
  }

  # read dimred
  if (!is.null(dimred_name) && dimred_name %in% anndata$obsm_keys()) {
    dimred <- anndata$obsm[[dimred_name]]
    rownames(dimred) <- obs_names
    colnames(dimred) <- paste0("comp_", seq_len(ncol(dimred)))

    traj <- traj %>% add_dimred(
      dimred = dimred
    )
  }

  traj
}
