#' Copy trajectory to h5ad
#'
#' @param anndata AnnData object
#' @param traj Trajectory object
#' @param gimp Gene Importance object
#' @param name_prefix A prefix to add to each created object in the anndata.
#'
#' @importFrom dplyr bind_cols select mutate
#' @importFrom purrr %>%
#' @export
to_h5ad <- function(
  traj,
  anndata = NULL,
  trajectory_prefix = "traj",
  dimred_prefix = trajectory_prefix,
  grouping_prefix = trajectory_prefix,
  fimp_prefix = trajectory_prefix
) {
  if (!is.null(anndata)) {
    obs_names <- anndata$obs_names
    var_names <- anndata$var_names
    X <- NULL # if anndata is not null, X is assumed to be filled in already
  } else {
    obs_names <- traj$cell_ids
    var_names <- traj$feature_ids
    X <- traj$expression
  }

  # process clustering
  dimred_segments <- bind_cols(
    traj$dimred_segment_progressions,
    as.data.frame(traj$dimred_segment_points)
  )
  dimred_ids <- colnames(traj$dimred)

  anndata$uns[[paste0(dimred_prefix, "dimred_ids")]] <- dimred_ids
  anndata$uns[[paste0(trajectory_prefix, "milestone_ids")]] <- traj$milestone_ids
  anndata$uns[[paste0(trajectory_prefix, "milestone_network")]] <- traj$milestone_network
  anndata$uns[[paste0(trajectory_prefix, "dimred_segments")]] <- dimred_segments
  anndata$uns[[paste0(trajectory_prefix, "dimred_milestones")]] <- traj$dimred_milestones

  # write to obs slots
  progr_obs <- traj$progressions %>% as.data.frame() %>% select(-cell_id) %>% mutate(
    from = match(from, traj$milestone_ids),
    to = match(to, traj$milestone_ids)
  )
  rownames(progr_obs) <- traj$progressions$cell_id
  progr_obs <- progr_obs[obs_names,,drop=FALSE]
  for (nam in names(progr_obs)) {
    anndata$obs[paste0(trajectory_prefix, nam)] <- progr_obs[[nam]]
  }

  # write to var slots
  fimp <- traj$feature_importance
  if (!is.null(fimp)) {
    fimp_var <- fimp %>% select(-feature_id) %>% as.data.frame()
    rownames(fimp_var) <- fimp$feature_id
    fimp_var <- fimp_var[var_names,,drop=FALSE]
    for (nam in names(fimp_var)) {
      anndata$var[paste0(fimp_prefix, nam)] <- fimp_var[[nam]]
    }
  }

  # write to obsm slots
  anndata$obsm[[paste0(dimred_prefix, "dimred")]] <- traj$dimred
}
