




#' Copy trajectory to h5ad
#'
#' @param anndata AnnData object
#' @param traj Trajectory object
#' @param anndata AnnData object
#' @param trajectory_prefix The prefix of a trajectory to write.
#' @param dimred_name A dimred to write, a matrix in `ad$obsm`.
#' @param grouping_name A grouping to write, a column in `ad$obs`.
#' @param fimp_name A feature importance to write, a column in `ad$var`.
#'
#' @importFrom dplyr bind_cols select mutate
#' @importFrom purrr %>%
#' @export
to_h5ad <- function(
  traj,
  anndata = NULL,
  trajectory_prefix = NULL,
  dimred_name = NULL,
  grouping_name = NULL,
  fimp_name = NULL
) {
  if (is.null(anndata)) {
    anndata <- anndata::AnnData(
      X = traj$expression,
      obs_names = traj$cell_ids,
      var_names = traj$feature_ids
    )
  }

  obs_names <- anndata$obs_names
  var_names <- anndata$var_names

  # write grouping
  if (dynwrap::is_wrapper_with_grouping(traj) && !is.null(grouping_name)) {
    group <- traj$grouping
    group_obs <- group %>% select(-cell_id) %>% as.data.frame()
    rownames(group_obs) <- fimp$cell_id
    anndata$obs[[grouping_name]] <- group_obs[,1]
  }

  # write dimred
  if (dynwrap::is_wrapper_with_dimred(traj) && !is.null(dimred_name)) {
    anndata$obsm[[dimred_name]] <- traj$dimred[obs_names,,drop=FALSE]
  }

  # write feature importance
  if (dynwrap::is_wrapper_with_feature_importance(traj) && !is.null(fimp_name)) {
    fimp <- traj$feature_importance
    fimp_var <- fimp %>% select(-feature_id) %>% as.data.frame()
    rownames(fimp_var) <- fimp$feature_id
    anndata$var[[fimp_name]] <- fimp_var[,1]
  }


  # process traj
  if (dynwrap::is_wrapper_with_trajectory(traj) && !is.null(trajectory_prefix)) {
    # process clustering
    dimred_segments <- bind_cols(
      traj$dimred_segment_progressions,
      as.data.frame(traj$dimred_segment_points)
    )
    dimred_ids <- colnames(traj$dimred)

    anndata$uns[[paste0(trajectory_prefix, "dimred_ids")]] <- dimred_ids
    anndata$uns[[paste0(trajectory_prefix, "milestone_ids")]] <- traj$milestone_ids
    anndata$uns[[paste0(trajectory_prefix, "milestone_network")]] <- traj$milestone_network
    anndata$uns[[paste0(trajectory_prefix, "dimred_segments")]] <- dimred_segments
    anndata$uns[[paste0(trajectory_prefix, "dimred_milestones")]] <- traj$dimred_milestones

    progr_obs <- traj$progressions %>% as.data.frame() %>% select(-cell_id) %>% mutate(
      from = match(from, traj$milestone_ids),
      to = match(to, traj$milestone_ids)
    )
    rownames(progr_obs) <- traj$progressions$cell_id
    progr_obs <- progr_obs[obs_names,,drop=FALSE]
    for (nam in names(progr_obs)) {
      anndata$obs[paste0(trajectory_prefix, nam)] <- progr_obs[[nam]]
    }

    anndata$obsm[[paste0(trajectory_prefix, "dimred")]] <- traj$dimred
  }



  anndata
}
