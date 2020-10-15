#' Copy trajectory to h5ad
#'
#' @param anndata AnnData object
#' @param traj Trajectory object
#' @param gimp Gene Importance object
#' @param name_prefix A prefix to add to each created object in the anndata.
#'
#' @importFrom dplyr bind_cols select mutate
#' @importFrom purrr %>%
#' @importFrom reticulate py_set_item
#' @export
to_h5ad <- function(anndata, traj, gimp = NULL, name_prefix = "") {
  obs_names <- python_builtins$list(anndata$obs_names)
  var_names <- python_builtins$list(anndata$var_names)

  # write to uns slots
  dimred_segments <- bind_cols(
    traj$dimred_segment_progressions,
    as.data.frame(traj$dimred_segment_points)
  )
  dimred_ids <- colnames(traj$dimred)
  reticulate::py_set_item(anndata$uns, paste0(name_prefix, "dimred_ids"), dimred_ids)
  reticulate::py_set_item(anndata$uns, paste0(name_prefix, "milestone_ids"), traj$milestone_ids)
  reticulate::py_set_item(anndata$uns, paste0(name_prefix, "milestone_network"), traj$milestone_network)
  reticulate::py_set_item(anndata$uns, paste0(name_prefix, "dimred_segments"), dimred_segments)
  reticulate::py_set_item(anndata$uns, paste0(name_prefix, "dimred_milestones"), traj$dimred_milestones)

  # write to obs slots
  obs <- traj$progressions %>% as.data.frame() %>% select(-cell_id) %>% mutate(
    from = match(from, traj$milestone_ids),
    to = match(to, traj$milestone_ids)
  )
  rownames(obs) <- traj$progressions$cell_id
  obs <- obs[obs_names,,drop=FALSE]
  for (nam in names(obs)) {
    anndata$obs[paste0(name_prefix, nam)] <- obs[[nam]]
  }

  # write to var slots
  if (!is.null(gimp)) {
    var <- gimp %>% select(-feature_id) %>% as.data.frame()
    rownames(var) <- gimp$feature_id
    var <- var[var_names,,drop=FALSE]
    for (nam in names(var)) {
      anndata$var[paste0(name_prefix, nam)] <- var[[nam]]
    }
  }

  # write to obsm slots
  reticulate::py_set_item(anndata$obsm, paste0(name_prefix, "dimred"), traj$dimred)
}
