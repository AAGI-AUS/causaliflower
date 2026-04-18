################################################################################
# Constants
################################################################################
.MODIS_UNSUPPORTED_INDICES <- c("NDWI", "NDMI", "MNDWI", "NBR", "NBR2", "BSI")
.DEFAULT_S2      <- 20
.DEFAULT_LANDSAT <- 30
.S2_START <- as.Date("2015-07-01")
.L8_START <- as.Date("2013-04-01")
.STAC_ENDPOINT <- "https://planetarycomputer.microsoft.com/api/stac/v1"

################################################################################
# Helpers
################################################################################

#' Null-coalesce operator
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Convert metres to decimal degrees
#' @export
metres_to_degrees <- function(latitude, distance = 1, pixel_size = 1){
  a <- 6378137.0
  e2 <- 0.00669437999014
  lat_c_rad <- latitude * pi / 180
  m_c <- a * (1 - e2) / (1 - e2 * sin(lat_c_rad)^2)^1.5
  half_lat <- (distance / 2) * 180 / (pi * m_c)
  lat_rad <- (latitude + c(-half_lat, 0, half_lat)) * pi / 180
  denom <- 1 - e2 * sin(lat_rad)^2
  n <- a / sqrt(denom)
  m <- a * (1 - e2) / denom^1.5
  w <- c(1, 4, 1)
  c(x_degrees = pixel_size * 180 / (pi * sum(n * cos(lat_rad) * w) / 6),
    y_degrees = pixel_size * 180 / (pi * sum(m * w) / 6))
}


#' Bounding box around a central point
#' @export
point_to_bbox <- function(longitude,
                          latitude,
                          x_metres,
                          y_metres
){
  x_half_width <- (x_metres / 2)
  y_half_width <- (y_metres / 2)

  point <- sf::st_sfc( sf::st_point(c(longitude, latitude)), crs = sf::st_crs(4326) )
  bbox_x <- point |> sf::st_buffer(dist = x_half_width) |> sf::st_bbox()
  bbox_y <- point |> sf::st_buffer(dist = y_half_width) |> sf::st_bbox()
  bbox <- c(bbox_x[1], bbox_y[2], bbox_x[3], bbox_y[4])
  sf::st_bbox(bbox, crs = sf::st_crs(4326))
}


#' Default crop_search_params values
#'
#' Returns the default parameter list consumed by \code{\link{crop_search}} and
#' the \code{crop_search_params} argument of \code{\link{get_rs_data}} /
#' \code{\link{get_rs_raster}}. Users can call this to inspect defaults, then
#' override selected fields before passing the list through.
#'
#' @return Named list of default crop-search parameters.
#' @examples
#' params <- crop_search_defaults()
#' params$VI_sensitivity <- 0.15
#' ## get_rs_data(..., crop_search_params = params)
#' @export
crop_search_defaults <- function(){
  list(VI_threshold             = 0,
       VI_sensitivity           = 0.10,
       expand_search_area       = FALSE,
       search_windows           = "all",
       enhanced_search          = NULL,
       initial_search_metres    = 300,
       max_search_metres        = 1000,
       search_start_date        = NULL,
       search_end_date          = NULL,
       search_max_cloud_cover   = NULL,
       search_scl_classes       = NULL,
       search_margin_behaviour  = "hard")
}


#' Merge user crop_search_params against defaults (internal)
#' @noRd
.mergeCropSearchParams <- function(user_params,
                                   start_date,
                                   end_date,
                                   max_cloud_cover,
                                   scl_classes,
                                   modis_interpolation
){
  defs <- crop_search_defaults()
  if (!is.list(user_params)) stop("crop_search_params must be a list.", call. = FALSE)
  unknown <- setdiff(names(user_params), names(defs))
  if (length(unknown)) stop("Unknown crop_search_params: ", paste(unknown, collapse = ", "), call. = FALSE)
  merged <- modifyList(defs, user_params)
  merged$search_start_date      <- merged$search_start_date      %||% start_date
  merged$search_end_date        <- merged$search_end_date        %||% end_date
  merged$search_max_cloud_cover <- merged$search_max_cloud_cover %||% max_cloud_cover
  merged$search_scl_classes     <- merged$search_scl_classes     %||% scl_classes
  if (is.null(merged$enhanced_search)){
    merged$enhanced_search <- !identical(modis_interpolation, "off")
  }
  if (!merged$search_margin_behaviour %in% c("hard", "soft")){
    stop("search_margin_behaviour must be \"hard\" or \"soft\".", call. = FALSE)
  }
  merged
}


#' Test a 3-pixel edge against a reference signature (internal)
#' @noRd
.edgeTest <- function(rows, cols, feat_arr, ref_sig, threshold, dist_cache,
                      vi_floor = 0
){
  alive <- logical(3)
  for (k in seq_len(3)){
    r <- rows[k]; c <- cols[k]
    if (r < 1 || r > nrow(dist_cache) || c < 1 || c > ncol(dist_cache)){
      alive[k] <- FALSE
      next
    }
    d <- dist_cache[r, c]
    if (is.na(d)){
      pixel_sig <- feat_arr[r, c, ]
      if (any(is.na(pixel_sig)) || all(pixel_sig < vi_floor)){
        d <- Inf
      } else {
        d <- min(abs(pixel_sig - ref_sig))
      }
      dist_cache[r, c] <- d
    }
    alive[k] <- d < threshold
  }
  list(alive = alive, dist_cache = dist_cache)
}


#' Expand search in one direction from a centre pixel (internal)
#'
#' The outer reserve (25\%) is respected via reserve_r_lo/hi, reserve_c_lo/hi —
#' boundary-pass logic activates only when the edge reaches the reserve
#' boundary (not the raster edge).
#' @noRd
.expandDirection <- function(centre_rc,
                             direction,
                             feat_arr,
                             ref_sig,
                             threshold,
                             dist_cache,
                             reserve_rows,
                             reserve_cols,
                             vi_floor = 0
){
  nr <- dim(feat_arr)[1]; nc <- dim(feat_arr)[2]
  horizontal <- direction %in% c("E", "W")
  step <- if (horizontal){
    if (direction == "E") 1L else -1L
  } else {
    if (direction == "N") -1L else 1L
  }

  if (horizontal){
    c_pos <- centre_rc[2] + step
    r_set <- centre_rc[1] + c(-1, 0, 1)
    c_set <- rep(c_pos, 3)
  } else {
    r_pos <- centre_rc[1] + step
    r_set <- rep(r_pos, 3)
    c_set <- centre_rc[2] + c(-1, 0, 1)
  }

  alive_rows        <- integer(0)
  alive_cols        <- integer(0)
  boundary_passes   <- 0L
  needs_widening    <- FALSE

  repeat {
    res <- .edgeTest(r_set, c_set, feat_arr, ref_sig, threshold, dist_cache,
                     vi_floor = vi_floor)
    dist_cache <- res$dist_cache
    alive      <- res$alive
    n_alive    <- sum(alive)

    if (n_alive > 0){
      alive_rows <- c(alive_rows, r_set[alive])
      alive_cols <- c(alive_cols, c_set[alive])
    }

    if (n_alive == 0) break

    at_reserve <- if (horizontal){
      c_set[1] <= reserve_cols[1] || c_set[1] >= reserve_cols[2]
    } else {
      r_set[1] <= reserve_rows[1] || r_set[1] >= reserve_rows[2]
    }

    if (horizontal){
      alive_rows_on_edge <- r_set[alive]
      next_c <- c_set[1] + step
      if (next_c < 1 || next_c > nc){
        if (at_reserve){
          boundary_passes <- boundary_passes + 1L
          if (boundary_passes >= 2L) needs_widening <- TRUE
        }
        break
      }
      pivot_r <- alive_rows_on_edge[which.min(abs(alive_rows_on_edge - r_set[2]))]
      r_set   <- pivot_r + c(-1, 0, 1)
      c_set   <- rep(next_c, 3)
    } else {
      alive_cols_on_edge <- c_set[alive]
      next_r <- r_set[1] + step
      if (next_r < 1 || next_r > nr){
        if (at_reserve){
          boundary_passes <- boundary_passes + 1L
          if (boundary_passes >= 2L) needs_widening <- TRUE
        }
        break
      }
      pivot_c <- alive_cols_on_edge[which.min(abs(alive_cols_on_edge - c_set[2]))]
      r_set   <- rep(next_r, 3)
      c_set   <- pivot_c + c(-1, 0, 1)
    }

    if (!at_reserve) boundary_passes <- 0L
  }

  list(alive_rc       = cbind(row = alive_rows, col = alive_cols),
       dist_cache     = dist_cache,
       needs_widening = needs_widening)
}


#' Non-parametric crop area search (internal)
#' @noRd
.nonparametricCropSearch <- function(fs,
                                     crop_search_m,
                                     VI_sensitivity = 0.10,
                                     VI_threshold = 0,
                                     search_window = "all",
                                     expand_search_area = FALSE,
                                     preserved_alive_xy = NULL
){
  log_msg <- function(...) message("    ", ...)

  crop_search_half_m <- crop_search_m / 2
  centre_xy          <- terra::geom(fs$centre_v)[1, c("x", "y")]
  crop_extent        <- terra::ext(centre_xy[["x"]] - crop_search_half_m,
                                   centre_xy[["x"]] + crop_search_half_m,
                                   centre_xy[["y"]] - crop_search_half_m,
                                   centre_xy[["y"]] + crop_search_half_m)
  search_rast        <- terra::crop(fs$feature_stack, crop_extent, snap = "out")
  n_layers_total     <- terra::nlyr(search_rast)
  nr <- terra::nrow(search_rast); nc <- terra::ncol(search_rast)

  window_idx <- if (identical(search_window, "all")){
    seq_len(n_layers_total)
  } else if (is.numeric(search_window) && length(search_window) >= 1){
    w <- as.integer(search_window)
    if (any(w < 1) || any(w > n_layers_total) || anyDuplicated(w)){
      return(list(success = FALSE, reason = "search_window out of range or duplicated"))
    }
    w
  } else {
    return(list(success = FALSE, reason = "search_window must be 'all' or integer(s)"))
  }
  active_rast <- search_rast[[window_idx]]

  centre_cell <- terra::cellFromXY(active_rast, matrix(centre_xy, nrow = 1))
  if (is.na(centre_cell)) return(list(success = FALSE, reason = "centre outside cropped grid"))
  centre_row  <- terra::rowFromCell(active_rast, centre_cell)
  centre_col  <- terra::colFromCell(active_rast, centre_cell)
  centre_rc   <- c(centre_row, centre_col)

  ## 75/25 reserve: boundary-pass trigger zone = inner 75%, outer 25% is reserve
  reserve_rows <- c(trunc(nr * 0.125) + 1L, nr - trunc(nr * 0.125))
  reserve_cols <- c(trunc(nc * 0.125) + 1L, nc - trunc(nc * 0.125))

  ## Preserve alive pixels from prior iteration (if any)
  preserved_rc <- NULL
  if (!is.null(preserved_alive_xy) && nrow(preserved_alive_xy) > 0){
    cells <- terra::cellFromXY(active_rast, preserved_alive_xy)
    valid <- !is.na(cells)
    if (any(valid)){
      preserved_rc <- cbind(row = terra::rowFromCell(active_rast, cells[valid]),
                            col = terra::colFromCell(active_rast, cells[valid]))
      log_msg(sprintf("Preserved %d alive pixels from prior iteration", nrow(preserved_rc)))
    }
  }

  union_rc       <- matrix(integer(0), ncol = 2, dimnames = list(NULL, c("row", "col")))
  any_widen      <- FALSE
  n_flagged      <- 0L
  n_parties      <- 4L * length(window_idx)

  for (wi in seq_along(window_idx)){
    win_label <- window_idx[wi]
    log_msg(sprintf("--- Window %d ---", win_label))

    feat_arr <- array(terra::as.matrix(active_rast[[wi]], wide = TRUE), c(nr, nc, 1L))

    win_res <- .searchOneWindow(feat_arr       = feat_arr,
                                centre_rc      = centre_rc,
                                VI_sensitivity = VI_sensitivity,
                                VI_threshold   = VI_threshold,
                                reserve_rows   = reserve_rows,
                                reserve_cols   = reserve_cols,
                                nr             = nr,
                                nc             = nc,
                                log_msg        = log_msg)

    n_flagged <- n_flagged + win_res$n_widen
    if (win_res$n_widen > 0) any_widen <- TRUE
    if (is.null(win_res$final_rc)){
      log_msg(sprintf("Window %d produced no patch, skipping", win_label))
      next
    }
    union_rc <- unique(rbind(union_rc, win_res$final_rc))
  }

  if (!is.null(preserved_rc)){
    union_rc <- unique(rbind(union_rc, preserved_rc))
  }

  ## Export UTM coordinates of every alive pixel for preservation across expansion
  alive_xy_out <- NULL
  if (nrow(union_rc) > 0){
    cells <- terra::cellFromRowCol(active_rast, union_rc[, "row"], union_rc[, "col"])
    alive_xy_out <- terra::xyFromCell(active_rast, cells)
  }

  if (nrow(union_rc) == 0){
    return(list(success = FALSE, reason = "no alive pixels found",
                needs_widening = any_widen, n_flagged = n_flagged, n_parties = n_parties,
                alive_xy = alive_xy_out))
  }

  patch_rast <- terra::rast(active_rast[[1]])
  terra::values(patch_rast) <- NA_real_
  cells <- terra::cellFromRowCol(patch_rast, union_rc[, "row"], union_rc[, "col"])
  patch_rast[cells] <- 1

  log_msg(sprintf("Union crop area patch: %d px", nrow(union_rc)))

  trimmed_patch     <- terra::trim(patch_rast)
  trimmed_patch_ext <- terra::ext(trimmed_patch)
  search_ext        <- terra::ext(active_rast)
  px                <- fs$pixel_size_metres

  x_metres <- as.numeric(trimmed_patch_ext$xmax - trimmed_patch_ext$xmin)
  y_metres <- as.numeric(trimmed_patch_ext$ymax - trimmed_patch_ext$ymin)

  truncated <- (trimmed_patch_ext$xmin <= search_ext$xmin + px) ||
    (trimmed_patch_ext$xmax >= search_ext$xmax - px) ||
    (trimmed_patch_ext$ymin <= search_ext$ymin + px) ||
    (trimmed_patch_ext$ymax >= search_ext$ymax - px)

  poly_utm <- terra::as.polygons(trimmed_patch, dissolve = TRUE)
  poly_wgs <- terra::project(poly_utm, "EPSG:4326")
  mask_sf  <- sf::st_as_sf(poly_wgs)

  cent_utm <- terra::centroids(poly_utm)
  cent_wgs <- terra::project(cent_utm, "EPSG:4326")
  xy       <- terra::geom(cent_wgs)[1, c("x", "y")]

  list(success        = TRUE,
       mask_sf        = mask_sf,
       x_metres       = x_metres,
       y_metres       = y_metres,
       centroid_lon   = as.numeric(xy["x"]),
       centroid_lat   = as.numeric(xy["y"]),
       n_pixels       = as.integer(nrow(union_rc)),
       truncated      = truncated,
       needs_widening = any_widen,
       n_flagged      = n_flagged,
       n_parties      = n_parties,
       alive_xy       = alive_xy_out)
}


#' Run Phase 1 + Phase 2 for a single window (internal)
#' @noRd
.searchOneWindow <- function(feat_arr,
                             centre_rc,
                             VI_sensitivity,
                             VI_threshold,
                             reserve_rows,
                             reserve_cols,
                             nr,
                             nc,
                             log_msg
){
  centre_row <- centre_rc[1]; centre_col <- centre_rc[2]
  threshold  <- VI_sensitivity

  centre_sig <- feat_arr[centre_row, centre_col, ]
  if (any(is.na(centre_sig))){
    log_msg("Centre pixel has no valid data in this window")
    return(list(final_rc = NULL, n_widen = 0L))
  }
  log_msg(sprintf("Initial centre signature: [%s]", paste(sprintf("%.3f", centre_sig), collapse = ", ")))

  dist_cache_init <- matrix(NA_real_, nr, nc)
  dist_cache_init[centre_row, centre_col] <- 0

  init_rc <- expand.grid(row = centre_row + -1:1, col = centre_col + -1:1,
                         KEEP.OUT.ATTRS = FALSE)
  init_rc <- as.matrix(init_rc[init_rc$row >= 1 & init_rc$row <= nr &
                                 init_rc$col >= 1 & init_rc$col <= nc, ])

  init_alive <- vapply(seq_len(nrow(init_rc)), function(i){
    r <- init_rc[i, "row"]; c <- init_rc[i, "col"]
    pixel_sig <- feat_arr[r, c, ]
    if (any(is.na(pixel_sig)) || all(pixel_sig < VI_threshold)) return(FALSE)
    d <- min(abs(pixel_sig - centre_sig))
    dist_cache_init[r, c] <<- d
    d < threshold
  }, logical(1))

  passing_rc <- init_rc[init_alive, , drop = FALSE]
  if (nrow(passing_rc) >= 2){
    passing_sigs <- t(apply(passing_rc, 1, function(rc) feat_arr[rc[1], rc[2], ]))
    centre_sig   <- apply(passing_sigs, 2, median, na.rm = TRUE)
    log_msg(sprintf("Refined centre signature (%d-pixel median): [%s]",
                    nrow(passing_rc),
                    paste(sprintf("%.3f", centre_sig), collapse = ", ")))
  } else {
    log_msg("False-start refinement skipped (< 2 passing pixels)")
  }

  dist_cache <- matrix(NA_real_, nr, nc)
  dist_cache[centre_row, centre_col] <- 0

  results <- list()
  n_widen <- 0L
  for (dir in c("N", "S", "E", "W")){
    res <- .expandDirection(centre_rc, dir, feat_arr, centre_sig,
                            threshold, dist_cache,
                            reserve_rows = reserve_rows, reserve_cols = reserve_cols,
                            vi_floor = VI_threshold)
    dist_cache <- res$dist_cache
    results[[dir]] <- res
    if (res$needs_widening) n_widen <- n_widen + 1L
  }

  n_success <- sum(vapply(results, function(r) nrow(r$alive_rc) > 0, logical(1)))
  if (n_success == 0){
    log_msg("No successful search direction in this window")
    return(list(final_rc = NULL, n_widen = n_widen))
  }

  ref_pixels_rc <- matrix(c(centre_row, centre_col), ncol = 2, dimnames = list(NULL, c("row", "col")))
  for (dir in names(results)){
    res <- results[[dir]]
    if (nrow(res$alive_rc) == 0) next
    d_centre   <- abs(res$alive_rc[, "row"] - centre_row) + abs(res$alive_rc[, "col"] - centre_col)
    first_edge <- res$alive_rc[d_centre == min(d_centre), , drop = FALSE]
    ref_pixels_rc <- rbind(ref_pixels_rc, first_edge)
  }
  ref_pixels_rc <- unique(ref_pixels_rc)

  ref_sigs    <- t(apply(ref_pixels_rc, 1, function(rc) feat_arr[rc[1], rc[2], ]))
  refined_sig <- apply(ref_sigs, 2, median, na.rm = TRUE)
  log_msg(sprintf("Refined signature: [%s] from %d pixels",
                  paste(sprintf("%.3f", refined_sig), collapse = ", "),
                  nrow(ref_pixels_rc)))

  all_alive <- do.call(rbind, lapply(results, `[[`, "alive_rc"))
  all_alive <- unique(rbind(all_alive, c(centre_row, centre_col)))
  colnames(all_alive) <- c("row", "col")

  local_outers <- list()
  for (dir in names(results)){
    res <- results[[dir]]
    if (nrow(res$alive_rc) == 0) next
    local_outers[[dir]] <- switch(dir,
                                  "N" = res$alive_rc[res$alive_rc[, "row"] == min(res$alive_rc[, "row"]), , drop = FALSE],
                                  "S" = res$alive_rc[res$alive_rc[, "row"] == max(res$alive_rc[, "row"]), , drop = FALSE],
                                  "E" = res$alive_rc[res$alive_rc[, "col"] == max(res$alive_rc[, "col"]), , drop = FALSE],
                                  "W" = res$alive_rc[res$alive_rc[, "col"] == min(res$alive_rc[, "col"]), , drop = FALSE])
  }
  local_rc <- do.call(rbind, local_outers)

  global_rc <- rbind(
    all_alive[all_alive[, "row"] == min(all_alive[, "row"]), , drop = FALSE][1, , drop = FALSE],
    all_alive[all_alive[, "row"] == max(all_alive[, "row"]), , drop = FALSE][1, , drop = FALSE],
    all_alive[all_alive[, "col"] == min(all_alive[, "col"]), , drop = FALSE][1, , drop = FALSE],
    all_alive[all_alive[, "col"] == max(all_alive[, "col"]), , drop = FALSE][1, , drop = FALSE])
  global_rc <- unique(global_rc)

  is_g_connected <- apply(local_rc, 1, function(p){
    any(p["row"] == global_rc[, "row"] | p["col"] == global_rc[, "col"])
  })
  g_separated <- local_rc[!is_g_connected, , drop = FALSE]
  g_connected <- local_rc[is_g_connected, , drop = FALSE]

  trace_line <- function(from_rc, to_rc){
    n_steps <- max(abs(to_rc[1] - from_rc[1]), abs(to_rc[2] - from_rc[2])) + 1
    rs <- round(seq(from_rc[1], to_rc[1], length.out = n_steps))
    cs <- round(seq(from_rc[2], to_rc[2], length.out = n_steps))
    cbind(row = rs, col = cs)
  }

  perimeter_rc <- matrix(integer(0), ncol = 2, dimnames = list(NULL, c("row", "col")))

  for (gi in seq_len(nrow(global_rc))){
    g <- global_rc[gi, ]
    for (li in seq_len(nrow(g_connected))){
      l <- g_connected[li, ]
      if (all(g == l)) next
      if (g["row"] == l["row"] || g["col"] == l["col"]){
        perimeter_rc <- rbind(perimeter_rc, trace_line(g, l))
      }
    }
  }
  if (nrow(g_separated) > 0 && nrow(g_connected) > 0){
    for (si in seq_len(nrow(g_separated))){
      s <- g_separated[si, ]
      dists  <- pmax(abs(g_connected[, "row"] - s["row"]),
                     abs(g_connected[, "col"] - s["col"]))
      target <- g_connected[which.max(dists), ]
      perimeter_rc <- rbind(perimeter_rc, trace_line(s, target))
    }
  }

  kept_perimeter <- matrix(integer(0), ncol = 2, dimnames = list(NULL, c("row", "col")))
  if (nrow(perimeter_rc) > 0){
    perimeter_rc <- unique(perimeter_rc)
    for (pi in seq_len(nrow(perimeter_rc))){
      r <- perimeter_rc[pi, "row"]; c <- perimeter_rc[pi, "col"]
      if (r < 1 || r > nr || c < 1 || c > nc) next
      d <- dist_cache[r, c]
      if (is.na(d)){
        pixel_sig <- feat_arr[r, c, ]
        d <- if (any(is.na(pixel_sig))) Inf else min(abs(pixel_sig - refined_sig))
        dist_cache[r, c] <- d
      }
      pixel_sig <- feat_arr[r, c, ]
      if (d < threshold && !any(is.na(pixel_sig)) && !all(pixel_sig < VI_threshold)){
        kept_perimeter <- rbind(kept_perimeter, c(r, c))
      }
    }
  }

  r_min <- min(global_rc[, "row"]); r_max <- max(global_rc[, "row"])
  c_min <- min(global_rc[, "col"]); c_max <- max(global_rc[, "col"])

  final_rc <- unique(rbind(all_alive, kept_perimeter))

  interior_grid <- expand.grid(row = r_min:r_max, col = c_min:c_max,
                               KEEP.OUT.ATTRS = FALSE)
  interior_rc <- as.matrix(interior_grid)
  interior_rc <- interior_rc[interior_rc[, "row"] >= 1 & interior_rc[, "row"] <= nr &
                               interior_rc[, "col"] >= 1 & interior_rc[, "col"] <= nc, , drop = FALSE]

  final_rc <- unique(rbind(final_rc, interior_rc))

  in_bbox  <- final_rc[, "row"] >= r_min & final_rc[, "row"] <= r_max &
    final_rc[, "col"] >= c_min & final_rc[, "col"] <= c_max
  is_perim <- apply(final_rc, 1, function(p){
    any(p["row"] == kept_perimeter[, "row"] & p["col"] == kept_perimeter[, "col"])
  })
  keep     <- in_bbox | is_perim
  final_rc <- final_rc[keep, , drop = FALSE]

  log_msg(sprintf("Window patch: %d px", nrow(final_rc)))
  list(final_rc = final_rc, n_widen = n_widen)
}


#' Detect field dimensions at a coordinate via non-parametric crop search
#'
#' Returns the bounding-box dimensions of the managed field at a given
#' coordinate. Iteratively widens the search area via area-based expansion
#' when Phase 1 detects boundary collisions, until clear buffer is found or
#' \code{max_search_metres} is reached.
#'
#' Intended for one-off use. Not for batch workflows — each call makes STAC
#' fetches per iteration.
#'
#' @details
#' The algorithm compares the signature of the centre pixel (refined via a
#' 3x3 false-start median) against neighbouring pixels in four cardinal
#' directions. Each selected window runs Phase 1 + Phase 2 independently and
#' results are unioned. Perimeter lines are traced between outer pixels to
#' close the bounding polygon.
#'
#' When \code{expand_search_area = TRUE}, each odd-numbered iteration may
#' trigger a new widened STAC fetch if boundary collisions were flagged. The
#' new fetch size is computed from the area-based multiplier
#' \code{sqrt(6 * n_flagged / n_parties)}, where flagged directions count
#' boundary second-passes across all selected windows. Alive pixels are
#' preserved across iterations via UTM coordinate translation.
#'
#' @section Crop search parameters:
#' \describe{
#'   \item{VI_threshold}{Minimum VI value. Pixels with all selected-window values below this are dead. Default 0.}
#'   \item{VI_sensitivity}{Per-window distance threshold (0-1 EVI2 units). Default 0.10.}
#'   \item{expand_search_area}{If TRUE, boundary-collision iterations trigger widened refetch. Default FALSE.}
#'   \item{search_windows}{"all" or integer vector of window indices. Default "all".}
#'   \item{enhanced_search}{Use MODIS-augmented feature stack. NULL auto-inherits from modis_interpolation != "off". Default NULL.}
#'   \item{initial_search_metres}{Starting search area. Default 300.}
#'   \item{max_search_metres}{Maximum search width. Default 1000.}
#'   \item{search_start_date, search_end_date}{Date range for detection feature stack. Default inherits from main start_date/end_date.}
#'   \item{search_max_cloud_cover}{Cloud cover for detection. Default inherits from max_cloud_cover.}
#'   \item{search_scl_classes}{SCL classes to mask during detection. Default inherits from scl_classes.}
#'   \item{search_margin_behaviour}{"hard" caps at max_search_metres. "soft" warns and continues at computed fetch. Default "hard".}
#' }
#'
#' @param longitude,latitude Centre coordinate in WGS84.
#' @param start_date,end_date Date range spanning the growing period.
#' @param crop_search_params Named list of parameters. See
#'   \code{\link{crop_search_defaults}} and the \strong{Crop search parameters}
#'   section for available fields. Empty list uses all defaults.
#' @param modis_interpolation Controls MODIS usage. One of \code{"off"}
#'   (default), \code{"merge"}, or \code{"sep"}. Only affects the default
#'   behaviour of \code{crop_search_params$enhanced_search} when user leaves
#'   it as NULL.
#' @param max_cloud_cover Default cloud cover for detection (overridden by
#'   \code{search_max_cloud_cover}). Default 50.
#' @param field_mask \code{FALSE} (default), \code{TRUE}, or an \code{sf} /
#'   \code{sfc} polygon for pixel-level masking.
#' @param scl_classes Default S2 SCL classes (overridden by
#'   \code{search_scl_classes}).
#' @return Named list with \code{x_metres}, \code{y_metres}, \code{longitude},
#'   \code{latitude}, \code{dim_source}, \code{detected}, \code{truncated},
#'   \code{reason}, and \code{mask_sf}.
#' @export
crop_search <- function(longitude,
                        latitude,
                        start_date,
                        end_date,
                        crop_search_params = list(),
                        modis_interpolation = "off",
                        max_cloud_cover = 50,
                        field_mask = FALSE,
                        scl_classes = NULL
){
  log_msg <- function(...) message("  ", ...)

  params <- .mergeCropSearchParams(crop_search_params,
                                   start_date, end_date,
                                   max_cloud_cover, scl_classes,
                                   modis_interpolation)

  fallback <- function(reason, metres){
    list(x_metres     = metres,
         y_metres     = metres,
         longitude    = longitude,
         latitude     = latitude,
         dim_source   = "default",
         detected     = FALSE,
         truncated    = NA,
         reason       = reason,
         mask_sf      = NULL)
  }

  iter           <- 0L
  crop_search_m  <- params$initial_search_metres
  fetch_metres   <- crop_search_m
  final_det      <- NULL
  final_fs       <- NULL
  preserved_xy   <- NULL
  last_reason    <- NA_character_
  coords         <- list(longitude = longitude, latitude = latitude)
  stack_fn       <- if (isTRUE(params$enhanced_search)) .enhancedFeatureStack else .fetchFeatureStack

  repeat {
    iter <- iter + 1L
    log_msg(sprintf("Iteration %d: fetching %.0f x %.0f m", iter, fetch_metres, fetch_metres))

    fs <- stack_fn(longitude       = longitude,
                   latitude        = latitude,
                   start_date      = params$search_start_date,
                   end_date        = params$search_end_date,
                   fetch_metres    = fetch_metres,
                   max_cloud_cover = params$search_max_cloud_cover,
                   scl_classes     = params$search_scl_classes)
    if (!fs$success){
      last_reason <- fs$reason
      log_msg("  Fetch failed: ", fs$reason)
      if (fetch_metres >= params$max_search_metres) return(fallback(last_reason, fetch_metres))
      fetch_metres <- min(fetch_metres * 2, params$max_search_metres)
      next
    }
    final_fs <- fs

    det <- .nonparametricCropSearch(fs                 = fs,
                                    crop_search_m      = fetch_metres,
                                    VI_sensitivity     = params$VI_sensitivity,
                                    VI_threshold       = params$VI_threshold,
                                    search_window      = params$search_windows,
                                    expand_search_area = params$expand_search_area,
                                    preserved_alive_xy = preserved_xy)

    if (!is.null(det$alive_xy)) preserved_xy <- det$alive_xy

    if (det$success && !det$truncated && !isTRUE(det$needs_widening)){
      final_det <- det
      coords    <- list(longitude = det$centroid_lon, latitude = det$centroid_lat)
      break
    }

    last_reason <- det$reason %||% NA_character_

    ## Odd-iteration expansion trigger
    trigger_expand <- isTRUE(params$expand_search_area) &&
      (iter %% 2L == 1L) &&
      isTRUE(det$needs_widening) &&
      det$n_flagged > 0

    at_cap <- fetch_metres >= params$max_search_metres

    if (at_cap){
      if (det$success){
        warning("Detected crop area touches bbox edge at max_search_metres = ",
                params$max_search_metres,
                ". Consider raising max_search_metres.", call. = FALSE)
        final_det <- det
        coords    <- list(longitude = det$centroid_lon, latitude = det$centroid_lat)
        break
      }
      return(fallback(last_reason, fetch_metres))
    }

    next_fetch <- if (trigger_expand){
      mult <- sqrt(6 * det$n_flagged / det$n_parties)
      trunc(fetch_metres * mult)
    } else {
      fetch_metres * 2L
    }

    if (next_fetch > params$max_search_metres){
      if (identical(params$search_margin_behaviour, "hard")){
        next_fetch <- params$max_search_metres
      } else {
        warning(sprintf("search_margin_behaviour = 'soft': proceeding at %d m (exceeds max_search_metres = %d).",
                        next_fetch, params$max_search_metres), call. = FALSE)
      }
    }
    fetch_metres <- next_fetch
  }

  px       <- final_fs$pixel_size_metres
  x_metres <- final_det$x_metres + 2 * px
  y_metres <- final_det$y_metres + 2 * px
  mask_sf  <- final_det$mask_sf
  if (inherits(field_mask, c("sf", "sfc"))){
    mask_sf <- if (inherits(field_mask, "sf")) field_mask else sf::st_sf(geometry = field_mask)
  }

  list(x_metres     = x_metres,
       y_metres     = y_metres,
       longitude    = coords$longitude,
       latitude     = coords$latitude,
       dim_source   = "detected",
       detected     = TRUE,
       truncated    = final_det$truncated,
       reason       = NA_character_,
       mask_sf      = mask_sf)
}


#' Fetch a multi-window EVI2 feature stack for phenology-based detection (internal)
#' @noRd
.fetchFeatureStack <- function(longitude,
                               latitude,
                               start_date,
                               end_date,
                               fetch_metres,
                               max_cloud_cover = 50,
                               num_windows = 3L,
                               scl_classes = NULL
){
  fetch_metres <- max(fetch_metres, (.DEFAULT_LANDSAT * 20))
  start_d <- suppressWarnings(as.Date(start_date))
  end_d   <- suppressWarnings(as.Date(end_date))
  if (is.na(start_d) || is.na(end_d)) return(list(success = FALSE, reason = "invalid dates"))
  if (as.numeric(end_d - start_d) < 2) return(list(success = FALSE, reason = "date range too short"))

  index_list   <- suppressMessages(.checkConfigs(start_date, "EVI2"))
  source       <- index_list[[1]]
  index_list   <- index_list[[2]]
  pixel_size_m <- index_list$res
  fetch_bbox   <- point_to_bbox(longitude, latitude, fetch_metres, fetch_metres)

  num_windows <- min(3L, trunc(as.numeric(end_d - start_d) / (pixel_size_m * 2) + 1))
  breaks  <- seq(start_d, end_d, length.out = num_windows + 1L)

  composites <- lapply(seq_len(num_windows), function(w){
    scenes <- tryCatch(
      .fetchStac(bbox            = fetch_bbox,
                 start_date      = as.character(breaks[w]),
                 end_date        = as.character(breaks[w + 1L]),
                 asset_names     = index_list$assets,
                 index_function  = index_list$fun,
                 source          = source,
                 max_cloud_cover = max_cloud_cover,
                 scl_classes     = scl_classes),
      error = function(e) NULL)
    if (is.null(scenes) || length(scenes) == 0) return(NULL)
    stack <- Reduce(c, lapply(scenes, `[[`, "rast"))
    if (terra::nlyr(stack) > 1) max(stack, na.rm = TRUE) else stack
  })

  if (any(vapply(composites, is.null, logical(1)))){
    return(list(success = FALSE, reason = "one or more windows empty"))
  }

  utm_epsg <- paste0("EPSG:",
                     if (latitude >= 0) "326" else "327",
                     sprintf("%02d", floor((longitude + 180) / 6) + 1))

  comps_utm <- lapply(composites, function(r){
    tryCatch(terra::project(r, utm_epsg, method = "bilinear"),
             error = function(e) NULL)
  })
  if (any(vapply(comps_utm, is.null, logical(1)))){
    return(list(success = FALSE, reason = "UTM projection failed"))
  }

  ref_r <- comps_utm[[1]]
  for (i in seq_along(comps_utm)){
    if (!terra::compareGeom(comps_utm[[i]], ref_r, stopOnError = FALSE)){
      comps_utm[[i]] <- terra::resample(comps_utm[[i]], ref_r, method = "bilinear")
    }
  }

  feature_stack <- Reduce(c, comps_utm)
  names(feature_stack) <- paste0("w", seq_len(num_windows))
  centre_v <- terra::project(terra::vect(cbind(longitude, latitude), crs = "EPSG:4326"), utm_epsg)

  list(success           = TRUE,
       feature_stack     = feature_stack,
       centre_v          = centre_v,
       utm_epsg          = utm_epsg,
       pixel_size_metres = as.numeric(terra::res(feature_stack)[1]),
       num_windows       = num_windows)
}


#' Check source and index name before running a query
#' @noRd
.checkConfigs <- function(start_date, index_name){
  d <- as.Date(start_date)
  if (is.na(d)) stop("start_date could not be parsed as a date: '", start_date, "'", call. = FALSE)

  if (d >= .S2_START){
    return(list("sentinel-2", .resolve_s2_index(index_name)))
  }
  if (d >= .L8_START){
    message("Sentinel-2 not available (start_date before ", format(.S2_START, "%Y-%m-%d"), "). Landsat 8 used.")
    return(list("landsat-8", .resolve_landsat_index(index_name)))
  }
  message("Sentinel-2 and Landsat 8 unavailable for start_date before ", format(.L8_START, "%Y-%m-%d"), "; Landsat 5 used.")
  list("landsat-5", .resolve_landsat_index(index_name))
}


#' @noRd
.resolve_s2_index <- function(index_name){
  il <- s2_index_list[[index_name]]
  if (is.null(il)) stop("Unknown index_name '", index_name, "'. Available: ", paste(names(s2_index_list), collapse = ", "), call. = FALSE)
  il
}

#' @noRd
.resolve_landsat_index <- function(index_name){
  il <- landsat_index_list[[index_name]]
  if (is.null(il)) stop("Unknown index_name '", index_name, "'. Available: ", paste(names(landsat_index_list), collapse = ", "), call. = FALSE)
  il
}

#' @noRd
.defaultMetres <- function(source){ if (source == "sentinel-2") .DEFAULT_S2 else .DEFAULT_LANDSAT }


#' Resolve plot dimensions for a single query (internal)
#' @noRd
.resolvePlotDims <- function(longitude,
                             latitude,
                             start_date,
                             end_date,
                             x_metres,
                             y_metres,
                             source,
                             field_mask,
                             modis_interpolation,
                             crop_search_params,
                             max_cloud_cover,
                             scl_classes
){
  mask_sf <- NULL

  if (!is.null(x_metres) && !is.null(y_metres)){
    dims <- list(x_metres = x_metres, y_metres = y_metres, dim_source = "user",
                 detected = NA, reason = NA_character_, mask_sf = NULL)

    if (isTRUE(field_mask) || inherits(field_mask, c("sf", "sfc"))){
      if (inherits(field_mask, c("sf", "sfc"))){
        mask_sf <- if (inherits(field_mask, "sf")) field_mask else sf::st_sf(geometry = field_mask)
      } else {
        params <- .mergeCropSearchParams(crop_search_params,
                                         start_date, end_date,
                                         max_cloud_cover, scl_classes,
                                         modis_interpolation)
        stack_fn <- if (isTRUE(params$enhanced_search)) .enhancedFeatureStack else .fetchFeatureStack
        fs <- stack_fn(longitude       = longitude,
                       latitude        = latitude,
                       start_date      = params$search_start_date,
                       end_date        = params$search_end_date,
                       fetch_metres    = max(x_metres, y_metres),
                       max_cloud_cover = params$search_max_cloud_cover,
                       scl_classes     = params$search_scl_classes)
        if (fs$success){
          det <- .nonparametricCropSearch(fs                 = fs,
                                          crop_search_m      = max(x_metres, y_metres),
                                          VI_sensitivity     = params$VI_sensitivity,
                                          VI_threshold       = params$VI_threshold,
                                          search_window      = params$search_windows,
                                          expand_search_area = params$expand_search_area)
          if (det$success) mask_sf <- det$mask_sf
          else message("  field mask failed: ", det$reason)
        }
      }
      dims$mask_sf <- mask_sf
    }
    return(dims)
  }

  det <- crop_search(longitude           = longitude,
                     latitude            = latitude,
                     start_date          = start_date,
                     end_date            = end_date,
                     crop_search_params  = crop_search_params,
                     modis_interpolation = modis_interpolation,
                     max_cloud_cover     = max_cloud_cover,
                     field_mask          = field_mask,
                     scl_classes         = scl_classes)

  list(x_metres = det$x_metres, y_metres = det$y_metres, dim_source = det$dim_source,
       detected = det$detected, reason = det$reason, mask_sf = det$mask_sf)
}


#' Fetch real + MODIS-interpolated scenes for time series (internal)
#' @noRd
.fetchFusedScenes <- function(bbox,
                              longitude,
                              latitude,
                              start_date,
                              end_date,
                              index_list,
                              source,
                              max_cloud_cover = 50,
                              scl_classes = NULL
){
  fine_scenes <- tryCatch(
    .fetchStac(bbox            = bbox,
               start_date      = start_date,
               end_date        = end_date,
               asset_names     = index_list$assets,
               index_function  = index_list$fun,
               source          = source,
               max_cloud_cover = max_cloud_cover,
               scl_classes     = scl_classes),
    error = function(e) NULL)
  n_fine <- length(fine_scenes %||% list())
  message(sprintf("    Real scenes fetched: %d", n_fine))

  modis_scenes <- .fetchMODIS(bbox            = bbox,
                              start_date      = start_date,
                              end_date        = end_date,
                              max_cloud_cover = max_cloud_cover)
  n_modis <- length(modis_scenes %||% list())
  message(sprintf("    MODIS scenes fetched: %d", n_modis))

  if (n_fine == 0 && n_modis == 0) return(NULL)
  if (n_fine == 0){
    message("    No real scenes available, fusion requires at least one anchor.")
    return(NULL)
  }

  utm_epsg <- paste0("EPSG:", if (latitude >= 0) "326" else "327",
                     sprintf("%02d", floor((longitude + 180) / 6) + 1))

  fine_utm <- lapply(fine_scenes, function(sc){
    r_utm <- tryCatch(terra::project(sc$rast, utm_epsg, method = "bilinear"),
                      error = function(e) NULL)
    if (is.null(r_utm)) return(NULL)
    list(date = sc$date, rast = r_utm, fused = FALSE,
         anchor_days = NA_integer_, residual = NULL, landsat_t1 = NULL)
  })
  fine_utm <- Filter(Negate(is.null), fine_utm)
  if (length(fine_utm) == 0) return(NULL)

  ref_r <- fine_utm[[1]]$rast
  for (i in seq_along(fine_utm)){
    if (!terra::compareGeom(fine_utm[[i]]$rast, ref_r, stopOnError = FALSE)){
      fine_utm[[i]]$rast <- terra::resample(fine_utm[[i]]$rast, ref_r, method = "bilinear")
    }
  }

  if (n_modis == 0){
    return(fine_utm[order(vapply(fine_utm, function(sc) as.numeric(sc$date), numeric(1)))])
  }

  fine_dates  <- vapply(fine_utm, function(sc) as.numeric(sc$date), numeric(1))
  modis_dates <- vapply(modis_scenes, function(sc) as.numeric(sc$date), numeric(1))
  tolerance_d <- 3L

  fused_scenes <- list()
  for (mi in seq_along(modis_scenes)){
    m_date <- modis_dates[mi]
    if (any(abs(fine_dates - m_date) <= tolerance_d)) next

    anchor <- fine_utm[[which.min(abs(fine_dates - m_date))]]
    m_t1   <- modis_scenes[[which.min(abs(modis_dates - as.numeric(anchor$date)))]]$rast
    fused  <- tryCatch(.fuseLandsatMODIS(anchor$rast, m_t1, modis_scenes[[mi]]$rast),
                       error = function(e) NULL)
    if (is.null(fused)) next

    fused_rast <- fused$rast
    if (!terra::compareGeom(fused_rast, ref_r, stopOnError = FALSE)){
      fused_rast <- terra::resample(fused_rast, ref_r, method = "bilinear")
    }
    resid_rast <- fused$residual
    anchor_t1  <- fused$landsat_t1
    if (!terra::compareGeom(resid_rast, ref_r, stopOnError = FALSE)){
      resid_rast <- terra::resample(resid_rast, ref_r, method = "bilinear")
      anchor_t1  <- terra::resample(anchor_t1,  ref_r, method = "bilinear")
    }
    fused_scenes[[length(fused_scenes) + 1]] <- list(
      date        = as.Date(m_date, origin = "1970-01-01"),
      rast        = fused_rast,
      fused       = TRUE,
      anchor_days = as.integer(abs(m_date - as.numeric(anchor$date))),
      residual    = resid_rast,
      landsat_t1  = anchor_t1)
  }

  message(sprintf("    Fused scenes created: %d", length(fused_scenes)))
  all_scenes <- c(fine_utm, fused_scenes)
  all_scenes[order(vapply(all_scenes, function(sc) as.numeric(sc$date), numeric(1)))]
}


#' Fetch MODIS EVI composites (internal)
#' @noRd
.fetchMODIS <- function(bbox,
                        start_date,
                        end_date,
                        max_cloud_cover = 50
){
  if (!requireNamespace("rstac", quietly = TRUE)){
    message("    rstac package required for MODIS fusion.")
    return(NULL)
  }
  aoi_bbox <- c(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]])

  items <- tryCatch({
    q <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
      rstac::stac_search(collections = "modis-13Q1-061",
                         bbox        = aoi_bbox,
                         datetime    = paste0(start_date, "/", end_date),
                         limit       = 100) |>
      rstac::post_request()
    if (is.null(q) || length(q$features) == 0) return(NULL)
    rstac::items_sign(q, rstac::sign_planetary_computer())
  }, error = function(e){
    message("    MODIS STAC query failed: ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(items) || length(items$features) == 0) return(NULL)

  features <- items$features
  is_terra <- vapply(features, function(f) grepl("^MOD", f$id), logical(1))
  if (any(is_terra)) features <- features[is_terra]

  dedup_key <- vapply(features, function(f){
    paste(substr(f$properties$start_datetime, 1, 10),
          f$properties[["modis:tile-id"]] %||% "",
          sep = "|")
  }, character(1))
  ids       <- vapply(features, `[[`, character(1), "id")
  proc_time <- sub(".*\\.(\\d+)$", "\\1", ids)
  features  <- features[order(proc_time, decreasing = TRUE)]
  dedup_key <- dedup_key[order(proc_time, decreasing = TRUE)]
  features  <- features[!duplicated(dedup_key)]

  scenes <- lapply(seq_along(features), function(f){
    evi_asset <- features[[f]]$assets[["250m_16_days_NDVI"]]
    if (is.null(evi_asset)) return(NULL)
    href <- evi_asset$href
    dt   <- tryCatch(as.Date(substr(features[[f]]$properties$start_datetime, 1, 10)),
                     error = function(e) NULL)
    if (length(dt) == 0) return(NULL)

    r <- tryCatch({
      evi_rast <- terra::rast(paste0("/vsicurl/", href))
      terra::scoff(evi_rast) <- NULL
      aoi_sfc     <- sf::st_as_sfc(bbox)
      aoi_in_crs  <- sf::st_transform(aoi_sfc, terra::crs(evi_rast))
      evi_rast    <- terra::crop(evi_rast, terra::ext(sf::st_bbox(aoi_in_crs)), snap = "out")
      evi_rast <- terra::ifel(evi_rast == -3000, NA, evi_rast)
      evi_rast <- terra::ifel(evi_rast < -2000 | evi_rast > 10000, NA, evi_rast)
      evi_rast <- evi_rast * 0.0001
      terra::project(evi_rast, "EPSG:4326", method = "bilinear")
    }, error = function(e) NULL)

    if (is.null(r)) return(NULL)
    list(date = dt, rast = r)
  })
  scenes <- Filter(Negate(is.null), scenes)
  if (length(scenes) == 0) return(NULL)
  scenes
}


#' Fuse a Landsat scene with MODIS temporal change (STARFM core, internal)
#' @noRd
.fuseLandsatMODIS <- function(landsat_t1, modis_t1, modis_t2){
  modis_t1_30m <- terra::project(modis_t1, landsat_t1, method = "bilinear")
  modis_t2_30m <- terra::project(modis_t2, landsat_t1, method = "bilinear")

  modis_delta <- modis_t2_30m - modis_t1_30m
  residual    <- abs(landsat_t1 - modis_t1_30m)
  max_resid   <- terra::global(residual, "max", na.rm = TRUE)[[1]]
  if (!is.finite(max_resid) || max_resid < 1e-6) max_resid <- 1

  weight    <- 1 - (residual / max_resid)
  weight    <- terra::ifel(weight < 0, 0, weight)
  predicted <- landsat_t1 + weight * modis_delta
  predicted <- terra::ifel(predicted < -0.2, -0.2, predicted)
  predicted <- terra::ifel(predicted >  1.0,  1.0, predicted)

  list(rast = predicted, residual = residual, landsat_t1 = landsat_t1)
}


#' Enhanced feature stack using STARFM fusion (internal)
#' @noRd
.enhancedFeatureStack <- function(longitude,
                                  latitude,
                                  start_date,
                                  end_date,
                                  fetch_metres,
                                  max_cloud_cover = 50,
                                  num_windows = 3L,
                                  scl_classes = NULL
){
  start_d <- suppressWarnings(as.Date(start_date))
  end_d   <- suppressWarnings(as.Date(end_date))
  if (is.na(start_d) || is.na(end_d)) return(list(success = FALSE, reason = "invalid dates"))
  if (as.numeric(end_d - start_d) < 2) return(list(success = FALSE, reason = "date range too short"))

  fetch_metres <- max(fetch_metres, (.DEFAULT_LANDSAT * 20))
  index_list   <- suppressMessages(.checkConfigs(start_date, "EVI2"))
  source       <- index_list[[1]]
  index_list   <- index_list[[2]]
  pixel_size_m <- index_list$res
  fetch_bbox   <- point_to_bbox(longitude, latitude, fetch_metres, fetch_metres)

  num_windows <- min(3L, trunc(as.numeric(end_d - start_d) / (pixel_size_m * 2) + 1))
  breaks      <- seq(start_d, end_d, length.out = num_windows + 1L)

  fallback <- function(){
    .fetchFeatureStack(longitude = longitude, latitude = latitude,
                       start_date = start_date, end_date = end_date,
                       fetch_metres = fetch_metres,
                       max_cloud_cover = max_cloud_cover,
                       num_windows = num_windows, scl_classes = scl_classes)
  }

  landsat_scenes <- tryCatch(
    .fetchStac(bbox            = fetch_bbox,
               start_date      = as.character(start_d),
               end_date        = as.character(end_d),
               asset_names     = index_list$assets,
               index_function  = index_list$fun,
               source          = source,
               max_cloud_cover = max_cloud_cover,
               scl_classes     = scl_classes),
    error = function(e) NULL)
  n_landsat <- length(landsat_scenes %||% list())

  modis_scenes <- .fetchMODIS(bbox            = fetch_bbox,
                              start_date      = as.character(start_d),
                              end_date        = as.character(end_d),
                              max_cloud_cover = max_cloud_cover)
  n_modis <- length(modis_scenes %||% list())

  if (n_modis == 0 || n_landsat == 0) return(fallback())

  utm_epsg <- paste0("EPSG:",
                     if (latitude >= 0) "326" else "327",
                     sprintf("%02d", floor((longitude + 180) / 6) + 1))

  ls_utm <- lapply(landsat_scenes, function(sc){
    r_utm <- tryCatch(terra::project(sc$rast, utm_epsg, method = "bilinear"),
                      error = function(e) NULL)
    if (is.null(r_utm)) return(NULL)
    list(date = sc$date, rast = r_utm)
  })
  ls_utm <- Filter(Negate(is.null), ls_utm)
  if (length(ls_utm) == 0) return(fallback())

  ref_ls <- ls_utm[[1]]$rast
  for (i in seq_along(ls_utm)){
    if (!terra::compareGeom(ls_utm[[i]]$rast, ref_ls, stopOnError = FALSE)){
      ls_utm[[i]]$rast <- terra::resample(ls_utm[[i]]$rast, ref_ls, method = "bilinear")
    }
  }

  ls_dates    <- vapply(ls_utm, function(sc) as.numeric(sc$date), numeric(1))
  modis_dates <- vapply(modis_scenes, function(sc) as.numeric(sc$date), numeric(1))

  composites <- lapply(seq_len(num_windows), function(w){
    w_start <- breaks[w]; w_end <- breaks[w + 1L]
    in_w <- ls_dates >= as.numeric(w_start) & ls_dates <= as.numeric(w_end)
    ls_in <- ls_utm[in_w]

    if (length(ls_in) > 0){
      stack <- Reduce(c, lapply(ls_in, `[[`, "rast"))
      return(if (terra::nlyr(stack) > 1) max(stack, na.rm = TRUE) else stack)
    }

    modis_in <- modis_scenes[modis_dates >= as.numeric(w_start) & modis_dates <= as.numeric(w_end)]
    if (length(modis_in) == 0) return(NULL)

    w_mid    <- as.numeric(w_start + (w_end - w_start) / 2)
    anchor   <- ls_utm[[which.min(abs(ls_dates - w_mid))]]
    modis_t1 <- modis_scenes[[which.min(abs(modis_dates - as.numeric(anchor$date)))]]$rast

    fused <- lapply(modis_in, function(m_sc){
      tryCatch(.fuseLandsatMODIS(anchor$rast, modis_t1, m_sc$rast)$rast,
               error = function(e) NULL)
    })
    fused <- Filter(Negate(is.null), fused)
    if (length(fused) == 0) return(NULL)
    stack <- Reduce(c, fused)
    if (terra::nlyr(stack) > 1) max(stack, na.rm = TRUE) else stack
  })

  if (any(vapply(composites, is.null, logical(1)))) return(fallback())

  ref_r <- composites[[1]]
  for (i in seq_along(composites)){
    if (!terra::compareGeom(composites[[i]], ref_r, stopOnError = FALSE)){
      composites[[i]] <- terra::resample(composites[[i]], ref_r, method = "bilinear")
    }
  }

  feature_stack <- Reduce(c, composites)
  names(feature_stack) <- paste0("w", seq_len(num_windows))
  centre_v <- terra::project(terra::vect(cbind(longitude, latitude), crs = "EPSG:4326"), utm_epsg)

  list(success           = TRUE,
       feature_stack     = feature_stack,
       centre_v          = centre_v,
       utm_epsg          = utm_epsg,
       pixel_size_metres = as.numeric(terra::res(feature_stack)[1]),
       num_windows       = num_windows)
}


#' Validate MODIS-Landsat fusion accuracy via leave-one-out holdout
#' @export
validate_data_fusion <- function(longitude,
                                 latitude,
                                 start_date,
                                 end_date,
                                 fetch_metres = 600,
                                 max_cloud_cover = 50
){
  fetch_bbox <- point_to_bbox(longitude, latitude, fetch_metres, fetch_metres)
  index_list <- suppressMessages(.checkConfigs(start_date, "EVI2"))
  source     <- index_list[[1]]
  idx        <- index_list[[2]]

  utm_epsg <- paste0("EPSG:",
                     if (latitude >= 0) "326" else "327",
                     sprintf("%02d", floor((longitude + 180) / 6) + 1))

  ls_raw <- tryCatch(
    .fetchStac(bbox            = fetch_bbox,
               start_date      = start_date,
               end_date        = end_date,
               asset_names     = idx$assets,
               index_function  = idx$fun,
               source          = source,
               max_cloud_cover = max_cloud_cover),
    error = function(e) NULL)

  if (is.null(ls_raw) || length(ls_raw) < 2) return(NULL)

  ls_scenes <- lapply(ls_raw, function(sc){
    r <- tryCatch(terra::project(sc$rast, utm_epsg, method = "bilinear"),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    list(date = sc$date, rast = r)
  })
  ls_scenes <- Filter(Negate(is.null), ls_scenes)
  if (length(ls_scenes) < 2) return(NULL)

  modis_scenes <- .fetchMODIS(bbox = fetch_bbox, start_date = start_date,
                              end_date = end_date, max_cloud_cover = max_cloud_cover)
  if (is.null(modis_scenes) || length(modis_scenes) == 0) return(NULL)

  ref_r <- ls_scenes[[1]]$rast
  for (i in seq_along(ls_scenes)){
    if (!terra::compareGeom(ls_scenes[[i]]$rast, ref_r, stopOnError = FALSE)){
      ls_scenes[[i]]$rast <- terra::resample(ls_scenes[[i]]$rast, ref_r, method = "bilinear")
    }
  }

  modis_dates <- vapply(modis_scenes, function(sc) as.numeric(sc$date), numeric(1))

  results <- lapply(seq_along(ls_scenes), function(holdout_i){
    target    <- ls_scenes[[holdout_i]]
    remaining <- ls_scenes[-holdout_i]
    rem_dates <- vapply(remaining, function(sc) as.numeric(sc$date), numeric(1))
    anchor    <- remaining[[which.min(abs(rem_dates - as.numeric(target$date)))]]

    m_t1 <- modis_scenes[[which.min(abs(modis_dates - as.numeric(anchor$date)))]]$rast
    m_t2 <- modis_scenes[[which.min(abs(modis_dates - as.numeric(target$date)))]]$rast

    predicted <- tryCatch(.fuseLandsatMODIS(anchor$rast, m_t1, m_t2),
                          error = function(e) NULL)
    if (is.null(predicted)) return(NULL)

    p <- predicted$rast
    if (!terra::compareGeom(p, target$rast, stopOnError = FALSE)){
      p <- terra::resample(p, target$rast, method = "bilinear")
    }

    obs   <- terra::values(target$rast, na.rm = FALSE)
    pred  <- terra::values(p, na.rm = FALSE)
    valid <- !is.na(obs) & !is.na(pred)
    if (sum(valid) < 4) return(NULL)

    rmse   <- sqrt(mean((obs[valid] - pred[valid])^2))
    ss_res <- sum((obs[valid] - pred[valid])^2)
    ss_tot <- sum((obs[valid] - mean(obs[valid]))^2)
    r_sq   <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

    data.frame(date = as.character(target$date), rmse = rmse, r_squared = r_sq,
               n_pixels = sum(valid), stringsAsFactors = FALSE)
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) return(NULL)

  do.call(rbind, results)
}


#' @noRd
.buildStacQuery <- function(bbox,
                            start_date,
                            end_date,
                            asset_names,
                            source,
                            max_cloud_cover = 50
){
  cfg <- .stacSourceConfig(source)
  list(
    endpoint        = cfg$endpoint,
    collection      = cfg$collection,
    bbox            = c(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]]),
    datetime        = paste0(start_date, "/", end_date),
    asset_names     = asset_names,
    cloud_key       = cfg$cloud_key,
    max_cloud_cover = max_cloud_cover,
    platforms       = cfg$platforms,
    scale_factor    = cfg$scale_factor,
    scale_offset    = cfg$scale_offset %||% 0,
    valid_range     = cfg$valid_range
  )
}


#' @noRd
.stacSourceConfig <- function(source){
  if (source == "sentinel-2") return(list(
    endpoint = .STAC_ENDPOINT, collection = "sentinel-2-l2a",
    cloud_key = "eo:cloud_cover", platforms = NULL,
    scale_factor = 0.0001, scale_offset = 0,
    valid_range = c(1, 10000), fill_value = 0))
  if (source == "landsat-8") return(list(
    endpoint = .STAC_ENDPOINT, collection = "landsat-c2-l2",
    cloud_key = "eo:cloud_cover", platforms = c("landsat-8", "landsat-9"),
    scale_factor = 0.0000275, scale_offset = -0.2,
    valid_range = c(7273, 43636), fill_value = 0))
  if (source == "landsat-5") return(list(
    endpoint = .STAC_ENDPOINT, collection = "landsat-c2-l2",
    cloud_key = "eo:cloud_cover", platforms = "landsat-5",
    scale_factor = 0.0000275, scale_offset = -0.2,
    valid_range = c(7273, 43636), fill_value = 0))
  stop("Unknown source: ", source, call. = FALSE)
}


#' @noRd
.executeStacQuery <- function(query){
  if (!requireNamespace("rstac", quietly = TRUE)){
    stop("rstac package required.", call. = FALSE)
  }
  items <- tryCatch({
    q <- rstac::stac(query$endpoint) |>
      rstac::stac_search(collections = query$collection,
                         bbox        = query$bbox,
                         datetime    = query$datetime,
                         limit       = 999) |>
      rstac::post_request()
    if (is.null(q) || length(q$features) == 0) return(NULL)
    rstac::items_sign(q, rstac::sign_planetary_computer())
  }, error = function(e){
    message("    STAC query failed: ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(items) || length(items$features) == 0) return(NULL)

  features <- items$features
  keep <- vapply(features, function(f){
    cc <- f$properties[[query$cloud_key]]
    if (is.null(cc) || cc > query$max_cloud_cover) return(FALSE)
    if (!is.null(query$platforms)){
      pf <- f$properties[["platform"]]
      if (is.null(pf) || !(pf %in% query$platforms)) return(FALSE)
    }
    TRUE
  }, logical(1))
  features <- features[keep]
  if (length(features) == 0) return(NULL)
  features
}


#' @noRd
.readStacScenes <- function(features,
                            asset_names,
                            index_function,
                            aoi_bbox,
                            source_cfg,
                            scl_classes = NULL
){
  aoi_sfc <- sf::st_as_sfc(aoi_bbox)
  mask_asset <- if (!is.null(scl_classes)){
    if (source_cfg$collection == "sentinel-2-l2a") "SCL" else "qa_pixel"
  } else NULL

  scenes <- lapply(features, function(f){
    dt <- tryCatch(as.Date(substr(f$properties$datetime %||%
                                    f$properties$start_datetime, 1, 10)),
                   error = function(e) NULL)
    if (is.null(dt) || is.na(dt)) return(NULL)

    read_band <- function(asset, scale = TRUE){
      href <- f$assets[[asset]]$href
      if (is.null(href)) return(NULL)
      b <- terra::rast(paste0("/vsicurl/", href))
      terra::scoff(b) <- NULL
      aoi_in_crs <- sf::st_transform(aoi_sfc, terra::crs(b))
      b <- terra::crop(b, terra::ext(sf::st_bbox(aoi_in_crs)), snap = "out")
      if (terra::ncell(b) == 0) return(NULL)
      if (scale){
        if (!is.null(source_cfg$valid_range)){
          b <- terra::ifel(b < source_cfg$valid_range[1] | b > source_cfg$valid_range[2], NA, b)
        }
        if (!is.null(source_cfg$scale_factor)){
          b <- b * source_cfg$scale_factor + (source_cfg$scale_offset %||% 0)
        }
      }
      names(b) <- asset
      b
    }

    r <- tryCatch({
      band_list <- lapply(asset_names, read_band, scale = TRUE)
      band_list <- Filter(Negate(is.null), band_list)
      if (length(band_list) == 0) return(NULL)
      bands <- do.call(c, band_list)

      if (!is.null(mask_asset)){
        mb <- tryCatch(read_band(mask_asset, scale = FALSE), error = function(e) NULL)
        if (!is.null(mb)){
          mb <- terra::resample(mb, bands, method = "near")
          cloud_mask <- Reduce(`|`, lapply(scl_classes, function(cls) mb == cls))
          bands <- terra::mask(bands, cloud_mask, maskvalue = TRUE)
        }
      }

      index_function(bands)
    }, error = function(e) NULL)

    if (is.null(r)) return(NULL)
    list(date = dt, rast = r)
  })
  Filter(Negate(is.null), scenes)
}


#' @noRd
.fetchStac <- function(bbox,
                       start_date,
                       end_date,
                       asset_names,
                       index_function,
                       source,
                       max_cloud_cover = 50,
                       scl_classes = NULL
){
  stopifnot("max_cloud_cover must be between 0 and 100" =
              max_cloud_cover >= 0 && max_cloud_cover <= 100)

  cfg   <- .stacSourceConfig(source)
  query <- .buildStacQuery(bbox            = bbox,
                           start_date      = start_date,
                           end_date        = end_date,
                           asset_names     = asset_names,
                           source          = source,
                           max_cloud_cover = max_cloud_cover)

  features <- .executeStacQuery(query)
  if (is.null(features)) return(NULL)

  scenes <- .readStacScenes(features       = features,
                            asset_names    = asset_names,
                            index_function = index_function,
                            aoi_bbox       = bbox,
                            source_cfg     = cfg,
                            scl_classes    = scl_classes)
  if (length(scenes) == 0) return(NULL)
  scenes
}


#' @noRd
.prepareFetch <- function(longitude,
                          latitude,
                          start_date,
                          end_date,
                          index_name,
                          x_metres,
                          y_metres,
                          field_mask,
                          modis_interpolation,
                          crop_search_params,
                          max_cloud_cover,
                          scl_classes
){
  index_list <- .checkConfigs(start_date, index_name)
  source     <- index_list[[1]]
  index_list <- index_list[[2]]

  dims <- .resolvePlotDims(longitude           = longitude,
                           latitude            = latitude,
                           start_date          = start_date,
                           end_date            = end_date,
                           x_metres            = x_metres,
                           y_metres            = y_metres,
                           source              = source,
                           field_mask          = field_mask,
                           modis_interpolation = modis_interpolation,
                           crop_search_params  = crop_search_params,
                           max_cloud_cover     = max_cloud_cover,
                           scl_classes         = scl_classes)

  bbox       <- point_to_bbox(longitude, latitude, dims$x_metres, dims$y_metres)
  fetch_x    <- max(dims$x_metres, .DEFAULT_LANDSAT * 10)
  fetch_y    <- max(dims$y_metres, .DEFAULT_LANDSAT * 10)
  fetch_bbox <- point_to_bbox(longitude, latitude, fetch_x, fetch_y)

  list(source     = source,
       index_list = index_list,
       bbox       = bbox,
       fetch_bbox = fetch_bbox,
       mask_sf    = dims$mask_sf)
}


#' @importFrom terra global mask project vect crs
#' @noRd
.aggregateStacValues <- function(scenes,
                                 bbox,
                                 mask_sf = NULL,
                                 field_mask,
                                 time_aggregate = "daily",
                                 aggregate_function = "mean",
                                 merge_mode = "off"
){
  if (is.null(scenes) || length(scenes) == 0){
    return(data.frame(date = character(0), value = numeric(0), stringsAsFactors = FALSE))
  }

  mask_v   <- if (!isFALSE(field_mask) && !is.null(mask_sf)) terra::vect(mask_sf) else NULL
  plot_sfc <- sf::st_as_sfc(bbox)

  has_fusion_meta <- !is.null(scenes[[1]]$fused)

  dates <- vapply(scenes, function(sc) as.character(sc$date), character(1))
  dates <- switch(time_aggregate,
                  "monthly" = , "month"  = substr(dates, 1, 7),
                  "yearly"  = , "year"   = , "annual" = substr(dates, 1, 4),
                  dates)

  scene_stats <- lapply(scenes, function(sc){
    r <- sc$rast
    plot_in_crs <- sf::st_transform(plot_sfc, terra::crs(r))
    r <- terra::crop(r, terra::ext(sf::st_bbox(plot_in_crs)), snap = "out")
    if (terra::ncell(r) == 0){
      return(list(value = NA_real_, mae = NA_real_, mape = NA_real_))
    }
    mask_proj <- if (!is.null(mask_v)) terra::project(mask_v, terra::crs(r)) else NULL
    if (!is.null(mask_proj)) r <- terra::mask(r, mask_proj)
    value <- terra::global(r, aggregate_function, na.rm = TRUE)[[1]]

    if (isTRUE(sc$fused) && !is.null(sc$residual)){
      resid <- terra::crop(sc$residual, terra::ext(sf::st_bbox(plot_in_crs)), snap = "out")
      lst1  <- terra::crop(sc$landsat_t1, terra::ext(sf::st_bbox(plot_in_crs)), snap = "out")
      if (!is.null(mask_proj)){
        resid <- terra::mask(resid, mask_proj)
        lst1  <- terra::mask(lst1,  mask_proj)
      }
      mae     <- as.numeric(terra::global(resid, "mean", na.rm = TRUE)[[1]])
      mean_ls <- as.numeric(terra::global(lst1,  "mean", na.rm = TRUE)[[1]])
      mape    <- if (is.finite(mean_ls) && abs(mean_ls) > 1e-6) mae / abs(mean_ls) else NA_real_
      return(list(value = value, mae = mae, mape = mape))
    }
    list(value = value, mae = NA_real_, mape = NA_real_)
  })

  values <- vapply(scene_stats, `[[`, numeric(1), "value")

  df <- if (has_fusion_meta){
    data.frame(
      date        = dates,
      value       = values,
      MODIS       = vapply(scenes, `[[`, logical(1), "fused"),
      anchor_days = vapply(scenes, `[[`, integer(1), "anchor_days"),
      MAE         = vapply(scene_stats, `[[`, numeric(1), "mae"),
      MAPE        = vapply(scene_stats, `[[`, numeric(1), "mape"),
      stringsAsFactors = FALSE)
  } else {
    data.frame(date = dates, value = values, stringsAsFactors = FALSE)
  }

  df <- df[!is.na(df$value) & !is.nan(df$value), , drop = FALSE]

  if (has_fusion_meta){
    df <- df[order(df$date, df$MODIS), , drop = FALSE]
    if (!identical(merge_mode, "sep")) df <- df[!duplicated(df$date), , drop = FALSE]
    if (!identical(time_aggregate, "daily")){
      df <- stats::aggregate(
        cbind(value, anchor_days, MAE, MAPE) ~ date + MODIS,
        data = df, FUN = function(x) mean(x, na.rm = TRUE), na.action = stats::na.pass)
    }
  } else {
    df <- stats::aggregate(value ~ date, data = df, FUN = aggregate_function)
  }

  df
}

################################################################################
# get_rs_xxxx functions
################################################################################

#' Collect remote sensing data
#'
#' Extracts spectral index values from Sentinel-2, Landsat 8, or Landsat 5 over
#' a user-specified location and date range. Optionally fills gaps with
#' MODIS-fused synthetic scenes, auto-detects field dimensions from imagery,
#' and masks pixels outside the detected field polygon.
#'
#' Source is chosen automatically from \code{start_date}: Sentinel-2 from
#' 2015-07-01, Landsat 8 from 2013-04-01, Landsat 5 before that.
#'
#' @details
#' ## Field dimensions
#'
#' If both \code{x_metres} and \code{y_metres} are \code{NULL} (default), field
#' dimensions are auto-detected via \code{\link{crop_search}}. Detection
#' behaviour is controlled through \code{crop_search_params}; see
#' \code{\link{crop_search_defaults}} for the full list.
#'
#' Supplying both \code{x_metres} and \code{y_metres} disables detection.
#'
#' ## MODIS interpolation
#'
#' \code{modis_interpolation} has three states:
#' \describe{
#'   \item{\code{"off"} (default)}{No MODIS used for output. Detection uses Landsat-only feature stack unless \code{crop_search_params$enhanced_search = TRUE}.}
#'   \item{\code{"merge"}}{MODIS used for detection (via auto-inherited \code{enhanced_search}) AND output, collapsing shared dates with real scenes taking precedence.}
#'   \item{\code{"sep"}}{MODIS used for detection AND output, keeping real and fused scenes as separate rows.}
#' }
#'
#' When \code{"merge"} or \code{"sep"}, the output includes \code{MODIS},
#' \code{anchor_days}, \code{MAE}, and \code{MAPE} columns. MAE/MAPE are
#' computed on the field-masked region when \code{field_mask = TRUE}, else on
#' the full fetch bbox.
#'
#' @param longitude,latitude Coordinate in WGS84.
#' @param x_metres,y_metres Bounding box dimensions in metres. If both
#'   \code{NULL}, field dimensions are auto-detected.
#' @param start_date,end_date Date range, "YYYY-MM-DD".
#' @param index_name Spectral index name from \code{\link{s2_index_list}}.
#' @param field_mask \code{FALSE}, \code{TRUE}, or an \code{sf}/\code{sfc} polygon.
#' @param crop_search_params Named list controlling detection behaviour. See
#'   \code{\link{crop_search_defaults}} and \code{\link{crop_search}}.
#' @param modis_interpolation One of \code{"off"} (default), \code{"merge"},
#'   or \code{"sep"}.
#' @param max_cloud_cover Maximum scene cloud cover (0-100). Default 50.
#' @param scl_classes Integer vector of S2 SCL classes to mask.
#' @param time_aggregate \code{"daily"}, \code{"monthly"}, or \code{"yearly"}.
#' @param aggregate_function Spatial aggregation function name. Default \code{"mean"}.
#' @return Data frame with \code{date} and \code{value}, plus fusion columns
#'   when \code{modis_interpolation} is \code{"merge"} or \code{"sep"}.
#' @examples
#' ## Auto-detected field, standard output
#' ndvi_df <- get_rs_data(longitude = 138.715,
#'                        latitude = -34.904,
#'                        start_date = "2023-01-01",
#'                        end_date = "2023-01-31",
#'                        index_name = "NDVI")
#'
#' ## Custom detection parameters
#' params <- crop_search_defaults()
#' params$VI_sensitivity <- 0.15
#' params$expand_search_area <- TRUE
#' ndvi_df <- get_rs_data(longitude = 138.715,
#'                        latitude = -34.904,
#'                        start_date = "2023-01-01",
#'                        end_date = "2023-01-31",
#'                        index_name = "NDVI",
#'                        crop_search_params = params)
#' @export
get_rs_data <- function(longitude,
                        latitude,
                        x_metres = NULL,
                        y_metres = NULL,
                        start_date,
                        end_date,
                        index_name,
                        field_mask = FALSE,
                        crop_search_params = list(),
                        modis_interpolation = "off",
                        max_cloud_cover = 50,
                        scl_classes = NULL,
                        time_aggregate = "daily",
                        aggregate_function = "mean"
){
  if (xor(is.null(x_metres), is.null(y_metres))){
    stop("x_metres and y_metres must both be set, or both left NULL.", call. = FALSE)
  }
  if (!identical(modis_interpolation, "off") && !identical(modis_interpolation, "merge") &&
      !identical(modis_interpolation, "sep")){
    stop("modis_interpolation must be \"off\", \"merge\", or \"sep\".", call. = FALSE)
  }
  if (!identical(modis_interpolation, "off") && index_name %in% .MODIS_UNSUPPORTED_INDICES){
    warning("MODIS interpolation not supported for index '", index_name,
            "'. Falling back to standard fetch.", call. = FALSE)
    modis_interpolation <- "off"
  }

  fetch_list <- .prepareFetch(longitude           = longitude,
                              latitude            = latitude,
                              start_date          = start_date,
                              end_date            = end_date,
                              index_name          = index_name,
                              x_metres            = x_metres,
                              y_metres            = y_metres,
                              field_mask          = field_mask,
                              modis_interpolation = modis_interpolation,
                              crop_search_params  = crop_search_params,
                              max_cloud_cover     = max_cloud_cover,
                              scl_classes         = scl_classes)
  if (is.null(fetch_list)){
    return(data.frame(date = character(0), value = numeric(0), stringsAsFactors = FALSE))
  }

  scenes <- if (identical(modis_interpolation, "merge") || identical(modis_interpolation, "sep")){
    .fetchFusedScenes(bbox            = fetch_list$fetch_bbox,
                      longitude       = longitude,
                      latitude        = latitude,
                      start_date      = start_date,
                      end_date        = end_date,
                      index_list      = fetch_list$index_list,
                      source          = fetch_list$source,
                      max_cloud_cover = max_cloud_cover,
                      scl_classes     = scl_classes)
  } else {
    .fetchStac(bbox            = fetch_list$fetch_bbox,
               start_date      = start_date,
               end_date        = end_date,
               asset_names     = fetch_list$index_list$assets,
               index_function  = fetch_list$index_list$fun,
               source          = fetch_list$source,
               max_cloud_cover = max_cloud_cover,
               scl_classes     = scl_classes)
  }

  .aggregateStacValues(scenes             = scenes,
                       bbox               = fetch_list$bbox,
                       mask_sf            = fetch_list$mask_sf,
                       field_mask         = field_mask,
                       time_aggregate     = time_aggregate,
                       aggregate_function = aggregate_function,
                       merge_mode         = modis_interpolation)
}


#' Collect remote sensing data as a composite raster
#'
#' @importFrom terra rast
#' @inheritParams get_rs_data
#' @param composite_function One of \code{"median"} (default), \code{"mean"},
#'   \code{"max"}, \code{"min"}, \code{"sum"}.
#' @return SpatRaster, or NULL if no scenes available.
#' @export
get_rs_raster <- function(longitude,
                          latitude,
                          x_metres = NULL,
                          y_metres = NULL,
                          start_date,
                          end_date,
                          index_name,
                          field_mask = FALSE,
                          crop_search_params = list(),
                          modis_interpolation = "off",
                          max_cloud_cover = 50,
                          scl_classes = NULL,
                          composite_function = "median"
){
  if (xor(is.null(x_metres), is.null(y_metres))){
    stop("x_metres and y_metres must both be set, or both left NULL.", call. = FALSE)
  }
  if (!identical(modis_interpolation, "off") && !identical(modis_interpolation, "merge") &&
      !identical(modis_interpolation, "sep")){
    stop("modis_interpolation must be \"off\", \"merge\", or \"sep\".", call. = FALSE)
  }
  if (!identical(modis_interpolation, "off") && index_name %in% .MODIS_UNSUPPORTED_INDICES){
    warning("MODIS interpolation not supported for index '", index_name,
            "'. Falling back to standard fetch.", call. = FALSE)
    modis_interpolation <- "off"
  }

  fetch_list <- .prepareFetch(longitude           = longitude,
                              latitude            = latitude,
                              start_date          = start_date,
                              end_date            = end_date,
                              index_name          = index_name,
                              x_metres            = x_metres,
                              y_metres            = y_metres,
                              field_mask          = field_mask,
                              modis_interpolation = modis_interpolation,
                              crop_search_params  = crop_search_params,
                              max_cloud_cover     = max_cloud_cover,
                              scl_classes         = scl_classes)
  if (is.null(fetch_list)) return(NULL)

  scenes <- if (identical(modis_interpolation, "merge") || identical(modis_interpolation, "sep")){
    .fetchFusedScenes(bbox            = fetch_list$fetch_bbox,
                      longitude       = longitude,
                      latitude        = latitude,
                      start_date      = start_date,
                      end_date        = end_date,
                      index_list      = fetch_list$index_list,
                      source          = fetch_list$source,
                      max_cloud_cover = max_cloud_cover,
                      scl_classes     = scl_classes)
  } else {
    .fetchStac(bbox            = fetch_list$fetch_bbox,
               start_date      = start_date,
               end_date        = end_date,
               asset_names     = fetch_list$index_list$assets,
               index_function  = fetch_list$index_list$fun,
               source          = fetch_list$source,
               max_cloud_cover = max_cloud_cover,
               scl_classes     = scl_classes)
  }
  if (is.null(scenes) || length(scenes) == 0) return(NULL)

  stack_r   <- Reduce(c, lapply(scenes, `[[`, "rast"))
  composite <- switch(composite_function,
                      "median" = terra::median(stack_r, na.rm = TRUE),
                      "mean"   = terra::mean(stack_r, na.rm = TRUE),
                      "max"    = max(stack_r, na.rm = TRUE),
                      "min"    = min(stack_r, na.rm = TRUE),
                      "sum"    = sum(stack_r, na.rm = TRUE),
                      terra::median(stack_r, na.rm = TRUE))

  plot_sfc    <- sf::st_as_sfc(fetch_list$bbox)
  plot_in_crs <- sf::st_transform(plot_sfc, terra::crs(composite))
  composite   <- terra::crop(composite, terra::ext(sf::st_bbox(plot_in_crs)), snap = "out")

  if (!is.null(fetch_list$mask_sf)){
    mask_v    <- terra::project(terra::vect(fetch_list$mask_sf), terra::crs(composite))
    composite <- terra::mask(composite, mask_v)
  }
  composite
}
