#' Sentinel-2 spectral indices list
#'
#' A named list of spectral indices compatible with the
#' \href{https://github.com/awesome-spectral-indices/awesome-spectral-indices}{Awesome
#' Spectral Indices} catalogue, for use with \code{\link{get_rs_data}} and
#' \code{\link{get_rs_raster}}.
#'
#' @format A named list where each element is a list with components:
#' \describe{
#'   \item{assets}{Character vector of Sentinel-2 band names.}
#'   \item{fun}{A function that takes a SpatRaster and returns a single-layer SpatRaster.}
#'   \item{res}{Recommended resolution (highest integer value used by assets).}
#' }
#' @examples
#' names(s2_index_list)
#'
#' ndvi_df <- get_rs_data(longitude = 138.715,
#'                       latitude = -34.904,
#'                       x_metres = 30,
#'                       y_metres = 30,
#'                       start_date = "2023-01-01",
#'                       end_date = "2023-01-10",
#'                       index_name = "NDVI")
#'
#' @export
s2_index_list <- list(
  NDVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04),
    res = 10 ),

  DVI = list(
    assets = c("B08", "B04"),
    fun = function(r) r$B08 - r$B04,
    res = 10 ),

  kNDVI = list(
    assets = c("B08", "B04"),
    fun = function(r) tanh( ((r$B08 - r$B04) / (r$B08 + r$B04))^2 ),
    res = 10 ),

  GNDVI = list(
    assets = c("B08", "B03"),
    fun = function(r) (r$B08 - r$B03) / (r$B08 + r$B03),
    res = 10 ),

  BNDVI = list(
    assets = c("B08", "B02"),
    fun = function(r) (r$B08 - r$B02) / (r$B08 + r$B02),
    res = 10 ),

  GBNDVI = list(
    assets = c("B08", "B03", "B02"),
    fun = function(r) (r$B08 - (r$B03 + r$B02)) / (r$B08 + (r$B03 + r$B02)),
    res = 10 ),

  GRNDVI = list(
    assets = c("B08", "B03", "B04"),
    fun = function(r) (r$B08 - (r$B03 + r$B04)) / (r$B08 + (r$B03 + r$B04)),
    res = 10 ),

  NIRv = list(
    assets = c("B08", "B04"),
    fun = function(r) ((r$B08 - r$B04) / (r$B08 + r$B04)) * r$B08,
    res = 10 ),

  NDWI = list(
    assets = c("B03", "B08"),
    fun = function(r) (r$B03 - r$B08) / (r$B03 + r$B08),
    res = 10 ),

  NDMI = list(
    assets = c("B08", "B11"),
    fun = function(r) (r$B08 - r$B11) / (r$B08 + r$B11),
    res = 20 ),

  EVI = list(
    assets = c("B08", "B04", "B02"),
    fun = function(r) 2.5 * (r$B08 - r$B04) / (r$B08 + 6 * r$B04 - 7.5 * r$B02 + 1),
    res = 10 ),

  SAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + 0.5) * 1.5,
    res = 10 ),

  MSAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (2 * r$B08 + 1 - sqrt((2 * r$B08 + 1)^2 - 8 * (r$B08 - r$B04))) / 2,
    res = 10 ),

  NDRE = list(
    assets = c("B08", "B05"),
    fun = function(r) (r$B08 - r$B05) / (r$B08 + r$B05),
    res = 20 ),

  MTCI = list(
    assets = c("B06", "B05", "B04"),
    fun = function(r) (r$B06 - r$B05) / (r$B05 - r$B04),
    res = 20 ),

  CIre_gitelson = list(
    assets = c("B07", "B05"),
    fun = function(r) (r$B07 / r$B05) - 1,
    res = 20 ),

  CIRE = list(
    assets = c("B08", "B05"),
    fun = function(r) (r$B08 / r$B05) - 1,
    res = 20 ),

  PSRI = list(
    assets = c("B06", "B04", "B02"),
    fun = function(r) (r$B04 - r$B02) / r$B06,
    res = 20 ),

  BSI = list(
    assets = c("B11", "B08", "B04", "B02"),
    fun = function(r) ((r$B11 + r$B04) - (r$B08 + r$B02)) / ((r$B11 + r$B04) + (r$B08 + r$B02)),
    res = 20 ),

  SR = list(
    assets = c("B08", "B04"),
    fun = function(r) r$B08 / r$B04,
    res = 10 ),

  EVI2 = list(
    assets = c("B08", "B04"),
    fun = function(r) 2.5 * (r$B08 - r$B04) / (r$B08 + 2.4 * r$B04 + 1),
    res = 10 ),

  OSAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + 0.16),
    res = 10 ),

  ARVI = list(
    assets = c("B08", "B04", "B02"),
    fun = function(r) (r$B08 - (2 * r$B04 - r$B02)) / (r$B08 + (2 * r$B04 - r$B02)),
    res = 10 ),

  GARI = list(
    assets = c("B08", "B03", "B02", "B04"),
    fun = function(r) (r$B08 - (r$B03 - (r$B02 - r$B04))) / (r$B08 + (r$B03 - (r$B02 - r$B04))),
    res = 10 ),

  CIG = list(
    assets = c("B08", "B03"),
    fun = function(r) (r$B08 / r$B03) - 1,
    res = 10 ),

  MNDWI = list(
    assets = c("B03", "B11"),
    fun = function(r) (r$B03 - r$B11) / (r$B03 + r$B11),
    res = 20 ),

  NBR = list(
    assets = c("B08", "B12"),
    fun = function(r) (r$B08 - r$B12) / (r$B08 + r$B12),
    res = 20 ),

  NBR2 = list(
    assets = c("B11", "B12"),
    fun = function(r) (r$B11 - r$B12) / (r$B11 + r$B12),
    res = 20 ),

  BAI = list(
    assets = c("B04", "B08"),
    fun = function(r) 1 / ((0.1 - r$B04)^2 + (0.06 - r$B08)^2),
    res = 10 ),

  MCARI = list(
    assets = c("B05", "B04", "B03"),
    fun = function(r) ((r$B05 - r$B04) - 0.2 * (r$B05 - r$B03)) * (r$B05 / r$B04),
    res = 20 ),

  IRECI = list(
    assets = c("B07", "B04", "B05", "B06"),
    fun = function(r) (r$B07 - r$B04) / (r$B05 / r$B06),
    res = 20 )
)

#' Landsat spectral indices list
#'
#' A named list of spectral indices compatible with Landsat 5 & 8 sources.
#'
#' @format A named list where each element is a list with components:
#' \describe{
#'   \item{assets}{Character vector of Sentinel-2 band names.}
#'   \item{fun}{A function that takes a SpatRaster and returns a single-layer SpatRaster.}
#'   \item{res}{Recommended resolution (highest integer value used by assets).}
#' }
#' @examples
#' names(landsat_index_list)
#'
#' ndvi_df <- get_rs_data(longitude = 138.715,
#'                       latitude = -34.904,
#'                       x_metres = 30,
#'                       y_metres = 30,
#'                       start_date = "2023-01-01",
#'                       end_date = "2023-01-10",
#'                       index_name = "NDVI")
#'
#' @export
landsat_index_list <- list(
  NDVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red),
    res = 30 ),

  DVI = list(
    assets = c("nir08", "red"),
    fun = function(r) r$nir08 - r$red,
    res = 30 ),

  kNDVI = list(
    assets = c("nir08", "red"),
    fun = function(r) tanh( ((r$nir08 - r$red) / (r$nir08 + r$red))^2 ),
    res = 30 ),

  GNDVI = list(
    assets = c("nir08", "green"),
    fun = function(r) (r$nir08 - r$green) / (r$nir08 + r$green),
    res = 30 ),

  BNDVI = list(
    assets = c("nir08", "blue"),
    fun = function(r) (r$nir08 - r$blue) / (r$nir08 + r$blue),
    res = 30 ),

  GBNDVI = list(
    assets = c("nir08", "green", "blue"),
    fun = function(r) (r$nir08 - (r$green + r$blue)) / (r$nir08 + (r$green + r$blue)),
    res = 30 ),

  GRNDVI = list(
    assets = c("nir08", "green", "red"),
    fun = function(r) (r$nir08 - (r$green + r$red)) / (r$nir08 + (r$green + r$red)),
    res = 30 ),

  NIRv = list(
    assets = c("nir08", "red"),
    fun = function(r) ((r$nir08 - r$red) / (r$nir08 + r$red)) * r$nir08,
    res = 30 ),

  NDWI = list(
    assets = c("green", "nir08"),
    fun = function(r) (r$green - r$nir08) / (r$green + r$nir08),
    res = 30 ),

  NDMI = list(
    assets = c("nir08", "swir16"),
    fun = function(r) (r$nir08 - r$swir16) / (r$nir08 + r$swir16),
    res = 30 ),

  EVI = list(
    assets = c("nir08", "red", "blue"),
    fun = function(r) 2.5 * (r$nir08 - r$red) / (r$nir08 + 6 * r$red - 7.5 * r$blue + 1),
    res = 30 ),

  SAVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red + 0.5) * 1.5,
    res = 30 ),

  MSAVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (2 * r$nir08 + 1 - sqrt((2 * r$nir08 + 1)^2 - 8 * (r$nir08 - r$red))) / 2,
    res = 30 ),

  BSI = list(
    assets = c("swir16", "nir08", "red", "blue"),
    fun = function(r) ((r$swir16 + r$red) - (r$nir08 + r$blue)) / ((r$swir16 + r$red) + (r$nir08 + r$blue)),
    res = 30 ),

  SR = list(
    assets = c("nir08", "red"),
    fun = function(r) r$nir08 / r$red,
    res = 30 ),

  EVI2 = list(
    assets = c("nir08", "red"),
    fun = function(r) 2.5 * (r$nir08 - r$red) / (r$nir08 + 2.4 * r$red + 1),
    res = 30 ),

  OSAVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red + 0.16),
    res = 30 ),

  ARVI = list(
    assets = c("nir08", "red", "blue"),
    fun = function(r) (r$nir08 - (2 * r$red - r$blue)) / (r$nir08 + (2 * r$red - r$blue)),
    res = 30 ),

  GARI = list(
    assets = c("nir08", "green", "blue", "red"),
    fun = function(r) (r$nir08 - (r$green - (r$blue - r$red))) / (r$nir08 + (r$green - (r$blue - r$red))),
    res = 30 ),

  CIG = list(
    assets = c("nir08", "green"),
    fun = function(r) (r$nir08 / r$green) - 1,
    res = 30 ),

  MNDWI = list(
    assets = c("green", "swir16"),
    fun = function(r) (r$green - r$swir16) / (r$green + r$swir16),
    res = 30 ),

  NBR = list(
    assets = c("nir08", "swir22"),
    fun = function(r) (r$nir08 - r$swir22) / (r$nir08 + r$swir22),
    res = 30 ),

  NBR2 = list(
    assets = c("swir16", "swir22"),
    fun = function(r) (r$swir16 - r$swir22) / (r$swir16 + r$swir22),
    res = 30 ),

  BAI = list(
    assets = c("red", "nir08"),
    fun = function(r) 1 / ((0.1 - r$red)^2 + (0.06 - r$nir08)^2),
    res = 30 )

)

################################################################################
# Constants
################################################################################
# Landsat native pixel resolution.
.PIXEL_SIZE_LANDSAT <- 30

# Default plot half-widths (metres) when user doesn't supply x_metres/y_metres or auto-detection is disabled/fails
.DEFAULT_HALF_S2      <- 10  # = 20 m × 20 m bbox
.DEFAULT_HALF_LANDSAT <- 15  # = 30 m × 30 m bbox

# Minimum bbox (half) width used for STAC downloads. Smaller = inefficient: GDAL fetches entire COG tile blocks.
.FETCH_HALF_WIDTH_MIN <- 150  # = 300 m × 300 m download window

# Operational start dates for source routing.
.S2_START <- as.Date("2015-07-01")
.L8_START <- as.Date("2013-04-01")

# Sentinel-2 -> Landsat (planetary_computer_v1) asset name map.
# Red-edge bands (B05/B06/B07) have no Landsat equivalent.
# B08 and B8A both map to nir08; see .translateToLandsat for collision handling.
.S2_TO_LANDSAT <- c(
  B02 = "blue",
  B03 = "green",
  B04 = "red",
  B08 = "nir08",
  B8A = "nir08",
  B11 = "swir16",
  B12 = "swir22"
)

.S2_PROVIDER <- "planetary_computer_v1"

################################################################################
# Helpers
################################################################################

#' Convert metres to decimal degrees
#'
#' @param latitude Latitude of centre point where the distance is being measured (WGS84).
#' @param distance Distance in metres. Default 1.
#' @param pixel_size Pixel resolution in metres. Default 1.
#' @return A vector containing x & y pixel sizes.
#' @export
metres_to_degrees <- function(latitude, distance = 1, pixel_size = 1){
  a <- 6378137.0
  e2 <- 0.00669437999014
  lat_c_rad <- latitude * pi / 180
  m_c <- a * (1 - e2) / (1 - e2 * sin(lat_c_rad)^2)^1.5
  half_lat <- (metres / 2) * 180 / (pi * m_c)
  lat_rad <- (latitude + c(-half_lat, 0, half_lat)) * pi / 180
  denom <- 1 - e2 * sin(lat_rad)^2
  n <- a / sqrt(denom)
  m <- a * (1 - e2) / denom^1.5
  w <- c(1, 4, 1)
  x_y_decimal_degrees <- c(x_degrees = pixel_size * 180 / (pi * sum(n * cos(lat_rad) * w) / 6),
                             y_degrees = pixel_size * 180 / (pi * sum(m * w) / 6))
  return(x_y_decimal_degrees)
}


#' Bounding box around a central point
#'
#' @param longitude Longitude in decimal degrees (WGS84).
#' @param latitude Latitude in decimal degrees (WGS84).
#' @param x_metres Bbox x width in metres.
#' @param y_metres Bbox y length in metres.
#' @return An sf bbox object.
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
  bbox <- sf::st_bbox(bbox, crs = sf::st_crs(4326))

  return(bbox)
}


#' Check source and index name before running a query
#'
#' Selects the STAC source based on start_date and resolves index name matches available indices.
#' @noRd
.checkConfigs <- function(start_date, index_name){
  d <- as.Date(start_date)
  if (is.na(d)) {
    stop("start_date could not be parsed as a date: '", start_date, "'", call. = FALSE)
  }
  if (d >= .S2_START){
    index_list <- .resolve_s2_index(index_name)
    return( list("sentinel-2", index_list) )
  }
  if (d >= .L8_START) {
    message("Sentinel-2 not available (start_date before ",
            format(.S2_START, "%Y-%m-%d"), "). Landsat 8 data obtained instead.")
    index_list <- .resolve_landsat_index(index_name)
    return( list("landsat-8", index_list) )
  }
  message("Sentinel-2 and Landsat 8 operations unavailable for start_date (prior to ",
          format(.L8_START, "%Y-%m-%d"), "); Landsat 5 data obtained instead.")
  index_list <- .resolve_landsat_index(index_name)
  return( list("landsat-5", index_list) )
}


#' Resolve an index_name to its s2_index_list entry
#' @noRd
.resolve_s2_index <- function(index_name){
  index_list <- s2_index_list[[index_name]]
  if (is.null(index_list)) {
    stop("Unknown index_name '", index_name, "'. Available indices: ",
         paste(names(s2_index_list), collapse = ", "), call. = FALSE)
  }
  return(index_list)
}


#' Resolve an index_name to its landsat_index_list entry
#' @noRd
.resolve_landsat_index <- function(index_name){
  index_list <- landsat_index_list[[index_name]]
  if (is.null(index_list)) {
    stop("Unknown index_name '", index_name, "'. Available indices: ",
         paste(names(landsat_index_list), collapse = ", "), call. = FALSE)
  }
  return(index_list)
}


#' Scene-level cloud cover filter (Sentinel-2)
#' @noRd
.cloudFilter <- function(max_cloud_cover){
  function(items, bbox, ...) {
    keep <- vapply(items[["features"]], function(item) {
      cc <- item[["properties"]][["eo:cloud_cover"]]
      !is.null(cc) && cc <= max_cloud_cover
    }, logical(1))
    items[["features"]] <- items[["features"]][keep]
    items
  }
}


#' Combined cloud + platform filter (Landsat)
#' @noRd
.landsatFilter <- function(max_cloud_cover, platforms){
  function(items, bbox, ...) {
    keep <- vapply(items[["features"]], function(item) {
      cc <- item[["properties"]][["eo:cloud_cover"]]
      pf <- item[["properties"]][["platform"]]
      (!is.null(cc) && cc <= max_cloud_cover) &&
        (!is.null(pf) && pf %in% platforms)
    }, logical(1))
    items[["features"]] <- items[["features"]][keep]
    items
  }
}


#' Sentinel-2 SCL pixel-level cloud/shadow mask
#' @noRd
.sclMask <- function(scl_classes){
  function(mask_rast) {
    result <- mask_rast == scl_classes[1]
    for (cls in scl_classes[-1]) {
      result <- result | mask_rast == cls
    }
    return(result)
  }
}


#' Fetch STAC data from Sentinel-2, Landsat 8, or Landsat 5
#'
#' @importFrom callr r
#' @importFrom sf st_as_sfc
#' @importFrom rsi get_stac_data landsat_mask_function
#' @noRd
.fetchStac <- function(bbox,
                       start_date,
                       end_date,
                       pixel_x_size,
                       pixel_y_size,
                       asset_names,
                       max_cloud_cover,
                       scl_classes = NULL,
                       composite_function = NULL,
                       source = "sentinel-2"
                       ){
  aoi <- sf::st_as_sfc(bbox)
  band_mapping <- rsi::landsat_band_mapping$planetary_computer_v1
  item_filter  <- .landsatFilter(max_cloud_cover, source)
  mask_band    <- NULL
  mask_fun     <- NULL
  scl_temp     <- scl_classes
  scl_classes  <- NULL
  if (source == "sentinel-2") {
    band_mapping <- rsi::sentinel2_band_mapping[[.S2_PROVIDER]]
    item_filter  <- .cloudFilter(max_cloud_cover)
    if (is.numeric(scl_temp)){
      scl_classes <- scl_temp
      mask_band <- "SCL"
      mask_fun  <- .sclMask(scl_classes)
    }
  }
  stac <- callr::r(
    function( aoi, start_date, end_date, pixel_x_size, pixel_y_size, stac_source, asset_names, collection, query_function,
              sign_function, item_filter, mask_band, mask_function, composite_function, output_filename, limit
      ){
      rsi::get_stac_data( aoi                  = aoi,
                          start_date           = start_date,
                          end_date             = end_date,
                          pixel_x_size         = pixel_x_size,
                          pixel_y_size         = pixel_y_size,
                          stac_source          = stac_source,
                          asset_names          = asset_names,
                          collection           = collection,
                          query_function       = query_function,
                          sign_function        = sign_function,
                          item_filter_function = item_filter,
                          mask_band            = mask_band,
                          mask_function        = mask_function,
                          composite_function   = composite_function,
                          output_filename      = output_filename,
                          limit                = limit )
      },
    args = list( aoi                = aoi,
                 start_date         = start_date,
                 end_date           = end_date,
                 pixel_x_size       = pixel_x_size,
                 pixel_y_size       = pixel_y_size,
                 stac_source        = attr(band_mapping, "stac_source"),
                 asset_names        = asset_names,
                 collection         = attr(band_mapping, "collection_name"),
                 query_function     = attr(band_mapping, "query_function"),
                 sign_function      = attr(band_mapping, "sign_function"),
                 item_filter        = item_filter,
                 mask_band          = mask_band,
                 mask_function      = mask_fun,
                 composite_function = composite_function,
                 output_filename    = tempfile(fileext = ".tif"),
                 limit              = 999 )
    )
  return(stac)
}


#' Extract per-scene mean index values from STAC files
#'
#' @importFrom terra rast global
#' @noRd
.extractStacMean <- function(stac_data,
                             index_function
                             ){
  bn <- basename(stac_data)
  dates <- regmatches(bn, regexpr("\\d{4}-\\d{2}-\\d{2}", bn))

  values <- vapply( stac_data, function(f){
    r <- terra::rast(f)
    terra::global(index_function(r), "mean", na.rm = TRUE)[[1]]
  }, numeric(1), USE.NAMES = FALSE )

  df <- data.frame(date = dates, value = values, stringsAsFactors = FALSE)
  df <- df[!is.nan(df$value), ]
  result <- aggregate(value ~ date, data = df, FUN = mean)
  return(result)
}


#' Validation + fetch wrapper
#' @noRd
.stacWrap <- function(bbox,
                      start_date,
                      end_date,
                      pixel_x_size,
                      pixel_y_size,
                      asset_names,
                      max_cloud_cover,
                      scl_classes,
                      composite_function,
                      source
                      ){
  stopifnot("max_cloud_cover must be between 0 and 100" = max_cloud_cover >= 0 && max_cloud_cover <= 100)

  stac <- .fetchStac(bbox = bbox,
                     start_date = start_date,
                     end_date = end_date,
                     pixel_x_size = pixel_x_size,
                     pixel_y_size = pixel_y_size,
                     asset_names = asset_names,
                     max_cloud_cover = max_cloud_cover,
                     scl_classes = scl_classes,
                     composite_function = composite_function,
                     source = source)
  return(stac)
}


################################################################################
# General functions
################################################################################

#' Collect remote sensing data (Sentinel-2, Landsat pre-2015)
#'
#' This is a wrapper for 'rsi' functions used to extract data from Sentinel-2 remote sensing spectral indices for a given location between two dates.
#'
#' Users simply provide coordinates (longitude & latitude), capture dimensions (x_metres and y_metres), period observed (start and end date), and an index name (must be included in \code{\link{s2_index_list}}).
#' Optional features include setting maximum cloud cover and cloud masking based on Sentinel-2 SCL classes.
#'
#' Note: Some locations may not return data for desired start and end dates. In addition, if \code{start_date} falls outside the Sentinel-2 operational window (before 2015-07-01), Landsat 8 (after 2013-04-01) or Landsat 5
#' (before 2013-04-01) STAC data is collected instead. Landsat data not compatible with red-edge bands (B05, B06, B07) will be skipped with a warning.
#'
#' @param longitude Longitude in decimal degrees (WGS84) of the centre point.
#' @param latitude Latitude in decimal degrees (WGS84) of the centre point.
#' @param x_metres Bounding box x width. Input \code{x_metres} represents the distance (metres) between xmin and xmax in a typical bbox object. It is converted to decimal degrees before it is used to extend the inputted coordinates to a bounding box of the dimensions \code{x_metres} by \code{y_metres}.
#' @param y_metres Bounding box y length. Input \code{y_metres} represents the distance (metres) between ymin and ymax in a typical bbox object. It is converted to decimal degrees and used to extend the inputted coordinates to create a bounding box of the dimensions \code{x_metres} by \code{y_metres}.
#' @param start_date First date to be observed, format "YYYY-MM-DD".
#' @param end_date Last date to be observed, format "YYYY-MM-DD".
#' @param index_name Name of a spectral index from \code{\link{s2_index_list}} (e.g. \code{"NDVI"}, \code{"EVI2"}, \code{"NBR"}). See \code{names(s2_index_list)}.
#' @param max_cloud_cover Maximum allowed cloud cover percentage (0–100). Default 50, e.g. scenes are dropped if cloud cover exceeds 50 per cent.
#' @param scl_classes Vector of integers corresponding to Sentinel-2 SCL classes to be masked (e.g. 9 = cloud high, 8 = cloud medium, 3 = cloud shadow, 10 = thin cirrus). Change from default \code{NULL} to enable. Not available for Landsat observations.
#' @returns A dataframe with columns \code{date} and \code{value}.
#' @examples
#' ndvi_df <- get_rs_data(longitude = 138.715,
#'                       latitude = -34.904,
#'                       x_metres = 30,
#'                       y_metres = 30,
#'                       start_date = "2023-01-01",
#'                       end_date = "2023-01-10",
#'                       index_name = "NDVI")
#'
#' # Pre-Sentinel-2 date: automatically falls back to Landsat
#' ndvi_df <- get_rs_data(longitude = 138.715,
#'                       latitude = -34.904,
#'                       x_metres = 30,
#'                       y_metres = 30,
#'                       start_date = "2000-01-01",
#'                       end_date = "2000-01-31",
#'                       index_name = "NDVI")
#'
#' @export
get_rs_data <- function(longitude,
                       latitude,
                       x_metres,
                       y_metres,
                       start_date,
                       end_date,
                       index_name,
                       max_cloud_cover = 50,
                       scl_classes = NULL
                       ){
  if (xor(is.null(x_metres), is.null(y_metres))) { stop("x_metres and y_metres must both be set, or both left NULL.", call. = FALSE) }

  index_list <- .checkConfigs(start_date, index_name)
  source <- index_list[[1]]
  index_list <- index_list[[2]]

  asset_names <- index_list$assets
  index_function <- index_list$fun
  pixel_size_m <- index_list$res

  pixel_size <- metres_to_degrees(latitude = latitude,
                                   metres = y_metres,
                                   pixel_size = pixel_size_m)

  bbox <- point_to_bbox(longitude,
                        latitude,
                        x_metres,
                        y_metres)

  stac_data <- .stacWrap(bbox = bbox,
                         start_date = start_date,
                         end_date = end_date,
                         pixel_x_size = pixel_size[1],
                         pixel_y_size = pixel_size[2],
                         asset_names = asset_names,
                         max_cloud_cover = max_cloud_cover,
                         scl_classes = scl_classes,
                         composite_function = NULL,
                         source = source)

  if (length(stac_data) == 0) { return(data.frame(date = character(0), value = numeric(0), stringsAsFactors = FALSE)) }

  df <- .extractStacMean(stac_data, index_function)

  return(df)
}


#' Collect remote sensing data (Sentinel-2, Landsat pre-2015)
#'
#' This is a wrapper for 'rsi' functions used to obtain a raster image from Sentinel-2 remote sensing spectral indices for a given location between two dates.
#'
#' Users simply provide coordinates (longitude & latitude), capture dimensions (x_metres and y_metres), period observed (start and end date), and an index name (must be included in \code{\link{s2_index_list}}).
#' Optional features include setting maximum cloud cover and cloud masking based on Sentinel-2 SCL classes.
#'
#' Note: Some locations may not return data for desired start and end dates. In addition, if \code{start_date} falls outside the Sentinel-2 operational window (before 2015-07-01), Landsat 8 (after 2013-04-01) or Landsat 5
#' (before 2013-04-01) STAC data is collected instead. Landsat data not compatible with red-edge bands (B05, B06, B07) will be skipped with a warning.
#'
#' @importFrom terra rast
#' @inheritParams get_rs_data
#' @param composite_function Passed to \code{rsi::get_stac_data}. Default \code{"median"}.
#'   Other options include \code{"mean"}, \code{"max"}, \code{"min"}, \code{"sum"}.
#' @returns A SpatRaster of composite index values, or \code{NULL} if no scenes were returned.
#' @examples
#' rast <- get_rs_raster(longitude = 138.715,
#'                      latitude = -34.904,
#'                      x_metres = 90,
#'                      y_metres = 90,
#'                      start_date = "2023-01-01",
#'                      end_date = "2023-01-10",
#'                      index_name = "NDVI")
#' terra::plot(rast)
#'
#' @export
get_rs_raster <- function(longitude,
                         latitude,
                         x_metres,
                         y_metres,
                         start_date,
                         end_date,
                         index_name,
                         max_cloud_cover = 50,
                         scl_classes = NULL,
                         composite_function = "median"
                         ){
  if (xor(is.null(x_metres), is.null(y_metres))) { stop("x_metres and y_metres must both be set, or both left NULL.", call. = FALSE) }

  index_list <- .checkConfigs(start_date, index_name)
  source <- index_list[[1]]
  index_list <- index_list[[2]]

  asset_names <- index_list$assets
  index_function <- index_list$fun
  pixel_size_m <- index_list$res

  pixel_size <- metres_to_degrees(latitude = latitude,
                                   metres = y_metres,
                                   pixel_size = pixel_size_m)

  bbox <- point_to_bbox(longitude,
                        latitude,
                        x_metres,
                        y_metres)

  stac_data <- .stacWrap(bbox = bbox,
                         start_date = start_date,
                         end_date = end_date,
                         pixel_x_size = pixel_size[1],
                         pixel_y_size = pixel_size[2],
                         asset_names = asset_names,
                         max_cloud_cover = max_cloud_cover,
                         scl_classes = scl_classes,
                         composite_function = composite_function,
                         source = source)

  if (length(stac_data) == 0) return(NULL)

  r <- terra::rast(stac_data)
  index_raster <- index_function(r)

  return(index_raster)
}
