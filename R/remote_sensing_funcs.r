#' Pre-defined Sentinel-2 spectral index definitions
#'
#' A named list of index definitions for use with \code{\link{getS2_data}} and
#' \code{\link{getS2_raster}}. Each entry contains \code{assets} (required band
#' names) and \code{fun} (the index function). The list name itself serves as
#' the index identifier.
#'
#' @format A named list where each element is a list with components:
#' \describe{
#'   \item{assets}{Character vector of Sentinel-2 band names.}
#'   \item{fun}{A function that takes a SpatRaster and returns a single-layer SpatRaster.}
#' }
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#'
#' getS2_data(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001,
#'            asset_names = s2_index$NDVI$assets,
#'            index_function = s2_index$NDVI$fun)
#'
#' names(s2_index)
#'
#' @export
s2_index <- list(
  NDVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04) ),

  NIRv = list(
    assets = c("B08", "B04"),
    fun = function(r) ((r$B08 - r$B04) / (r$B08 + r$B04)) * r$B08 ),

  NDWI = list(
    assets = c("B03", "B08"),
    fun = function(r) (r$B03 - r$B08) / (r$B03 + r$B08) ),

  NDMI = list(
    assets = c("B8A", "B11"),
    fun = function(r) (r$B8A - r$B11) / (r$B8A + r$B11) ),

  EVI = list(
    assets = c("B08", "B04", "B02"),
    fun = function(r) 2.5 * (r$B08 - r$B04) / (r$B08 + 6 * r$B04 - 7.5 * r$B02 + 1) ),

  SAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + 0.5) * 1.5 ),

  MSAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (2 * r$B08 + 1 - sqrt((2 * r$B08 + 1)^2 - 8 * (r$B08 - r$B04))) / 2 ),

  GNDVI = list(
    assets = c("B08", "B03"),
    fun = function(r) (r$B08 - r$B03) / (r$B08 + r$B03) ),

  NDRE = list(
    assets = c("B08", "B05"),
    fun = function(r) (r$B08 - r$B05) / (r$B08 + r$B05) ),

  MTCI = list(
    assets = c("B06", "B05", "B04"),
    fun = function(r) (r$B06 - r$B05) / (r$B05 - r$B04) ),

  CIre = list(
    assets = c("B07", "B05"),
    fun = function(r) (r$B07 / r$B05) - 1 ),

  PSRI = list(
    assets = c("B06", "B04", "B02"),
    fun = function(r) (r$B04 - r$B02) / r$B06 ),

  BSI = list(
    assets = c("B11", "B08", "B04", "B02"),
    fun = function(r) ((r$B11 + r$B04) - (r$B08 + r$B02)) / ((r$B11 + r$B04) + (r$B08 + r$B02)) )
)


################################################################################
# Internal helpers
################################################################################

#' Create a scene-level cloud cover filter
#'
#' @noRd
.cloud_filter <- function(max_cloud_cover){

  function(items, bbox, ...) {
    keep <- vapply(items[["features"]], function(item) {
      cc <- item[["properties"]][["eo:cloud_cover"]]
      !is.null(cc) && cc <= max_cloud_cover
    }, logical(1))
    items[["features"]] <- items[["features"]][keep]
    items
  }
}


#' SCL pixel-level cloud/shadow mask
#'
#' @param scl_classes Integer vector of SCL values to mask
#' @return A function suitable for mask_function in rsi::get_stac_data
#' @noRd
.scl_mask <- function(scl_classes){
  function(mask_rast) {
    result <- mask_rast == scl_classes[1]
    for (cls in scl_classes[-1]) {
      result <- result | mask_rast == cls
    }
    return(result)
  }
}


#' Fetch Sentinel-2 STAC data with cloud filtering
#'
#' @importFrom sf st_as_sfc
#' @importFrom rsi get_stac_data
#' @param bbox bbox object
#' @param start_date,end_date Date strings
#' @param pixel_x_size,pixel_y_size Resolution in decimal degrees
#' @param asset_names Character vector of band names
#' @param max_cloud_cover Scene-level threshold (0–100)
#' @param composite_function NULL for per-scene, or "median" for composites
#' @param provider Name of the STAC provider in rsi::sentinel2_band_mapping
#' @param scl_classes Integer vector of SCL classes to mask, or NULL for no masking
#' @noRd
.fetch_s2 <- function(bbox,
                      start_date,
                      end_date,
                      pixel_x_size,
                      pixel_y_size,
                      asset_names,
                      max_cloud_cover,
                      composite_function = NULL,
                      scl_classes = NULL,
                      provider = "planetary_computer_v1"
){

  band_mapping <- rsi::sentinel2_band_mapping[[provider]]
  if (is.null(band_mapping)) {
    stop(
      "Unknown provider '", provider, "'. Available providers: ",
      paste(names(rsi::sentinel2_band_mapping), collapse = ", "),
      call. = FALSE
    )
  }

  aoi <- sf::st_as_sfc(bbox)

  if (is.null(scl_classes)) {
    mask_band <- NULL
    mask_function <- NULL
  } else {
    mask_band <- "SCL"
    mask_function <- .scl_mask(scl_classes)
  }

  s2_stac <- rsi::get_stac_data(
    aoi = aoi,
    start_date = start_date,
    end_date = end_date,
    pixel_x_size = pixel_x_size,
    pixel_y_size = pixel_y_size,
    stac_source = attr(band_mapping, "stac_source"),
    asset_names = asset_names,
    collection = attr(band_mapping, "collection_name"),
    query_function = attr(band_mapping, "query_function"),
    sign_function = attr(band_mapping, "sign_function"),
    item_filter_function = .cloud_filter(max_cloud_cover),
    mask_band = mask_band,
    mask_function = mask_function,
    composite_function = composite_function,
    output_filename = tempfile(fileext = ".tif"),
    limit = 999
  )
  return(s2_stac)
}


#' Compute a per-scene mean index time series from STAC results
#'
#' @importFrom terra rast global
#' @param stac_data Character vector of file paths from .fetch_s2()
#' @param index_function A function(rast) -> single-layer SpatRaster of index values
#' @noRd
.compute_timeseries <- function(stac_data,
                                index_function
){

  num_features <- length(stac_data)
  df <- data.frame(date = character(num_features),
                   value = numeric(num_features),
                   stringsAsFactors = FALSE)

  for (i in seq_len(num_features)) {
    df[i, "date"] <- regmatches(basename(stac_data[i]),
                                regexpr("\\d{4}-\\d{2}-\\d{2}", basename(stac_data[i])))

    rast <- terra::rast(stac_data[i])
    df[i, "value"] <- terra::global(index_function(rast), "mean", na.rm = TRUE)[[1]]
  }

  df <- df[!is.nan(df$value), ]
  result <- aggregate(value ~ date, data = df, FUN = mean)

  return(result)
}


################################################################################
# General functions
################################################################################

#' Custom vegetation index time series using Sentinel-2 data
#'
#' getS2_data() calculates a user-defined vegetation index from Sentinel-2
#' satellite data, returning a time series of per-scene mean values. This
#' allows computation of any index expressible as a function of Sentinel-2
#' band assets.
#'
#' @param bbox A 'bbox' object (WGS84). Create with sf::st_bbox(..., crs = sf::st_crs(4326)).
#' @param start_date Start of the date range, "YYYY-MM-DD".
#' @param end_date Final date to be observed, "YYYY-MM-DD".
#' @param pixel_x_size Longitude resolution in decimal degrees.
#' @param pixel_y_size Latitude resolution in decimal degrees.
#' @param asset_names Character vector of Sentinel-2 band names required by the index (e.g. \code{c("B08", "B04")}).
#' @param index_function A function that takes a \code{SpatRaster} and returns a single-layer \code{SpatRaster} of index values. Band layers are accessible by name (e.g. \code{r$B08}).
#' @param max_cloud_cover Maximum scene-level cloud cover percentage (0–100). Default 50.
#' @param provider Name of the STAC provider to use. Must be a valid entry in
#'   \code{rsi::sentinel2_band_mapping}. Default \code{"planetary_computer_v1"}.
#' @param scl_classes Integer vector of SCL (Scene Classification Layer) classes
#'   to mask. Common classes: 9 = cloud high probability, 8 = cloud medium probability,
#'   3 = cloud shadow, 10 = thin cirrus, 1 = saturated/defective.
#'
#'   Set to \code{NULL} (default) to disable pixel-level masking.
#' @returns A dataframe with columns for date and value.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906), crs = sf::st_crs(4326))
#'
#' # Using s2_index
#' ndvi_df <- getS2_data(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001,
#'                       asset_names = s2_index$NDVI$assets,
#'                       index_function = s2_index$NDVI$fun)
#'
#' # Custom index with cloud masking
#' evi_df <- getS2_data(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001,
#'                      asset_names = c("B08", "B04", "B02"),
#'                      index_function = function(r) {
#'                        2.5 * (r$B08 - r$B04) / (r$B08 + 6 * r$B04 - 7.5 * r$B02 + 1)
#'                      },
#'                      scl_classes = c(8, 9, 10))
#'
#' @export
getS2_data <- function(bbox,
                       start_date,
                       end_date,
                       pixel_x_size,
                       pixel_y_size,
                       asset_names,
                       index_function,
                       max_cloud_cover = 50,
                       scl_classes = NULL,
                       provider = "planetary_computer_v1"
){

  stopifnot(
    "asset_names must be a character vector" = is.character(asset_names),
    "index_function must be a function" = is.function(index_function),
    "max_cloud_cover must be between 0 and 100" = max_cloud_cover >= 0 && max_cloud_cover <= 100,
    "scl_classes must be NULL or an integer vector" = is.null(scl_classes) || is.numeric(scl_classes)
  )

  stac_data <- .fetch_s2(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                         asset_names = asset_names,
                         max_cloud_cover = max_cloud_cover,
                         provider = provider,
                         scl_classes = scl_classes)

  if (length(stac_data) == 0) {
    warning("No scenes returned. Try increasing max_cloud_cover or widening the date range.",
            call. = FALSE)
    return(data.frame(date = character(0), value = numeric(0), stringsAsFactors = FALSE))
  }

  result <- .compute_timeseries(stac_data, index_function)

  return(result)
}


#' Custom vegetation index raster using Sentinel-2 data
#'
#' getS2_raster() produces a median composite raster of a user-defined
#' vegetation index from Sentinel-2 satellite data.
#'
#' @importFrom terra rast
#' @inheritParams getS2_data
#' @returns A SpatRaster of median composite index values.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#'
#' ndvi_rast <- getS2_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001,
#'                           asset_names = s2_index$NDVI$assets,
#'                           index_function = s2_index$NDVI$fun)
#' terra::plot(ndvi_rast)
#'
#' @export
getS2_raster <- function(bbox,
                         start_date,
                         end_date,
                         pixel_x_size,
                         pixel_y_size,
                         asset_names,
                         index_function,
                         max_cloud_cover = 50,
                         scl_classes = NULL,
                         provider = "planetary_computer_v1"
){

  stopifnot(
    "asset_names must be a character vector" = is.character(asset_names),
    "index_function must be a function" = is.function(index_function),
    "max_cloud_cover must be between 0 and 100" = max_cloud_cover >= 0 && max_cloud_cover <= 100,
    "scl_classes must be NULL or an integer vector" = is.null(scl_classes) || is.numeric(scl_classes)
  )

  comp <- .fetch_s2(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                    asset_names = asset_names,
                    max_cloud_cover = max_cloud_cover,
                    composite_function = "median",
                    provider = provider,
                    scl_classes = scl_classes)

  r <- terra::rast(comp)
  index_raster <- index_function(r)

  return(index_raster)
}


################################################################################
# Convenience functions
################################################################################

#' @title NDVI using Sentinel-2 data
#' @description Normalised Difference Vegetation Index (NDVI) time series.
#' @inheritParams getS2_data
#' @returns A dataframe of mean NDVI per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndvi_df <- getNDVI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getNDVI <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B04"),
                   index_function = function(r) (r$B08 - r$B04) / (r$B08 + r$B04),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title NDVI raster using Sentinel-2 data
#' @description Normalised Difference Vegetation Index (NDVI) median composite raster.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median NDVI.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndvi_rast <- getNDVI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(ndvi_rast)
#' @export
getNDVI_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  NDVI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B08", "B04"),
                              index_function = function(r) (r$B08 - r$B04) / (r$B08 + r$B04),
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(NDVI_raster)
}


#' @title NIRv using Sentinel-2 data
#' @description Near-Infrared Reflectance of Vegetation (NIRv) time series.
#' @inheritParams getS2_data
#' @returns A dataframe of mean NIRv per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' nirv_df <- getNIRv(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getNIRv <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B04"),
                   index_function = function(r) ((r$B08 - r$B04) / (r$B08 + r$B04)) * r$B08,
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title NIRv raster using Sentinel-2 data
#' @description Near-Infrared Reflectance of Vegetation (NIRv) median composite raster.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median NIRv.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' nirv_rast <- getNIRv_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(nirv_rast)
#' @export
getNIRv_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  NIRv_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B08", "B04"),
                              index_function = function(r) ((r$B08 - r$B04) / (r$B08 + r$B04)) * r$B08,
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(NIRv_raster)
}


#' @title NDWI using Sentinel-2 data
#' @description Normalised Difference Water Index (NDWI) time series. Detects and monitors water bodies using green and near-infrared reflectance.
#' @inheritParams getS2_data
#' @returns A dataframe of mean NDWI per date.
#' @references McFeeters, S. K. (1996). The use of the Normalized Difference Water
#'   Index (NDWI) in the delineation of open water features. International Journal
#'   of Remote Sensing, 17(7), 1425-1432.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndwi_df <- getNDWI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getNDWI <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B03", "B08"),
                   index_function = function(r) (r$B03 - r$B08) / (r$B03 + r$B08),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title NDWI raster using Sentinel-2 data
#' @description Normalised Difference Water Index (NDWI) median composite raster. Detects and monitors water bodies using green and near-infrared reflectance.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median NDWI.
#' @references McFeeters, S. K. (1996). The use of the Normalized Difference Water
#'   Index (NDWI) in the delineation of open water features. International Journal
#'   of Remote Sensing, 17(7), 1425-1432.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndwi_rast <- getNDWI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(ndwi_rast)
#' @export
getNDWI_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  NDWI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B03", "B08"),
                              index_function = function(r) (r$B03 - r$B08) / (r$B03 + r$B08),
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(NDWI_raster)
}


#' @title NDMI using Sentinel-2 data
#' @description Normalised Difference Moisture Index (NDMI) time series. Monitors vegetation water content using near-infrared and shortwave infrared reflectance.
#' @inheritParams getS2_data
#' @returns A dataframe of mean NDMI per date.
#' @references Gao, B.-C. (1996). NDWI - A normalized difference water index for
#'   remote sensing of vegetation liquid water from space. Remote Sensing of
#'   Environment, 58(3), 257-266.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndmi_df <- getNDMI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getNDMI <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B8A", "B11"),
                   index_function = function(r) (r$B8A - r$B11) / (r$B8A + r$B11),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title NDMI raster using Sentinel-2 data
#' @description Normalised Difference Moisture Index (NDMI) median composite raster. Monitors vegetation water content using near-infrared and shortwave infrared reflectance.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median NDMI.
#' @references Gao, B.-C. (1996). NDWI - A normalized difference water index for
#'   remote sensing of vegetation liquid water from space. Remote Sensing of
#'   Environment, 58(3), 257-266.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndmi_rast <- getNDMI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(ndmi_rast)
#' @export
getNDMI_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  NDMI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B8A", "B11"),
                              index_function = function(r) (r$B8A - r$B11) / (r$B8A + r$B11),
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(NDMI_raster)
}


#' @title PSRI using Sentinel-2 data
#' @description Plant Senescence Reflectance Index (PSRI) time series.
#' @inheritParams getS2_data
#' @returns A dataframe of mean PSRI per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' psri_df <- getPSRI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getPSRI <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B06", "B04", "B02"),
                   index_function = function(r) (r$B04 - r$B02) / r$B06,
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title PSRI raster using Sentinel-2 data
#' @description Plant Senescence Reflectance Index (PSRI) median composite raster.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median PSRI.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' psri_rast <- getPSRI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(psri_rast)
#' @export
getPSRI_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  PSRI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B06", "B04", "B02"),
                              index_function = function(r) (r$B04 - r$B02) / r$B06,
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(PSRI_raster)
}


#' @title CIre using Sentinel-2 data
#' @description Chlorophyll Index Red-Edge (CIre) time series.
#' @inheritParams getS2_data
#' @returns A dataframe of mean CIre per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' cire_df <- getCIre(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getCIre <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B07", "B05"),
                   index_function = function(r) (r$B07 / r$B05) - 1,
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title CIre raster using Sentinel-2 data
#' @description Chlorophyll Index Red-Edge (CIre) median composite raster.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median CIre.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' cire_rast <- getCIre_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(cire_rast)
#' @export
getCIre_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  CIre_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B07", "B05"),
                              index_function = function(r) (r$B07 / r$B05) - 1,
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(CIre_raster)
}


#' @title MTCI using Sentinel-2 data
#' @description MERIS Terrestrial Chlorophyll Index (MTCI) time series.
#' @inheritParams getS2_data
#' @returns A dataframe of mean MTCI per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' mtci_df <- getMTCI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getMTCI <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B06", "B05", "B04"),
                   index_function = function(r) (r$B06 - r$B05) / (r$B05 - r$B04),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title MTCI raster using Sentinel-2 data
#' @description MERIS Terrestrial Chlorophyll Index (MTCI) median composite raster.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median MTCI.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' mtci_rast <- getMTCI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(mtci_rast)
#' @export
getMTCI_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  MTCI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B06", "B05", "B04"),
                              index_function = function(r) (r$B06 - r$B05) / (r$B05 - r$B04),
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(MTCI_raster)
}


#' @title NDRE using Sentinel-2 data
#' @description Normalised Difference Red-Edge Index (NDRE) time series.
#' @inheritParams getS2_data
#' @returns A dataframe of mean NDRE per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndre_df <- getNDRE(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getNDRE <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B05"),
                   index_function = function(r) (r$B08 - r$B05) / (r$B08 + r$B05),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title NDRE raster using Sentinel-2 data
#' @description Normalised Difference Red-Edge Index (NDRE) median composite raster.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median NDRE.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' ndre_rast <- getNDRE_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(ndre_rast)
#' @export
getNDRE_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1"
){

  NDRE_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B08", "B05"),
                              index_function = function(r) (r$B08 - r$B05) / (r$B08 + r$B05),
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(NDRE_raster)
}


#' @title EVI using Sentinel-2 data
#' @description Enhanced Vegetation Index (EVI) time series. Corrects for atmospheric and canopy background effects, and is more sensitive than NDVI in dense vegetation.
#' @inheritParams getS2_data
#' @returns A dataframe of mean EVI per date.
#' @references Huete, A. et al. (2002). Overview of the radiometric and biophysical
#'   performance of the MODIS vegetation indices. Remote Sensing of Environment, 83(1-2), 195-213.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' evi_df <- getEVI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getEVI <- function(bbox,
                   start_date,
                   end_date,
                   pixel_x_size,
                   pixel_y_size,
                   max_cloud_cover = 50,
                   provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B04", "B02"),
                   index_function = function(r) {2.5 * (r$B08 - r$B04) / (r$B08 + 6 * r$B04 - 7.5 * r$B02 + 1)},
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title EVI raster using Sentinel-2 data
#' @description Enhanced Vegetation Index (EVI) median composite raster. Corrects for atmospheric and canopy background effects, and is more sensitive than NDVI in dense vegetation.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median EVI.
#' @references Huete, A. et al. (2002). Overview of the radiometric and biophysical
#'   performance of the MODIS vegetation indices. Remote Sensing of Environment, 83(1-2), 195-213.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' evi_rast <- getEVI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(evi_rast)
#' @export
getEVI_raster <- function(bbox,
                          start_date,
                          end_date,
                          pixel_x_size,
                          pixel_y_size,
                          max_cloud_cover = 50,
                          provider = "planetary_computer_v1"
){

  EVI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                             asset_names = c("B08", "B04", "B02"),
                             index_function = function(r) {2.5 * (r$B08 - r$B04) / (r$B08 + 6 * r$B04 - 7.5 * r$B02 + 1)},
                             max_cloud_cover = max_cloud_cover,
                             provider = provider)
  return(EVI_raster)
}


#' @title SAVI using Sentinel-2 data
#' @description Soil Adjusted Vegetation Index (SAVI) time series. Minimises soil brightness effects in areas with sparse vegetation cover.
#' @inheritParams getS2_data
#' @param L Soil brightness correction factor (0–1). Default 0.5. Use lower values for high vegetation cover, higher values for low vegetation cover.
#' @returns A dataframe of mean SAVI per date.
#' @references Huete, A. R. (1988). A soil-adjusted vegetation index (SAVI).
#'   Remote Sensing of Environment, 25(3), 295-309.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' savi_df <- getSAVI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getSAVI <- function(bbox,
                    start_date,
                    end_date,
                    pixel_x_size,
                    pixel_y_size,
                    max_cloud_cover = 50,
                    provider = "planetary_computer_v1",
                    L = 0.5
){

  stopifnot("L must be between 0 and 1" = L >= 0 && L <= 1)
  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B04"),
                   index_function = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + L) * (1 + L),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title SAVI raster using Sentinel-2 data
#' @description Soil Adjusted Vegetation Index (SAVI) median composite raster. Minimises soil brightness effects in areas with sparse vegetation cover.
#' @inheritParams getSAVI
#' @returns A SpatRaster of median SAVI.
#' @references Huete, A. R. (1988). A soil-adjusted vegetation index (SAVI).
#'   Remote Sensing of Environment, 25(3), 295-309.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' savi_rast <- getSAVI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(savi_rast)
#' @export
getSAVI_raster <- function(bbox,
                           start_date,
                           end_date,
                           pixel_x_size,
                           pixel_y_size,
                           max_cloud_cover = 50,
                           provider = "planetary_computer_v1",
                           L = 0.5
){

  stopifnot("L must be between 0 and 1" = L >= 0 && L <= 1)
  SAVI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                              asset_names = c("B08", "B04"),
                              index_function = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + L) * (1 + L),
                              max_cloud_cover = max_cloud_cover,
                              provider = provider)
  return(SAVI_raster)
}


#' @title MSAVI using Sentinel-2 data
#' @description Modified Soil Adjusted Vegetation Index (MSAVI) time series. Self-adjusting improvement over SAVI that does not require a user-defined soil brightness parameter.
#' @inheritParams getS2_data
#' @returns A dataframe of mean MSAVI per date.
#' @references Qi, J. et al. (1994). A modified soil adjusted vegetation index.
#'   Remote Sensing of Environment, 48(2), 119-126.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' msavi_df <- getMSAVI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getMSAVI <- function(bbox,
                     start_date,
                     end_date,
                     pixel_x_size,
                     pixel_y_size,
                     max_cloud_cover = 50,
                     provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B04"),
                   index_function = function(r) {(2 * r$B08 + 1 - sqrt((2 * r$B08 + 1)^2 - 8 * (r$B08 - r$B04))) / 2},
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title MSAVI raster using Sentinel-2 data
#' @description Modified Soil Adjusted Vegetation Index (MSAVI) median composite raster. Self-adjusting improvement over SAVI that does not require a user-defined soil brightness parameter.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median MSAVI.
#' @references Qi, J. et al. (1994). A modified soil adjusted vegetation index.
#'   Remote Sensing of Environment, 48(2), 119-126.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' msavi_rast <- getMSAVI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(msavi_rast)
#' @export
getMSAVI_raster <- function(bbox,
                            start_date,
                            end_date,
                            pixel_x_size,
                            pixel_y_size,
                            max_cloud_cover = 50,
                            provider = "planetary_computer_v1"
){

  MSAVI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                               asset_names = c("B08", "B04"),
                               index_function = function(r) {(2 * r$B08 + 1 - sqrt((2 * r$B08 + 1)^2 - 8 * (r$B08 - r$B04))) / 2},
                               max_cloud_cover = max_cloud_cover,
                               provider = provider)
  return(MSAVI_raster)
}


#' @title GNDVI using Sentinel-2 data
#' @description Green Normalised Difference Vegetation Index (GNDVI) time series. More sensitive than NDVI to variation in chlorophyll and nitrogen content.
#' @inheritParams getS2_data
#' @returns A dataframe of mean GNDVI per date.
#' @references Gitelson, A. A. et al. (1996). Use of a green channel in remote sensing
#'   of global vegetation from EOS-MODIS. Remote Sensing of Environment, 58(3), 289-298.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' gndvi_df <- getGNDVI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getGNDVI <- function(bbox,
                     start_date,
                     end_date,
                     pixel_x_size,
                     pixel_y_size,
                     max_cloud_cover = 50,
                     provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B08", "B03"),
                   index_function = function(r) (r$B08 - r$B03) / (r$B08 + r$B03),
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title GNDVI raster using Sentinel-2 data
#' @description Green Normalised Difference Vegetation Index (GNDVI) median composite raster. More sensitive than NDVI to variation in chlorophyll and nitrogen content.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median GNDVI.
#' @references Gitelson, A. A. et al. (1996). Use of a green channel in remote sensing
#'   of global vegetation from EOS-MODIS. Remote Sensing of Environment, 58(3), 289-298.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' gndvi_rast <- getGNDVI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(gndvi_rast)
#' @export
getGNDVI_raster <- function(bbox,
                            start_date,
                            end_date,
                            pixel_x_size,
                            pixel_y_size,
                            max_cloud_cover = 50,
                            provider = "planetary_computer_v1"
){

  GNDVI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                               asset_names = c("B08", "B03"),
                               index_function = function(r) (r$B08 - r$B03) / (r$B08 + r$B03),
                               max_cloud_cover = max_cloud_cover,
                               provider = provider)
  return(GNDVI_raster)
}


#' @title BSI using Sentinel-2 data
#' @description Bare Soil Index (BSI) time series. Highlights exposed soil, useful for monitoring crop emergence, tillage, and erosion.
#' @inheritParams getS2_data
#' @returns A dataframe of mean BSI per date.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' bsi_df <- getBSI(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' @export
getBSI <- function(bbox,
                   start_date,
                   end_date,
                   pixel_x_size,
                   pixel_y_size,
                   max_cloud_cover = 50,
                   provider = "planetary_computer_v1"
){

  df <- getS2_data(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                   asset_names = c("B11", "B08", "B04", "B02"),
                   index_function = function(r) {((r$B11 + r$B04) - (r$B08 + r$B02)) / ((r$B11 + r$B04) + (r$B08 + r$B02))},
                   max_cloud_cover = max_cloud_cover,
                   provider = provider)
  return(df)
}


#' @title BSI raster using Sentinel-2 data
#' @description Bare Soil Index (BSI) median composite raster. Highlights exposed soil, useful for monitoring crop emergence, tillage, and erosion.
#' @inheritParams getS2_data
#' @returns A SpatRaster of median BSI.
#' @examples
#' bbox <- sf::st_bbox(c(xmin = 138.712, xmax = 138.718, ymin = -34.902, ymax = -34.906),
#'                     crs = sf::st_crs(4326))
#' bsi_rast <- getBSI_raster(bbox, "2023-01-01", "2023-01-10", 0.0001, 0.0001)
#' terra::plot(bsi_rast)
#' @export
getBSI_raster <- function(bbox,
                          start_date,
                          end_date,
                          pixel_x_size,
                          pixel_y_size,
                          max_cloud_cover = 50,
                          provider = "planetary_computer_v1"
){

  BSI_raster <- getS2_raster(bbox, start_date, end_date, pixel_x_size, pixel_y_size,
                             asset_names = c("B11", "B08", "B04", "B02"),
                             index_function = function(r) {((r$B11 + r$B04) - (r$B08 + r$B02)) / ((r$B11 + r$B04) + (r$B08 + r$B02))},
                             max_cloud_cover = max_cloud_cover,
                             provider = provider)
  return(BSI_raster)
}
