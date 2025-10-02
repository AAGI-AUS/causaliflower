
################################################################################
# Functions:
# getNDVI()
# getNDVIraster()
################################################################################

#' Provides NDVI at different time points using satellite data (Sentinel-2)
#'
#' getNDVI() calculates mean NDVI in an area at a user-specified resolution for Sentinel-2 satellite image data in the supplied range. This function is essentially a wrapper for the 'rsi' functions rsi_query_api() and get_stac_data() that calculates NDVI from Sentinel-2 observations.
#'
#' @importFrom magrittr %>%
#' @importFrom terra rast
#' @importFrom sf st_as_sfc
#' @importFrom rsi get_stac_data
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr count
#' @param bbox A 'bbox' object indicating the area to be observed. This can be created using sf::st_bbox() with the argument crs = sf::st_crs(4326).
#' @param start_date Start of the date range, "YYYY-MM-DD".
#' @param end_date Final date to be observed, "YYYY-MM-DD".
#' @param pixel_x_size Longitude resolution using decimal degree notation.
#' @param pixel_y_size Latitude resolution using decimal degree notation.
#' @returns A dataframe containing columns for NDVI, and observation dates. NDVI is calculated from the mean B04 and B08 asset values at each time point, taken from pixels in the area defined by the 'bbox' object. Mean NDVI is used for duplicate dates.
#' @examples
#' st_bbox <- sf::st_bbox(c(xmin = 150, xmax = 150.5, ymin = -30, ymax = -30.2), crs = sf::st_crs(4326))
#' start_date <- "2023-01-01"
#' end_date <- "2023-01-10"
#' pixel_y_size <- 0.01
#' pixel_x_size <- 0.01
#'
#' ndvi_df <- getNDVI(st_bbox, start_date, end_date, pixel_x_size, pixel_y_size)
#' ndvi_df
#'
#' @export
getNDVI <- function(bbox, start_date, end_date, pixel_x_size, pixel_y_size){ # Function to produce a composite NDVI raster using B08 and B04 assets

  aoi <- sf::st_as_sfc(bbox)
  stac_source <- "https://planetarycomputer.microsoft.com/api/stac/v1/"
  collection <- "sentinel-2-l2a"
  asset_names <- c("B08", "B04")

  query_function <- attr(rsi::sentinel2_band_mapping[["planetary_computer_v1"]], "query_function")

  s2_stac_data <- rsi::get_stac_data(
    aoi = aoi,
    start_date= start_date,
    end_date = end_date,
    pixel_x_size = pixel_x_size,
    pixel_y_size = pixel_y_size,
    stac_source = stac_source,
    asset_names = asset_names,
    collection = collection,
    query_function = query_function,
    sign_function = NULL,
    item_filter_function = NULL,
    mask_band = NULL,
    mask_function = NULL,
    composite_function = NULL,
    output_filename = tempfile(fileext = ".tif"),
    limit = 999
  )

  num_features <- length(s2_stac_data)

  ndvi_df <- data.frame(matrix(vector(), 0, 2, dimnames=list(c(), c("date", "ndvi"))),
                        stringsAsFactors=F)

  feature <- 1

  for(feature in 1:num_features) {
    date <- sub(".*file(.*?)\\..*", "\\1", s2_stac_data[feature]) # date values are extracted from the file name (could be improved)
    ndvi_df[feature, 1] <- sub(".*_(.*?)\\T.*", "\\1", date)

    rast <- terra::rast(s2_stac_data[feature])
    ndvi_df[feature, 2] <- as.numeric(sum(na.omit(as.data.frame((rast$B08 - rast$B04) / (rast$B08 + rast$B04)))) / dplyr::count(na.omit(as.data.frame(rast))))

    feature <- feature + 1
  }

  ndvi_df <- ndvi_df %>% dplyr::group_by(date) %>% dplyr::summarise(ndvi = mean(as.numeric(ndvi)))

  return(ndvi_df)
}


#' Produces an NDVI image using satellite data (Sentinel-2)
#'
#' getNDVIraster() produces a raster image of the median NDVI values for the supplied date range, area and resolution. This function is essentially a wrapper for the 'rsi' functions and , which collect and transform satellite image data. Sentinel-2 STAC are queried using rsi_query_api(), and a composite function from rsi::get_stac_data() combines asset layers by taking the median value of pixels for each image collected. The median NDVI is calculated for each pixel in the supplied area and resolution.
#'
#' @importFrom magrittr %>%
#' @importFrom lubridate as_date
#' @importFrom terra rast
#' @importFrom sf st_as_sfc
#' @importFrom rsi rsi_query_api
#' @importFrom rsi get_stac_data
#' @param bbox A 'bbox' object indicating the area to be observed. This can be created using sf::st_bbox() with the argument crs = sf::st_crs(4326).
#' @param start_date Start of the date range, "YYYY-MM-DD".
#' @param end_date Final date to be observed, "YYYY-MM-DD".
#' @param pixel_x_size Longitude resolution using decimal degree notation.
#' @param pixel_y_size Latitude resolution using decimal degree notation.
#' @returns A SpatRaster object of NDVI values, calculated from band 4 and band 8 assets. A median composite function is used in rsi::get_stac_date() to combine STAC images, collected from sentinel 2 satellite observation data using rsi::rsi_query_api().
#' @examples
#' st_bbox <- sf::st_bbox(c(xmin = 150, xmax = 150.5, ymin = -30, ymax = -30.2), crs = sf::st_crs(4326))
#' start_date <- "2023-01-01"
#' end_date <- "2023-01-10"
#' pixel_y_size <- 0.01
#' pixel_x_size <- 0.01
#'
#' ndvi_rast <- getNDVIraster(st_bbox, start_date, end_date, pixel_x_size, pixel_y_size)
#' terra::plot(ndvi_rast)
#'
#' @export
getNDVIraster <- function(bbox, start_date, end_date, pixel_x_size, pixel_y_size){ # Function to produce a composite NDVI raster using B08 and B04 assets

  aoi <- sf::st_as_sfc(bbox)
  stac_source <- "https://planetarycomputer.microsoft.com/api/stac/v1/"
  collection <- "sentinel-2-l2a"
  asset_names <- c("B08", "B04")

  query_function <- attr(rsi::sentinel2_band_mapping[["planetary_computer_v1"]], "query_function")
  # Band 8 median composite using rsi::get_stac_data()
  median_composite <- rsi::get_stac_data(
    aoi = aoi,
    start_date= start_date,
    end_date = end_date,
    pixel_x_size = pixel_x_size,
    pixel_y_size = pixel_y_size,
    asset_names = asset_names,
    stac_source = stac_source,
    collection = collection,
    query_function = query_function,
    sign_function = NULL,
    item_filter_function = NULL,
    mask_band = NULL,
    mask_function = NULL,
    composite_function = "median",
    output_filename = tempfile(fileext = ".tif"),
    limit = 999
  )

  rast <- terra::rast(median_composite)
  rast_B08 <- rast$B08
  rast_B04 <- rast$B04

  ndvi_raster <- (rast_B08 - rast_B04) / (rast_B08 + rast_B04)
  return(ndvi_raster)
}
