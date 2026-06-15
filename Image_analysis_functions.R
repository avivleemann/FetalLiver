# ==============================================================================
# FILE: Image_analysis_funcs_2026_03_12.R
# ==============================================================================

poly_area <- function(x, y) {
  abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1]))) / 2
}

vecsplit = function(strvec, del, i) {
#	 unlist(lapply(sapply(strvec, strsplit, del), "[[", i))
    apply(sapply(i, function(j) unlist(lapply(sapply(strvec, strsplit, del), "[[", j))),
        1, paste0, collapse = del)
}

date_str <- function(year, month, day, format = "%Y_%m_%d") {
  format(as.Date(sprintf("%04d-%02d-%02d", year, month, day)), format)
}

measure_intensity_q95 <- function(pixels) {
  # Fast 95th percentile calculation
  return(collapse::fquantile(pixels, probs = 0.95, names = FALSE))
}

measure_intensity_q90 <- function(pixels) {
  # Fast 90th percentile calculation
  return(collapse::fquantile(pixels, probs = 0.90, names = FALSE))
}


# ==============================================================================
# LIBRARY MANAGEMENT - Load all required libraries for the project, ensuring they are installed and available
# ==============================================================================
load_project_libraries <- function() {
  pkgs <- c("tidyverse","reshape2", "plyr", "plotrix",
    "magick", "scales", "deldir", "stringr", "dplyr", "grid", 
    "Matrix", "foreach", "doParallel", "ggplot2", "progress", 
    "spatstat", "dbscan", "SPIAT", "patchwork", "logger", "grid",
    "purrr", "sp","collapse","Rfast","ggpubr","ggrepel","tidyr",
    "data.table","grid","EBImage"
  )
  if (!requireNamespace("pacman", quietly = FALSE)) install.packages("pacman")
  pacman::p_load(char = pkgs, character.only = TRUE, install = FALSE)
}



# ==============================================================================
# Amirs function - a collection of small helper functions that I have found useful across multiple projects.
# ==============================================================================
venn = function(A, B) {
    both = union(A, B)
    X = cbind(A = both %in% A, B = both %in% B)
    table(apply(X, 1, function(x) paste0(colnames(X)[x], collapse = ",")))
}

# This function calculates the area of a polygon defined by its vertices (x, y) using the shoelace formula. 
poly_area <- function(x, y) {
  abs(sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1]))) / 2
}

measure_intensity_q95 <- function(pixels) {
  # Fast 95th percentile calculation
  return(collapse::fquantile(pixels, probs = 0.95, names = FALSE))
}

# This function ensures that when you convert a vector to a factor,
# you only include levels that are actually present in the data.
# This prevents issues with empty factor levels and ensures that your analyses and plots are based on the actual data.
full.factor = function(fac) {
     f.fac = factor(fac, levels = names(which(table(fac) > 0)))
     names(f.fac) = names(fac)
     f.fac
}

# This function takes a vector of strings, splits each string by a specified delimiter,
# and then extracts the i-th element from each split.
vecsplit = function(strvec, del, i) {
#	 unlist(lapply(sapply(strvec, strsplit, del), "[[", i))
    apply(sapply(i, function(j) unlist(lapply(sapply(strvec, strsplit, del), "[[", j))),
        1, paste0, collapse = del)
}

# This function creates an empty plot with specified x and y limits, and no axes or labels.
plot.empty = function(xlim = c(0,1), ylim = c(0,1)) {
     plot(1,1, type="n", axes=F, xlab="", ylab="")
}

# This function takes a table (or matrix) and adds a "total" column that sums each row,
# and then adds a "total" row that sums each column (including the new total column).
summarize.table = function(X) {
    X2 = cbind(X, total = rowSums(X))
    rbind(X2, total = colSums(X2))
}

# ==============================================================================
# date formatting - for consistent naming of output files and folders
# ==============================================================================
date_str <- function(year, month, day, format = "%Y_%m_%d") {
  format(as.Date(sprintf("%04d-%02d-%02d", year, month, day)), format)
}


# ==============================================================================
# PURE IMAGE LOADER - Replaces all magick loading routines EBimage package
# ==============================================================================
# one channel image files can sometimes be saved as multi-plane images with one plane containing the actual data and the others being empty.
read_image_to_matrix_ebimage <- function(file_path) {
  if (!file.exists(file_path)) stop("[ERROR] Target file not found: ", file_path)
  # Read natively into normalized float memory [0, 1]
  img <- EBImage::readImage(file_path)
  img_dims <- dim(img)
  if (length(img_dims) > 2) {
    n_planes <- img_dims[3]
    plane_variance <- numeric(n_planes)
    # Calculate the standard deviation of each plane to find real structural data
    for (i in 1:n_planes) {
      plane_variance[i] <- sd(img[, , i])
    }
    # Extract the plane with the highest variance (true signal)
    best_plane <- which.max(plane_variance)
    mat_array <- img[, , best_plane]
  } else {
    mat_array <- as.array(img)
  }
  # Transpose to align with R matrix specs and scale back to [0, 255]
  return(mat_array * 255)
}


# ==============================================================================
# filter by morphological shape and area (CD41 MEGAKARYOCYTES)
# ==============================================================================
## adding mask and features to return for downstream analysis and visualization
filter_cd41_morphological <- function(
  img_matrix,
  bg_brush_size = 255,
  median_size = 5,
  noise_floor = 0.03,
  min_cell_area = 2800,
  max_irregularity = 0.65,
  min_radius_mean = 40,
  return_mask = TRUE,
  otsu_multiplier = 0.50,
  max_intensity_image_pixel = 255,
  intensity_threshold_for_glare = 160,
  epsilon = 1e-5
) {
  zero_var  <-  0
  img_temp <- EBImage::Image(t(img_matrix) / max_intensity_image_pixel, colormode = "Grayscale")
  bg_kern  <- EBImage::makeBrush(size = bg_brush_size, shape = "disc")
  top_hat  <- EBImage::whiteTopHat(img_temp, bg_kern)
  smoothed  <- EBImage::medianFilter(top_hat, size = median_size)

  proc_mat <- t(as.array(smoothed))

  otsu_thresh <- EBImage::otsu(smoothed)
  mask_thresh <- max(noise_floor, otsu_thresh * otsu_multiplier)

  binary_mask <- proc_mat > mask_thresh
  filled_mask  <- EBImage::fillHull(binary_mask)

  labeled_components <- EBImage::bwlabel(filled_mask)
  object_features    <- EBImage::computeFeatures.shape(labeled_components)

  final_clean_mat   <- proc_mat * zero_var
  valid_pixels_mask <- matrix(FALSE, nrow = nrow(proc_mat), ncol = ncol(proc_mat))
  valid_cell_ids    <- integer(0)

  if (!is.null(object_features) && is.matrix(object_features) && nrow(object_features) > 0) {
    areas    <- object_features[, "s.area"]
    rad_sd   <- object_features[, "s.radius.sd"]
    rad_mean <- object_features[, "s.radius.mean"]
    shape_irregularity <- rad_sd / (rad_mean + epsilon)

    valid_cell_ids <- which(
      areas >= min_cell_area &
      shape_irregularity < max_irregularity &
      rad_mean >= min_radius_mean
    )

    if (length(valid_cell_ids) > zero_var) {
      valid_pixels_mask <- labeled_components %in% valid_cell_ids
      final_clean_mat[valid_pixels_mask] <- proc_mat[valid_pixels_mask]

      raw_glare_mask <- (img_matrix > intensity_threshold_for_glare) & valid_pixels_mask
      if (any(raw_glare_mask)) {
        final_clean_mat[raw_glare_mask] <- img_matrix[raw_glare_mask] / max_intensity_image_pixel
      }
    }
  }

  out <- list(
    clean_image = final_clean_mat * max_intensity_image_pixel,
    mk_mask = valid_pixels_mask,
    labels = labeled_components,
    features = object_features,
    kept_ids = valid_cell_ids
  )

  if (return_mask) out else out$clean_image
}

# ==============================================================================
# 2. Ly6G (Granulocytes) Robust Global Filtering - currently not in use
# ==============================================================================
filter_matrix_ly6g_global_robust <- function(img_matrix, k_sigma = 3) {
  if (is.null(img_matrix) || !is.matrix(img_matrix)) {
    stop("Input must be a valid numeric matrix.")
  }
  # Convert to column-major EBImage array scaled to [0, 1]
  img <- EBImage::Image(t(img_matrix) / 255, colormode = "Grayscale")
  # Fixed 5x5 Median Filter to crush sharp antibody aggregates
  despeckled <- EBImage::medianFilter(img, size = 5)
  mat <- as.array(despeckled)
  # Model the whole-image background noise statistics dynamically
  img_mean <- mean(mat)
  img_sd   <- sd(mat)
  # Calculate a dynamic cutoff threshold based on image variance, k_sigma = 3.5 means a pixel must stand out 3.5 standard deviations above baseline
  cutoff <- img_mean + (k_sigma * img_sd)
  # Execute a global intensity truncation mask
  clean_mat <- mat
  clean_mat[clean_mat < cutoff] <- 0
  # Rebuild EBImage array and smooth cell perimeters gently
  img_out <- EBImage::Image(clean_mat, colormode = "Grayscale")
  smoothed <- EBImage::gblur(img_out, sigma = 1)
  # 2. Return back to standard pipeline matrix layout [0, 255]
  return(t(EBImage::imageData(smoothed)) * 255)
}

# ==============================================================================
# 3. F4_80 - Gaussian Smoothing Filter
# ==============================================================================
filter_matrix_gaussian <- function(img_matrix, sigma = 1.5) {
  if (is.null(img_matrix) || !is.matrix(img_matrix)) {
    stop("Error in Gaussian: Input is NOT a matrix.")
  }
  # Convert, blur, and return uniformly matching the [0, 255] pipeline specs
  img <- EBImage::Image(t(img_matrix) / 255, colormode = "Grayscale")
  filtered <- EBImage::gblur(img, sigma = sigma)
  return(t(as.array(filtered)) * 255)
}


# ==============================================================================
# DISPATCHER: General Image Filter Routing
# ==============================================================================
apply_marker_filter <- function(img_matrix, marker_name, config = CONFIG) {
  
  # Standardize marker name: remove special characters (e.g., F4/80 -> F480)
  clean_marker <- stringr::str_remove_all(marker_name, "[^A-Za-z0-9]")
  
  # 1. No Filter (Dapi, Il1b)
  if (stringr::str_detect(clean_marker, stringr::regex("^(DAPI|IL1B)$", ignore_case = TRUE))) {
    cat(sprintf("     -> [%s] Routing to NO FILTER (Returning raw matrix)\n", marker_name))
    return(img_matrix)
    
  # 2. Gaussian Filter (F4/80, Clec4f, S100A8)
  } else if (stringr::str_detect(clean_marker, stringr::regex("^(F480|CLEC4F|S100A8)$", ignore_case = TRUE))) {
    cat(sprintf("     -> [%s] Routing to GAUSSIAN FILTER\n", marker_name))
    # Pass Gaussian-specific parameters here:
    return(filter_matrix_gaussian(img_matrix, sigma = config$SIGMA_NOISE))
    
  # 3. Morphological Gate (CD41)
  } else if (stringr::str_detect(clean_marker, stringr::regex("^CD41$", ignore_case = TRUE))) {
    cat(sprintf("     -> [%s] Routing to MORPHOLOGICAL GATE\n", marker_name))
    # Pass CD41-specific parameters here:
    return(filter_cd41_morphological(img_matrix, return_mask = FALSE))
    
  # Fallback for unrecognized markers
  } else if(stringr::str_detect(clean_marker, stringr::regex("^Ly6G$", ignore_case = TRUE))) {
    cat(sprintf("     -> [%s] Routing to LY6G FILTER\n", marker_name))
    # Pass Ly6G-specific parameters here:
    return(filter_matrix_ly6g_global_robust(img_matrix))

  }else {
    warning(sprintf("     -> [WARNING] Marker '%s' not recognized. Returning raw matrix.", marker_name))
    return(img_matrix)
  }
}

# ==============================================================================
# PROCESSING single image - magick
# ==============================================================================
process_raw_image <- function(file_path) {
  # 1. Read single image
  img <- magick::image_read(file_path)
# Use on.exit to ensure cleanup happens even if the code errors out
  # Use on.exit to ensure cleanup happens (both R pointer AND C-level RAM)
  on.exit({
    rm(img)
    gc(verbose = FALSE) # Forces the MagickCore destructor to run immediately
  }, add = TRUE)

  data_array <- magick::image_data(img)
  n_rows <- dim(data_array)[2]
  
  # 2. Extract and collapse channels
  if (dim(data_array)[1] > 1) {
    c1 <- strtoi(paste0("0x", data_array[1, , ]))
    c2 <- strtoi(paste0("0x", data_array[2, , ]))
    c3 <- strtoi(paste0("0x", data_array[3, , ]))
    mat <- matrix(pmax(pmax(c1, c2), c3), nrow = n_rows)
  } else {
    mat <- matrix(strtoi(paste0("0x", data_array[1, , ])), nrow = n_rows)
  }
  return(mat)
}
# ==============================================================================
# PROCESSING single image with integers - magick
# ==============================================================================
process_raw_image_int <- function(file_path) {
  # 1. Read single image
  img <- magick::image_read(file_path)
  # Use on.exit to ensure cleanup happens (both R pointer AND C-level RAM)
  on.exit({
    rm(img)
    gc(verbose = FALSE) # Forces the MagickCore destructor to run immediately
  }, add = TRUE)

  data_array <- magick::image_data(img)
  n_rows <- dim(data_array)[2]
  # 2. Extract and collapse channels (Direct integer conversion for speed)
  if (dim(data_array)[1] > 1) {
    c1 <- as.integer(data_array[1, , ])
    c2 <- as.integer(data_array[2, , ])
    c3 <- as.integer(data_array[3, , ])
    mat <- matrix(pmax(pmax(c1, c2), c3), nrow = n_rows)
  } else {
    mat <- matrix(as.integer(data_array[1, , ]), nrow = n_rows)
  }
  return(mat)
}


# ==============================================================================
# Plot tesselation 
# ==============================================================================
plot_tessellation <- function(nuclei_sub,
                              data_roi,
                              output,
                              tilesRDS,
                              mfrow_num = 3,
                              mfcol_num = 2,
                              height_out = 2000,
                              width_out = 2000) {
      library(grid) # Explicit load
      col_ord <- names(data_roi)
      tiles <- readRDS(tilesRDS)
      x_range <- seq_len(nrow(data_roi[[1]]))
      y_range <- seq_len(ncol(data_roi[[1]]))
      png(output, height = height_out, width = width_out)
      # Apply the requested 2x2 layout
      par(mfrow = c(mfrow_num, mfcol_num))
      # Plot Channels
      sapply(col_ord, function(x) {
            image(x_range, y_range,
                  data_roi[[x]],
                  main = x, col = colorRampPalette(c("white", "chocolate", "black"))(1000), zlim = c(0, 255)
            )
            grid(col = "black")
            plot(tiles, pch = 19, cex = 0, add = TRUE, border = "blue")
            with(nuclei_sub, text(Location_Center_X, Location_Center_Y, Number_Object_Number))
      })
      # Summary Plot
      plot(1, 1, type = "n", xlim = quantile(x_range, c(0, 1)), ylim = quantile(y_range, c(0, 1)))
      plot(tiles, pch = 19, cex = 0, border = "blue", add = TRUE)
      dev.off()
}

# ==============================================================================
# get tile geometry with only polygon constraints  - for calculation intensity 
# ==============================================================================

get_tile_within <- function(tile, img_rows, img_cols) {
  center_x <- tile$pt[1]
  center_y <- tile$pt[2]
  pol <- cbind(tile$x, tile$y)
  pol <- rbind(pol, pol) # Close the polygon
  # B. Bounding Box (Optimized with Max Radius)
  box <- pmax(c(apply(pol, 2, range)), 1)
  # Restrict box by biological radius
  box[1] <- max(box[1], center_x )
  box[2] <- min(box[2], center_x )
  box[3] <- max(box[3], center_y )
  box[4] <- min(box[4], center_y )
    # C. Define Search Ranges (Clamped to image size)
  x_range <- seq(max(1, floor(box[1])), min(ceiling(box[2]), img_rows))
  y_range <- seq(max(1, floor(box[3])), min(ceiling(box[4]), img_cols))
  if (length(x_range) == 0 || length(y_range) == 0) return(NULL)

  # D. Create Grid & Mask
  x_grid <- rep(x_range, each = length(y_range))
  y_grid <- rep(y_range, length(x_range))
  # Point-in-Polygon Math
  n <- nrow(pol) / 2
  all_pos <- apply(matrix(1:n), 1, function(i) {
    sign((pol[i + 1, 1] - pol[i, 1]) * (y_grid - pol[i, 2]) -
      (pol[i + 1, 2] - pol[i, 2]) * (x_grid - pol[i, 1]))
  })
  
  return(list(x_range = x_range, y_range = y_range, mask = rowSums(all_pos < 0) == 0))
}

# ==============================================================================
# Get cell geometry with both polygon and radius constraints 
# ==============================================================================
get_cell_geometry <- function(tile, max_radius, img_rows, img_cols) {
  center_x <- tile$pt[1]
  center_y <- tile$pt[2]
  # A. Setup Polygon
  pol <- cbind(tile$x, tile$y)
  pol <- rbind(pol, pol) # Close the polygon
  # B. Bounding Box (Optimized with Max Radius)
  box <- pmax(c(apply(pol, 2, range)), 1)
  # Restrict box by biological radius
  box[1] <- max(box[1], center_x - max_radius)
  box[2] <- min(box[2], center_x + max_radius)
  box[3] <- max(box[3], center_y - max_radius)
  box[4] <- min(box[4], center_y + max_radius)

  # C. Define Search Ranges (Clamped to image size)
  x_range <- seq(max(1, floor(box[1])), min(ceiling(box[2]), img_rows))
  y_range <- seq(max(1, floor(box[3])), min(ceiling(box[4]), img_cols))
  if (length(x_range) == 0 || length(y_range) == 0) return(NULL)

  # D. Create Grid & Mask
  x_grid <- rep(x_range, each = length(y_range))
  y_grid <- rep(y_range, length(x_range))
  # Point-in-Polygon Math
  n <- nrow(pol) / 2
  all_pos <- apply(matrix(1:n), 1, function(i) {
    sign((pol[i + 1, 1] - pol[i, 1]) * (y_grid - pol[i, 2]) -
      (pol[i + 1, 2] - pol[i, 2]) * (x_grid - pol[i, 1]))
  })
  inside_poly <- rowSums(all_pos < 0) == 0  # Radius Math
  dist_sq <- (x_grid - center_x)^2 + (y_grid - center_y)^2
  inside_radius <- dist_sq < max_radius^2
  # Final Mask
  mask <- inside_poly & inside_radius
  if (sum(mask) == 0) {
    return(NULL)
  }

  return(list(x_range = x_range, y_range = y_range, mask = mask))
}


# ==============================================================================
# Get cell geometry with both polygon and radius constraints 
# ==============================================================================

get_dual_geometry <- function(tile, max_radius, img_rows, img_cols) {
  center_x <- tile$pt[1]
  center_y <- tile$pt[2]
  pol_x <- tile$x
  pol_y <- tile$y

  # Bounding Box covering the ENTIRE Voronoi polygon
  x_min <- max(1, floor(min(pol_x)))
  x_max <- min(img_rows, ceiling(max(pol_x)))
  y_min <- max(1, floor(min(pol_y)))
  y_max <- min(img_cols, ceiling(max(pol_y)))

  if (x_min > x_max || y_min > y_max) return(NULL)

  x_grid <- rep(x_min:x_max, times = length(y_min:y_max))
  y_grid <- rep(y_min:y_max, each = length(x_min:x_max))

  # 1. IL1b Footprint: The full Voronoi polygon using fast C code
  mask_voronoi <- sp::point.in.polygon(x_grid, y_grid, pol_x, pol_y) > 0

  # 2. Labels Footprint: The tight circle radius (clipped by the Voronoi borders)
  inside_radius <- (x_grid - center_x)^2 + (y_grid - center_y)^2 <= max_radius^2
  mask_radius <- mask_voronoi & inside_radius

  return(list(
    x_range = x_min:x_max, 
    y_range = y_min:y_max, 
    x_grid = x_grid,
    y_grid = y_grid,
    mask_labels = mask_radius,   #  Use for Phenotypes
    mask_quantiles = mask_voronoi #  Use for IL1b
  ))
}

# ==============================================================================
# Classification - Classify single cells based on intensity thresholds and co-expression patterns, with robust handling of edge cases and noise
# ==============================================================================



classify_single_cell_rfast <- function(pixel_matrix, int_thresh,  min_pixels = 20) {
  # Safety check: Count the number of valid pixels (rows)
  if (nrow(pixel_matrix) < min_pixels) return("none") 

  # Background Correction (Fast Median subtraction)
  vec_n <- pmax(pixel_matrix - Rfast::rowMedians(pixel_matrix, parallel = FALSE), 0)
  
  #  Thresholding (Fast 95th percentile)
  q95s <- apply(pixel_matrix, 2, function(v) collapse::fquantile(v, probs = 0.95, names = FALSE))
  high_channels <- names(which(q95s > int_thresh))
  
  # Decision Tree
  if (length(high_channels) == 0) return("none")
  if (length(high_channels) == 1) return(high_channels)
  
  # 4. Complex Logic
  sub_mat <- vec_n[, high_channels, drop = FALSE]
  
  # Safety: If any channel has zero variance, default to the one with highest q95
  if (any(Rfast::colVars(sub_mat, std = TRUE) == 0)) {
    return(high_channels[which.max(q95s[high_channels])])
  }
  
  
  # Find the dominant channel based on normalized intensity
  lead_channel <- high_channels[which.max(Rfast::colsums(sub_mat))]
  
  # Return STRICTLY the lead channel (Biology constraint: no double-labels)
  return(lead_channel)
}

# ==============================================================================
# plot all channels with labels (Updated to use color coding for labels and improved layout)
# ==============================================================================
plot_all_channels <- function(nuclei_sub, data_sub, tiles, cell_types, output_file, image_height = 2000, image_width = 2000, mfrow_num = 2, mfcol_num = 2) {
    # Cells classified as "CD41" get Red, "F4_80" get Green & "none" get gray
    library(grid) 
    type_colors <- c(
        "CD41" = "red", # Red
        "F4_80" = "green", # Green
        "CD41,F4_80" = "navy", # Yellow
        "F4_80,CD41" = "navy", # Yellow (order safety)
        "none" = "gray" # gray for unclassified cells
    )

    png(output_file, height = image_height, width = image_width)
    # Apply the requested 2x2 layout
    par(mfrow = c(mfrow_num, mfcol_num))
    # Loop through EVERY image (Dapi, Il1b, etc.)
    sapply(names(data_sub), function(img_name) {
        # 1. Draw Image
        image(seq_len(nrow(data_sub[[img_name]])), seq_len(ncol(data_sub[[img_name]])),
            data_sub[[img_name]],
            main = img_name, # Title is the channel name
            col = colorRampPalette(c("white", "chocolate", "black"))(1000),
            zlim = c(0, 255), axes = FALSE
        )

        grid(col = "black", lwd = 1)

        # 2. Draw Tiles (ALL channels) - This happens for Dapi, Il1b, CD41, and F4_80
        plot(tiles, pch = 19, cex = 0, add = TRUE, border = "blue", lwd = 1)

        if (grepl("Dapi|Il1b", img_name, ignore.case = TRUE)) {
            # --- ROW 1: PLOT OBJECT NUMBERS ---
            text(
                x = nuclei_sub$location_center_x,
                y = nuclei_sub$location_center_y,
                labels = nuclei_sub$number_object_number,
                col = "black", # Black is very visible on light backgrounds
                cex = 1.3, font = 1
            )
        } else {
            # --- ROW 2: PLOT CLASSIFICATION LABELS ---
            # (Matches anything NOT Dapi/Il1b, e.g. CD41, F4_80)

            txt_cols <- type_colors[cell_types]
            has_label <- !is.na(txt_cols)
            if (any(has_label)) {
                text(
                    x = nuclei_sub$location_center_x[has_label],
                    y = nuclei_sub$location_center_y[has_label],
                    labels = cell_types[has_label],
                    col = txt_cols[has_label],
                    cex = 1, font = 2
                )
            }
        }
    })
    dev.off()
}

 

# ==============================================================================
# Get Fucsions for Gating FACS Plots
# ==============================================================================

gate_facs <- function(facs, xlab, ylab, gate = TRUE, colgate = NULL, roof = 1e4, rect = FALSE, log = "", polygon = NULL, n = 1000) {
  if (is.null(polygon)) {
    disp_facs <- facs[gate, ]
    if (roof < dim(disp_facs)[1]) {
      ind <- sample(dim(disp_facs)[1], roof)
    } else {
      ind <- 1:dim(disp_facs)[1]
    }

    x <- disp_facs[ind, xlab]
    y <- disp_facs[ind, ylab]

    if (!is.null(colgate)) {
      colgate <- colgate[gate, ]
      colgate <- colgate[ind, ]
      colors <- 
        (apply(matrix(1:dim(colgate)[1]), MARGIN = 1, FUN = function(x) {
        paste(colgate[x, ] * 1, collapse = "")
      }))
      # smoothScatter(x, y)
      # points(x, y, col = colors, pch=21, cex = 0.4)
      plot(x, y, pch = 21, cex = 0.4, log = log, col = colors + 1)
    } else {
      plot(x, y, pch = 21, cex = 0.4, log = log)
    }

    coords <- locator(n, type = "l") # add lines
    C <- unlist(coords)
    if (!rect) {
      n <- length(coords$x)
      c <- rep(0, n * 4)
      dim(c) <- c(n * 2, 2)
      c[1:n, ] <- C
      c[(n + 1):(n * 2), ] <- C
    } else {
      n <- 4
      c <- matrix(0, n * 2, 2)
      c[1, ] <- c(C[1], C[3])
      c[2, ] <- c(C[1], C[4])
      c[3, ] <- c(C[2], C[4])
      c[4, ] <- c(C[2], C[3])
      c[5:8, ] <- c[1:4, ]
    }
    lines(c, col = "red")
  } else {
    c <- polygon
    n <- nrow(c) / 2
  }

  x <- facs[, xlab]
  y <- facs[, ylab]

  all_pos <- apply(matrix(1:n),
    MARGIN = 1, FUN =
      function(i) {
        sign((c[i + 1, 1] - c[i, 1]) * (y - c[i, 2]) - (c[i + 1, 2] - c[i, 2]) * (x - c[i, 1]))
      }
  )
  w <- unlist(all_pos)
  within <- apply(w < 0, 1, prod) == 1
  return(list(gate = within, polygon = c))
}


map_relative_clusters <- function(nuclei_df, cluster_df, current_series, tiles, max_area = 4000) {
  
  # Initialize a standalone mapping dataframe using the primary key
  mapping_df <- data.frame(
    number_object_number = nuclei_df$number_object_number,
    mapped_cluster_id = NA,
    is_target = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Calculate areas and find the valid tissue polygons
  areas <- sapply(tiles, function(t) poly_area(t$x, t$y))
  small_idx <- which(areas < max_area) 
  
  # Filter clusters for this specific image
  img_clusters <- cluster_df[cluster_df$image_batch == current_series, ]
  
  if(nrow(img_clusters) > 0) {
    for(i in seq_len(nrow(img_clusters))) {
      cid <- img_clusters$cluster_id[i] # Fixed 'cluster_idq' typo here
      
      # Extract relative cells and clean the string
      rel_str <- img_clusters$cluster_cells_relative[i]
      rel_str_clean <- gsub("[^0-9,]", "", rel_str) 
      
      # Map it if valid
      if(nchar(rel_str_clean) > 0) {
        rel_idx <- as.numeric(unlist(strsplit(rel_str_clean, ",")))
        full_idx <- small_idx[rel_idx] # Translate to full image index
        
        # Assign safely to the isolated mapping dataframe
        mapping_df$mapped_cluster_id[full_idx] <- cid
        mapping_df$is_target[full_idx] <- TRUE
      }
    }
  }
  
  return(mapping_df)
}





plot_clec4f_channels <- function(nuclei_sub, data_sub, tiles, cell_types, output_file, image_height = 2000, image_width = 2000, mfrow_num = 2, mfcol_num = 2) {
    # Cells classified as "CD41" get Red, "F4_80" get Green & "none" get gray
    library(grid) 
    type_colors <- c(
        "cd41" = "red", # Red
        "clec4f" = "green", # Green
        "cd41,clec4f" = "navy", # Yellow
        "clec4f,cd41" = "navy", # Yellow (order safety)
        "none" = "gray" # gray for unclassified cells
    )

    png(output_file, height = image_height, width = image_width)
    # Apply the requested 2x2 layout
    par(mfrow = c(mfrow_num, mfcol_num))
    # Loop through EVERY image (Dapi, Il1b, etc.)
    sapply(names(data_sub), function(img_name) {
        # 1. Draw Image

        img_low <- tolower(img_name)
        channel_col <- case_when(
            grepl("dapi", img_low)   ~ "blue",
            grepl("clec4f", img_low) ~ "#00FFFF", # Cyan
            grepl("cd41", img_low)   ~ "#FF00FF", # Magenta
            TRUE                     ~ "white"    # Default for IL1b/others
        )
        
        curr_palette <- colorRampPalette(c("black", channel_col))(256)
        image(seq_len(nrow(data_sub[[img_name]])), seq_len(ncol(data_sub[[img_name]])),
            data_sub[[img_name]],
            main = img_name, # Title is the channel name
            col = curr_palette,
            zlim = c(0, 255), axes = FALSE,
        )

        grid(col = "gray20", lty = 3)

        # 2. Draw Tiles (ALL channels) - This happens for Dapi, Il1b, CD41, and clec4f
        plot(tiles, pch = 20, cex = 0, add = TRUE, border = "blue", lwd = 1.5)

        if (grepl("dapi", img_name, ignore.case = TRUE)) {
            # --- ROW 1: PLOT OBJECT NUMBERS ---
            text(
                x = nuclei_sub$location_center_x,
                y = nuclei_sub$location_center_y,
                labels = nuclei_sub$number_object_number,
                col = "yellow", # Black is very visible on light backgrounds
                cex = 0.7, font = 2
            )
        } else {
            # --- ROW 2: PLOT CLASSIFICATION LABELS ---
            # (Matches anything NOT Dapi/Il1b, e.g. CD41, clec4f)

            txt_cols <- type_colors[cell_types]
            has_label <- !is.na(txt_cols)
            if (any(has_label)) {
                text(
                    x = nuclei_sub$location_center_x[has_label],
                    y = nuclei_sub$location_center_y[has_label],
                    labels = cell_types[has_label],
                    col = txt_cols[has_label],
                    cex = 0.8, font = 2
                )
            }
        }
    })
    dev.off()
}


