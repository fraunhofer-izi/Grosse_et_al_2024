# small set of recurrently needed plotting functions

# custom basic plotting function, modelled after Michael Rades plotting functions using STutility
# now using semla instead + it distributes the color scale arguments properly, so you can e.g. just say transform = log1p and be done
# currently does not support colouring more than 1 section or blended images, these are just passed through to MapFeatures
# set quantile limits (with squishing) by using the min and max_cutoff args to MapFeatures. To set absoluet cutoffs, pass the 'limits' argument to the palettte, like c(1000, NA) (potentially with oob=scales::squish)
if (!(exists(x = "default_plot_focus"))){warning("'feature_overlay(se)' with default parameters assumes that you defined 'default_plot_focus' as a valid 'crop_area' argument for the semla::MapFeatures function.")}
feature_overlay = function(se, features, titles = NULL,
                           pt.size = 1.5, pt.stroke = 0.25,    # changed defaults
                           scalename = "Expr", assay = NULL, # default assay
                           gg_guide_colorbar = guide_colorbar(title = scalename, draw.ulim = T, draw.llim = T),
                           colorbar_width = unit(0.6, "lines"), colorbar_height = unit(2., "lines"),
                           slot = NULL, layer = NULL, # slot and layer are the same
                           crop_area = default_plot_focus,   # default plot cropping
                           image_use = "transformed",
                           blend = F, section_number = 1, # repeated defaults from MapFeatures (we wrap plots ourselves)
                           ncol = NULL, arrange_features = c("col", "row"),
                           base_theme = mytheme_histoplot(),
                           palette = "romaO", palette_end = 0.8, palette_begin = 0., # our new default scico palette
                           na.value = "#000000", palette_reversed = T, palette_midpoint = NA,
                           limits = NULL, oob = scales::censor, transform = "identity", # arguments from ggplot2::continuous_scale which can be used to modify the color scale and breaks
                           breaks = ggplot2::waiver(), minor_breaks = ggplot2::waiver(), labels = ggplot2::waiver(),
                           ...) {                            # forward all other args
  if (is.null(titles)) {
    titles = features
  }
  direction = if (palette_reversed) {-1} else {1}
  # slot and layer is the same, Seurat v5 uses only layer and has deprecated
  # slot, but semla still uses slot currently, so translate if necessary
  slot_use = if (is.null(slot)) {
    if (is.null(layer)) {
      "data"
    } else {
      layer
    }
  } else {
    if (is.null(layer)) {
      slot
    } else {
      if (slot != layer) {stop("Slot and layer are equivalent, can't set both to different values.")}
      slot
    }
  }

  # since we do not return the seurat object, we can simply set the assay
  # before plotting without worrying about changing it back
  if (!(is.null(assay))) {
    DefaultAssay(se) <- assay
  }
  # fixed below so that scalebars are properly scaled
  p = MapFeatures(se, features = features, section_number = section_number,
                pt_size = pt.size, pt_stroke = pt.stroke, slot = slot_use,
                image_use = image_use, crop_area = crop_area,
                blend = blend, ncol = ncol, arrange_features = arrange_features,
                return_plot_list = !blend,
                ...)
  if (blend | (section_number != 1)) {
    # we let MapFeatures do it's things normally
    return(p)

  } else {
    # change color scheme and titles
    # p is a nested list of plots here, mapply feature names into this
    # p[[1]] is just the first nesting level corresponding to sections
    # if there are more than 1 sections, we let semla do it's normal things above
    # this function here does not exactly deal with that case
    # the following code mostly replicates semla:::.arrange_plots,
    # but we can ignore blending here and change the titles instead
    p <- Reduce(c, p)
    # each "plot" consists of a 2 patch patchwork:
    # [[1]]: the overlay plot <- this needs to get title and color scale changes
    # [[2]]: the histo picture
    p = mapply(p, titles, FUN = function(plt, nm){
      plt[[1]] = plt[[1]] + ggtitle(nm, subtitle = NULL)
      return(plt)
    })
    if (is.null(ncol)) {
      ncol = ceiling(sqrt(length(p)))
    }
    plts = wrap_plots(p, ncol = ncol)
    plts = plts &
      scale_fill_scico(palette = palette, end = palette_end, begin  = palette_begin,
                       direction = direction, na.value = na.value, midpoint = palette_midpoint,
                       limits = limits, oob = oob, transform = transform,
                       breaks = breaks, minor_breaks = minor_breaks, labels = labels) &
      guides(fill = gg_guide_colorbar) &
      base_theme &
      theme(
        legend.key.width = colorbar_width,
        legend.key.height = colorbar_height
      )
    return(plts)
  }
}

# TODO fix color schemes
feature_overlay_discrete = function(se, features, titles = NULL,
                           pt.size = 1.5, pt.stroke = 0.25,    # changed defaults
                           scalename = "",
                           crop_area = default_plot_focus,   # default plot cropping
                           image_use = "transformed",
                           section_number = 1, # repeated defaults from MapFeatures (we wrap plots ourselves)
                           ncol = NULL, scalebar_position = c(0.8, 0.7),
                           base_theme = mytheme_histoplot(), split_labels = F,
                           # palette = "romaO", palette_end = 0.8, palette_begin = 0., # our new default scico palette
                           # na.value = "#000000", palette_reversed = T, palette_midpoint = NA,
                           # limits = NULL, oob = scales::censor, transform = "identity", # arguments from ggplot2::continuous_scale which can be used to modify the color scale and breaks
                           # breaks = ggplot2::waiver(), minor_breaks = ggplot2::waiver(), labels = ggplot2::waiver(),
                           ...) {                            # forward all other args
  if (is.null(titles)) {
    titles = features
  }
  #direction = if (palette_reversed) {-1} else {1}

  # split labels works properly only with 1 feature, let the normal function handle this
  if (split_labels) {
    p = MapLabels(object = se, column_name = features, section_number = section_number,
                  pt_size = pt.size, pt_stroke = pt.stroke,
                  image_use = image_use, crop_area = crop_area,
                  scalebar_position = scalebar_position, split_labels = T,
                  return_plot_list = T, ...)
  } else {
    if (length(crop_area)==4 && !is.list(crop_area)) {
      # this is a normal crop spec, to get mapply to work, wrap this in a list
      crop_area = list(crop_area)
    }
    # same with scalebar_position
    if (length(scalebar_position)==2 && !is.list(scalebar_position)) {
      scalebar_position = list(scalebar_position)
    }
    p = mapply(FUN = .feature_overlay_discrete_single, SIMPLIFY = F,
               USE.NAMES = F, MoreArgs = list(object = se, ...),
               column_name = features, section_number = section_number,
               pt_size = pt.size, pt_stroke = pt.stroke,
               image_use = image_use, crop_area = crop_area,
               scalebar_position = scalebar_position,
               return_plot_list = T)
    p <- Reduce(c, p)
  }

  # change color scheme and titles
  # p is a nested list of plots here, mapply feature names into this
  # p[[1]] is just the first nesting level corresponding to sections
  # if there are more than 1 sections, we let semla do it's normal things above
  # this function here does not exactly deal with that case
  # the following code mostly replicates semla:::.arrange_plots,
  # but we can ignore blending here and change the titles instead
  # each "plot" consists of a 2 patch patchwork:
  # [[1]]: the overlay plot <- this needs to get title and color scale changes
  # [[2]]: the histo picture
  p = mapply(p, titles, FUN = function(plt, nm){
    plt[[1]] = plt[[1]] + ggtitle(nm, subtitle = NULL)
    return(plt)
  })
  if (is.null(ncol)) {
    ncol = ceiling(sqrt(length(p)))
  }
  plts = wrap_plots(p, ncol = ncol)
  plts = plts &
    # scale_fill_scico(palette = palette, end = palette_end, begin  = palette_begin,
    #                  direction = direction, na.value = na.value, midpoint = palette_midpoint,
    #                  limits = limits, oob = oob, transform = transform,
    #                  breaks = breaks, minor_breaks = minor_breaks, labels = labels) &
    # guides(fill = guide_colorbar(title = scalename, draw.ulim = T, draw.llim = T)) &
    base_theme &
    theme(
      legend.key.width = unit(0.6, "lines"),
      legend.key.height = unit(2., "lines")
    )
  return(plts)
}

# a helper for above
.feature_overlay_discrete_single = function(...)MapLabels(...)

# hist with mad and/or quantile annotations
# adapted from Michaels earlier code
annotatedHist <- function(data, bins = 30,
                             show_quantiles = FALSE, which_quantiles = c(0.05, 0.95), mads = 1:3,
                             show_zero = FALSE) {

  med    <- median(data)
  MAD    <- mad(data, center = med, na.rm = TRUE)

  outs <- lapply(mads, function(n) {
    lower  <- med - n * MAD
    higher <- med + n * MAD
    n_low  <- sum(data < lower)
    n_high <- sum(data > higher)

    list(n = n, lower = lower, higher = higher,
         n_low = n_low, n_high = n_high)
  })

  gg <- ggplot(data.frame(dat=data), aes(x=dat)) +
    geom_histogram(bins = bins, fill = "#BBBBBB") +
    geom_vline(xintercept = med, colour = "#004488") +
    annotate("text", x = med, y = 0, label = "Median", colour = "#004488",
             angle = 90, hjust = -0.1, vjust = -0.5, size = theme_get()$text[["size"]]/2.5) +
    annotate("text", x = med, y = Inf, label = round(med, 2), colour = "#004488",
             angle = 90, hjust = 1.25, vjust = -0.5, size = theme_get()$text[["size"]]/2.5) +
    theme(legend.position = "none")

  cols <- scales::brewer_pal(palette = "Reds",
                             direction = -1)(length(mads) + 1)
  if(length(mads) == 1) {cols = "#CB181D"}

  for (i in seq_along(outs)) {
    out <- outs[[i]]

    if (show_zero) {
      show_low  <- TRUE
      show_high <- TRUE
    } else {
      show_low  <- out$n_low > 0
      show_high <- out$n_high > 0
    }

    if (show_low) {
      gg <- gg +
        geom_vline(xintercept = out$lower, colour = cols[i]) +
        annotate("text", x = out$lower, y = 0,
                 label = paste0("-", out$n, " MADs"), colour = cols[i],
                 angle = 90, hjust = -0.1, vjust = -0.5, size = theme_get()$text[["size"]]/2.5) +
        annotate("text", x = out$lower, y = Inf,
                 label = paste0(round(out$lower, 2),
                                " (", out$n_low, " lower)"),
                 colour = cols[i], angle = 90, hjust = 1.05, vjust = -0.5,
                 size = theme_get()$text[["size"]]/2.5)
    }

    if (show_high) {
      gg <- gg +
        geom_vline(xintercept = out$higher, colour = cols[i]) +
        annotate("text", x = out$higher, y = 0,
                 label = paste0("+", out$n, " MADs"), colour = cols[i],
                 angle = 90, hjust = -0.1, vjust = -0.5, size = theme_get()$text[["size"]]/2.5) +
        annotate("text", x = out$higher, y = Inf,
                 label = paste0(round(out$higher, 2),
                                " (", out$n_high, " higher)"),
                 colour = cols[i], angle = 90, hjust = 1.05, vjust = -0.5
                 , size = theme_get()$text[["size"]]/2.5)
    }
  }
  if (show_quantiles) {
    quantile_cols <- scales::brewer_pal(palette = "Greens",
                                        direction = -1)(length(which_quantiles)+3) # +3 is pretty much darkness offset
    quantiles <- quantile(data, probs = which_quantiles)
    for (q_i in seq_along(which_quantiles)) {
      q_percent <- which_quantiles[q_i]*100
      position <- quantiles[q_i]
      n_higher <- sum(data > position)
      gg <- gg +
        geom_vline(xintercept = position, colour = quantile_cols[q_i]) +
        annotate("text", x = position, y = 0,
                 label = paste0(q_percent, " % quantile"), colour = quantile_cols[q_i],
                 angle = 90, hjust = -0.1, vjust = -0.5, size = theme_get()$text[["size"]]/2.5) +
        annotate("text", x = position, y = Inf,
                 label = paste0(round(position, 2),
                                " (", n_higher, " higher)"),
                 colour = quantile_cols[q_i], angle = 90, hjust = 1.05, vjust = -0.5
                 , size = theme_get()$text[["size"]]/2.5)
    }
  }

  return(gg)
}

# estimate a good-looking point size based on the size of the crop window
estim_point_size <- function(sample_manifest){
  crop = sample_manifest$crop_area
  hsize = crop[3] - crop[1]
  vsize = crop[4] - crop[2]
  return(1. / max(hsize, vsize))
}

library(semla)

fixed_MapFeatures <- function (
    object,
    crop_area = NULL,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    shape = "point",
    spot_side = NULL,
    scale_alpha = FALSE,
    section_number = NULL,
    label_by = NULL,
    ncol = NULL,
    colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
    center_zero = FALSE,
    scale = c("shared", "free"),
    arrange_features = c("col", "row"),
    dims = NULL,
    coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    return_plot_list = FALSE,
    drop_na = FALSE,
    blend = FALSE,
    blend_order = 1:3,
    add_scalebar = FALSE,
    scalebar_gg = NULL,
    scalebar_height = 0.05,
    scalebar_position = c(0.8, 0.8),
    ...
) {

  # Set global variables to NULL
  barcode <- sampleID <- pxl_col_in_fullres <- pxl_row_in_fullres <- image_use <- NULL

  # Check data
  .prep_data_for_plotting(object, colors, label_by, scale, arrange_features, coords_columns)

  # Expand colors if length is 1
  if (length(colors) == 1) {
    colors <- c("lightgray", colors)
  }

  # get features
  features <- setdiff(colnames(object), c("barcode", coords_columns, "sampleID", label_by, "pxl_col_in_fullres", "pxl_row_in_fullres"))

  # Split data by sampleID
  data <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # Check section number and subset data
  if (!is.null(section_number)) {
    if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for",
                                                " {cli::col_br_green('section_number')}, expected 'numeric'"))
    if (!section_number %in% seq_along(data)) abort(glue("'section_number' = {section_number} is out of range. ",
                                                         "Select one of {paste(seq_along(data), collapse = ', ')}"))
    data <- data[section_number]
  }

  # get image dimensions
  dims <- .get_dims(dims)

  # get feature limits
  feature_limits <- .get_feature_limits(data, coords_columns, scale = ifelse(blend, "shared", scale))

  # add blend colors if blend=TRUE
  if (blend) {
    if (!requireNamespace("farver", quietly = TRUE)) {
      abort(glue("Package {cli::col_br_magenta('farver')} is required. Please install it with: \n",
                 "install.packages('farver')"))
    }
    data <- .color_blender(data, features, blend_order, feature_limits, scale_alpha)
    extreme_colors <- farver::encode_colour(diag(ncol = 3, nrow = 3)*255, from = "rgb")
    extreme_colors <- extreme_colors[blend_order[1:length(features)]]
    #features <- features[blend_order[1:length(features)]]
  }

  # Check crop_area
  if (!is.null(crop_area)) {
    if (!is.numeric(crop_area)) abort(glue("Invalid class '{class(crop_area)}' for 'crop_area', expected 'numeric'"))
    if (length(crop_area) != 4) abort(glue("Invalid length for 'crop_area', expected a 'numeric' vector of length 4"))
    if (!all(between(x = crop_area, left = 0, right = 1))) abort("'crop_area' can only take values between 0-1")
    if (!crop_area[1] < crop_area[3]) abort("'left' value needs to be lower that 'right' value")
    if (!crop_area[2] < crop_area[4]) abort("'top' value needs to be lower that 'bottom' value")
  }

  # Edit dims of a crop area if provided
  if (!is.null(crop_area) & shape == "point") {
    c(dims, data) %<-% .crop_dims(dims, crop_area, data, coords_columns)
  } else if (!is.null(crop_area) & shape != "point") {
    c(dims, data) %<-% .crop_array(dims, crop_area, data, coords_columns)
  }

  # Plot features on spatial coordinates for each sample
  sample_plots <- setNames(lapply(names(data), function(nm) {

    # Get data for plotting
    gg <- data[[nm]]

    # Create an appropriate plot title
    if (!is.null(label_by)) {
      cur_label <- unique(gg |> pull(all_of(label_by))) |> as.character()
    } else {
      cur_label <- paste0("section ", nm)
    }

    # Evaluate plotting options and return plots
    switch(shape,
           "point" = if (!blend) { # Default plotting for each feature when blend = FALSE
             feature_plots <- lapply(features, function(ftr) {
               .spatial_feature_plot(
                 gg = gg,
                 nm = nm,
                 ftr = ftr,
                 feature_limits = feature_limits,
                 colors = colors,
                 dims = dims,
                 pt_size = pt_size,
                 pt_alpha = pt_alpha,
                 pt_stroke = pt_stroke,
                 scale_alpha = scale_alpha,
                 coords_columns = coords_columns,
                 cur_label = cur_label,
                 drop_na = drop_na,
                 center_zero = center_zero
               )
             })
             # Add names to feature_plots
             p <- setNames(feature_plots, nm = features)
           } else {
             p <- .spatial_feature_plot(
               gg = gg,
               nm = nm,
               colors = colors,
               dims = dims,
               all_features = features,
               extreme_colors = extreme_colors,
               pt_size = pt_size,
               pt_alpha = pt_alpha,
               pt_stroke = pt_stroke,
               scale_alpha = scale_alpha,
               cur_label = cur_label,
               coords_columns = coords_columns,
               drop_na = drop_na,
               center_zero = center_zero
             )
           },
           "tile" = , # fall through to next scenario
           "raster" = if (!blend) {   # Default plotting for each feature when blend = FALSE
             feature_plots <- lapply(features, function(ftr) {
               .spatial_feature_grid_plot(
                 gg = gg,
                 nm = nm,
                 ftr = ftr,
                 feature_limits = feature_limits,
                 colors = colors,
                 dims = dims,
                 image_use = image_use,
                 shape = shape,
                 spot_side = spot_side[[nm]],
                 pt_alpha = pt_alpha,
                 scale_alpha = scale_alpha,
                 coords_columns = coords_columns,
                 cur_label = cur_label,
                 drop_na = drop_na,
                 center_zero = center_zero
               )
             })
             # Add names to feature_plots
             p <- setNames(feature_plots, nm = features)
           } else {
             p <- .spatial_feature_grid_plot(
               gg = gg,
               nm = nm,
               colors = colors,
               dims = dims,
               image_use = image_use,
               all_features = features,
               extreme_colors = extreme_colors,
               shape = shape,
               spot_side = spot_side[[nm]],
               pt_alpha = pt_alpha,
               scale_alpha = scale_alpha,
               coords_columns = coords_columns,
               cur_label = cur_label,
               drop_na = drop_na,
               center_zero = center_zero
             )
           }
    )

    return(p)
  }), nm = names(data))

  # Add scalebar
  if (add_scalebar & all(coords_columns != c("x", "y"))) {
    if (!requireNamespace("dbscan")) {
      abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with: \n",
                 "install.packages('dbscan')"))
    }
    scalebar_width <- scalebar_gg$labels$scalebar_width
    sample_plots <- lapply(names(sample_plots), function(nm) {
      gg <- data[[nm]]
      nn_dist <- dbscan::kNN(gg |> select(all_of(coords_columns)), k = 1)$dist[, 1] |> min()
      plots <- sample_plots[[nm]]
      d <- dims |> filter(sampleID == nm)
      sf <- scalebar_width/100
      prop_width <- (nn_dist*sf)/(d$full_width - d$x_start)
      scalebar_pos <- scalebar_position %||% c(0.8, 0.8)
      scalebar_pos[1] <- ifelse((1 - prop_width) > scalebar_pos[1], scalebar_pos[1], (1 - prop_width))
      scalebar_pos[2] <- ifelse((1 - scalebar_height) > scalebar_pos[2], scalebar_pos[2], (1 - scalebar_height))
      if (blend) {
        plots <- plots + # changed full -> plot here
          inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "plot",
                        right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
      } else {
        plots <- lapply(names(plots), function(ftr_nm) {
          plots[[ftr_nm]] + # and here
            inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "plot",
                          right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
        }) |> setNames(nm = names(plots))
      }
    }) |> setNames(nm = names(sample_plots))
  } else if (add_scalebar & all(coords_columns == c("x", "y"))) {
    abort(glue("For shape {col_br_green({shape})} without HE, no scalebar can be produced. Try setting {col_br_magenta('image_use')} to {col_br_green({'raw'})} or {col_br_green({'transformed'})}"))
  }

  # Create final patchwork
  if (!return_plot_list) {
    wrapped_plots <- .arrange_plots(sample_plots, features, blend, arrange_features, ncol)
  } else {
    # return list of ggplot objects if return_plot_list = TRUE
    wrapped_plots <- sample_plots
  }

  return(wrapped_plots)
}

environment(fixed_MapFeatures) <- asNamespace("semla")
assignInNamespace("MapFeatures.default", fixed_MapFeatures, asNamespace("semla"))

fixed_MapLabels <- function (
    object,
    spot_side = NULL,
    crop_area = NULL,
    pt_size = 1,
    pt_alpha = 1,
    pt_stroke = 0,
    shape = "point",
    section_number = NULL,
    label_by = NULL,
    split_labels = FALSE,
    ncol = NULL,
    colors = NULL,
    dims = NULL,
    coords_columns = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    return_plot_list = FALSE,
    drop_na = FALSE,
    add_scalebar = FALSE,
    scalebar_gg = NULL,
    scalebar_height = 0.05,
    scalebar_position = c(0.8, 0.8),
    ...
) {

  # Set global variables to NULL
  barcode <- sampleID <- NULL

  # Check data
  .prep_data_for_plotting(object = object, label_by = label_by, coords_columns = coords_columns)

  # get label column
  label <- setdiff(colnames(object), c("barcode", coords_columns, "sampleID", label_by, "pxl_col_in_fullres", "pxl_row_in_fullres"))

  # Check if label column only contains NA values
  if (object |> pull(all_of(label)) |> is.na() |> all()) abort(glue("Selected feature only contains NA values."))

  # Convert label column to factor
  object <- object |>
    mutate(across(all_of(label), ~ factor(.x)))

  # Split data by sampleID
  data <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(nm = unique(object$sampleID))

  # Check section number and subset data
  if (!is.null(section_number)) {
    if (!is.numeric(section_number)) abort(glue("Invalid class '{class(section_number)}' for",
                                                " {cli::col_br_green('section_number')}, expected 'numeric'"))
    if (!section_number %in% seq_along(data)) abort(glue("'section_number' = {section_number} is out of range. ",
                                                         "Select one of {paste(seq_along(data), collapse = ', ')}"))
  }
  if (split_labels) {
    section_number <- section_number %||% {
      warn("No section_number selected. Selecting section 1.")
      1L
    }
    data <- data[section_number]
    c(data, dims, colors) %<-% .split_data_by_label(data, dims[dims$sampleID == section_number, ], label, colors, drop_na)
  }
  if (!is.null(section_number) & !split_labels) {
    data <- data[section_number]
  }

  # Obtain colors
  colors <- colors %||% .gg_color_hue(length(levels(object |> select(all_of(label)) |> pull(all_of(label)))))

  # get image dimensions
  dims <- .get_dims(dims)

  # Check crop_area
  if (!is.null(crop_area)) {
    if (!is.numeric(crop_area)) abort(glue("Invalid class '{class(crop_area)}' for 'crop_area', expected 'numeric'"))
    if (length(crop_area) != 4) abort(glue("Invalid length for 'crop_area', expected a 'numeric' vector of length 4"))
    if (!all(between(x = crop_area, left = 0, right = 1))) abort("'crop_area' can only take values between 0-1")
    if (!crop_area[1] < crop_area[3]) abort("'left' value needs to be lower that 'right' value")
    if (!crop_area[2] < crop_area[4]) abort("'top' value needs to be lower that 'bottom' value")
  }

  # Edit dims of a crop area is provided
  if (!is.null(crop_area) & shape == "point") {
    c(dims, data) %<-% .crop_dims(dims, crop_area, data, coords_columns)
  } else if (!is.null(crop_area) & shape != "point") {
    c(dims, data) %<-% .crop_array(dims, crop_area, data, coords_columns)
  }

  # Plot features on spatial coordinates for each sample
  sample_plots <- setNames(lapply(names(data), function(nm) {

    # Get data for plotting
    gg <- data[[nm]]

    # Set label if available
    if (!is.null(label_by)) {
      cur_label <- unique(gg |> pull(all_of(label_by)) |> as.character())
    } else {
      cur_label <- paste0("section ", nm)
    }

    # Overwrite label if split_labels = TRUE
    if (split_labels) {
      cur_label <- paste0("label: ", nm)
    }

    # Evaluate options and return plots
    p <- switch(shape,
                "point" = .spatial_label_plot(
                  gg = gg,
                  nm = nm,
                  lbl = label,
                  colors = colors,
                  dims = dims,
                  pt_size = pt_size,
                  pt_alpha = pt_alpha,
                  pt_stroke = pt_stroke,
                  coords_columns = coords_columns,
                  cur_label = cur_label,
                  drop_na = drop_na
                ),
                "tile" = , # fall through to next scenario
                "raster" = .spatial_label_grid_plot(
                  gg = gg,
                  nm = nm,
                  lbl = label,
                  shape = shape,
                  spot_side = spot_side[[nm]],
                  colors = colors,
                  dims = dims,
                  pt_alpha = pt_alpha,
                  coords_columns = coords_columns,
                  cur_label = cur_label,
                  drop_na = drop_na
                ))

    return(p)
  }), nm = names(data))

  # Add scalebar
  if (add_scalebar & all(coords_columns != c("x", "y"))) {
    if (!requireNamespace("dbscan")) {
      abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with: \n",
                 "install.packages('dbscan')"))
    }

    scalebar_width <- scalebar_gg$labels$scalebar_width
    sample_plots <- lapply(names(sample_plots), function(nm) {
      gg <- data[[nm]]
      nn_dist <- dbscan::kNN(gg |> select(all_of(coords_columns)), k = 1)$dist[, 1] |> min()
      plots <- sample_plots[[nm]]
      d <- dims |> filter(sampleID == nm)
      sf <- scalebar_width/100
      prop_width <- (nn_dist*sf)/(d$full_width - d$x_start)
      scalebar_pos <- scalebar_position %||% c(0.8, 0.8)
      scalebar_pos[1] <- ifelse((1 - prop_width) > scalebar_pos[1], scalebar_pos[1], (1 - prop_width))
      scalebar_pos[2] <- ifelse((1 - scalebar_height) > scalebar_pos[2], scalebar_pos[2], (1 - scalebar_height))
      plots <- plots +
        inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "plot",
                      right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
    }) |> setNames(nm = names(sample_plots))
  } else if (add_scalebar & all(coords_columns == c("x", "y"))) {
    abort(glue("For shape {col_br_green({shape})} without HE, no scalebar can be produced. Try setting {col_br_magenta('image_use')} to {col_br_green({'raw'})} or {col_br_green({'transformed'})}"))
  }

  if (add_scalebar) {
    if (!requireNamespace("dbscan")) {
      abort(glue("Package {cli::col_br_magenta('dbscan')} is required. Please install it with: \n",
                 "install.packages('dbscan')"))
    }
    scalebar_width <- scalebar_gg$labels$scalebar_width
    sample_plots <- lapply(names(sample_plots), function(nm) {
      gg <- data[[nm]]
      nn_dist <- dbscan::kNN(gg |> select(all_of(coords_columns)), k = 1)$dist[, 1] |> min()
      plots <- sample_plots[[nm]]
      d <- dims |> filter(sampleID == nm)
      sf <- scalebar_width/100
      prop_width <- (nn_dist*sf)/(d$full_width - d$x_start)
      scalebar_pos <- scalebar_position %||% c(0.8, 0.8)
      scalebar_pos[1] <- ifelse((1 - prop_width) > scalebar_pos[1], scalebar_pos[1], (1 - prop_width))
      scalebar_pos[2] <- ifelse((1 - scalebar_height) > scalebar_pos[2], scalebar_pos[2], (1 - scalebar_height))
      plots <- plots +
        inset_element(p = scalebar_gg, left = scalebar_pos[1], bottom = scalebar_pos[2], align_to = "plot",
                      right = scalebar_pos[1] + prop_width, top = scalebar_pos[2] + scalebar_height, on_top = TRUE)
    }) |> setNames(nm = names(sample_plots))
  }

  if (!return_plot_list) {
    # Create final patchwork
    ncol <- ncol %||% ceiling(sqrt(length(data)))
    wrapped_plots <- wrap_plots(sample_plots, ncol = ncol)

    return(wrapped_plots)
  } else {
    return(sample_plots)
  }
}

environment(fixed_MapLabels) <- asNamespace("semla")
assignInNamespace("MapLabels.default", fixed_MapLabels, asNamespace("semla"))

# fix scalebar scaling
fixed_scalebar = function (
    x = 500,
    breaks = 6,
    highlight_breaks = c(1, 6),
    breakheight = 0.5,
    title_position = c("top", "bottom"),
    flip_bar = FALSE,
    text_height = 2,
    ...
) {

  # load ggfittext
  if (!requireNamespace("ggfittext")) {
    abort(glue("Package {cli::col_br_magenta('ggfittext')} is required. Please install it with: \n",
               "install.packages('ggfittext')"))
  }

  # Set global variables to NULL
  #ord <- min_y <- max_y <- NULL

  # Check input
  if (!inherits(x = x, what = c("numeric", "integer"))) {
    abort(glue("Invalid class '{class(x)}' for {col_br_magenta('x')}. Expected a numeric or integer of length 1"))
  }
  if (!between(x = x, left = 100, right = 1.1e4)) {
    abort(glue("Invalid value for {col_br_magenta('x')}. Expected a numeric or integer of length 1"))
  }
  if (!inherits(breaks, what = c("numeric", "integer"))) {
    abort(glue("Invalid class '{class(breaks)}' for {col_br_magenta('breaks')}. Expected a numeric or integer of length 1"))
  }
  if (length(breaks) > 1) {
    abort(glue("Invalid length '{length(breaks)}' for {col_br_magenta('breaks')}. Expected a numeric or integer of length 1"))
  }
  stopifnot(inherits(highlight_breaks, what = c("numeric", "integer")),
            inherits(title_position, what = "character"),
            inherits(flip_bar, what = "logical"),
            length(flip_bar) == 1)

  # Match title position args
  title_position <- match.arg(title_position, choices = c("top", "bottom"))

  # Create title
  title = ifelse(x >= 1e3, paste0(x/1e3, "mm"), paste0(x, "\u00B5m"))

  sb_breaks <- tibble(min_y = rep(-breakheight, breaks),
                      max_y = rep(breakheight, breaks)) |>
    mutate(ord = 1:n()) |>
    mutate(min_y = case_when(ord %in% highlight_breaks ~ min_y*1.5,
                             TRUE ~ min_y),
           max_y = case_when(ord %in% highlight_breaks ~ max_y*1.5,
                             TRUE ~ max_y))

  p <- ggplot() +
    geom_segment(aes(x = 1, xend = max(sb_breaks$ord), y = 0, yend = 0), ...) +
    geom_segment(data = sb_breaks, aes(x = ord, xend = ord, y = min_y, yend = max_y), ...) +
    scale_x_continuous(expand = c(0., 0.)) +
    theme_void()
  if (title_position == "top") {
    p <- p +
      ggfittext::geom_fit_text(aes(xmin = 1, xmax = breaks, ymin = max(sb_breaks$max_y), ymax = max(sb_breaks$max_y) + text_height, label = title),
                               grow = TRUE, angle = ifelse(flip_bar, 90, 0), min.size = 0.1, padding.x = grid::unit(0, "mm"), ...) +
      scale_y_continuous(limits = c(min(sb_breaks$min_y), max(sb_breaks$max_y) + text_height))
  } else {
    p <- p +
      ggfittext::geom_fit_text(aes(xmin = 1, xmax = breaks, ymin = min(sb_breaks$min_y) - text_height, ymax = min(sb_breaks$min_y), label = title),
                               grow = TRUE, angle = ifelse(flip_bar, 90, 0), min.size = 0.1, padding.x = grid::unit(0, "mm"), ...) +
      scale_y_continuous(limits = c(min(sb_breaks$min_y) - text_height, max(sb_breaks$max_y)))
  }

  # Flip bar?
  if (flip_bar) {
    p <- p + coord_flip()
  }

  p$labels$scalebar_width <- x
  return(p)
}
environment(fixed_scalebar) <- asNamespace("semla")
assignInNamespace("scalebar", fixed_scalebar, asNamespace("semla"))
