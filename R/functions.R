#' @useDynLib rawsewage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import Rcpp
NULL

#' Get coefficients from a Waters _HEADER.TXT file
#' @description
#' Automatically extracts vector of calibration coefficients from Waters imaging data
#' @param filename A folder path to a .RAW directory (with or without trailing slash), or direct path to _HEADER.TXT file
#' @return Numeric vector.
#' @export
#' @examples
#' getCoefs("C:/My Folder/image.raw")
getCoefs <- function(filename = NULL){
  if (fs::is_dir(filename) == TRUE) {
    if (substr(filename, nchar(filename), nchar(filename)) != 
        "/") {
      filename <- paste0(filename, "/")
    }
    filename <- paste0(filename, "_HEADER.TXT")
  }
  
  header <- utf8::as_utf8(readChar(filename, nchars = 1e+07))
  
  if(length(strsplit(strsplit(header, "Cal Function 1: ")[[1]][2], ",T0")[[1]]) > 1){
    warning("Calibration function may be of wrong type!")
  }
  
  coefs <- as.numeric(unlist(strsplit((strsplit(strsplit(header, "Cal Function 1: ")[[1]][2], ",T1")[[1]][1]), ",")))
  return(coefs)
}

#' Get memory addresses of scans from a Waters .IDX metadata file
#' @description
#' Automatically extracts vector of memory addresses indicating the start of each scan in an associated .DAT file 
#' @param idx_path A file path to a .IDX file
#' @return Numeric vector.
#' @export
#' @examples
#' getAdds("C:/My Folder/image.raw/_FUNC001.IDX")
getAdds <- function(idx_path) {
  byte_idx <- readBin(file(idx_path, "rb"), what = "raw", n = file.size(idx_path))
  address_bytes <- 23:30
  record_size <- 30
  n_records <- floor(length(byte_idx) / record_size)
  # extract addresses without built in integer overflow nonsense
  adds <- numeric(n_records)
  for (i in seq_along(address_bytes)) {
    raw_component <- byte_idx[seq(address_bytes[i], length(byte_idx), by = record_size)]
    # convert to unsigned integer (0–255) and shift left, should ensure address is always correct
    # regardless of how long the address bytes are
    byte_vals <- as.numeric(bitwAnd(as.integer(raw_component), 0xFF))
    adds <- adds + byte_vals * 2^((i - 1) * 8)
  }
  return(adds)
}

#' Get metadata from a Waters .IDX metadata file
#' @description
#' Automatically extracts vector of memory addresses indicating the start of each scan in an associated .DAT file 
#' @param idx_path A file path to a .IDX file
#' @param calib_coefs A vector of calibration coefficients computed by getCoefs
#' @return Data frame
#' @export
#' @examples
#' getMeta("C:/My Folder/image.raw/_FUNC001.IDX")
getMeta <- function(idx_path, calib_coefs = NULL) {
  con <- file(idx_path, "rb")
  byte <- readBin(con, what = "raw", n = file.size(idx_path))
  # round number in case of file length issues
  num_records <- as.integer(length(byte) / 30)
  
  numpeaks <- integer(num_records)
  bp_mz <- numeric(num_records)
  bp_int <- numeric(num_records)
  tic <- numeric(num_records)
  
  # Extract positions and truncate to length of num_records
  pos_5 <- head(seq(5, length(byte), by = 30), num_records)
  pos_6 <- head(seq(6, length(byte), by = 30), num_records)
  pos_7 <- head(seq(7, length(byte), by = 30), num_records)
  
  pos_9 <- head(seq(9, length(byte), by = 30), num_records)
  pos_10 <- head(seq(10, length(byte), by = 30), num_records)
  pos_11 <- head(seq(11, length(byte), by = 30), num_records)
  pos_12 <- head(seq(12, length(byte), by = 30), num_records)
  
  pos_13 <- head(seq(13, length(byte), by = 30), num_records)
  pos_14 <- head(seq(14, length(byte), by = 30), num_records)
  pos_15 <- head(seq(15, length(byte), by = 30), num_records)
  pos_16 <- head(seq(16, length(byte), by = 30), num_records)
  
  pos_17 <- head(seq(17, length(byte), by = 30), num_records)
  pos_18 <- head(seq(18, length(byte), by = 30), num_records)
  
  pos_19 <- head(seq(19, length(byte), by = 30), num_records)
  pos_20 <- head(seq(20, length(byte), by = 30), num_records)
  pos_21 <- head(seq(21, length(byte), by = 30), num_records)
  pos_22 <- head(seq(22, length(byte), by = 30), num_records)
  
  # Process numpeaks (24-bit integer)
  numpeaks <- as.integer(as.hexmode(paste0(
    sprintf("%02x", as.integer(byte[pos_7])),
    sprintf("%02x", as.integer(byte[pos_6])),
    sprintf("%02x", as.integer(byte[pos_5]))
  )))
  
  # Process bp_int (16-bit integer, little-endian)
  bp_int <- readBin(as.raw(rbind(byte[pos_17], byte[pos_18])), what = "integer", size = 2, endian = "little", n = num_records)
  
  addr <- sapply(byte[pos_19], function(x) {
    bits <- rawToBits(x)
    # extract first 4 bits
    flags <- as.integer(bits[1:4])
    # sum adder value flags
    addr_value <- sum(flags * c(2, 4, 8, 16))
    return(addr_value)
  })
  
  # base peak intensity = integer part *2^(flags)
  bp_int <- bp_int * 2^addr
  
  bits_list <- lapply(byte[pos_19], rawToBits)
  
  # scaler values = (8 - integer) * -1
  scalr <- sapply(bits_list, function(bits) {
    scalr_value <- (8 - bitsToInt(rev(bits[5:8]))) * -1
    return(scalr_value)
  })
  
  # base peak mz raised by 2^scalr
  int_mz <- as.integer(byte[pos_22])
  int_mz <- as.numeric(int_mz)
  
  frac_mz <- as.integer(byte[pos_21])
  frac_mz <- as.numeric((1/256)*frac_mz)  
  
  adj_mz <- strtoi(substr(paste0((rawToBits(as.raw(byte[pos_20])))[1], (rawToBits(as.raw(byte[pos_20])))[2], (rawToBits(as.raw(byte[pos_20])))[3], (rawToBits(as.raw(byte[pos_20])))[4]), 1, 7), base = 2)
  
  adj_mz <- as.numeric((1/65536)*(adj_mz*15))  
  
  bp_mz <- data.frame(int = int_mz, frac = frac_mz, adj = adj_mz)
  bp_mz <- rowSums(bp_mz)
  
  bp_mz <- bp_mz*2^scalr
  
  # recalibrate m/z values if coefficients are provided
  if (is.null(calib_coefs) == FALSE) {
    powers <- outer(sqrt(bp_mz), 0:(length(calib_coefs) - 1), `^`)
    bp_mz <- (powers %*% calib_coefs)^2
  }
  
  # tic (32-bit float, little-endian)
  tic <- readBin(as.raw(rbind(byte[pos_9], byte[pos_10], byte[pos_11], byte[pos_12])), what = "numeric", size = 4, endian = "little", n = num_records)
  
  # scan_time, in minutes (32-bit float, little-endian)
  scan_time <- readBin(as.raw(rbind(byte[pos_13], byte[pos_14], byte[pos_15], byte[pos_16])), what = "numeric", size = 4, endian = "little", n = num_records)
  
  # 
  scan_meta <- data.frame(
    numpeaks = numpeaks,
    bp = bp_mz,
    bp_int = bp_int,
    scan_time = scan_time,
    tic = tic
  )
  on.exit(close(con))
  return(scan_meta)
}

#' Read Waters imaging data directly into an R session
#' @description
#' Believe it or not, reads Waters imaging data directly into an R session
#' @param file_path A file path to a .DAT file
#' @param adds A vector of memory addresses computed by getAdds
#' @param indices Indices of specific scans to be processed
#' @param batch_size The number of scans to process per thread
#' @param calib_coefs A vector of calibration coefficients computed by getCoefs
#' @return A list of lists. Each scan generates a list containing lists for mz and intensity values
#' @export
#' @examples
#' readScans("C:/My Folder/image.raw/_FUNC001.IDX", adds, batch_size = 100, calib_coefs = coefs)
readScans <- function(file_path, adds, indices = c(1:length(adds)), batch_size = 250, calib_coefs = NULL) {

  
  scan_starts <- adds[indices]
  scan_ends <- sapply(indices, function(i) {
    if (i == length(adds)) {
      file.size(file_path)
    } else {
      adds[i + 1]
    }
  })
  
  scan_lengths <- scan_ends - scan_starts
  
  n_blocks <- length(indices)
  batches <- split(seq_len(n_blocks), ceiling(seq_len(n_blocks) / batch_size))
  
  get_uint <- function(x) bitwAnd(as.integer(x), 0xFF)
  
  process_batch <- function(batch_idx, batches, p) {
    batch_ids <- batches[[batch_idx]]
    batch_indices <- indices[batch_ids]
    
    batch_start <- min(adds[batch_indices])
    batch_end <- max(sapply(batch_indices, function(i) {
      if (i == length(adds)) file.size(file_path) else adds[i + 1]
    }))
    block_size <- batch_end - batch_start
    
    con <- file(file_path, "rb")
    on.exit(close(con))
    seek(con, where = batch_start, origin = "start")
    raw_batch <- readBin(con, what = "raw", n = block_size)
    
    result <- lapply(batch_indices, function(i) {
      local_start <- adds[i] - batch_start + 1
      local_end <- if (i == length(adds)) {
        file.size(file_path) - batch_start
      } else {
        adds[i + 1] - batch_start
      }
      byte <- raw_batch[local_start:local_end]
      
      # m/z byte indices
      base_index  <- seq(8, length(byte), by = 8)
      add_index   <- base_index - 1
      frac_index  <- base_index - 2
      micro_index <- base_index - 3
      # intensity byte indices
      int4 <- base_index - 4
      int3 <- base_index - 5
      int2 <- base_index - 6
      int1 <- base_index - 7
      
      byte_ints  <- get_uint(byte[base_index])
      masked     <- bitwAnd(byte_ints, 0x7F)
      num4bit    <- bitwShiftR(masked, 3)
      num3bit    <- bitwAnd(masked, 0x07)
      
      base_mz    <- 8 * (2^(num4bit - 6)) * num3bit
      add_ints   <- get_uint(byte[add_index]) * 2^(num4bit - 11)
      frac_ints  <- get_uint(byte[frac_index]) * (1 / 2^(19 - num4bit))
      micro_ints <- get_uint(byte[micro_index]) * (1 / 2^(27 - num4bit))
      
      mz <- as.numeric(base_mz + add_ints + frac_ints + micro_ints)
      
      if (!is.null(calib_coefs)) {
        powers <- outer(sqrt(mz), 0:(length(calib_coefs) - 1), `^`)
        mz <- as.numeric((powers %*% calib_coefs)^2)
      }
      
      b1 <- get_uint(byte[int1])
      b2 <- get_uint(byte[int2])
      b3 <- get_uint(byte[int3])
      b4 <- get_uint(byte[int4])
      
      frac1 <- sapply(0:7, function(i) as.integer(bitwAnd(b1, bitwShiftL(1, 7 - i)) > 0) * 2^-(14 + i))
      frac1 <- rowSums(matrix(frac1, ncol = 8))
      
      frac2 <- sapply(0:7, function(i) as.integer(bitwAnd(b2, bitwShiftL(1, 7 - i)) > 0) * 2^-(6 + i))
      frac2 <- rowSums(matrix(frac2, ncol = 8))
      
      mult <- sapply(0:7, function(i) as.integer(bitwAnd(b3, bitwShiftL(1, 7 - i)) > 0))[, 1:2]
      mult[,1] <- mult[,1]*4
      mult[,2] <- mult[,2]*2
      mult[which(mult == 0)] <- 1
      mult <- mult[,1] * mult[,2]
      subtract1 <- bitwAnd(b3, 32) > 0
      
      recip_bits <- sapply(4:8, function(bit_pos) {
        as.integer(bitwAnd(b3, bitwShiftL(1, 8 - bit_pos)) > 0) * 2^-(bit_pos - 3)
      })
      frac3 <- rowSums(matrix(recip_bits, ncol = 5))
      
      tempsum <- frac1 + frac2 + frac3
      tempsum[subtract1] <- tempsum[subtract1] - 1
      
      exp_factor <- rep(0, length(b4))
      exp_factor <- exp_factor + ifelse(bitwAnd(b4, 8) > 0, 32, 0)
      exp_factor <- exp_factor + ifelse(bitwAnd(b4, 4) > 0, 16, 0)
      exp_factor <- exp_factor + ifelse(bitwAnd(b4, 2) > 0, 8, 0)
      exp_factor <- exp_factor + ifelse(bitwAnd(b4, 1) > 0, 4, 0)
      
      intensity <- tempsum * mult * 2^exp_factor
      
      return(list(mz = mz, intensity = intensity))
    })
    
    p()
    return(result)
  }
  
  p <- progressr::progressor(steps = length(batches))
  
  result <- future_lapply(seq_along(batches), function(batch_idx) {
    process_batch(batch_idx, batches, p)
  }, future.seed = TRUE)
  
  return(unlist(result, recursive = FALSE))
}

#' Read Waters imaging data directly into an R session, using a wildly experimental C++ function
#' @description
#' Believe it or not, reads Waters imaging data directly into an R session
#' @param file_path A file path to a .DAT file
#' @param adds A vector of memory addresses computed by getAdds
#' @param indices Indices of specific scans to be processed
#' @param batch_size The number of scans to process per thread
#' @param calib_coefs A vector of calibration coefficients computed by getCoefs
#' @param use_integer_intensity Whether to return integer intensities
#' @param use_float_mz Whether to return 32 bit floats. No logic yet to call float so still taking up 64 bits in memory
#' @return A list of lists. Each scan generates a list containing lists for mz and intensity values
#' @export
#' @examples
#' readScans_cpp("C:/My Folder/image.raw/_FUNC001.IDX", adds, batch_size = 100, calib_coefs = coefs)
readScans_cpp <- function(file_path, adds, indices = seq_along(adds),
                          batch_size = 250, calib_coefs = NULL, use_integer_intensity = FALSE, use_float_mz = FALSE) {
  
  n_blocks <- length(indices)
  batches <- split(seq_len(n_blocks), ceiling(seq_len(n_blocks) / batch_size))
  
  p <- progressr::progressor(steps = length(batches))
  
  # Use future_lapply for parallel processing
  result <- future.apply::future_lapply(seq_along(batches), future.seed = NULL, function(batch_idx) {
    
    batch_ids <- batches[[batch_idx]]
    batch_indices <- indices[batch_ids]
    
    start_pos <- min(adds[batch_indices])
    end_pos <- if (max(batch_indices) == length(adds)) {
      file.size(file_path)
    } else {
      max(adds[batch_indices + 1])
    }
    
    con <- file(file_path, "rb")
    seek(con, start_pos)
    raw_data <- readBin(con, "raw", n = end_pos - start_pos)
    close(con)
    
    p()  # progress update
    
    # Call your Rcpp function for this batch
    rcpp_process_batch(raw_data, adds, batch_indices, start_pos, calib_coefs)
  })
  
  # future_lapply returns a list of batch results (each a list of scans)
  # flatten to a single list of scans matching indices
  do.call(c, result)
}


#' Convert raw binary to integer
#' @description
#' Convenience function used within other functions, ignore
#' @param x Raw binary data
#' @return Integer
#' @examples
#' bitsToInt(as.raw(1))
bitsToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

#' Convert Waters .raw file to .mzML via proteowizard
#' @description
#' Convenience function interfacing to proteowizard so the vignette can look relatively clean
#' @param filename Path to a .raw folder of imaging data
#' @param outpath Path for .mzML output to be written to
#' @param msconvert Path to msconvert.exe
#' @return Returns nothing in the session, outputs a .mzML file somewhere on your hard drive, probably
#' @export
#' @examples
#' convertR(filename = "C:/My Folder/image.raw/", outpath = "C:/My Folder/", msconvert = "C:/Users/whoever/AppData/Local/Apps/ProteoWizard 3.0.24094.d2966db 64-bit/msconvert.exe")
convertR <- function(
    filename = NULL,
    outpath = NULL,
    msconvert = msconexe){
  args = c("--mzML",
           "--64",
           "--zlib")
  # strip trailing / if it exists
  if(substr(filename, nchar(filename), nchar(filename)) == "/"){
    filename <- substr(filename, 1, (nchar(filename)-1))
  }
  # if no output folder provided, set to default wd
  if(is.null(outpath) == TRUE){
    args <- c(args, "-o ", shQuote(getwd(), type = "cmd"))
  } else {
    args <- c(args, "-o ", shQuote(outpath, type = "cmd"))
  }
  tempout <- NULL
  system2(command = msconvert, args = c(shQuote(filename, type = "cmd"), args))
}

#' Estimate dimensions of a Waters imaging file from the _extern.inf file
#' @description
#' Convenience function to attempt to calculate dimensions of a Waters imaging file, with a warning for obviously wrong values
#' @param filename Path to a .raw folder of imaging data
#' @return Vector of x,y pixel values
#' @export
#' @examples
#' estDims(filename = "C:/My Folder/image.raw/")
estDims <- function(filename = NULL) {
  if (fs::is_dir(filename) == TRUE) {
    if (substr(filename, nchar(filename), nchar(filename)) != 
        "/") {
      filename <- paste0(filename, "/")
    }
    filename <- paste0(filename, "_extern.inf")
  }
  file_in <- readChar(filename, nchars = 1e+07)
  xlen <- as.numeric(strsplit(strsplit(utf8::as_utf8(file_in), 
                                       "DesiXLength\t\t")[[1]][2], "\r")[[1]][1])
  ylen <- as.numeric(strsplit(strsplit(utf8::as_utf8(file_in), 
                                       "DesiYLength\t\t")[[1]][2], "\r")[[1]][1])
  xstep <- as.numeric(strsplit(strsplit(utf8::as_utf8(file_in), 
                                        "DesiXStep\t\t")[[1]][2], "\r")[[1]][1])
  ystep <- as.numeric(strsplit(strsplit(utf8::as_utf8(file_in), 
                                        "DesiYStep\t\t")[[1]][2], "\r")[[1]][1])
  dim1 <- round(xlen/xstep)
  dim2 <- round(ylen/ystep)
  if (xstep != ystep) {
    warning("Pixel step values do not match, probably returning a wacky dimension")
  }
  return(c(dim1, dim2))
}

#' Estimate potential dimensions of a Waters imaging file as factor pairs of the number of scans
#' @description
#' Convenience function to attempt to calculate dimensions of a Waters imaging file
#' @param filename Path to a .raw folder of imaging data
#' @return List of potential x,y pixel values
#' @export
#' @examples
#' estDims2(filename = "C:/My Folder/image.raw/")
estDims2 <- function(n){
  factors <- function(num) {
    result <- c()
    for (i in 1:sqrt(num)) {
      if (num%%i == 0) {
        result <- c(result, i)
        if (i != num/i) {
          result <- c(result, num/i)
        }
      }
    }
    sort(result)
  }
  factor_list <- factors(n)
  combinations <- list()
  for (i in seq_along(factor_list)) {
    for (j in i:length(factor_list)) {
      if (factor_list[i] * factor_list[j] == n) {
        combinations <- c(combinations, list(c(factor_list[i], 
                                               factor_list[j])))
      }
    }
  }
  combinations <- combinations[-1]
  combinations <- c(combinations, lapply(combinations, rev))
  return(combinations)
}

#' Plot images from an estDims2 list
#' @description
#' Plots data according to potential x,y coordinates in a list generated via estDims2 list
#' @param dims_list List of potential x,y coordinates generated by estDims2
#' @param input_data Matrix of data to be plotted, for example a TIC value from getMeta
#' @return ggplot2 plots
#' @export
#' @examples
#' plotDims(estDims2(filename = "C:/My Folder/image.raw/", input_data = as.matrix(getMeta("C:/My Folder/image.raw/_FUNC001.IDX")$tic)
plotDims <- function(dims_list, input_data){
  plotlist <- list()
  for (i in 1:length(dims_list)) {
    cimgplot <- matrix(data = NA, nrow = dims_list[[i]][1], 
                       ncol = dims_list[[i]][2])
    cimgplot <- sapply(1:nrow(input_data), function(i) {
      sum(input_data[i, ])
    })
    cimgplot <- imager::as.cimg(cimgplot, x = dims_list[[i]][1], 
                                y = dims_list[[i]][2])
    p <- ggplot(data = as.data.frame(cimgplot)) + geom_raster(aes(x = x, 
                                                                  y = rev(y), fill = log(value))) + scale_fill_viridis_c(option = "magma", 
                                                                                                                         name = "Intensity (A.U)") + coord_fixed() + scale_x_continuous(expand = c(0, 
                                                                                                                                                                                                   0)) + scale_y_continuous(expand = c(0, 0), trans = scales::reverse_trans()) + 
      xlab(NULL) + ylab(NULL) + ggtitle(paste0(dims_list[[i]][1], 
                                               "x", dims_list[[i]][2])) + theme(plot.title = element_text(hjust = 0.5))
    plotlist[[i]] <- p
  }
  return(plotlist)
}

#' Estimate potential image dimensions by scan time disjunction
#' @description
#' Searches for sudden jumps in scan acquisition time to calculate raster length
#' @param input_vec A vector of a subset of linearly increasing scan times
#' @param test_vec The full length vector of scan times for the file in question
#' @return Vector of scan time disjunction indixes
#' @export
#' @examples
#' timeTest(scan_meta$scan_time[1:50], scan_meta$scan_time)
#' @details
#' timeTest simply looks at the difference between sequential scan times in a given selection of scans vs the difference over the entire course of the run. Where acquisition stops for the stage to rasterise, this results in a sudden increase in scan time at the rasterisation position, which should be greater than the standard deviation.
timeTest <- function(input_vec, test_vec) {
  # check if the input series has at least 5 elements
  if (length(input_vec) < 5) {
    stop("Input must have at least 5 observations")
  }
  input_mean <- mean(input_vec)
  threshold <- 2*sd(input_vec)
  test_differences <- abs(diff(test_vec))
  disjunct_indices <- which(test_differences > threshold)
  return(disjunct_indices)
}

#' Get dataframe of x/y positions reported by maldichrom.exe
#' @description
#' Calculates x/y positions via maldichrom.exe
#' @param folder_path 
#' @param exe Location of installed maldichrom.exe, defaults to assuming you have one distributed with HDImaging installed in the default location
#' @return Dataframe of x/y positions
#' @export
#' @examples
#' getPos("C:/My Folder/image.raw/")
#' @details
#' This will interface to maldichrom.exe, perform peak picking on your data (looking only for 1 peak per scan), output the data, then find the x and y coordinates associated with each pixel. That is to say, it is extremely slow and wasteful. In rare instances, it appears that HDImaging performs some kind of transformation, this will help understand what's going on.
getPos <- function(folder_path, exe = "C:/HDImaging/lib/maldichrom.exe") {
  text <- as.character(paste0("[MALDICHROM PROCESS]
Type=0
Mass1=0
Mass2=0
Window=0.02
MS Resolution=20000
Number of Peaks=1
Output File=", normalizePath(folder_path), "\\dims.txt",
                              "

[MASSMEASURE MALDICHROM PROCESS]
Do Subtract=0
Do Smooth=0
Do Tof Accurate Mass=0

[BACKSUB MALDICHROM PROCESS]
Percent Below=35
Polynomial Order=5

[SMOOTH MALDICHROM PROCESS]
Smooth Type=0
Smooth Width=3
Number Of Smooths=2

[CENTER MALDICHROM PROCESS]
Top Percent=80.0
Min Peak Width Channels=4"))
  
  fileConn<-file(paste0(folder_path, "test.ini"))
  writeLines(text, fileConn)
  
  
  file1 <- normalizePath(folder_path)
  
  system2(command = exe, args = paste0("-p ", paste0("\"", normalizePath(folder_path), "\\test.ini\"") ," -d ", paste0("\"", normalizePath(file1), "\"")), stderr = "TRUE")
  
  temptsv <- read.csv(paste0(folder_path, "/dims.txt"), sep = "\t", row.names = NULL, skip = 4, header = FALSE)
  xpos <- temptsv[3:nrow(temptsv), 2]
  ypos <- temptsv[3:nrow(temptsv), 3]
  pos_df <- data.frame(xpos = xpos, ypos = ypos)
  return(pos_df)
}

#' Median/mean raster normalisation
#' @description
#' Normalises values by the median/mean of a given raster to mitigate "false" spatial variance over an acquisition
#' @param input_df Your input data, which is best given actually as a matrix not a data frame, with rows as scans and columns as features
#' @param image_dims Vector of x/y dimensions of length 2
#' @param method Whether to normalise by "median" or "mean"
#' @param offset_to_zero Whether to rescale data such that the minimum value is 0 (useful if you need to subsequently log something)
#' @return Matrix of normalised values
#' @export
#' @examples
#' normalize_raster(as.matrix(cbind(log(scan_meta$tic), scan_meta$numpeaks)), imgDims, method = "median")
normRaster <- function(input_df, image_dims, method = c("mean", "median"), offset_to_zero = FALSE) {
  method <- match.arg(method)
  # Preallocate output matrix
  normalized_df <- matrix(data = NA, nrow = nrow(input_df), ncol = ncol(input_df))
  for (i in 1:ncol(input_df)) {
    #cat("\r", "Processing column", i, "of", ncol(input_df))
    # Convert column to matrix representing an image
    image <- matrix(data = input_df[, i], nrow = image_dims[1], ncol = image_dims[2])
    # Preallocate vectors to store means/medians for each raster (column)
    raster_values <- numeric(image_dims[2])
    for (z in 1:image_dims[2]) {
      # Compute the mean or median for all pixels in the current raster
      if (method == "mean") {
        raster_values[z] <- mean(image[, z])
      } else {
        raster_values[z] <- median(image[, z])
      }
    }
    # Normalize each raster by subtracting its mean/median
    for (z in 1:image_dims[2]) {
      image[, z] <- image[, z] - raster_values[z]
    }
    # Optionally offset values so they are all non-negative
    if (offset_to_zero) {
      min_value <- min(image)
      image <- image + abs(min_value)
    }
    # Flatten the normalized image back into a column and store it
    normalized_df[, i] <- as.vector(image)
  }
  return(normalized_df)
}

#' FFT deband an image
#' @description
#' Removes horizontal artifacting (ie persistent rasters in REIMS, false spatial variance due to signal loss over time) via scaling down directional Fourier coefficients
#' @param input_img Your input data, which is best given actually as a matrix not a data frame, with rows as scans and columns as features
#' @param image_dims Vector of x/y dimensions of length 2
#' @return Matrix of normalised values, with dimensions of image. You may need to flatten these to a long dataframe
#' @export
#' @examples
#' as.numeric(deband(as.matrix(scan_meta$tic), numScans))
deband <- function(input_img, image_dims){
  fft_image <- matrix(data = input_img, nrow = image_dims[1], ncol = image_dims[2])
  original_image <- fft_image
  fft_image <- fft(fft_image)
  # Normalize FFT image
  fft_image[1, ] <- fft_image[1, ] / (var(fft_image[1, ]) / mean(var(fft_image[-1, ])))
  fft_image <- Re(fft(fft_image, inverse = TRUE)) / length(fft_image)
  # Offset to reset to zero
  if(abs(min(fft_image)) < 0){
    offset <- abs(min(fft_image))
    fft_image <- fft_image + offset
  }
  return(fft_image)
}

#' Freehand select a region of interest using points
#' @description
#' Allows you to click on your pretty picture and determine the pixel locations of everything within a region. Requires packages X11 and gatepoints. Left click to set points around your ROI, right click and select stop to return ROI.
#' @param input_df Your input data, which is best given actually as a matrix not a data frame, with rows as scans and columns as features
#' @param image_dims Vector of x/y dimensions of length 2
#' @param log Whether to log transform your values for visualisation
#' @param flipy Whether to flip your Y dimension, as different packages will use data in forward or reverse order
#' @param invert Whether to invert pixel values for visualisation
#' @return Vector of scan indices
#' @export
#' @examples
#' roiSelect(as.matrix(scan_meta$tic), numScans)
roiSelect <- function(image_df, image_dims, log = FALSE, flipy = FALSE, invert = FALSE) {
  cimgplot <- matrix(data = NA, nrow = image_dims[1], ncol = image_dims[2])
  tics <- rowSums(image_df)
  for(i in 1:prod(image_dims)){
    cimgplot[i] <- tics[i]
  }
  cimgplot <- as.cimg(cimgplot)
  if(flipy == TRUE){
    cimgplot <- mirror(cimgplot, axis = "y")
  }
  X11()
  plot.new()
  if(invert == TRUE){
    cimgplot <- cimgplot*-1
  }
  # standard
  plot(cimgplot)
  # log tic
  if(log == TRUE){
    plot(log(cimgplot+1))
  }
  selectedPoints <- gatepoints::fhs(as.data.frame(cimgplot), mark = TRUE)
  return(as.numeric(selectedPoints))
}