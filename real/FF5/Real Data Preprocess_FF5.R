################################################################################
# Real Data Preprocess (FF5)
# ------------------------------------------------------------------------------
# Builds the industry-excess-return and Fama-French five-factor panel from the
# two official Kenneth French archives.
#
# Workflow:
#   1) Download the two zip archives into raw/ (skipped if already present).
#   2) Parse the first table of each archive, trim the header names, replace the
#      missing codes -99.99 and -999, and convert percent to decimal returns.
#   3) Inner join, truncate the sample, and form the excess return
#        Y_t = R^{industry}_t - R_{f,t}.
#   4) Write the model data.
#
# The sample ends at 2020-12-31, a calendar boundary for year, half,
# quarter, bimonth, and month blocks.
#
# Outputs:
#   - raw/49_Industry_Portfolios_Daily_CSV.zip
#   - raw/F-F_Research_Data_5_Factors_2x3_daily_CSV.zip
#   - financial_banks_ff5_model_data.csv
################################################################################

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
  out_dir <- dirname(script_path)
} else {
  out_dir <- getwd()
}

raw_dir  <- file.path(out_dir, "raw")
dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)

INDUSTRY   <- "Banks"
DATE_START <- as.Date("2000-01-01")
DATE_END   <- as.Date("2020-12-31")

industry_url <- paste0(
  "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/",
  "49_Industry_Portfolios_Daily_CSV.zip"
)
factor_url <- paste0(
  "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/",
  "F-F_Research_Data_5_Factors_2x3_daily_CSV.zip"
)

industry_zip <- file.path(raw_dir, "49_Industry_Portfolios_Daily_CSV.zip")
factor_zip   <- file.path(raw_dir, "F-F_Research_Data_5_Factors_2x3_daily_CSV.zip")

download_if_needed <- function(url, destination) {
  if (!file.exists(destination) || file.info(destination)$size == 0L) {
    message("Downloading: ", url)
    download.file(url, destination, mode = "wb", quiet = FALSE)
  } else {
    message("Using existing raw file: ", destination)
  }
}

download_if_needed(industry_url, industry_zip)
download_if_needed(factor_url,   factor_zip)

read_zip_lines <- function(zip_path) {
  members <- unzip(zip_path, list = TRUE)$Name
  csv_member <- members[grepl("\\.csv$", members, ignore.case = TRUE)][1L]
  if (is.na(csv_member)) stop("No CSV file found in: ", zip_path)
  temporary_dir <- tempfile("ff_zip_")
  dir.create(temporary_dir)
  on.exit(unlink(temporary_dir, recursive = TRUE, force = TRUE), add = TRUE)
  extracted <- unzip(zip_path, files = csv_member, exdir = temporary_dir)
  readLines(extracted[1L], warn = FALSE, encoding = "latin1")
}

extract_csv_section <- function(lines, header_pattern) {
  header_idx <- grep(header_pattern, lines)[1L]
  if (is.na(header_idx)) stop("Header not found: ", header_pattern)

  after_header <- lines[(header_idx + 1L):length(lines)]
  blank_rel <- which(!nzchar(trimws(after_header)))[1L]
  if (is.na(blank_rel)) {
    end_idx <- length(lines)
  } else {
    end_idx <- header_idx + blank_rel - 1L
  }

  text <- paste(lines[header_idx:end_idx], collapse = "\n")
  out <- read.csv(text = text, check.names = FALSE, strip.white = TRUE)

  names(out) <- trimws(names(out))
  out
}

parse_ff_date <- function(x) {
  as.Date(as.character(x), format = "%Y%m%d")
}

industry_lines <- read_zip_lines(industry_zip)
factor_lines   <- read_zip_lines(factor_zip)

industry <- extract_csv_section(industry_lines, "^,Agric,")
names(industry)[1L] <- "date"
industry$date <- parse_ff_date(industry$date)

if (!INDUSTRY %in% names(industry)) {
  stop(INDUSTRY, " column not found in industry data. Available: ",
       paste(setdiff(names(industry), "date"), collapse = ", "))
}
industry[[INDUSTRY]] <- as.numeric(industry[[INDUSTRY]])

industry[[INDUSTRY]][industry[[INDUSTRY]] <= -99] <- NA_real_
industry[[INDUSTRY]] <- industry[[INDUSTRY]] / 100
industry <- industry[, c("date", INDUSTRY)]

factors <- extract_csv_section(factor_lines, "^,Mkt-RF,SMB,HML,RMW,CMA,RF")
names(factors)[1L] <- "date"
factors$date <- parse_ff_date(factors$date)

factor_names <- c("Mkt-RF", "SMB", "HML", "RMW", "CMA", "RF")
missing_factor_names <- setdiff(factor_names, names(factors))
if (length(missing_factor_names) > 0L) {
  stop("Missing factor columns: ", paste(missing_factor_names, collapse = ", "))
}
for (nm in factor_names) {
  factors[[nm]] <- as.numeric(factors[[nm]])
  factors[[nm]][factors[[nm]] <= -99] <- NA_real_
  factors[[nm]] <- factors[[nm]] / 100
}
factors <- factors[, c("date", factor_names)]

dat <- merge(industry, factors, by = "date", all = FALSE, sort = TRUE)
dat <- dat[dat$date >= DATE_START & dat$date <= DATE_END, ]

n_before <- nrow(dat)
dat <- dat[stats::complete.cases(dat), ]
n_dropped <- n_before - nrow(dat)
if (n_dropped > 0L) {
  message("Dropped ", n_dropped, " rows with missing values inside the sample period.")
}

excess_name <- paste0(INDUSTRY, "_excess")
dat[[excess_name]] <- dat[[INDUSTRY]] - dat$RF

dat <- dat[, c("date", excess_name, "Mkt-RF", "SMB", "HML", "RMW", "CMA", "RF")]
row.names(dat) <- NULL

if (anyDuplicated(dat$date)) stop("Duplicate dates found after merging.")
if (anyNA(dat))              stop("Missing values remain in the final model data.")
if (!all(diff(dat$date) > 0)) stop("Dates are not strictly increasing.")
if (nrow(dat) == 0L)          stop("The final model data are empty.")

output_csv <- file.path(out_dir, "financial_banks_ff5_model_data.csv")
write.csv(dat, output_csv, row.names = FALSE, quote = FALSE, na = "")

message("Created: ", output_csv)
message("Rows: ", nrow(dat), "; dates: ", min(dat$date), " to ", max(dat$date))
