#' Create rodent species table
#'
#' Processes rodent capture data so it can be used for LDA analysis
#'
#' @param periods vector of trapping periods to be included
#' @param plots vector of plot numbers to be included
#' @param species vector of species codes to be included
#'
#' @return Table of species counts per period
#'
#' @examples
#' create_rodent_table()
#'
#' @export 
#'
create_rodent_table <- function(periods = 1:436, 
                                plots = c(2, 4, 8, 11, 12, 14, 17, 22),
                                species = c("BA", "DM", "DO", "DS", "NA", 
                                            "OL", "OT", "PB", "PE", "PF", 
                                            "PH", "PI", "PL", "PM", "PP", 
                                            "RF", "RM", "RO", "SF", "SH", 
                                            "SO")){

  base_url <- "raw.githubusercontent.com/weecology/PortalData/master/Rodents/"
  rodent_file <- "Portal_rodent.csv"
  moon_file <- "moon_dates.csv"
  trapping_file <- "Portal_rodent_trapping.csv"
  rodent_url <- paste("https://", base_url, rodent_file, sep = "")
  moon_url <- paste("https://", base_url, moon_file, sep = "")
  trapping_url <- paste("https://", base_url, trapping_file, sep = "")
  rodent_csv <- RCurl::getURL(rodent_url)
  moon_csv <- RCurl::getURL(moon_url)
  trapping_csv <- RCurl::getURL(trapping_url)
  rodent_data <- read.csv(text = rodent_csv, stringsAsFactors = FALSE,
                   na.strings = c(""), colClasses = c("tag" = "character"))
  moon_data <- read.csv(text = moon_csv, stringsAsFactors = FALSE)
  trapping_data <- read.csv(text = trapping_csv)

  grouping <- rlang::quos(period, species)
  rodents <- rodent_data %>%
             dplyr::filter(period %in% periods, plot %in% plots) %>%
             dplyr::filter(species %in% species) %>%
             dplyr::count(!!!grouping) %>%
             tidyr::complete(!!!grouping, fill = list(n = as.integer(0))) %>%
             tidyr::spread(species, n) %>%
             dplyr::select(period, species)

  trapping <- trapping_data %>%
              dplyr::filter(period %in% periods, plot %in% plots) %>%
              dplyr::group_by(period) %>%
              dplyr::summarise(nplots = sum(sampled), ntraps = sum(effort))

  moons <- moon_data %>%
           dplyr::select(newmoonnumber, newmoondate, period) %>%
           na.omit() %>%
           tibble::as_tibble()

  out <- moons %>%
         dplyr::right_join(trapping, "period") %>%
         dplyr::right_join(rodents, "period") %>%
         dplyr::select(-period)

  return(out)
}