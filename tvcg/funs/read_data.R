#' Read and preprocess German dataset
#'
#' @param file 
#'
#' @return matrix
read_german <- function(file = "data/socio_economics_germany_full.rds") {
  german_orig <- readRDS(file)
  german <- german_orig[, 14:45]
  #german_factors <- german_orig[, 3:5]
  scale(as.matrix(german))
}


#' Read and preprocess air quality dataset
#'
#' @return data.frame
read_aq <- function(file = "data/AirQualityUCI.csv") {
  AQ_raw <- read.csv(file, header=TRUE, sep=";", dec=",", stringsAsFactors = F)
  all_na_rows = is.na(AQ_raw[,"CO.GT."])
  all_na_cols = sapply(AQ_raw, function(col) all(is.na(col)))
  AQ = AQ_raw[!all_na_rows, !all_na_cols]
  for(i in 3:15) 
    AQ[AQ[,i]==-200, i] <- NA
  AQ$timestamp = as.POSIXct(paste(AQ$Date, AQ$Time), format = "%d/%m/%Y %H.%M.%S")
  AQ = AQ[!is.na(AQ$timestamp),] # Two dates are NA for some reason.
  j_numeric = sapply(AQ, is.numeric)
  AQ[, j_numeric] = scale(AQ[, j_numeric])
  j = c("timestamp", names(AQ)[!names(AQ) %in% c("Date", "Time", "timestamp")])
  AQ[ , j]
}

#' Read SMEAR datasets
#'
#' @param data_file RDS file of particle concentration data over time.
#' The first column is a timestamp, and the rest are are particle concentrations
#' (d100e1 means Particle concentration dNdlogDp Dp=1.00nm)
#' @param event_data_file RDS file of event classification data
#'
#' @return list of two data frames
read_smear <- function(data_file = "data/smear_2010.rds", event_data_file = "data/smear_2010_events.rds") {
  smear_data <- readRDS(data_file) 
  colnames(smear_data) = gsub("HYY_DMPS.", "", colnames(smear_data))
  smear_events <- readRDS(event_data_file) # This dataset is not public.
  
  list(data = smear_data, 
       events = smear_events, 
       i_non_event = smear_events %in% "non.event", 
       i_event = smear_events %in% c("event.Ia", "event.Ib"))
}

#' Read SMEAR datasets
#'
#' @param data_file RDS file of particle concentration data over time.
#' The first column is a timestamp, and the rest are are particle concentrations
#' (d100e1 means Particle concentration dNdlogDp Dp=1.00nm)
#' @param event_data_file txt file of event classification data
#'
#' @return list of two data frames
read_smear_all <- function(data_file = "../ideas-for-smear-data/data/banana-data.rds", 
                           event_data_file = "../ideas-for-smear-data/data/DMPS_Event_Classification.txt") {
  if (file.exists(data_file)) {
    smear_data <- readRDS(data_file) 
  } else {
    cat("Downloading data from SMEAR API...\n")
    smear_data <- get_smear_dmps(from="1996-01-20", to="2018-12-31") # 225 MB
    cat("Saving data to", data_file)
    saveRDS(smear_data, data_file) 
  }
  smear_data = smear_data[!is.na(smear_data$timestamp), ]
  smear_data = smear_data[apply(smear_data[,-1], 1, function(x) !all(is.na(x))), ]
  colnames(smear_data) = gsub("HYY_DMPS.", "", colnames(smear_data))
  events = read_smear_events(event_data_file)
  non_event_days = strftime(events[events$type == "non.event",]$timestamp, "%Y-%m-%d")
  i_non_event = strftime(smear_data$timestamp, "%Y-%m-%d") %in% non_event_days
  event_days = strftime(events[events$type == c("event.Ia", "event.Ib"),]$timestamp, "%Y-%m-%d")
  i_event = strftime(smear_data$timestamp, "%Y-%m-%d") %in% event_days
  
  list(data = smear_data, 
       events = events, 
       i_non_event = i_non_event, 
       i_event = i_event)
}

#' Read and clean event classification
#'
#' @param path_to_events path to DMPS_Event_Classification.txt
#'
#' @return data.frame
read_smear_events <- function(path_to_events="data/DMPS_Event_Classification.txt") {
  events <- read.csv(path_to_events, sep=" ", comment.char = "%", 
                     col.names=c("datenum", "event.Ia", "event.Ib", "event.II", "event.apple", "eventbump", "event.rain", "event.featureless", "non.event", "undefined", "bad.data", "partly.bad.data", "checksum"))
  
  events$timestamp <- as.POSIXct((events$datenum - 719529)*86400, origin = "1970-01-01", tz = "UTC") # ?as.POSIXct() matlab datenum conversion
  events$type <- NA #initialize
  types.in.dmps.data <- c("event.Ia", "event.Ib", "event.II", "non.event", "undefined", "bad.data", "partly.bad.data")
  for (ev.type in types.in.dmps.data) 
    events$type[as.logical(events[ , ev.type])] <- ev.type
  events[c("timestamp", "type")]
}

#' Get data for bananaplot from SMEAR API
#' 
#' Download measurements from the SMEAR API https://avaa.tdata.fi/web/smart/smear/api
#'
#' @param from,to dates as %Y-%m-%d
#'
#' @return
#' @export
#'
#' @examples
#' d = get_smear_dmps("2010-01-01", "2010-01-03")
get_smear_dmps <- function(from, to) {
  format_date <- function(d) paste0(format(as.POSIXct(d, origin="1970-01-01"), 
                                           format = "%Y-%m-%d"), "%20", 
                                    format(as.POSIXct(d, origin="1970-01-01"), 
                                           format = "%H:%M:%S"))
  varnames = c("d100e1", "d112e1",
               "d126e1", "d141e1", "d158e1", "d178e1", "d200e1", "d224e1", "d251e1",
               "d282e1", "d316e1", "d355e1", "d398e1", "d447e1", "d501e1", "d562e1",
               "d631e1", "d708e1", "d794e1", "d891e1", "d100e2", "d112e2", "d126e2",
               "d141e2", "d158e2", "d178e2", "d200e2", "d224e2", "d251e2", "d282e2",
               "d316e2", "d355e2", "d398e2", "d447e2", "d501e2", "d562e2", "d631e2",
               "d708e2", "d794e2", "d891e2", "d100e3", "d112e3", "d126e3", "d141e3",
               "d158e3", "d178e3", "d200e3", "d224e3", "d251e3", "d282e3", "d316e3",
               "d355e3", "d398e3", "d447e3", "d501e3", "d562e3", "d631e3", "d708e3",
               "d794e3", "d891e3", "d100e4")
  url <- paste0("https://avaa.tdata.fi/smear-services/smeardata.jsp?",
                "table=HYY_DMPS",
                "&variables=", paste(varnames, collapse=","),
                "&from=", format_date(from),
                "&to=", format_date(to),
                "&quality=any",
                "&averaging=none",
                "&type=arithmetic",
                "&format=csv")
  data = read.csv(url, stringsAsFactors = F)
  timestamp <- with(data, as.POSIXct(paste(Year, Month, Day, Hour, Minute, Second), 
                                     format = "%Y %m %d %H %M %OS"))
  vars_keep = setdiff(names(data), 
                      c(names(data)[sapply(data, function(x) all(is.na(x)))], 
                        "Year", "Month", "Day", "Hour", "Minute", "Second"))
  cbind(timestamp=timestamp, data[ , vars_keep])
}
