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

#' Download and prepare smear data
#'
#' @param file_smear 
#' @param file_smear_raw 
#' @param file_events 
#'
#' @return
#' @export
#'
#' @examples
check_smear = function(file_smear="data/smear.rds", 
                       file_smear_raw="data/smear_raw.rds", 
                       file_events="data/DMPS_Event_Classification.txt") {
  if (!file.exists(file_smear)) {
    if (!file.exists(file_smear_raw))  {
      cat("Downloading data from SMEAR API (225MB)...\n")
      # Apparently there is some new unspecified limit for how much data to download per request.
      dates = seq(as.Date("1996-01-20"), as.Date("2018-12-31"), by=50)
      D = list() 
      for (i in 1:(length(dates)-1)) {
        print(sprintf("%d/%d from %s to %s", i, (length(dates)-1), dates[i], dates[i+1]))
        D[[i]] = get_smear_dmps(dates[i], dates[i+1])
      }
      raw = Reduce(rbind, D, init=data.frame())
      cat("Saving to", file_smear_raw, "\n")
      saveRDS(raw, file_smear_raw) 
    }
    cat("Preparing raw SMEAR data...\n")
    d = prepare_smear(file_smear_raw, file_events)
    cat("Saving to", file_smear, "\n")
    saveRDS(d, file_smear)
  }
}

#' Prepare SMEAR data
#'
#' @param file_smear_raw RDS file of particle concentration data over time.
#' The first column is a timestamp, and the rest are are particle concentrations
#' (d100e1 means Particle concentration dNdlogDp Dp=1.00nm)
#' @param file_events txt file of event classification data
#'
#' @return list of two data frames
prepare_smear <- function(file_smear_raw = "data/smear_raw.rds", 
                          file_events = "data/DMPS_Event_Classification.txt") {
  stopifnot(file.exists(file_smear_raw), file.exists(file_events))
  d = readRDS(file_smear_raw)
  d = d[,sapply(d, function(x) !all(is.na(x)))]
  events = read_smear_events(file_events)
  d = d[!is.na(d$timestamp), ]
  d = d[apply(d[,-1], 1, function(x) !all(is.na(x))), ]
  colnames(d) = gsub("HYY_DMPS.", "", colnames(d))
  non_event_days = strftime(events[events$type == "non.event",]$timestamp, "%Y-%m-%d")
  i_non_event = strftime(d$timestamp, "%Y-%m-%d") %in% non_event_days
  event_days = strftime(events[events$type == c("event.Ia", "event.Ib"),]$timestamp, "%Y-%m-%d")
  i_event = strftime(d$timestamp, "%Y-%m-%d") %in% event_days
  
  list(data = d, 
       events = events, 
       i_non_event = i_non_event, 
       i_event = i_event)
}

#' Read and clean SMEAR event classification
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
#' Download measurements from the SMEAR API:  
#' https://smear-backend.rahtiapp.fi/swagger-ui/#/default/get_search_timeseries_csv
#'
#' @param from,to dates as %Y-%m-%d
#'
#' @return
#' @export
#'
#' @examples
#' d = get_smear_dmps("2010-01-01", "2010-12-31")
get_smear_dmps <- function(from, to) {
  # Timestamp in ISO 8601 format 
  format_ts = function(x) gsub(":", "%3A", strftime(x, format="%FT%T.000"))
  varnames = c("d100e1", "d112e1", "d126e1", "d141e1", "d158e1", "d178e1",
               "d200e1", "d224e1", "d251e1", "d282e1", "d316e1", "d355e1", "d398e1",
               "d447e1", "d501e1", "d562e1", "d631e1", "d708e1", "d794e1", "d891e1",
               "d100e2", "d112e2", "d126e2", "d141e2", "d158e2", "d178e2", "d200e2",
               "d224e2", "d251e2", "d282e2", "d316e2", "d355e2", "d398e2", "d447e2",
               "d501e2", "d562e2", "d631e2", "d708e2", "d794e2", "d891e2", "d100e3",
               "d112e3", "d126e3", "d141e3", "d158e3", "d178e3", "d200e3", "d224e3",
               "d251e3", "d282e3", "d316e3", "d355e3", "d398e3", "d447e3", "d501e3",
               "d562e3", "d631e3", "d708e3", "d794e3", "d891e3", "d100e4")
  var_query = paste0("&tablevariable=HYY_DMPS.", varnames, collapse="")
  url = paste0("https://smear-backend.rahtiapp.fi/search/timeseries/csv?",
               "aggregation=ARITHMETIC", 
               "&from=", format_ts(from),
               "&to=", format_ts(to),
               "&interval=10", 
               "&quality=ANY",
               var_query)
  data = read.csv(url, stringsAsFactors = FALSE)
  timestamp <- with(data, as.POSIXct(paste(Year, Month, Day, Hour, Minute, Second), 
                                     format = "%Y %m %d %H %M %OS"))
  vars_keep = setdiff(names(data), 
                      #c(names(data)[sapply(data, function(x) all(is.na(x)))], 
                      c("Year", "Month", "Day", "Hour", "Minute", "Second"))
  d = cbind(timestamp=timestamp, data[ , vars_keep])
  j_ordered = names(d)[-1][order(as.numeric(gsub("HYY_DMPS.d", "", names(d)[-1], fixed=T)))]
  d[,c("timestamp", j_ordered)]
}
