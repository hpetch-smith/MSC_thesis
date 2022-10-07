# Takes date and time string
get_date <-function(text){
  # Split into list of chars 
  text <- strsplit(text,'')[[1]]
  #blank list to store parse
  date <- list()
  # loop through string
  for (i in 1:length(text)){
    # isolate each char
    char <-  text[[i]]
    # Check for blank spaces 
    if (char == " "){
      date <- toString(date)
      date <- gsub(',','',date)
      date <- gsub(' ','',date)
      return(date)
    }
    # append char to new string
    date <- append(date, char)
  }
}

# Takes date and time string
get_time <-function(text){
  # Split into list of chars 
  text <- strsplit(text,'')[[1]]
  #blank list to store parse
  time <- list()
  # loop through string
  for (i in length(text):1){
    # isolate each char
    char <-  text[[i]]
    # Check for blank spaces 
    if (char == " "){
      # reverse the date 
      time <- rev(time)
      # return formmatted string
      time <- toString(time)
      time <- gsub(',','',time)
      time <- gsub(' ','',time)
      time <- gsub('-',':',time)
      return(time)
    }
    # append char to new string
    time <- append(time, char)
  }
}

# load tidyverse for filtering 
library(tidyverse)
# function to retireve days
get_day <- function(query){
  # Filter selected day 
    filt <- filter(final_data, final_data$date == query)
    return(filt)
}
