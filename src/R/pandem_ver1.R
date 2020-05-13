# setClass('myDate')
# setAs("character","myDate", function(from) as.Date(from, format="%m-%d-%Y") )

library(ggplot2)

path.dir <- "C:/Users/rouna/Documents/GitHub/PanDem"
data.dir <- file.path(path.dir,"src/temp")
data.USA <- read.table(file.path(data.dir,"USARegionStats.csv"),sep=",",
                       colClasses = c("character","NULL","integer","integer","integer"),
                       header = TRUE)

data.USA$Date <- as.Date(data.USA$Date,format = "%m-%d-%y")
data.confirmed = ts(data.USA$Confirmed)
data.deaths = ts(data.USA$Deaths)
data.recovered = ts(data.USA$Recovered)

plot(data.confirmed,type='h')
plot(data.deaths,type='h')
plot(data.recovered,type='h')

log.confirmed <- log(data.confirmed)

ggplot(data = data.USA, aes(x = Date, y = Confirmed)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Total confirmed cases in USA",
       subtitle = "Spring 2020",
       x = "Date", y = "Daily Confirmed CoViD-19 cases")