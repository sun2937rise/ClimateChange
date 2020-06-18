# cleaning the environment
rm(list=ls())

# calling package ncdf4 for working with netcdf4 format of the data
library(ncdf4)

# Declare data frame
df <- NULL

#Open all hourly files
files <- list.files("C:/Users/13127/OneDrive/Documents/consulting/group project/HDF", pattern='*.nc4', full.names=TRUE)

# Loop over files to open all files
for(i in seq_along(files)) {
  
  nc <- nc_open(files[i])
  
  if (is.null(df)){
    #getting longitude and latitude
    ylat <- nc$dim$nlat$vals
    xlon <- nc$dim$nlon$vals
    
    # getting precipitation
    rain_matrix <- ncvar_get(nc, "precipitation")
    
    # create dataframe -- reshape data
    # matrix (nlon*nlat rows by 2 columns) of lons and lats
    lonlat <- as.matrix(expand.grid(xlon, ylat))
    
    # transposes rain_matrix to match up with grid and then creates a vector of rain values
    rain_vec <- as.vector(t(rain_matrix))
    
    # create dataframe, add names
    rain_df01 <- data.frame(cbind(lonlat, rain_vec))
    names(rain_df01) <- c("lon","lat", "precipitation")
    df <- rbind(df, data.frame(rain_df01))
  }
  
  else {
     
    rain_matrix <- ncvar_get(nc, "precipitation")
    rain_vec <- as.vector(t(rain_matrix))
    rain_df02 <- data.frame(cbind(lonlat, rain_vec))
    names(rain_df02) <- c("lon","lat", "precipitation")
    
    #Add the precipitation from each file to a single data.frame
    df <- cbind(df, rain_df02$precipitation) 
  }
  
  # closing each file after using it(preventing error)
  print(i)
  nc_close(nc)
}

# checking for the missing data
sum(is.na(df)) # 11174 NA 
sum(is.na(df[,6058]))  # 8875 of them from 01/27/2016 9pm(column 6058)

# imputing NA's with seasonal adjustment+interpolation
library(imputeTS)
dfdf=na_seadec(df, algorithm = "interpolation")
sum(is.na(dfdf))

# calculating daily precipitation by adding up hourly precipitation
df1 = NULL

for(i in 1:(length(files)/8)){
  t=8*(i-1)+3
  s=8*i+2
  df1=cbind(df1,rowSums(dfdf[,t:s]))
  print(i)
}

# combining lon and lat with daily precipitation
df2 = cbind(dfdf[,c(1,2)],df1)

setwd("C:/Users/13127/OneDrive/Documents/consulting/group project/daily_data")
fileNames = Sys.glob("*.nc4")
# extracting date from the name of the daily files
col_names = substring(fileNames,12,19)

# adding date as the name of the columns
colnames(df2)[3:1828] = col_names

# getting the name of the states for each longitude and latitude
library(sp)
library(maps)
library(maptools)

latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  #IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  IDs <- states$names
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

testPoints <- data.frame(df2$lon, df2$lat)
## testPoints is the dataframe of longitude and latitude columns

# add name of the states to data frame
df2$states <- latlong2state(testPoints)

# checking for missing data
sum(is.na(df2$states))

# remove the rows with N/A states which are out of the U.S
df3 = na.omit(df2)

# rearrange the column and bring states to the third column
library(dplyr)
df4 = select(df3, lon, lat, states, "20140101":"20181231")

# writing exel file of daily data
setwd("C:/Users/13127/OneDrive/Documents/consulting/group project")
write.csv(df4, file = "dailydata.csv")

#############################################################################

# cleaning the environment
rm(list=ls())

# read the daily file
df5 = read.csv("dailydata.csv")

# converting data from wide to long
library(tidyr)
df6 = gather(df5, date, precipitation, -c("lon","lat","states"))

# get a subset of daily data with only consider wet days
df7 = subset(df6, precipitation > 0.254)

# convert data from long to wide again
library(tidyr)
df8= spread(df7, date, precipitation)

# calculating average of daily precipitation from 2014 to 2018
df8$mean_daily_precip = rowMeans(df8[,-c(1:3)], na.rm = TRUE)
write.csv(df8[,c(1,2,3,1830)], file = "Avg_dailydata.csv")

# drawing heat map of U.S
library(ggplot2)
library(maps)
library(mapdata)

states <- map_data("state")
US_map = ggplot(data = states) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  coord_fixed(1.3) 
  
US_map + geom_tile(aes(lon, lat, fill=mean_daily_precip), data=df8, alpha=0.8)+
  geom_polygon(data=states,aes(x = long, y = lat, group = group),color = "black", fill = NA)+
  scale_fill_gradientn(colours = rev(rainbow(7)))+
  scale_x_continuous("Longitude",breaks = seq(-120,-70, by = 10)) + 
  scale_y_continuous("Latitude",breaks = seq(25,50, by = 5)) +
  theme_bw()+
  labs(title = "Average daily precipitation(mm/day) during 2014-2018 in the U.S") +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 18, face = "bold"),
        rect = element_blank(),panel.grid.major = element_line(),
        legend.justification = "top")

###############  separate data for each season ##################

# filter on data to get only winter season
library(dplyr)
df_winter = df7 %>% filter(substring(date,6,7) %in% c("12","01","02"))

# reshape data and make it wide again
df_winter = spread(df_winter, date, precipitation)

# calculating means of the daily precipitation for winter season from 2014 to 2018
df_winter$mean_winter = rowMeans(df_winter[,-c(1,2,3,4)], na.rm = TRUE)

# drawing heatmap of U.S for winter season
US_map + geom_tile(aes(lon, lat, fill=mean_winter), data=df_winter, alpha=0.8)+
  geom_polygon(data=states,aes(x = long, y = lat, group = group),color = "black", fill = NA)+
  scale_fill_gradientn(colours = rev(rainbow(7)))+
  scale_x_continuous("Longitude",breaks = seq(-120,-70, by = 10)) + 
  scale_y_continuous("Latitude",breaks = seq(25,50, by = 5)) +
  theme_bw()+
  labs(title = "Average daily precipitation(mm/day)in winter during 2014-2018") +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 18, face = "bold"),
        rect = element_blank(),panel.grid.major = element_line(),
        legend.justification = "top")

# filter on data to get only spring season
library(dplyr)
df_spring = df7 %>% filter(substring(date,6,7) %in% c("03","04","05"))

# reshape data and make it wide again
df_spring = spread(df_spring, date, precipitation)

# calculating means of the daily precipitation for spring season from 2014 to 2018
df_spring$mean_spring = rowMeans(df_spring[,-c(1,2,3,4)], na.rm = TRUE)

# drawing heatmap of U.S for spring season
US_map + geom_tile(aes(lon, lat, fill=mean_spring), data=df_spring, alpha=0.8)+
  geom_polygon(data=states,aes(x = long, y = lat, group = group),color = "black", fill = NA)+
  scale_fill_gradientn(colours = rev(rainbow(7)))+
  scale_x_continuous("Longitude",breaks = seq(-120,-70, by = 10)) + 
  scale_y_continuous("Latitude",breaks = seq(25,50, by = 5)) +
  theme_bw()+
  labs(title = "Average daily precipitation(mm/day) in spring during 2014-2018 in the U.S") +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 18, face = "bold"),
        rect = element_blank(),panel.grid.major = element_line(),
        legend.justification = "top")

# filter on data to get only summer season
df_summer = df7 %>% filter(substring(date,6,7) %in% c("06","07","08"))
df_summer= spread(df_summer, date, precipitation)
df_summer$mean_summer = rowMeans(df_summer[,-c(1,2,3,4)], na.rm = TRUE)

# drawing heatmap of U.S for summer season
US_map + geom_tile(aes(lon, lat, fill=mean_summer), data=df_summer, alpha=0.8)+
  geom_polygon(data=states,aes(x = long, y = lat, group = group),color = "black", fill = NA)+
  scale_fill_gradientn(colours = rev(rainbow(7)))+
  scale_x_continuous("Longitude",breaks = seq(-120,-70, by = 10)) + 
  scale_y_continuous("Latitude",breaks = seq(25,50, by = 5)) +
  theme_bw()+
  labs(title = "Average daily precipitation (mm/day) in summer during 2014-2018 in the U.S") +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 18, face = "bold"),
        rect = element_blank(),panel.grid.major = element_line(),
        legend.justification = "top")
  

# filter on data to get only fall season
df_fall = df7 %>% filter(substring(date,6,7) %in% c("09","10","11"))
df_fall= spread(df_fall, date, precipitation)
df_fall$mean_fall = rowMeans(df_fall[,-c(1,2,3,4)], na.rm = TRUE)

# drawing heatmap of U.S for fall season
US_map + geom_tile(aes(lon, lat, fill=mean_fall), data=df_fall, alpha=0.8)+
  geom_polygon(data=states,aes(x = long, y = lat, group = group),color = "black", fill = NA)+
  scale_fill_gradientn(colours = rev(rainbow(7)))+
  scale_x_continuous("Longitude",breaks = seq(-120,-70, by = 10)) + 
  scale_y_continuous("Latitude",breaks = seq(25,50, by = 5)) +
  theme_bw()+
  labs(title = "Average daily precipitation (mm/day) in fall during 2014-2018 in the U.S") +
  theme(panel.grid = element_blank(), 
        plot.title = element_text(size = 18, face = "bold"),
        rect = element_blank(),panel.grid.major = element_line(),
        legend.justification = "top")
#########################################################################################
# creating monthly data

# reshape data frame from wide to long
library(tidyr)
df9 = gather(df8, date, precipitation, -c("lon","lat","states","mean_daily_precip"))


library(dplyr)
df_jan = df9 %>% filter(substring(date,6,7) %in% "01")
df_jan= spread(df_jan, date, precipitation)
df_jan$mean_jan = rowSums(df_jan[,-c(1,2,3,4)], na.rm = TRUE)/5
df_jan = select(df_jan, states, mean_jan)
colnames(df_jan)[2] = "Jan"

df_feb = df9 %>% filter(substring(date,6,7) %in% "02")
df_feb= spread(df_feb, date, precipitation)
df_feb$mean_feb = rowSums(df_feb[,-c(1,2,3,4)], na.rm = TRUE)/5
df_feb = select(df_feb, states, mean_feb)
colnames(df_feb)[2] = "Feb"

df_mar = df9 %>% filter(substring(date,6,7) %in% "03")
df_mar= spread(df_mar, date, precipitation)
df_mar$mean_mar = rowSums(df_mar[,-c(1,2,3,4)], na.rm = TRUE)/5
df_mar = select(df_mar, states, mean_mar)
colnames(df_mar)[2] = "Mar"

df_Apr = df9 %>% filter(substring(date,6,7) %in% "04")
df_Apr= spread(df_Apr, date, precipitation)
df_Apr$mean_apr = rowSums(df_Apr[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Apr = select(df_Apr, states, mean_apr)
colnames(df_Apr)[2] = "Apr"

df_May = df9 %>% filter(substring(date,6,7) %in% "05")
df_May= spread(df_May, date, precipitation)
df_May$mean_may = rowSums(df_May[,-c(1,2,3,4)], na.rm = TRUE)/5
df_May = select(df_May, states, mean_may)
colnames(df_May)[2] = "May"

df_Jun = df9 %>% filter(substring(date,6,7) %in% "06")
df_Jun= spread(df_Jun, date, precipitation)
df_Jun$mean_jun = rowSums(df_Jun[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Jun = select(df_Jun, states, mean_jun)
colnames(df_Jun)[2] = "Jun"

df_Jul = df9 %>% filter(substring(date,6,7) %in% "07")
df_Jul= spread(df_Jul, date, precipitation)
df_Jul$mean_jul = rowSums(df_Jul[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Jul = select(df_Jul, states, mean_jul)
colnames(df_Jul)[2] = "Jul"

df_Aug = df9 %>% filter(substring(date,6,7) %in% "08")
df_Aug= spread(df_Aug, date, precipitation)
df_Aug$mean_aug = rowSums(df_Aug[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Aug = select(df_Aug, states, mean_aug)
colnames(df_Aug)[2] = "Aug"

df_Sep = df9 %>% filter(substring(date,6,7) %in% "09")
df_Sep= spread(df_Sep, date, precipitation)
df_Sep$mean_sep = rowSums(df_Sep[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Sep = select(df_Sep, states, mean_sep)
colnames(df_Sep)[2] = "Sep"

df_Oct = df9 %>% filter(substring(date,6,7) %in% "10")
df_Oct= spread(df_Oct, date, precipitation)
df_Oct$mean_oct = rowSums(df_Oct[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Oct = select(df_Oct, states, mean_oct)
colnames(df_Oct)[2] = "Oct"

df_Nov = df9 %>% filter(substring(date,6,7) %in% "11")
df_Nov= spread(df_Nov, date, precipitation)
df_Nov$mean_nov = rowSums(df_Nov[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Nov = select(df_Nov, states, mean_nov)
colnames(df_Nov)[2] = "Nov"

df_Dec = df9 %>% filter(substring(date,6,7) %in% "12")
df_Dec= spread(df_Dec, date, precipitation)
df_Dec$mean_dec = rowSums(df_Dec[,-c(1,2,3,4)], na.rm = TRUE)/5
df_Dec = select(df_Dec, states, mean_dec)
colnames(df_Dec)[2] = "Dec"

monthly_data = cbind(df_jan, df_feb[,-1], df_mar[,-1], df_Apr[,-1], df_May[,-1], df_Jun[,-1], df_Jul[,-1], df_Aug[,-1], df_Sep[,-1],
      df_Oct[,-1], df_Nov[,-1], df_Dec[,-1])

names(monthly_data)[3:13] = c("Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# converting date variable to a date format
#library(lubridate)
#df10$date = ymd(df10$date)
# separating date and time
#library(tidyr)
#df2 = separate(df1, col = date_time, into = c("date", "time"), sep = "[^[:alnum:]]+")

# filtering to get only midwest states
midwest_data <- subset(monthly_data, states %in% c("illinois", "iowa", "indiana", "minnesota", "missouri", "michigan:north", "michigan:south", "wisconsin", "ohio"))

# calculating the average of the monthly precip for each state(from whole pixel in each state)
month_sums <- NULL
statelist <- unique(midwest_data$states)
for (i in statelist){
  single_state_m <- midwest_data[midwest_data$states == i,]
  state_sum_m <- colMeans(single_state_m[,2:13])
  month_sums <- cbind(month_sums, state_sum_m)
}

monthly_total <- data.frame("state" = rep(c("MN", "IA", "MO", "WI", "IL", "MIn", "IN", "MIs", "OH"), each = 12), 
                            "month" = c(1:12), 
                            "rain" = c(month_sums[1:12,]))

# ploting monthly precipitation for midwest
library(ggplot2)
ggplot(data = monthly_total) + 
  # making background white
  theme_bw() + 
  geom_line(aes(x = month, y=rain, linetype = state, colour = state), size = 1.5) +
  # giving lables to the plot axis and title
  labs(x="Month", y="Precipitation in mm",
       title="Average Monthly Total Precipitation 2014-2018", subtitle = "Monthly Trend", size = 1) +
  scale_x_continuous(breaks = seq(1, 12, by = 1), expand = c(0,0) , 
                     labels = c("J","F","M","A","M","J","J","A","S","O","N","D"))+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  # set legend position
  theme(legend.position="top", legend.direction="horizontal", 
        legend.text = element_text(size = 16), legend.title = element_text(size = 20),
        axis.text = element_text(size = 14), axis.title = element_text(size = 20),
        panel.grid = element_blank(), 
        plot.title = element_text(size = 22), plot.subtitle = element_text(size = 18))
 
# do the same thing for other states
other_states <- subset(monthly_data, states %in% c("washington:main", "florida", "texas", "california"))

month_sums <- NULL
statelist <- unique(other_states$states)
for (i in statelist){
  single_state_m <- other_states[other_states$states == i,]
  state_sum_m <- colMeans(single_state_m[,2:13])
  month_sums <- cbind(month_sums, state_sum_m)
}

monthly_total <- data.frame("state" = rep(c("WA", "CA", "TX", "FL"), each = 12), 
                            "month" = c(1:12), 
                            "rain" = c(month_sums[1:12,]))

library(ggplot2)
ggplot(data = monthly_total) + 
  # making background white
  theme_bw() + 
  geom_line(aes(x = month, y=rain, linetype = state, colour = state), size = 1.5) +
  # giving lables to the plot axis and title
  labs(x="Month", y="Precipitation in mm",
       title="Average Monthly Total Precipitation 2014-2018", subtitle = "Monthly Trend", size = 1) +
  scale_x_continuous(breaks = seq(1, 12, by = 1), expand = c(0,0) , 
                     labels = c("J","F","M","A","M","J","J","A","S","O","N","D"))+
  scale_y_continuous(limits = c(0,70), expand = c(0,0))+
  # set legend position
  theme(legend.position="top", legend.direction="horizontal", 
        legend.text = element_text(size = 16), legend.title = element_text(size = 20),
        axis.text = element_text(size = 14), axis.title = element_text(size = 20),
        panel.grid = element_blank(), 
        plot.title = element_text(size = 22), plot.subtitle = element_text(size = 18)) 


#######################################################################################

# Pdf plot of daily precip for each state
df11 = read.csv("Avg_dailydata.csv")
df12 = df11[,c(3,4)]

midwest_data <- subset(df12, states %in% c("illinois", "iowa", "indiana", "minnesota", "missouri", "michigan:north", "michigan:south", "wisconsin", "ohio"))

day_sums <- NULL
statelist <- unique(midwest_data$states)
for (i in statelist){
  single_state_m <- midwest_data[midwest_data$states == i,]
  state_sum_m <- colMeans(single_state_m[2])
  day_sums <- cbind(day_sums, state_sum_m)
}

daily_total <- data.frame("state" = c("MN", "IA", "MO", "WI", "IL", "MIu", "IN", "MIl", "OH"), 
                            "rain" = c(day_sums[1:9]))

library(ggplot2)
ggplot(daily_total, aes(x= state, y = rain)) + 
  geom_bar(aes(fill = state),stat = "identity") +
  labs(y = "Average daily preciptation(mm)",x="States",title = " Average daily total precipitation (mm/day) in Midwest", size = 1) +
  scale_y_continuous(limits = c(0,5), expand = c(0,0))+
  theme(plot.title=element_text(size=22,face="bold"),legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),axis.text = element_text(size = 14), 
        axis.title = element_text(size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


other_states = subset(df12, states %in% c("washington:main", "florida", "texas", "california","arkansas","louisiana","mississippi"))

day_sums <- NULL
statelist <- unique(other_states$states)
for (i in statelist){
  single_state_m <- other_states[other_states$states == i,]
  state_sum_m <- colMeans(single_state_m[2])
  day_sums <- cbind(day_sums, state_sum_m)
}

daily_total <- data.frame("state" = c("WA", "CA", "TX", "AR","LA","MS","FL"), 
                          "rain" = c(day_sums[1:7]))

library(ggplot2)
ggplot(daily_total, aes(x= state, y = rain)) + 
  geom_bar(aes(fill = state),stat = "identity") +
  labs(y = "Average daily preciptation(mm)",x="States",title = " Average daily total precipitation (mm/day) in other states", size = 1) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0))+
  theme(plot.title=element_text(size=22,face="bold"),legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),axis.text = element_text(size = 14), 
        axis.title = element_text(size = 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
############################################################################################
