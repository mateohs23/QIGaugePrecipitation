library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(tidyr)
library(zoo)
library(data.table)
library(ggplot2)
library(hydrostats)
library(FAdist)
library(nortsTest)
library(ggrepel) 

# Define path
path <- "I:/Unidades compartidas/WadiLab_General/Study_cases/Mateo/Aricanduva/Quality_info/dados_plu"
files_list <- list.files(path, full.names = TRUE)
path_save <- "I:/Unidades compartidas/WadiLab_General/Study_cases/Mateo/Aricanduva/Quality_info/quality_info_code/Results/"

# Load files into a list of dataframes
dataframe_list_0 <- lapply(files_list, fread)

## Process dataframe ##
# Arrange the dataframes

# Create list of gauges #
gauges_PLU <- list()
gauges_NPLU <- list()

for (i in 1:length(dataframe_list_0)){
  prova <-  dataframe_list_0[[i]]
  # Identification of reapeted column names
  if (any(duplicated(names(prova)))){
    # Identify duplicate column names
    col_names <- names(prova)
    duplicate_indices <- which(duplicated(col_names))
    # Rename duplicate columns
    for (n in duplicate_indices) {
      # Create a new name by appending a suffix
      new_name <- paste0(col_names[n], "_", sum(col_names[1:n] == col_names[n]))
      col_names[n] <- new_name
    }
    # Assign the new names back to the data frame
    names(prova) <- col_names
  }
  # Establish the DATA format's date 
  prova <-  prova%>%
      mutate(DATA = ymd_hm(DATA)) 
  # Organize the data 
  if (Reduce(`|`, colnames(prova)%in%"PLU(mm)")) {
    if (Reduce(`|`, colnames(prova)%in%"V3") && mean(prova$"V3", na.rm = TRUE) < 100) {
      if (Reduce(`|`, colnames(prova)%in%"Q(m3/s)")) {
        prova <- prova %>%
          rename(bat = `Q(m3/s)`,
                 `Q(m3/s)` = `FLU(m)`,
                 `FLU(m)` = `PLU(mm)`,
                 `PLU(mm)` = `V3`)
      } else if (Reduce(`|`, colnames(prova)%in%"Vazao (m2s)")) {
        prova <- prova %>%
          rename(bat = `Vazao (m2s)`,
                 `Q(m3/s)` = `FLU(m)`,
                 `FLU(m)` = `PLU(mm)`,
                 `PLU(mm)` = `V3`)
      } else if (Reduce(`|`, colnames(prova)%in%"FLU(m)")){
        prova <- prova %>%
          rename(bat = `FLU(m)`,
                 `FLU(m)` = `PLU(mm)`,
                 `PLU(mm)` = `V3`)
      } else if (Reduce(`|`, colnames(prova)%in%"FLU(m).1")){ 
        prova <- prova %>%
          rename(bat = `FLU(m).1`,
                 `FLU(m)_2` = `FLU(m)`,
                 `FLU(m)` = `PLU(mm)`,
                 `PLU(mm)` = `V3`)
      } else {
        prova <- prova %>%
          rename(`FLU(m)` = `PLU(mm)`,
                 `PLU(mm)` = `V3`)
      }
    } else if (Reduce(`|`, colnames(prova)%in%"Vazao (m2s)")) {
      prova <- prova %>% rename(`Q(m3/s)` = `Vazao (m2s)`)
    }
    gauges_PLU <- append(gauges_PLU,list(prova))
  } else {
    gauges_NPLU <- append(gauges_NPLU,list(prova))
  }
}

# Create list for normality test
gauge_normality_test <- list()

for(j in 1:length(gauges_PLU)){
  ## Calculate intensities in 10 min ##
  # Create a null dataframe with all data period with intervals of 10 min
  prova2 <- gauges_PLU[[j]] 
  
  first_r <- prova2 %>% filter(minute(DATA) %% 10 == 0) %>% slice(1) %>% pull(DATA)
  last_r <- max(prova2$DATA)
  intervals <- seq(from = first_r, to = last_r, by = "10 mins")
  new_df <- tibble(DATA = intervals)

  # Fill the data  
  if (Reduce(`|`, colnames(prova2)%in%"Q(m3/s)")) {
    new_df <- left_join(new_df, prova2 %>% select(DATA, Posto, `PLU(mm)`, `FLU(m)`, `Q(m3/s)`), by = "DATA")
    } else if (Reduce(`|`, colnames(prova2)%in%"FLU(m)")) {
      new_df <- left_join(new_df, prova2 %>% select(DATA, Posto, `PLU(mm)`, `FLU(m)`), by = "DATA")%>%
      mutate(`Q(m3/s)` = NA)
    } else {
      new_df <- left_join(new_df, df %>% select(DATA, Posto, `PLU(mm)`), by = "DATA")%>% 
      mutate(`FLU(m)` = NA, `Q(m3/s)` = NA)
    }

  # Intensities calculation
  #new_df <- new_df %>% 
  #  mutate(`Precipitation(mm/h)` = pmax((`PLU(mm)` - lag(`PLU(mm)`)) * 6, 0, na.rm = TRUE))

  # Instant PLU (mm) calculation
  new_df <- new_df %>% 
    mutate(`Instant PLU(mm)` = pmax((`PLU(mm)` - lag(`PLU(mm)`)),0, na.rm = FALSE))

  ### Statistics analysis ###
  new_df3 <- subset(new_df, select = -c(`PLU(mm)`,`FLU(m)`,`Q(m3/s)`) )
  new_df3 <- drop_na(new_df3)
  ## Normal Test ##
  ADTest <- normal.test(new_df3$`Instant PLU(mm)`, normality = "ad", alpha =0.05)
  ## Save the test ##
  gauge_normality_test[[j]] <- list()
  names(gauge_normality_test)[[j]] <- mean(summary(new_df3$Posto))
  gauge_normality_test[[j]] <- ADTest
  

  ### Quality analysis ###
  ## First analysis ##
  # Found the values of Instant Precipitation "Instant PLU(mm)") > 40 mm
  prova3 <- new_df %>%
    filter(`Instant PLU(mm)` > 40)

  new_df2 <- new_df %>%
    left_join(select(prova3,DATA,`Instant PLU(mm)`), by = "DATA")%>% 
    mutate (Match = ifelse(!is.na(`Instant PLU(mm).y`),2,1))

  #clean the dataframe
  new_df2 <- new_df2 %>% 
    select(-`Instant PLU(mm).y`) %>%
    rename( `Instant PLU(mm)` = `Instant PLU(mm).x`)

  # Replace values greater than 40
  prova4 <- new_df %>%
    mutate(`Instant PLU(mm)` = case_when(
    `Instant PLU(mm)` > 40 ~ 40,
    TRUE ~  `Instant PLU(mm)`))

  ## Second analysis##
  # Accumulate the Instant PLU(mm) in daily data
  prova5 <- prova4 %>%
    mutate(date = as.Date(DATA)) %>%
    group_by(date) %>%
    reframe(AcumDailyPre = sum(`Instant PLU(mm)`),MeanDailyPre = mean(`Instant PLU(mm)`))

  # Identify sequences greater than 60 days without Data 
  prova6 <- prova5 %>%
    mutate(missing = is.na(AcumDailyPre)) %>%
    arrange(date) %>%
    mutate(consecutive_missing = cumsum(missing != lag(missing, default = first(missing)))) %>%
    group_by(consecutive_missing, missing) %>%
    reframe(start_date = min(date), end_date = max(date),length = n(), .groups = "drop" ) %>%
    filter(missing == TRUE & length >= 60)

  #Save the periods of station's missing Data
  if(nrow(prova6)!=0){
    for (i in 1:nrow(prova6)) {
      first_r <- ymd_hms(paste(prova6$start_date[i], "00:00:00"))
      last_r <- ymd_hms(paste(prova6$end_date[i], "23:50:00"))
      intervals <- seq(from = first_r, to = last_r, by = "10 mins")
      pivot_df <- tibble(DATA = intervals)
      new_df2 <- new_df2 %>%
        mutate(Match = ifelse(DATA %in% pivot_df$DATA, 0, Match))
    }
  }
  ## Third analysis ##
  #Replace NA by 0 for rolling average analysis
  prova7 <- prova5 %>%
    mutate(AcumDailyPre = ifelse(is.na(AcumDailyPre), 0, AcumDailyPre)) 

  #Calculate a rolling average of 35 days of accumulation daily Precipitation
  RA35DayAcum <- rollmean(prova7$AcumDailyPre,k=35,fill=NA, align = "center")
  prova7 <- cbind(prova7,RA35DayAcum)

  #Identify if RA35DayAcum is equal to 0.2 (+-10%) mm (gauge's minimum measure) 
  prova8 <- prova7 %>%
    arrange(date) %>%
    mutate(Pre_minValue = RA35DayAcum >= 0.18 & RA35DayAcum <= 0.22) %>%
    group_by(Pre_minValue) %>%
    reframe(start_date = date - days(35), end_date = date + days(35), .groups = "drop" ) %>%
    filter(Pre_minValue == TRUE)

  #Save the periods of station's sequences of RA35DayAcum near to gauge's minimum measure
  if(nrow(prova8)!=0){
    for (i in 1:nrow(prova8)) {
      first_r <- ymd_hms(paste(prova8$start_date[i], "00:00:00"))
      last_r <- ymd_hms(paste(prova8$end_date[i], "23:50:00"))
      intervals <- seq(from = first_r, to = last_r, by = "10 mins")
      pivot_df <- tibble(DATA = intervals)
      new_df2 <- new_df2 %>%
        mutate(Match = ifelse(DATA %in% pivot_df$DATA, 0, Match))
    }
  }
  ## Fourth analysis##
  # Identify the sequences of 200 days or greater of RA35DayAcum is equal to zero  
  prova9 <- prova7 %>%
    arrange(date) %>%
    mutate(Pre_0mm = RA35DayAcum <= 0.01) %>%
    mutate(consecutive_0mm = cumsum(Pre_0mm != lag(Pre_0mm, default = first(Pre_0mm)))) %>%
    group_by(consecutive_0mm, Pre_0mm) %>%
    reframe(start_date = min(date), end_date = max(date),length = n(), .groups = "drop" ) %>%
    filter(Pre_0mm == TRUE & length >= 200)

  #Save the periods of station's sequences of RA35DayAcum equal to zero
  if(nrow(prova9)!=0){
    for (i in 1:nrow(prova9)) {
      first_r <- ymd_hms(paste(prova9$start_date[i], "00:00:00"))
      last_r <- ymd_hms(paste(prova9$end_date[i], "23:50:00"))
      intervals <- seq(from = first_r, to = last_r, by = "10 mins")
      pivot_df <- tibble(DATA = intervals)
      new_df2 <- new_df2 %>%
        mutate(Match = ifelse(DATA %in% pivot_df$DATA, 0, Match))
    }
  }  
    gauges_PLU[[j]] <- new_df2
}

##Plot##
# Specify the path for the directory
period_data_path <- paste(path_save, "Period data",sep = "")

# Check if the directory exists
if (!dir.exists(period_data_path)) {
  # If it doesn't exist, create the directory
  dir.create(period_data_path)
  cat("Directory created:", period_data_path, "\n")
} else {
  cat("Directory already exists:", period_data_path, "\n")
}

for (i in 1:length(gauges_PLU)){
  ## Dataframe arrange ##
  new_df2 <- gauges_PLU[[i]]%>%
    mutate(year=year(DATA),month=(month(DATA)))#%>%
    #mutate(hydro_year = ifelse (month>= 10, year, year-1))
  gauge_idx <- new_df2$Posto[1]
  new_df2 <- subset(new_df2, select = -c(`Posto`,`PLU(mm)`,`FLU(m)`,`Q(m3/s)`,`Instant PLU(mm)`) )
  
  #create a column with the 3 classifications
  new_df2$Category<-ifelse(new_df2$Match == 0, "Invalid Data",ifelse(new_df2$Match == 1,"Valid Data","Maximum-change"))

  # Convert the Category column to a factor with specified levels to control the order
  new_df2$Category <-factor(new_df2$Category, levels = c("Invalid Data", "Valid Data", "Maximum-change"))

  # Plotting histogram of valid data every hydrological year
  df_pivot<-new_df2%>%
    group_by(year,month,Category)%>%
    summarise(n=n(),.groups = "drop")

  df_pivot$Category <- droplevels(df_pivot$Category)# Remove unused levels if necessary

  fig_pivot <- ggplot(df_pivot, aes(x = as.factor(month),y = n,fill= Category))+
    geom_bar(stat = "identity", position ="fill") +
    labs(title = paste("Valid / Invalid Data - Gauge:",gauge_idx),
         x = "Month",
         y = "Portion") +
    facet_wrap(~year)+
    scale_fill_manual(values = c("Invalid Data" = "red", 
                                 "Valid Data" = "blue",
                                 "Maximum-change" = "green"), 
                      labels = c("Invalid Data", "Valid Data", "Maximum-change")) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7.5), # Rotate x-axis labels for better readability
      axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
      axis.ticks.length = unit(0.1, "cm"))

  print(fig_pivot)
  
  # Saving plot
  fig_name <- paste(period_data_path,"/","Valid_Data",gauge_idx,".jpg",sep = "")
  ggsave(filename = fig_name, 
         plot = fig_pivot, 
         width = 10,  # Width in inches
         height = 6,  # Height in inches
         dpi = 300)   # Resolution in dots per inch (DPI)

}

##############################################################################################################
### Double-mass analysis###
double_mass_path <- paste(path_save,"Double_mass",sep = "") #path to save the figures
# Check if the directory exists
if (!dir.exists(double_mass_path)) {
  # If it doesn't exist, create the directory
  dir.create(double_mass_path)
  cat("Directory created:", double_mass_path, "\n")
} else {
  cat("Directory already exists:", double_mass_path, "\n")
}

## Double-mass calculation ##
double_mass_table <- data.frame(Year=0) #creating the double mass table

for (i in 1:length(gauges_PLU)){
  ## Dataframe arrange ##
  prova10 <- gauges_PLU[[i]]%>%
    mutate(year=year(DATA),month=(month(DATA)))
  prova10 <- subset(prova10, select = -c(`PLU(mm)`,`FLU(m)`,`Q(m3/s)`,Match) )%>%
    mutate(hydro_year = ifelse (month>= 10, year, year-1))
  gauge <- prova10$Posto[1]

  ## calculating the accumulation by year for each gauge ##
  years <- unique(prova10$hydro_year[duplicated(prova10$hydro_year)]) #record period#
  double_mass_table$gauge <- NA
  colnames(double_mass_table)[colnames(double_mass_table) == "gauge"] <- as.character(gauge)
  yearly_sum_pre <- c(gauge) #vector that will be fill with yearly accumulation precipitation values#
  for (j in  1:length(years)){
    if (!(years[j] %in% double_mass_table$Year)){
      # If the year does not exist, add it to the dataframe
      double_mass_table <- bind_rows(double_mass_table, data.frame( Year = years[j]))
    }
    df_pivot <- prova10 %>%
      filter(prova10$hydro_year == years[j]) #filter the year#
    sum_pre <- sum(df_pivot$`Instant PLU(mm)`, na.rm = TRUE) #accumulation of the precipitation#
    row_index <- which(double_mass_table$Year == years[j])
    col_name <- as.character(gauge)
    double_mass_table[row_index, col_name] <- sum_pre #save the result#
  }
}

## Dataframe rearrange ##
double_mass_table <- double_mass_table %>% 
  arrange(Year) %>%
  filter(Year != 0)%>%
  mutate(across(everything(), ~ if_else(. == 0, NA_real_, .)))%>%
  drop_na()

## Calculating the mean per year and add to the table ##
row_means <- rowMeans(double_mass_table[, -1], na.rm = TRUE) #mean per year#
double_mass_table$Mean <- row_means #adding to table#

## Calculating the cumulative precipitation ##
for (i in 2:length(colnames(double_mass_table))){
  cum_sum <- cumsum(ifelse(is.na(double_mass_table[,i]), 0, double_mass_table[,i])) #cumulative sum per year
  new_col_name <- paste('Cum_Sum',colnames(double_mass_table)[i]) #new column name
  double_mass_table[[new_col_name]] <- cum_sum #adding the result in the double mass table
}

## Plot the double mass curve by gauge ##
mean_idx <- which(colnames(double_mass_table) == 'Mean') #position of mean column
double_mass_simp<-double_mass_table[,-c(2:mean_idx)]

# Remove the prefix "Cum_Sum" from all relevant columns
colnames(double_mass_simp) <- gsub("^Cum_Sum\\s*", "", colnames(double_mass_simp))

# Figure for each gauge
keygauge_idx <- as.character(colnames(double_mass_simp[3:length(double_mass_simp)-1]))

for (i in 1:length(keygauge_idx)){
  gauge_idx <- keygauge_idx[i]
  # Plot creation
  fig_pivot <- ggplot(double_mass_simp, aes(x = Mean, y = .data[[gauge_idx]])) +
    geom_line(color = "blue") + # Line for the first dataset
    geom_point(color = "red") + # Points for visibility
    geom_text_repel(aes(label = Year), size = 5) +
    labs(title = paste("Double-Mass Curve", "Gauge:",gauge_idx),
        x = "Cumulative Mean Precipitation (mm)",
        y = "Cumulative Precipitation Gauge (mm)") +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.2),  # Customize both x and y axis lines
      axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
      axis.ticks.length = unit(0.1, "cm"))

  print(fig_pivot)
  
  # Saving plot
  fig_name <- paste(double_mass_path,"/","dmc_",gauge_idx,".jpg",sep = "")
  ggsave(filename = fig_name, 
        plot = fig_pivot, 
        width = 10,  # Width in inches
        height = 6,  # Height in inches
        dpi = 300)   # Resolution in dots per inch (DPI)
}

## Plot between two stations (e.g. 1000858 and 1000490) ##
gauge_idx1 <- "1000858"
gauge_idx2 <- "1000490"
fig_pivot <- ggplot(double_mass_simp) +
  geom_line(aes(x = Mean, y = .data[[gauge_idx1]], color = "Gauge 1"),linewidth = 0.8) + # Line for the first dataset
  geom_line(aes(x = Mean, y = .data[[gauge_idx2]], color = "Gauge 2"),linewidth = 0.8) + # Line for the second dataset
  geom_point(aes(x = Mean, y = .data[[gauge_idx1]]),color = "red") + # Points for visibility idx1
  geom_point(aes(x = Mean, y = .data[[gauge_idx2]]),color = "red") + # Points for visibility idx2
  geom_text_repel(aes(x = Mean, y = .data[[gauge_idx1]],label = Year), size = 5) +
  geom_text_repel(aes(x = Mean, y = .data[[gauge_idx2]],label = Year), size = 5) +
  labs(title = paste("Double-Mass Curve", "Gauge:",gauge_idx1, "and", gauge_idx2),
       x = "Cumulative Mean Precipitation (mm)",
       y = "Cumulative Precipitation Gauge (mm)",
       color = "Gauge") +
  scale_color_manual(values = c("Gauge 1" = "blue", "Gauge 2" = "green"),
                     labels = c(gauge_idx1, gauge_idx2)) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.2),  # Customize both x and y axis lines
    axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
    axis.ticks.length = unit(0.1, "cm"))

print(fig_pivot)

# Saving plot
fig_name <- paste(double_mass_path,"/","dmc_",gauge_idx1,"_",gauge_idx2,".jpg",sep = "")
ggsave(filename = fig_name, 
       plot = fig_pivot, 
       width = 10,  # Width in inches
       height = 6,  # Height in inches
       dpi = 300)   # Resolution in dots per inch (DPI)

##Plotting the General Double Mass Curve##
# Get the index of the 'Mean' column
mean_idx <- which(colnames(double_mass_table) == 'Mean')
double_mass_long<-double_mass_table[,-c(2:mean_idx)]
# Remove the prefix "Cum_Sum" from all relevant columns
colnames(double_mass_long) <- gsub("^Cum_Sum\\s*", "", colnames(double_mass_long))

# Create a vector of gauge names (assuming they are in columns after 'Mean')
gauge_names <- as.character(colnames(double_mass_long)[3:ncol(double_mass_long)-1])

# Initialize a long format dataframe for ggplot
double_mass_long <- double_mass_long %>%
  select(`Mean`, Year, all_of(gauge_names)) %>%
  pivot_longer(cols = all_of(gauge_names), 
               names_to = "Gauge", 
               values_to = "Cumulative_Precipitation")

# Create a color palette and line types for each gauge
colors <- c("black", rainbow(length(gauge_names)-1))  # Generate distinct colors
line_types <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length.out = length(gauge_names)) #became as a long table
double_mass_long <- double_mass_long %>%
  mutate(Year = ifelse(Gauge == as.character(gauge_names[1]), Year, NA)) # Only label for one gauge

general_dmc <- ggplot(double_mass_long, aes(x = Mean, y = Cumulative_Precipitation, color = Gauge)) +
  geom_line() +         # Add lines for each gauge
  geom_point(aes(color = Gauge)) +                # Add points for visibility
  geom_text_repel(aes(label = Year), size = 5, na.rm = TRUE) + # Add labels
  labs(title = "General Double-Mass Curve",
       x = "Cumulative Mean Precipitation (mm)",
       y = "Cumulative Precipitation by Gauge (mm)") +
  #theme_minimal()+
  theme(
    axis.line = element_line(color = "black", linewidth  = 0.2),  # Customize both x and y axis lines
    axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
    axis.ticks.length = unit(0.1, "cm"))+                   # Length of tick
  scale_color_manual(values = colors)+          # Assign colors to gauges
  scale_linetype_manual(values = line_types)    # Assign line types to gauges

print(general_dmc)

# Save the plot with specified dimensions and resolution
fig_name <- paste(double_mass_path,"/","general_dmc.jpg",sep = "")
ggsave(filename = fig_name, 
       plot = general_dmc, 
       width = 10,  # Width in inches
       height = 6,  # Height in inches
       dpi = 300)   # Resolution in dots per inch (DPI)

############################################################################################################
### Accumulation time analysis###
Accum_time_path <- paste(path_save,"Accum_time",sep = "") #path to save the figures
# Check if the directory exists
if (!dir.exists(Accum_time_path)) {
  # If it doesn't exist, create the directory
  dir.create(Accum_time_path)
  cat("Directory created:", Accum_time_path, "\n")
} else {
  cat("Directory already exists:", Accum_time_path, "\n")
}

## Accumulation time calculation ##
accum_time_table <- data.frame(Day="1900-01-01") #creating the double mass table

for (i in 1:length(gauges_PLU)){
  ## Dataframe arrange ##
  prova11 <- gauges_PLU[[i]]%>%
    mutate(year=year(DATA),month=month(DATA),day=day(DATA))
  prova11 <- subset(prova11, select = -c(`PLU(mm)`,`FLU(m)`,`Q(m3/s)`,Match) )
  gauge_idx <- prova11$Posto[1]
  
  ## Routine to accumulate daily values each gauge ##
  prova11 <- prova11%>%
    group_by(year, month, day) %>%
    summarise(daily_total = sum(`Instant PLU(mm)`), .groups = "drop") %>%
    unite(.,"Day",year,month,day,sep = "-")
  
  
  ## Add the result in the table ##
  record_days <- unique(prova11$Day) #record period#
  for (j in  1:length(record_days)){
    if (!(record_days[j] %in% accum_time_table$Day)){
      # If the year does not exist, add it to the dataframe
      accum_time_table <- bind_rows(accum_time_table, data.frame(Day = record_days[j]))
    }
  }
  accum_time_table <- accum_time_table %>%
    left_join(prova11, by = "Day")
  colnames(accum_time_table)[colnames(accum_time_table) == "daily_total"] <- as.character(gauge_idx)
}

## Dataframe rearrange ##
accum_time_table <- accum_time_table %>% 
  filter(Day != "1900-01-01")%>%
  mutate(Day = ymd(Day))%>%
  arrange(Day)
first_valid_row <- which(complete.cases(accum_time_table)) # Identify the first row without any NAs
if (length(first_valid_row) > 0) {
  first_valid_index <- first_valid_row[1] # Get the index of the first valid row
  accum_time_table <- accum_time_table[first_valid_index:nrow(accum_time_table), ] # Filter the dataframe to keep only rows from the first valid row onward
}else{
  print("No valid rows found. Check data")
}

## Calculating the mean daily precipitation ##
row_means <- rowMeans(accum_time_table[, -1], na.rm = TRUE) #mean per day#
accum_time_table$Mean <- row_means #adding to table#

## Calculating the cumulative daily precipitation ##
for (i in 2:length(colnames(accum_time_table))){
  cum_sum <- cumsum(ifelse(is.na(accum_time_table[,i]), 0, accum_time_table[,i])) #cumulative daily precipitation sum 
  new_col_name <- paste('Cum_Sum',colnames(accum_time_table)[i]) #new column name
  accum_time_table[[new_col_name]] <- cum_sum #adding the result in the accum time table
}

## Plot the Mass-Time Accumulation Time by gauge ##
mean_idx <- which(colnames(accum_time_table) == 'Mean') #position of mean column
accum_time_simp<-accum_time_table[,-c(2:mean_idx)]
accum_time_simp<-accum_time_simp %>%
  mutate(Day_Count = row_number())

# Remove the prefix "Cum_Sum" from all relevant columns
colnames(accum_time_simp) <- gsub("^Cum_Sum\\s*", "", colnames(accum_time_simp))

# Figure for each gauge
keygauge_idx <- as.character(colnames(accum_time_simp[3:length(accum_time_simp)-1]))

for (i in 1:length(keygauge_idx)){
  gauge_idx <- keygauge_idx[i]
  # Plot creation
  fig_pivot <- ggplot(accum_time_simp) +
    geom_line( aes(x = Day_Count, y = .data[[gauge_idx]],color = "Gauge")) + # Line for the gauge
    geom_line( aes(x = Day_Count, y = .data$Mean, color = "Mean")) + # Line for the mean
    labs(title = paste("Mass Accumulation", "Gauge:",gauge_idx),
         x = "Cumulative Time (days)",
         y = "Cumulative Precipitation Gauge (mm)") +
    scale_color_manual(values = c("Gauge" = "blue", "Mean" = "red"),
                        labels = c(gauge_idx, "Mean")) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.2),  # Customize both x and y axis lines
      axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
      axis.ticks.length = unit(0.1, "cm"))
  
  print(fig_pivot)
  
  # Saving plot
  fig_name <- paste(Accum_time_path,"/","mat_",gauge_idx,".jpg",sep = "")
  ggsave(filename = fig_name, 
         plot = fig_pivot, 
         width = 10,  # Width in inches
         height = 6,  # Height in inches
         dpi = 300)   # Resolution in dots per inch (DPI)
}

## Plot between two stations (e.g. 1000858 and 1000490) ##
gauge_idx1 <- "1000858"
gauge_idx2 <- "1000490"
fig_pivot <- ggplot(accum_time_simp) +
  geom_line( aes(x = Day_Count, y = .data[[gauge_idx1]],color = "Gauge 1"), linewidth = 1) + # Line for the gauge 1
  geom_line( aes(x = Day_Count, y = .data[[gauge_idx2]],color = "Gauge 2"), linewidth = 1) + # Line for the gauge 2
  geom_line( aes(x = Day_Count, y = .data$Mean, color = "Mean"), linewidth = 1) + # Line for the mean
  labs(title = paste("Mass Accumulation", "Gauge:",gauge_idx1, "and", gauge_idx2),
       x = "Cumulative Time (days)",
       y = "Cumulative Precipitation Gauge (mm)",
       color = "Gauge") +
  scale_color_manual(values = c("Gauge 1" = "blue", "Gauge 2" = "green","Mean" = "red"),
                     labels = c(gauge_idx1, gauge_idx2, "Mean")) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.2),  # Customize both x and y axis lines
    axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
    axis.ticks.length = unit(0.1, "cm"))

print(fig_pivot)

# Saving plot
fig_name <- paste(Accum_time_path,"/","mat_",gauge_idx1,"_",gauge_idx2,".jpg",sep = "")
ggsave(filename = fig_name, 
       plot = fig_pivot, 
       width = 10,  # Width in inches
       height = 6,  # Height in inches
       dpi = 300)   # Resolution in dots per inch (DPI)

##Plotting the General Mass Accumulation Time Curve##
# Get the index of the 'Mean' column
mean_idx <- which(colnames(accum_time_table) == 'Mean')
accum_time_long<-accum_time_table[,-c(2:mean_idx)]
# Remove the prefix "Cum_Sum" from all relevant columns
colnames(accum_time_long) <- gsub("^Cum_Sum\\s*", "", colnames(accum_time_long))
accum_time_long<-accum_time_long %>%
  mutate(Day_Count = row_number())

# Create a vector of gauge names (assuming they are in columns after 'Mean')
gauge_names <- as.character(colnames(accum_time_long)[3:ncol(accum_time_long)-1])

# Initialize a long format dataframe for ggplot
accum_time_long <- accum_time_long %>%
  select(Day_Count,Day, all_of(gauge_names)) %>%
  pivot_longer(cols = all_of(gauge_names), 
               names_to = "Gauge", 
               values_to = "Cumulative_Daily_Precipitation")

# Create a color palette and line types for each gauge
colors <- c(rainbow(length(gauge_names)-1),"black")  # Generate distinct colors
line_types <- rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length.out = length(gauge_names)) #became as a long table
accum_time_long <- accum_time_long %>%
  mutate(Day = ifelse(Gauge == as.character(gauge_names[1]), ymd(Day), NA)) # Only label for one gauge

general_mat <- ggplot(accum_time_long, aes(x = Day_Count, y = Cumulative_Daily_Precipitation, color = Gauge)) +
  geom_line() +         # Add lines for each gauge
  labs(title = "General Mass Accumulation Curve",
       x = "Cumulative Time (days)",
       y = "Cumulative Precipitation Gauge (mm)") +
  theme_minimal()+
  theme(
    axis.line = element_line(color = "black", linewidth = 0.2),  # Customize both x and y axis lines
    axis.ticks = element_line(color = "black", linewidth = 0.2), # Customize ticks
    axis.ticks.length = unit(0.1, "cm"))+                   # Length of tick
  scale_color_manual(values = colors) +          # Assign colors to gauges
  scale_linetype_manual(values = line_types)    # Assign line types to gauges

print(general_mat)

# Save the plot with specified dimensions and resolution
fig_name <- paste(Accum_time_path,"/","general_mat.jpg",sep = "")
ggsave(filename = fig_name, 
       plot = general_mat, 
       width = 10,  # Width in inches
       height = 6,  # Height in inches
       dpi = 300)   # Resolution in dots per inch (DPI)
