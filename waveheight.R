# Calculate the wave heights for the AWH data logger.
# The calculation uses the zero-crossing method to determine wave heights and
# wave periods.
# Greg Nishihara
# 2020 Sep 30
################################################################################
library(tidyverse)
library(lubridate)
library(oceanwaves)

# If the time is messed up when being read in, then maybe it needs to be setup correctly.
Sys.setlocale("LC_TIME", "en_US.UTF-8") 

# Functions --------------------------------------------------------------------

# This function is to read the data files for the AWH datalogger.
read_awh = function(fnames) {
  X = read_lines(fnames, n_max = 100)
  interval = str_subset(X, "^Interval") %>% str_extract("[0-9]+") %>% as.numeric()
  startfrom = str_which(X, "\\[Item\\]")
  
  read_csv(fnames, skip = startfrom, locale = locale(encoding = "CP932")) %>% 
    select(datetime = 日時,
           mpa = contains("MPa"),
           depth = matches("\\[m\\]")) %>% 
    mutate(intervals = seq(0, by = 0.1, length = nrow(.))) %>% 
    mutate(intervals = as.integer(10*(intervals %% 1))) %>% 
    mutate(datetime = str_glue("{datetime}.{intervals}")) %>% 
    mutate(datetime = parse_date_time(datetime, orders = "Ymd HMOS")) %>% 
    select(datetime, mpa, depth)
}

# The original version of this function is from the oceanwaves package. 
# However, this version will not work on the data logger data. This is because,
# the data logger data is intermittent.
# I removed some code and added modifications.
# However, there will be errors in the results since not all of the data
# segments provide enough information to calculate the wave heights.
# To circumvent this, use possibly().
waveStatsZC_gnn = function (data, Fs, threshold = NULL) 
{
  data <- -data
  detrended <- oceanwaves::detrendHeight(data)
  data = detrended[["pt"]][]
  d0 <- data[!oceanwaves:::almostZero(data, 0)]
  back0 <- 1:length(data)
  back0 <- back0[!oceanwaves:::almostZero(data, 0)]
  f <- which(d0[1:(length(d0) - 1)] * d0[2:length(d0)] < 0)
  crossing <- back0[f]
  if (data[1] > 0) {
    crossing <- crossing[-1]
  }
  
  crossing <- crossing[seq(1, length(crossing), by = 2)]
  wave <- matrix(NA, nrow = length(crossing) - 1, ncol = 4)
  for (i in 1:(length(crossing) - 1)) {
    wave[i, 2] <- max(data[crossing[i]:crossing[(i + 1)]])
    wave[i, 3] <- -min(data[crossing[i]:crossing[(i + 1)]])
  }
  if (nrow(wave) > 0) {
    wave[, 4] <- diff(crossing)/Fs
    if (is.null(threshold)) {
      threshold <- 0.01 * max(wave[, 2] + wave[, 3])
    }
    else if (threshold < 0) {
      stop("threshold must be greater than 0")
    }
    i <- 0
    while (i < nrow(wave)) {
      i <- i + 1
      if (wave[i, 2] < threshold) {
        if (i != 1) {
          wave[i - 1, 2:4] <- c(max(wave[(i - 1):i, 2]), 
                                max(wave[(i - 1):i, 3]), 
                                sum(wave[(i - 1):i, 4]))
        }
        wave <- wave[-i, ]
      }
      else if (wave[i, 3] < threshold) {
        if (i != nrow(wave)) {
          wave[i, 2:4] <- c(max(wave[i:(i + 1), 2]), 
                            max(wave[i:(i + 1), 3]), 
                            sum(wave[i:(i + 1), 4]))
          wave <- wave[-(i + 1), ]
        }
        else {
          wave <- wave[-i, ]
        }
      }
    }
    wave[, 1] <- wave[, 2] + wave[, 3]
  }
  nb <- nrow(wave)
  wave <- wave[order(wave[, 1], decreasing = TRUE), ]
  wave = matrix(wave, ncol = 4)
  tibble(H = wave[,1], Tau = wave[, 4])
}

# possible() will return NULL if the function fails. This is the safe way to
# run the code with map().
ZC_try = possibly(waveStatsZC_gnn, NULL)

run_ZC = function(df, Fs = 10) {
  ZC_try(df$depth, Fs = Fs)
}

# This function is used to calculate the mean significant wave height or
# the mean of the significant wave periods, which are calculated from the top third. 
sig = function(x) {
  x = x[order(x, decreasing = TRUE)]
  mean(x[1:round(length(x) / 3)])
}

# This function is used to calculate the mean of the top 10% of the wave heights
# or wave periods.
top10 = function(x) {
  x = x[order(x, decreasing = TRUE)]
  mean(x[1:round(length(x) / 10)])
}

################################################################################

# Read the wave logger data. ---------------------------------------------------
data = tibble(fnames = dir("~/Lab_Data/wave_logger/", full = TRUE))

# The data is taken at 10Hz for 10 seconds every 10 minutes.
data = data %>% mutate(data = map(fnames, read_awh))

data = data %>% mutate(location = str_extract(fnames, "garamo|inside")) %>% 
  select(location, data) %>% unnest(data)

HZ = 10 # 10 Hz or 10 Samples/second

# Group the data according to the sampling setup, which is 100 samples at
# 10 Hz for 1 minute.
data = data %>% 
  mutate(sampling_group = floor_date(datetime, "minute")) %>% 
  group_nest(location, sampling_group)


# Run the zero-crossing function on the data groups.
# This will take a few seconds to run. 
data = data %>% mutate(data = map(data, run_ZC, Fs = HZ))


# Remove the data that could not be processed, then expand the nested data.
data = data %>% 
  mutate(date = as_date(sampling_group)) %>% 
  select(location, date, data) %>% 
  drop_na() %>% 
  unnest(data)

# Calculate the mean wave height and period, significant mean height and period,
# top 10 % of the wave height and periods for each day.

data2 = data %>% group_by(date, location) %>% 
  summarise(N = length(H),
            Hmean = mean(H),
            Hsig = sig(H),
            H10 = top10(H),
            Tmean = mean(Tau),
            Tsig = sig(Tau),
            T10 = top10(Tau))

# Just a plot to take a look at the results. ------------------------------------
 
ylabel  = "Daily mean wave heights (m)"
p1 = ggplot(data2) + 
  geom_point(aes(x = date, y = Hmean, color = location)) +
  geom_line(aes(x = date, y = Hmean, color = location))  +
  scale_x_date() +
  scale_y_continuous(ylabel, limits = c(0, 0.5)) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.background = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank())

ylabel  = "Daily mean signficant wave heights (m)"
p2 = ggplot(data2) + 
  geom_point(aes(x = date, y = Hsig, color = location)) +
  geom_line(aes(x = date, y = Hsig, color = location))  +
  scale_x_date() +
  scale_y_continuous(ylabel, limits = c(0, 0.5)) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.background = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank())

ylabel  = "Daily mean top 10% wave heights (m)"
p3 = ggplot(data2) + 
  geom_point(aes(x = date, y = H10, color = location)) +
  geom_line(aes(x = date, y = H10, color = location))  +
  scale_x_date() +
  scale_y_continuous(ylabel, limits = c(0, 0.5)) +
  scale_color_brewer(palette = "Set2") +
  theme(legend.background = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title = element_blank())

ggpubr::ggarrange(p1,p2,p3, ncol = 3)









