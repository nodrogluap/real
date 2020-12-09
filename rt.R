library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(grid)
library(stringr)
library(zoo)

# Jurisdiction-specific data
data <- read.csv("ab_cases.csv", row.names=1)
df <- read.csv("ab_demographics.csv", header=FALSE); 
age_group2pop <- setNames(df[,2], df[,1])
testing <- read.csv("ab_testing.csv")
testing$Unknown <- testing[,2] - rowSums(testing[,3:ncol(testing)])

# SARS-CoV-2 specific settings for Rt calculation
rolling_window_size <- 15
viral_shedding_days <- 2:10
viral_shedding_proportions <- c(0.05,0.2,0.2,0.2,0.15,0.1,0.05,0.03,0.02)

data[,1] <- as.Date(data[,1])
data[,2] <- as.character(data[,2])
# ddz, i.e. data_by_date_and_zone
data %>% group_by_at(c(1,2)) %>% tally %>% mutate_if(is.character, str_replace_all, pattern = ' ', replacement = '.') %>% collect() -> ddz

Rt_per_zone <- function(zones, case_dates, cases_per_date, lead_days, day_weights, testing){
        # Sanity checks
        if(length(cases_per_date) != length(case_dates)){
                stop("'cases' and 'dates' arrays have different lengths")
        }
        if(length(lead_days) != length(day_weights)){
                stop("'lead_days' and 'day_weights' arrays have different lengths")
        }
        zones <- as.character(zones)
        cases_per_date <- as.integer(cases_per_date)
        lead_days <- as.integer(lead_days)
        case_dates <- as.Date(case_dates)
        day_weights <- as.numeric(day_weights)
        testing$Date.reported <- as.Date(testing$Date.reported)
        case_but_no_testing_dates = case_dates[which(!(case_dates %in% testing$Date.reported))]
        if(length(case_but_no_testing_dates) != 0){
                stop(paste0("cases were reported on date(s) without any testing statistics: ", case_but_no_testing_dates, "\n"))
        }
        # Use R's built-in hash table implementation
        zone_and_date2Rt <- new.env(hash=TRUE)
        for (zone in unique(zones)){
                case_but_no_testing_date_for_zone <- is.na(testing[which(case_dates %in% testing$Date.reported),zone]) || testing[which(case_dates %in% testing$Date.reported),zone] < 1
                if(!is.na(case_but_no_testing_date_for_zone) & sum(case_but_no_testing_date_for_zone) > 0){
                        stop(paste0("cases were reported in zone '", zone, "' on date(s) without any testing data: ", case_dates[case_but_no_testing_dates]))
                }
                zone_case_dates <- case_dates[zones == zone]
                zone_cases_per_date <- cases_per_date[zones == zone]
                for (i in 1:length(zone_case_dates)){
                        case_date <- zone_case_dates[i]
                        zone_and_date_key <- paste0(zone,as.character(case_date))
                        cases_this_date <- zone_cases_per_date[i]
                        if(cases_this_date < 1 || case_date + max(lead_days) > max(zone_case_dates)){ # No cases, or not enough dates downstream available, means no Rt estimate
                                zone_and_date2Rt[[zone_and_date_key]] <- NA
                                next
                        }
                        tests_this_date <- testing[(which(testing$Date.reported == case_date))[[1]],zone]
			if(length(tests_this_date) == 0){
				stop(paste0("Number of cases reported (", cases_this_date, ") in zone '", zone, "' on date ", case_date, " but no tests were performed that day"))
			}
                        if(cases_this_date > tests_this_date){
                                stop(paste0("Number of cases reported (", cases_this_date, ") in zone '", zone, "' on date ", case_date, " exceeds the number of tests performed (", tests_this_date, ")"))
                        }
                        valid_lead_day_num_cases <- c()
                        valid_lead_day_indices <- c()
                        valid_testing_index <- c()
                        for (ld in lead_days){
                                downstream_infection_date <- case_date + ld
                                which_testing_index <- which(testing$Date.reported == downstream_infection_date)
                                if(length(which_testing_index) == 0){
                                        # There were no tests performed on this date
                                        next;
                                }
                                # See if there were any cases on that future day
                                downstream_infection_date_cases_index <- which(zone_case_dates == downstream_infection_date)
                                # If no cases reported (but we know we had testing), assume cases count is 0 for this date
                                ld_cases <- ifelse(length(downstream_infection_date_cases_index) == 0, 0, zone_cases_per_date[downstream_infection_date_cases_index])
                                if(ld > testing[which_testing_index[[1]],zone]){
                                        stop(paste0("Number of cases reported (", ld_cases, ") in zone '", zone, "' on date ", downstream_infection_date,
                                                    " exceeds the number of tests performed (", testing[which_testing_index,zone], ")"))
                                }
                                valid_lead_day_indices <- append(valid_lead_day_indices, which(lead_days == ld))
                                valid_lead_day_num_cases <- append(valid_lead_day_num_cases, ld_cases)
                                valid_testing_index <- append(valid_testing_index, which_testing_index[[1]])
                        }
                        # Normalize the day_weights to add up to one based on what's available for future dates
                        valid_day_weights <- day_weights[valid_lead_day_indices]/sum(day_weights[valid_lead_day_indices])
                        # Use geometric math
                        attributed_downstream_cases <- 0
                        for (j in 1:length(valid_lead_day_indices)){
                                testing_weight <- tests_this_date/testing[valid_testing_index[j],zone]
                                attributed_downstream_cases <- attributed_downstream_cases + valid_lead_day_num_cases[j]*valid_day_weights[j]*testing_weight
                        }
                        # Rt is the derivative term (slope) from current new cases -> attributable new cases downstream (as defined by lead days and weights, tempered by testing numbers)
                        zone_and_date2Rt[[zone_and_date_key]] <- attributed_downstream_cases/cases_this_date
                }
        }
        vapply(paste0(zones,as.character(case_dates)), function(x){zone_and_date2Rt[[x]]}, numeric(1))
}
ddz$Rt <- Rt_per_zone(pull(ddz, colnames(ddz)[2]), pull(ddz, colnames(ddz)[1]), pull(ddz, colnames(ddz)[3]), viral_shedding_days, viral_shedding_proportions, testing)
dmin <- min(pull(ddz, colnames(ddz)[1]))
dmax <- max(pull(ddz, colnames(ddz)[1]))
ddz %>% arrange(Alberta.Health.Services.Zone) %>% group_by(Alberta.Health.Services.Zone) %>% complete(Date.reported = seq.Date(dmin, dmax, by="day")) %>% fill(colnames(ddz)[2]) %>% collect() -> ddz

ddz %>% group_by(Alberta.Health.Services.Zone) %>% mutate(rolling_mean = rollapply(Rt, rolling_window_size, mean, na.rm=TRUE, na.pad=TRUE)) %>% mutate(rolling_sd = rollapply(Rt, rolling_window_size, sd, na.rm=TRUE, na.pad=TRUE)) %>% collect() -> ddz
ddz$shapiro.p <- sapply(1:length(ddz$Rt), function(i){if(i+rolling_window_size>length(ddz$Rt)) return(NA); subset <- ddz$Rt[i:(i+rolling_window_size)]; ifelse(sum(!is.na(subset)) < 3, NA,unlist(shapiro.test(subset)[2]))})
# Pick between Gaussian and Poisson processes for CI estimation
ddz$ci_low <- ifelse(ddz$shapiro.p > 0.05, ddz$rolling_mean-2*ddz$rolling_sd, ddz$rolling_mean*(1-1.96/sqrt(rolling_window_size)))
# Obviously can't be negative
ddz$ci_low[ddz$ci_low < 0] <- 0
ddz$ci_high <- ifelse(ddz$shapiro.p > 0.05, ddz$rolling_mean+2*ddz$rolling_sd, ddz$rolling_mean*(1+1.96/sqrt(rolling_window_size)))
# Pick the more optimistic one
ddz$ci_high <- ifelse(ddz$ci_high > ddz$rolling_mean*(1+1.96/sqrt(rolling_window_size)), ddz$rolling_mean*(1+1.96/sqrt(rolling_window_size)), ddz$ci_high)

# Show the age range incidence levels as a heatmap
data %>% group_by_at(c(1,4)) %>% tally %>% mutate_if(is.character, str_replace_all, pattern = ' ', replacement = '.') %>% collect() -> age_group_tallies
# Order by first number in the age group label, or NA for Unknown and 0 for Under 1
age_group_tallies$order <- as.numeric(sapply(age_group_tallies$Age.group, function(x){o <- str_split(x, "[- +]")[[1]][1]; ifelse(o == "Under", 0, o)}))
age_group_tallies$Age.group.ordered <- fct_reorder(age_group_tallies$Age.group, age_group_tallies$order)

X11.options(type="dbcairo")

cimax <- max(ddz$ci_high, na.rm=TRUE)
# Display points beyond the averaging window as hints for future time (and info on very start)
ddz$Rt_only <- ifelse(is.na(ddz$ci_high) & ddz$n > 0, ddz$Rt, NA)

# Define the drawing order (last factor levels are drawn last and so are the most visible)
ddz$Alberta.Health.Services.Zone <- as.factor(ddz$Alberta.Health.Services.Zone) 
ddz$Alberta.Health.Services.Zone <- ordered(ddz$Alberta.Health.Services.Zone, levels = c("Unknown", "North.Zone", "Central.Zone", "South.Zone", "Edmonton.Zone", "Calgary.Zone"))
Rlevel <- levels(ddz$Alberta.Health.Services.Zone)
Rlabels <- sapply(Rlevel, function(zone){Rts_for_zone <- ddz[ddz$Alberta.Health.Services.Zone == zone, "rolling_mean"]; Rts_for_zone <- Rts_for_zone[!is.na(Rts_for_zone$rolling_mean),]; round(unlist(Rts_for_zone[dim(Rts_for_zone)[1],"rolling_mean"]), digits=2)})
Rlabels <- paste0(levels(ddz$Alberta.Health.Services.Zone), " (",Rlabels,")")

cip <- ggplot(ddz, aes(x=Date.reported, y=Rt_only, color=Alberta.Health.Services.Zone)) + geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill=Alberta.Health.Services.Zone), color=NA, alpha=0.2, show.legend = F) + geom_line(aes(y=rolling_mean), size=2) + theme(axis.text.x = element_text(angle=90)) + scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") + scale_y_continuous(position = "right", limits=c(0,cimax)) + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + geom_hline(yintercept=1, linetype="dotted", color="black") + labs(color="Virus reprod. number\ninstantaneous estimate\n(latest value)") +theme(axis.title.y=element_blank(), axis.ticks.y=element_blank())+geom_point()+ scale_colour_viridis_d(labels=Rlabels)+annotate("text", x=max(ddz$Date.reported),y=1,label = "Steady #\nof cases", hjust=0, vjust=0.5)
gA <- ggplotGrob(cip)

#ddz$n <- ifelse(is.na(ddz$n), 0,ddz$n)
np <- ggplot(ddz, aes(x=Date.reported, y=Alberta.Health.Services.Zone, fill=n))+geom_tile()+theme(panel.background=element_rect(fill="white", colour="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(axis.text.x = element_text(angle=90)) + scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y") + scale_y_discrete(position = "right") + coord_fixed(ratio=4)+labs(fill="New cases by zone")+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) + scale_fill_gradient(na.value = 'white', low = "grey95", high="darkslateblue")
gB <- ggplotGrob(np)

agep <- ggplot(age_group_tallies, aes(x=Date.reported, y=Age.group.ordered, fill=100*age_group_tallies$n/age_group2pop[levels(age_group_tallies$Age.group.ordered)[as.numeric(age_group_tallies$Age.group.ordered)]] ))+geom_tile()+theme(panel.background=element_rect(fill="white", colour="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_discrete(position = "right") + coord_fixed(ratio=4)+labs(fill="New cases as % of\npop. in age category")+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(), plot.margin = unit(c(0, 1, 1, 1), "cm")) + scale_fill_gradient(low = "grey95", high="darkred")
gC <- ggplotGrob(agep);

# Ensures that the three stacked graphs have the same width for the plot area, (i.e. x-axis dates line up) legend, etc.
g = rbind(gA, gB, gC, size = "last")
grid.draw(g)
