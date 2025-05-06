library(tidyverse)
library(data.table)
library(paletteer)
library(survival)
library(lubridate)
library(ggfortify)
library(survminer)
library(ggforestplot)
library(ggalluvial)
library(stringr)
library(ggrepel)
library(ggforce)
library(stringi)
library(patchwork)
library(lmtest)
library(pROC)
library(svglite)
library(conover.test)
library(epitools)
library(ggbeeswarm)
library(scales)

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

progression_colours <- c(
    "OAC" =  "#7f0000",
    "HGD, IMC" = "#ab1111",
    "LGD" = "#fc8d59",
    "Crypt dysplasia" = "#c6a6f5",
    "IND" = "#7dd5f0",
    "BE" = "#f0e4ce",
    "0" = "#c6c6c6", 
    "1" = "#ce1267a7",
    "2" = "#1280ce91",
    "3" = "#6ace1291",
    "Low risk" = "#fdf07d", 
    "Moderate risk" = "#fbcc76", 
    "High risk" = "#fd6159", 
    "High Confidence Positive" = "#cc83fc", 
    "Low Confidence Positive" = "#e8aeff", 
    "Negative" = "#acb1fb",
    "Hidden" = "#ffffff00"
)
three_scale <- c(
    "OAC" =  "#51004f",
    "HGD, IMC" = "#8856a7",
    "HGD, IMC, OAC" = "#8856a7",
    "LGD" = "#8c96c6",
    "Crypt dysplasia" = "#e6edff",
    "IND" = "#a6dff0",
    "BE" = "#f0e4ce",
    "0" = "#c6c6c6", 
    "1" = "#ffa600",
    "2" = "#d2749b",
    "3" = "#7d7195",
    "Low risk" = "#86c9b3", 
    "Moderate risk" ="#61a47e" , 
    "High risk" = "#4c6e5a", 
    "High Confidence Positive" = "#7d7195", 
    "Low Confidence Positive" = "#a16f9a", 
    "Negative" = "#c986a2",
    "Hidden" = "#ffffff00"
)
three_scale_v2 <- c(
    "OAC" =  "#51004f",
    "HGD/IMC" = "#8856a7",
    "HGD" = "#9070a5", 
    "IMC" = "#8856a7",
    "HGD, IMC, OAC" = "#8856a7",
    "LGD" = "#8c96c6",
    "LGD/Crypt dysplasia" = "#8c96c6",
    "Crypt dysplasia" = "#e7d9ff",
    "IND" = "#a6dff0",
    "BE" = "#e3e6e9",
    "0" = "#c8c9c7", 
    "1" = "#a5149a", #ffd700",
    "2" = "#e693dc",
    "3" = "#ffd43b",#e6f3df"
    "Low risk" = "#c0ddc0", 
    "Moderate risk" = "#90d692", 
    "High risk" = "#18a96d", 
    "High Confidence Positive" = "#a5149a",#8365c8"
    "Low Confidence Positive" = "#e693dc", #c453ba", 
    "Negative" ="#f5ecf4",
    "Hidden" = "#ffffff00", 
    "C" = "#ffb509", 
    "B" = "#ffd43b", 
    "A" = "#ffd43b", 
    " " = "#ffffff00", # single space 
    "  " = "#ffffff00" # double space 
)
three_scale_v3 <- c(
    "OAC" =  "#51004f",
    "HGD, IMC" = "#8856a7",
    "HGD" = "#9070a5", 
    "IMC" = "#8856a7",
    "HGD, IMC, OAC" = "#8856a7",
    "LGD" = "#8c96c6",
    "LGD/Crypt dysplasia" = "#8c96c6",
    "Crypt dysplasia" = "#e7d9ff",
    "IND" = "#a6dff0",
    "BE" = "#e3e6e9",
    "0" = "#c8c9c7", 
    "1" = "#a5149a", #ffd700",
    "2" = "#e693dc",
    "3" = "#ebb800",#e6f3df"
    "Low risk" = "#525252", 
    "Moderate risk" = "#90d692", 
    "High risk" = "#18a96d", 
    "High Confidence Positive" = "#a5149a",#8365c8"
    "Low Confidence Positive" = "#e693dc", #c453ba", 
    "Negative" ="#f5ecf4",
    "Hidden" = "#ffffff00", 
    "C" = "#ffb509", 
    "B" = "#ffd43b", 
    "A" = "#ffd43b", 
    " " = "#ffffff00", # single space 
    "  " = "#ffffff00" # double space 
)
progression_colours <- three_scale_v3
df <- fread(
    "/home/somers01/OptimalScreening/Data/DELTA_NHSE_cohort_TS_01Apr2025.csv", 
    na.strings = "", 
    stringsAsFactors = FALSE
)
kf <- fread(
    "/home/somers01/OptimalScreening/DataCleaning/Keith+Greta.csv", 
    na.string = ""
)
kf <- kf[, 
    c(
        "Age", 
        "Gender", 
        "Study_Number",
        "Date_of_Endo0",
        "Hx_dysplasia"
    )
]
kf2 <- fread(
    "/home/somers01/OptimalScreening/DataCleaning/DELTA sponge master RCF.csv",
     na.strings = ""
)
new_colnames = c(
    "Site", 
    "DELTA_ID",
    "Co_consented", 
    "Age", 
    "Sex", 
    "Date_sponge", 
    "Pot_number", 
    "CYT_ID", 
    "Date_endo0", 
    "Date_endo1", 
    "Histopathology_endo1", 
    "Date_endo2", 
    "Histopathology_endo2", 
    "Date_endo3", 
    "Histopathology_endo3", 
    "Date_endo4", 
    "Histopathology_endo4", 
    "Sponge_risk_strat", 
    "Comments", 
    "Prague_C0", 
    "Prague_M0", 
    "Atypia", 
    "P53", 
    "Greta_calls", 
    "Greta_comments"
)
colnames(df) <- new_colnames
df[df == "NA"] <- NA
df$Age <- df$Age %>% as.numeric
df <- df %>% filter(
    !is.na(Date_sponge)
)
StratifySpongeRisk <- function(Prague_C0, Prague_M0, Age, Sex, Atypia, P53, ...){
    C <- Prague_C0
    M <- Prague_M0 
    age <- Age
    sex <- Sex
    atypia <- Atypia %in% c("Y", "AUS")
    p53 <- P53 %in% c("Y", "E", "e")

    clin_mod_crit1 <- C > 6 || M > 10
    clin_mod_crit2 <- (M > 5 || C >=3) && (sex == "Male" || age > 60)

    if(p53 || atypia){
        return("High risk")
    } else if(clin_mod_crit1 || clin_mod_crit2){
        return("Moderate risk")
    } else{
        return("Low risk")
    }
}
StratifyClinRisk <- function(Prague_C0, Prague_M0, Age, Sex, ...){
    C <- Prague_C0
    M <- Prague_M0 
    age <- Age
    sex <- Sex
    
    clin_mod_crit1 <- C > 6 || M > 10
    clin_mod_crit2 <- (M > 5 || C >=3) && (sex == "Male" || age > 60)

    if(clin_mod_crit1 || clin_mod_crit2){
        return("Moderate risk")
    } else{
        return("Low risk")
    }
}
GetTimeToEvent <- function(
    Date_endo0, 
    Date_sponge, 
    Date_endo1, 
    Histopathology_endo1, 
    Date_endo2, 
    Histopathology_endo2, 
    Date_endo3, 
    Histopathology_endo3, 
    Date_endo4, 
    Histopathology_endo4, 
    ...){
    const <- 12/365.25
    events <- c(
        Histopathology_endo1,
        Histopathology_endo2, 
        Histopathology_endo3, 
        Histopathology_endo4
    )
    times <- c(
        Date_endo1,
        Date_endo2, 
        Date_endo3, 
        Date_endo4
    )
    outcomes <- c(
        "OAC", 
        "HGD", 
        "IMC", 
        "LGD", 
        "Crypt dysplasia"
    )
    endo0 <- as.Date(Date_endo0, format = "%d/%m/%Y")
    baseline <- as.Date(Date_sponge, format="%d/%m/%Y")
    time_to_sponge <- as.numeric((baseline - endo0)*const)
    gen_hist <- NA
    time_to_gen_hist <- NA

    # time to general histopathology ###########################################
    # if(Histopathology_endo1 == "IND"){
    #     stop_crit <- c(
    #         "NDBE",
    #         "BE", 
    #         "Crypt dysplasia", 
    #         "HGD", 
    #         "IMC", 
    #         "OAC", 
    #         "LGD"
    #     )
    #     for(i in seq_along(events[-1])){
    #         if(events[i+1] %in% stop_crit){
    #             gen_hist <- events[i+1]
    #             date_gen_hist <- as.Date(times[i+1], format = "%d/%m/%Y")
    #             time_to_gen_hist <- as.numeric(date_gen_hist - baseline)*const
    #             break
    #         }
    #     }
    # }else{
    #     gen_hist <- Histopathology_endo1 
    #     date_gen_hist <- as.Date(Date_endo1, format = "%d/%m/%Y")
    #     time_to_gen_hist <- as.numeric((date_gen_hist - baseline)*const)
    # }
    # if(is.na(gen_hist)){ # if no event occurs, get time as last endoscopy 
    #     last_endo_index <- max(which(!is.na(times)))
    #     date_gen_hist <- as.Date(times[last_endo_index], format="%d/%m/%Y")
    #     time_to_gen_hist <- as.numeric((date_gen_hist - baseline)*const)
    #     gen_hist <- events[last_endo_index]
    # }
    # if(gen_hist %in% outcomes){
    #     event <- 1
    # }else{
    #     event <- 0
    # }
    # time_to <- time_to_gen_hist

    #time to first diagnosis of dysplaisa ############################################

    for(i in seq_along(events)){
        if(events[i] %in% outcomes){
            gen_hist <- events[i] 
            date_endo <- as.Date(times[i], format = "%d/%m/%Y")
            time_to_gen_hist <- as.numeric((date_endo - baseline)*const)
            event <- 1 
            break
        }
    }
    if(is.na(gen_hist)){
        times <- times %>% rev 
        events <- events %>% rev
        for(i in seq_along(events)){
            if(!is.na(events[i])){
                gen_hist <- events[i] 
                date_endo <- as.Date(times[i], format = "%d/%m/%Y")
                time_to_gen_hist <- as.numeric((date_endo - baseline)*const)
                event <- 0 
                break
            }
        }
    }

    time_to <- time_to_gen_hist 
    
    tibble::tibble(time_to_sponge, time_to, event)
}
SurvivalStrat <- function(Sponge_risk_strat, Atypia, P53, ...){
    # return a new column with patients classified by sponge risk strat 
    # and if high risk, sub-stratified by biomarker status 
    if(Sponge_risk_strat == "High risk"){
        a <- toupper(Atypia)
        p <- toupper(P53)
        
        out <- paste0(a, p) 
        if(out == "YY"){
            return("High conf positive")
        }else{
            return("Low conf positive")
        }
    }else{
        return(Sponge_risk_strat)
    }
}
SurvivalStratAnyAtypia<- function(Sponge_risk_strat, Atypia, P53, ...){
    if(Sponge_risk_strat %in% c("Low risk", "Moderate risk")){
        return( Sponge_risk_strat) 
    }else{
        a <- toupper(Atypia) 
        p <- toupper(P53)
        if(a == "Y"){
            out <- "High Atypia"
        }else{
            out <- "High other"
        }
        return(out)
    }
}
GetBiomarkStatus <- function(Atypia, P53, ...){
    atypia <- toupper(Atypia)
    p53 <- toupper(P53) 

    if(atypia == "Y" && p53 == "Y"){
        return("High Confidence Positive")
    }else if(atypia == "Y" || p53 == "Y"){
        return("Low Confidence Positive")
    }else if(atypia == "AUS" || p53 == "E"){
        return("Low Confidence Positive")
    }else{
        return("Negative")
    }
}
GetFinalHist <- function(Histopathology_endo1, Histopathology_endo2, Histopathology_endo3, Histopathology_endo4, ...){
    paths <- c(
        Histopathology_endo4, 
        Histopathology_endo3, 
        Histopathology_endo2, 
        Histopathology_endo1
    )
    for(path in paths){
        if(is.na(path)){
            next
        }else{
            return(path)
        }
    }
}
GetAnyDysplasia <- function(Histopathology_endo1, Histopathology_endo2, Histopathology_endo3, Histopathology_endo4, ...){
    hists <- c(
        Histopathology_endo4, 
        Histopathology_endo3, 
        Histopathology_endo2, 
        Histopathology_endo1
    )
    value <- "F"
    for(hist in hists){
        if(hist %in% c("IND", "LGD", "HGD", "Crypt dysplasia", "OAC", "IMC")){
            value <- "T"
        }
        stop
    }
    return(value)
}
GetTimeE0E1 <- function(Date_endo0, Date_endo1, ...){
    const <- 12/365.25
    endo0 <- as.Date(Date_endo0, format = "%d/%m/%Y")
    endo1 <- as.Date(Date_endo1, format = "%d/%m/%Y")
    time <- as.numeric(endo1-endo0)*const 
    return(time)
}
GetTimeSpongeE1 <- function(Date_sponge, Date_endo1, ...){
    const <- 12/365.25
    sponge <- as.Date(Date_sponge, format = "%d/%m/%Y")
    endo1 <- as.Date(Date_endo1, format = "%d/%m/%Y")
    time <- as.numeric(endo1 - sponge)*const 
    return(time)
}
GetBSG <- function(Prague_M0, ...){
    if(Prague_M0 >= 3){
        return("Moderate risk")
    }else{
        return("Low risk")
    }
}
GetHistGeneral <- function(Histopathology_endo1, Histopathology_endo2, Histopathology_endo3, Histopathology_endo4, ...){
    if(Histopathology_endo1 == "IND"){
        hists = c(
            Histopathology_endo2, 
            Histopathology_endo3, 
            Histopathology_endo4
        )
        stop_crit <- c(
            "BE", 
            "Crypt dysplasia", 
            "HGD", 
            "IMC", 
            "OAC", 
            "LGD"
        )
        next_hist <- "IND"
        for(hist in hists){
            if(hist %in% stop_crit){
                next_hist <- hist
            }
        }
        return(next_hist)
    }else{
        return(Histopathology_endo1)
    }
}
df$Sponge_risk_strat <- pmap(df, StratifySpongeRisk) %>% as.character # applies pmap per row, single output 
df$Clin_risk_strat <- pmap(df, StratifyClinRisk) %>% as.character
df$BSG_risk_strat <- pmap(df, GetBSG) %>% as.character
df$Biomark_status <- pmap(df, GetBiomarkStatus) %>% as.character
df$Final_histopathology <- pmap(df, GetFinalHist) %>% as.character
df$General_histopathology <- pmap(df, GetHistGeneral) %>% as.character
df <- pmap_dfr(df, GetTimeToEvent) %>% bind_cols(df, .) # create multi-column df with pmap, appended to original df 
df$SurvStrat <- pmap(df, SurvivalStrat) %>% as.character %>% factor(levels = c(
    "Low risk", 
    "Moderate risk", 
    "Low conf positive",
    "High conf positive"
))
df$SurvStratAnyAtypia <- pmap(df, SurvivalStratAnyAtypia) %>% as.character %>% factor(levels = c(
    "Low risk", 
    "Moderate risk", 
    "High other", 
    "High Atypia"
))
df$Any_dysplasia <- pmap(df, GetAnyDysplasia) %>% as.character
df$time_to_E0E1 <- pmap(df, GetTimeE0E1) %>% as.numeric
df$time_to_spongeE1 <- pmap(df, GetTimeSpongeE1) %>% as.numeric
new_id_df <- data.frame()
for(risk in unique(df$Sponge_risk_strat)){
    new_rows <- df %>% 
        filter(Sponge_risk_strat == risk)%>% 
        arrange(time_to) %>%  
        group_by(time_to, DELTA_ID, CYT_ID) %>%
        mutate(ID = cur_group_id()) %>% 
        ungroup() 
    if(risk == "Low risk"){
        r <- "LR"
    }else if(risk == "Moderate risk"){
        r <- "MR"
    }else{
        r <- "HR"
    }
    new_rows$ID <- sapply(new_rows$ID, function(x) paste0(r, x))

    new_id_df <- rbind(new_id_df, new_rows)
}
df <- new_id_df
#fwrite(df, "Data/DELTA_working_df.csv")
#write_csv(df, "/home/somers01/OptimalScreening/Data/DELTA_df_tim.csv")
# ================= Data Cleanup for Sharing =============================
# run all the above code first 
df_clean <- df
df_clean$P53 <- toupper(df$P53) 
cleanup_colnames <- c(
    "Site", "DELTA_ID", "Co_consented",
    "Age", "Sex", "Date_sponge",
    "Pot_number", "CYT_ID", "Date_endo0",
    "Date_endo1", "Histopathology_endo1", "Date_endo2",
    "Histopathology_endo2", "Date_endo3", "Histopathology_endo3",
    "Date_endo4", "Histopathology_endo4", "Sponge_risk_strat",
    "Comments", "Prague_C0", "Prague_M0",
    "Atypia", "P53", "AI_calls",
    "AI_comments", "Clinical_risk_strat", "BSG_risk_strat",
    "Biomarker_status", "Final_histopathology", "Resolved_histopathology",
    "time_endo0_to_sponge", "time_sponge_to_resolved", "event",
    "SurvStrat", "SurvStratAnyAtypia", "Any_dysplasia",
    "time_endo0_to_endo1", "time_sponge_to_endo1", "Risk_strat_ID"
)
colnames(df_clean) <- cleanup_colnames
cleanup_colorder <- c(
    "Site", "DELTA_ID", "Co_consented",
    "Age", "Sex", "Date_sponge",
    "Pot_number", "CYT_ID", "Date_endo0",
    "Date_endo1", "Histopathology_endo1", "Date_endo2",
    "Histopathology_endo2", "Date_endo3", "Histopathology_endo3",
    "Date_endo4", "Histopathology_endo4", "Sponge_risk_strat",
    "Clinical_risk_strat", "BSG_risk_strat",
    "Comments", "Prague_C0", "Prague_M0",
    "Atypia", "P53", "Final_histopathology", 
    "Resolved_histopathology", "Biomarker_status", 
    "time_endo0_to_sponge", "time_endo0_to_endo1", "time_sponge_to_endo1", 
    "time_sponge_to_resolved", "event",
    "SurvStrat", "SurvStratAnyAtypia", "Any_dysplasia",
    "Risk_strat_ID", "AI_calls", "AI_comments"
)
df_clean[, cleanup_colorder]
SeparateCoConsent <- function(Co_consented, ...){
    joined <- toupper(Co_consented)
    split <- strsplit(joined, " ")[[1]]
    AHM <- NA 
    BEST <- NA
    for(string in split){
        string <- gsub("(?<=AHM)\\s", "", string, perl=T)
        if(grepl("AHM|AD", string)){
            AHM <- trimws(string)
        }else if(grepl("BEST", string)){
            BEST <- trimws(string)
        } 
    }
    tibble::tibble(AHM, BEST)
}
pmap_dfr(df_clean, SeparateCoConsent) %>% bind_cols(., df) %>% View
fwrite(df_clean, "Data/DELTA_NHSE_clean.csv", na = "")
# ================= Rate Calculations ====================================
GetRates <- function(df, primary=FALSE){
    if(primary){
        outcomes <- c(
            "OAC", 
            "HGD", 
            "IMC"
        )
    }else{
        outcomes <- c(
            "OAC", 
            "HGD", 
            "IMC", 
            "LGD", 
            "Crypt dysplasia"
        )
    }

    n_total <- nrow(df)
    time_total <- sum(df$time_to_E0E1, na.rm=T)
    n_outcomes <- sum(df$General_histopathology %in% outcomes)
    prevalence = n_outcomes/n_total 
    ci = prop.test(n_outcomes, n_total)$conf.int[1:2]
    ci_low = ci[1] 
    ci_hi = ci[2]

    const <- 1200 # convert events/month to events/100-years
    incidence <- (n_outcomes/time_total)*const
    ci <- poisson.test(n_outcomes, time_total)$conf.int[1:2]
    ci_low2 <- ci[1]*const
    ci_hi2 <- ci[2]*const
     

    return(c(n_total, n_outcomes, prevalence, ci_low, ci_hi, incidence, ci_low2, ci_hi2))
}

# Any dysplasia 
rate_df <- rbind(
    GetRates(df), 
    GetRates(df[df$Clin_risk_strat == "Low risk", ]),
    GetRates(df[df$Clin_risk_strat == "Moderate risk", ]),
    GetRates(df[df$Sponge_risk_strat == "Low risk", ]),
    GetRates(df[df$Sponge_risk_strat == "Moderate risk", ]),
    GetRates(df[df$Sponge_risk_strat == "High risk", ]), 
    GetRates(df[df$SurvStrat == "Low conf positive", ]), 
    GetRates(df[df$SurvStrat == "High conf positive", ])
)
col <- c("N_group", "N_outcome", "Prev", "Prev_ci_low", "Prev_ci_hi", "Inc", "Inc_ci_low", "Inc_ci_hi")
colnames(rate_df) <- col

#HGD+
rate_df <- rbind(
    GetRates(df, primary=TRUE),
    GetRates(df[df$Clin_risk_strat == "Low risk", ], primary=TRUE),
    GetRates(df[df$Clin_risk_strat == "Moderate risk", ], primary=TRUE),
    GetRates(df[df$Sponge_risk_strat == "Low risk", ], primary=TRUE),
    GetRates(df[df$Sponge_risk_strat == "Moderate risk", ], primary=TRUE),
    GetRates(df[df$Sponge_risk_strat == "High risk", ], primary=TRUE), 
    GetRates(df[df$SurvStrat == "Low conf positive", ], primary=TRUE), 
    GetRates(df[df$SurvStrat == "High conf positive", ], primary=TRUE)
)
col <- c("N_group", "N_outcome", "Prev", "Prev_ci_low", "Prev_ci_hi", "Inc", "Inc_ci_low", "Inc_ci_hi")
colnames(rate_df) <- col

year_to_E1_df <- rbind(
    GetRates(df[df$Sponge_risk_strat == "Moderate risk" & df$year_to_E1 == T,]),
    GetRates(df[df$Sponge_risk_strat == "Moderate risk" & df$year_to_E1 == F,])
)
col <- c("N_group", "N_outcome", "Prev", "Prev_ci_low", "Prev_ci_hi", "Inc", "Inc_ci_low", "Inc_ci_hi")
colnames(year_to_E1_df) <- col
View(year_to_E1_df)

#================ Patient Demographics ========================================
GetDemographics <- function(df){
    Age <- median(df$Age, na.rm=T)
    Age_sd <- sd(df$Age, na.rm=T)
    Age_min <- quantile(df$Age, 1/4, na.rm=T)
    Age_max <- quantile(df$Age, 3/4, na.rm=T)

    M_freq <- sum(df$Sex == "Male", na.rm=T)
    M_prop <- M_freq/nrow(df)
    F_freq <- sum(df$Sex == "Female", na.rm=T)
    F_prop <- F_freq/nrow(df) 

    C <- median(df$Prague_C0, na.rm=T) 
    C_IQR_upper <- quantile(df$Prague_C0, 3/4)
    C_IQR_lower <- quantile(df$Prague_C0, 1/4)
    # C_min <- min(df$Prague_C0, na.rm=T)
    # C_max <- max(df$Prague_C0, na.rm=T)

    M <- median(df$Prague_M0, na.rm=T) 
    # M_min<- min(df$Prague_M0, na.rm=T) 
    # M_max <- max(df$Prague_M0, na.rm=T)
    M_IQR_upper <- quantile(df$Prague_M0, 3/4)
    M_IQR_lower <- quantile(df$Prague_M0, 1/4)

    NDBE_freq <- sum(df$General_histopathology == "BE", na.rm=T)
    NDBE_prop <- NDBE_freq/nrow(df)
    IND_freq <- sum(df$General_histopathology == "IND", na.rm=T)
    IND_prop <- IND_freq/nrow(df)
    Crypt_freq <- sum(df$General_histopathology == "Crypt dysplasia", na.rm=T)
    Crypt_prop <- Crypt_freq/nrow(df)
    LGD_freq <- sum(df$General_histopathology =="LGD", na.rm=T)
    LGD_prop <- LGD_freq/nrow(df)
    HGD_freq <- sum(df$General_histopathology %in% c("HGD", "IMC"), na.rm=T)
    HGD_prop <- HGD_freq/nrow(df)
    OAC_freq <- sum(df$General_histopathology == "OAC", na.rm=T)
    OAC_prop <- OAC_freq/nrow(df)

    time_to_mean <- median(df$time_to_spongeE1, na.rm=TRUE)
    time_to_min <- quantile(df$time_to_spongeE1, 1/4)
    time_to_max <- quantile(df$time_to_spongeE1, 3/4)

    time_to_sponge_median <- median(df$time_to_E0Sponge, na.rm=TRUE)
    time_to_sponge_min <- quantile(df$time_to_E0Sponge, 1/4, na.rm=T)
    time_to_sponge_max <- quantile(df$time_to_E0Sponge, 3/4, na.rm=T)

    c(
        Age, 
        Age_min, 
        Age_max, 
        M_freq, 
        M_prop, 
        F_freq, 
        F_prop, 
        C, 
        C_IQR_lower, 
        C_IQR_upper, 
        M, 
        M_IQR_lower, 
        M_IQR_upper, 
        NDBE_freq, 
        NDBE_prop, 
        IND_freq, 
        IND_prop, 
        Crypt_freq, 
        Crypt_prop, 
        LGD_freq, 
        LGD_prop,
        HGD_freq, 
        HGD_prop, 
        OAC_freq, 
        OAC_prop,
        time_to_mean, 
        time_to_min, 
        time_to_max,
        time_to_sponge_median, 
        time_to_sponge_min, 
        time_to_sponge_max
    )
}
demographics <- data.frame(
    rbind(
        GetDemographics(df[df$Sponge_risk_strat == "Low risk", ]), 
        GetDemographics(df[df$Sponge_risk_strat == "Moderate risk", ]), 
        GetDemographics(df[df$Sponge_risk_strat == "High risk", ])
    )
) %>% round(digits = 3)
colnames(demographics) <- c(
        "Age", 
        "Age_IQR_low", 
        "Age_IQR_high",
        "Male_freq", 
        "Male_prop", 
        "Female_freq", 
        "Female_prop", 
        "C_median", 
        "C_IQR_low", 
        "C_IQR_high", 
        "M_median", 
        "M_IQR_low", 
        "M_IQR_high", 
        "NDBE_freq", 
        "NDBE_prop", 
        "IND_freq", 
        "IND_prop", 
        "Crypt_freq", 
        "Crypt_prop", 
        "LGD_freq", 
        "LGD_prop", 
        "HGD_freq", 
        "HGD_prop", 
        "OAC_freq", 
        "OAC_prop",
        "Time_to_endo1_mean", 
        "Time_to_endo1_min", 
        "Time_to_endo1_max", 
        "Time_to_sponge_median", 
        "Time_to_sponge_IQR_low", 
        "Time_to_sponge_IQR_high"
)
demographics <- demographics %>% t()
# Krukal wallis test followed by Conover-Iman 
ggplot(aes(x = Sponge_risk_strat, y = time_to_spongeE1), data = df)+ 
    geom_violin()
ggsave("time_to_spongeE1.png", width = 50, limitsize = F)

a <- kruskal.test(time_to ~ Sponge_risk_strat, data = df)
conover.test(df$time_to_E0E1 - df$time_to_sponge, df$Sponge_risk_strat, method = "bh")
kruskal.test(Age ~ Sponge_risk_strat, data = df)
conover.test(df$Age, df$Sponge_risk_strat, method = "bh")
kruskal.test(Prague_C0 ~ Sponge_risk_strat, data = df) 
conover.test(df$Prague_C0, df$Sponge_risk_strat, method = 'bh')
kruskal.test(Prague_M0 ~ Sponge_risk_strat, data = df) 
conover.test(df$Prague_M0, df$Sponge_risk_strat, method = 'bh')
# chi-squared test 
chisq_table <- xtabs(~Sex + Sponge_risk_strat, data = df)
chisq.test(chisq_table, correct=FALSE)
# PPV's 
dysplasia <- c("Crypt dysplasia", "LGD", "OAC", "HGD", "IMC")
dysplasia <- c("HGD", "OAC", 'IMC')
RG <- df[df$Clin_risk_strat == "Low risk", ]
n_nondys <- sum(!(RG$Histopathology_endo1 %in% dysplasia)) 
PPV1 <- n_nondys/nrow(RG)
prop.test(n_nondys, nrow(RG))
n_dys <- sum((RG$General_histopathology %in% dysplasia)) 
PPV2 <- n_dys/nrow(RG)
prop.test(n_dys, nrow(RG))
RG <- df[df$Sponge_risk_strat == "Moderate risk", ]
n_nondys <- sum(!(RG$General_histopathology  %in% dysplasia)) 
n_dys <- sum((RG$General_histopathology  %in% dysplasia)) 
PPV1 <- n_nondys/nrow(RG)
prop.test(n_nondys, nrow(RG))
PPV2 <- n_dys/nrow(RG)
prop.test(n_dys, nrow(RG))
RG <- df[df$Sponge_risk_strat == "High risk", ]
n_dys <- sum((RG$General_histopathology  %in% dysplasia)) 
PPV1 <- n_dys/nrow(RG)
prop.test(n_dys, nrow(RG))
# senstivity & specificity 
n_dys <- sum(df$General_histopathology %in% dysplasia)
pred_dys <- sum(df$General_histopathology %in% dysplasia & df$Sponge_risk_strat %in% c("Moderate risk", "High risk"))
sens <- pred_dys/n_dys 
prop.test(pred_dys, n_dys)
n_ndbe <- sum(df$General_histopathology %in% c("NDBE", "BE"))
pred_ndbe <- sum(df$General_histopathology %in% c("NDBE", "BE") & df$Sponge_risk_strat == "Low risk")
spec <- pred_ndbe/n_ndbe
prop.test(pred_ndbe, n_ndbe)

bar_df <- df 
bar_df$new_hist <- sapply(bar_df$Histopathology_endo1, function(x){
    if(x %in% c("HGD", "IMC")){
        return("HGD, IMC")
    }else{
        return(x)
    }
})
bar_1 <- bar_df %>% group_by(new_hist, Sponge_risk_strat) %>% count %>% rename(risk_strat = Sponge_risk_strat)
bar_1 <- bar_1 %>% 
    group_by(risk_strat) %>% 
    mutate(risk_strat_label = paste0(risk_strat, "\n n=", sum(n))) 
bar_1$Strat <- "Capsule Sponge & Clinical Risk Groups"
bar_2 <- bar_df %>% group_by(new_hist, Clin_risk_strat) %>% count %>% rename(risk_strat = Clin_risk_strat)
bar_2 <- bar_2 %>% 
    group_by(risk_strat) %>% 
    mutate(risk_strat_label = paste0(risk_strat, "\n n=", sum(n))) 
bar_2$Strat <- "Clinical Risk Groups"
bar_3 <- bar_df %>% group_by(new_hist, BSG_risk_strat) %>% count %>% rename(risk_strat = BSG_risk_strat)
bar_3 <- bar_3 %>% 
    group_by(risk_strat) %>% 
    mutate(risk_strat_label = paste0(risk_strat, "\n n=", sum(n))) 
bar_3$Strat <- "BSG Stratification"
barplot_df <- rbind(bar_1, bar_2)
barplot_df$Histopathology_endo1 <- barplot_df$new_hist %>% factor(levels = c(
        "BE", 
        "IND",
        "Crypt dysplasia",
        "LGD/Crypt dysplasia",  
        "LGD", 
        "HGD",
        "HGD, IMC",
        "IMC",
        "OAC" 
    ) %>% rev 
)
barplot_df$risk_strat <- barplot_df$risk_strat %>% factor(levels = c(
    "Low risk", 
    "Moderate risk", 
    "High risk"
))
barplot_df$risk_strat_label <- barplot_df$risk_strat_label %>% factor(levels = c(
    c(unique(barplot_df$risk_strat_label)[-1], unique(barplot_df$risk_strat_label[1]))
))
bars <- ggplot(aes(x = risk_strat_label, fill = Histopathology_endo1, y = n), data = barplot_df) + 
    theme_minimal() + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(name = "Histopathology at \npost-sponge endoscopy", values = progression_colours, labels = c(
        "OAC", 
        "HGD, IMC", 
        "LGD",
        "Crypt dysplasia", 
        "Indefinite", 
        "NDBO"
    )) +
    facet_grid(~Strat %>% factor(levels = c("Clinical Risk Groups", "Capsule Sponge & Clinical Risk Groups")), scales = "free", space = "free")+ 
    ylab("Number of patients") + 
    xlab("Risk Groups") + 
    theme(
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16), 
        strip.text = element_text(size = 16), 
        legend.text = element_text(size = 19), 
        legend.title = element_text(size = 20), 
        panel.spacing = unit(6, "lines")) 
print(bars) 
ggsave("Plots/demographic_barplots_stack.pdf", bars, dpi = 600, width = 15, height = 9.6)
#================= Sankey Diagram ==============================================
sankey_df <- data.frame()
clin_status <- unique(df$Clin_risk_strat) %>% factor(
    levels = c("Moderate risk", "Low risk")
)
biomark_status <- unique(df$Biomark_status)
sponge_status <- unique(df$Sponge_risk_strat)
final_status <- unique(df$Histopathology_endo1)

for(clin in clin_status){
    for(biomark in biomark_status){
        for(sponge in sponge_status){
            for(final in final_status){
                highlight <- 0 
                # if(clin == "Low risk" & biomark != "Negative"){
                #     highlight <- "A"
                # }else if(
                #     sponge == "Low risk" & 
                #     final %in% c("IND", "LGD", "HGD", "IMC", "OAC")){
                #         highlight <- "B"
                # }else if(
                #     clin == "Moderate risk" & 
                #     biomark != "Negative"){
                #         highlight <- "C"
                # }

                freq <- sum(
                    df$Clin_risk_strat == clin & 
                    df$Biomark_status == biomark & 
                    df$Sponge_risk_strat == sponge & 
                    df$Histopathology_endo1 == final
                )
                if(final %in% c("HGD", "IMC", "OAC")){
                    final <- "HGD, IMC, OAC"
                }else if(final %in% c("Crypt dysplasia", "LGD")){
                    final <- "LGD/Crypt dysplasia"
                }
                row <- c(
                    "Clinical Risk Group" = clin, 
                    "Biomarker Status" = biomark, 
                    "Sponge Risk Group" = sponge, 
                    "Post-Sponge Endoscopy" = final, 
                    "Highlight" = highlight, 
                    "Freq" = freq
                )
                sankey_df <- rbind(sankey_df, data.frame(t(row)))
            }
        }
    }
}

# sankey_df <- sankey_df %>% filter(Freq > 0)
sankey_df$Clinical.Risk.Group <- sankey_df$Clinical.Risk.Group %>% factor(levels=c(
    "High risk", 
    "Moderate risk", 
    "Low risk"))
sankey_df$Biomarker.Status <- sankey_df$Biomarker.Status %>% factor(levels = c(
    "Negative", 
    "Low Confidence Positive", 
    "High Confidence Positive"
) %>% rev)
sankey_df$Sponge.Risk.Group <- sankey_df$Sponge.Risk.Group %>% factor(levels = c(
    "Low risk", 
    "Moderate risk", 
    "High risk"
) %>% rev )
sankey_df$Post.Sponge.Endoscopy <- sankey_df$Post.Sponge.Endoscopy %>% factor(levels = c(
    "BE", 
    "Crypt dysplasia", 
    "IND", 
    "LGD/Crypt dysplasia", 
    "HGD, IMC, OAC"
) %>% rev)
sankey_df$Freq <- as.numeric(sankey_df$Freq)

ggplot( data = sankey_df, 
        aes(
            axis1 = Clinical.Risk.Group, 
            axis2 = Biomarker.Status, 
            axis3 = Sponge.Risk.Group, 
            axis4 = Post.Sponge.Endoscopy, 
            y = Freq))+
    theme_minimal() +
    scale_x_discrete(
        limits = c(
            "Clinical Risk", 
            "Biomarker Status", 
            "Sponge Risk", 
            "Histopathology at Post-sponge Endoscopy") %>% str_wrap(width=15),
            expand = c(.2, .05)
    ) +
    geom_flow(aes(fill=Highlight)) + 
    geom_stratum(aes(fill = after_stat(stratum)), color=alpha("black", 0.3))+ 
    geom_text(
        stat = "stratum", 
        aes(label = str_wrap(after_stat(stratum), width = 8)), 
        size = 3
    ) + 
    scale_fill_manual(values = progression_colours) + 
    ylab("Number of Patients")

#initialise lode-form data frame
lode_df <- to_lodes_form(sankey_df, axes=1:4)
# highlight partial paths in lode-form df 
for(alluvial in lode_df$alluvium){
    flow_rows <- which(lode_df$alluvium == alluvial)
    flow_df <- lode_df[lode_df$alluvium == alluvial, ]
    # if(flow_df$stratum[3] == "Low risk" && 
    #     flow_df$stratum[4] %in% c("HGD, IMC, OAC", "LGD/Crypt dysplasia", "IND")){ 
    #         # highlight if sponge risk is low and goes to dysplasia
    #         lode_df$Highlight[flow_rows[-c(1,2)]] <- "A"
    # }else if(flow_df$stratum[3] == "Moderate risk" && 
    #         flow_df$stratum[4] %in% c("HGD, IMC, OAC", "LGD/Crypt dysplasia", "IND")){ 
    #             #highlight if sponge risk is mod and goes to dysplasia 
    #             lode_df$Highlight[flow_rows[-c(1,2)]] <- "B"
    # }
    if(flow_df$stratum[2] == "High Confidence Positive"){
        lode_df$Highlight[flow_rows[-1]] <- "1"
    }else if(flow_df$stratum[2] == "Low Confidence Positive"){
        lode_df$Highlight[flow_rows[-1]] <- "2"
    }else if(flow_df$stratum[3] == "Low risk" ){
        lode_df$Highlight[flow_rows[-c(1,2)]] <- "3"
    }
}
lode_plot <- ggplot(
    data = lode_df %>% filter(Freq > 0), 
    aes(x = x, 
        stratum=stratum, 
        alluvium=alluvium, 
        fill=Highlight, 
        y = Freq))+ 
    theme_minimal()+
    geom_flow(stat="alluvium", color=alpha("darkgrey", 0.1))+ 
    scale_fill_manual(values = progression_colours, guide="none") +
    ggnewscale::new_scale_fill() + 
    geom_stratum(
        data = lode_df %>% filter(Freq>0),
        aes(fill = stratum), 
        color = alpha("darkgrey", alpha=0.9)) + 
    geom_text(
        stat = "stratum", 
        aes(
            label = ifelse(after_stat(prop) > 0.1, paste0(round(100*after_stat(prop), digits = 1) %>% format(nsmall=1), "%"), ""), 
            angle = 90
        ), 
        size = 5,
        fontface = "bold"
    ) +
    geom_text(
        stat = "stratum", 
        aes(
            label = ifelse(after_stat(prop) < 0.1 & round(100*after_stat(prop), digits=1) != 3.0, paste0(round(100*after_stat(prop), digits = 1) %>% format(nsmall=1), "%"), "")
        ), 
        nudge_x = 0.35,
        size = 5, 
        fontface = "bold"
    ) +
    geom_text(
        stat = "stratum", 
        aes(
            label = ifelse(round(100*after_stat(prop), digits=1) == 3.0 , paste0(round(100*after_stat(prop), digits =1)%>% format(nsmall=1), "%"), "")
        ), 
        colour = "white",
        size = 5,
        fontface = "bold"
    ) +
    geom_text(
        stat = "stratum", 
        aes(
            label = ifelse(round(100*after_stat(prop), digits=1) %in% c(62.1, 54.4) , paste0(round(100*after_stat(prop), digits =1)%>% format(nsmall=1), "%"), ""),
            angle=90
        ), 
        colour = "white",
        size = 5,
        fontface = "bold"
    ) +
    scale_fill_manual(
        name = "Strata", 
        limits= c(
            "High risk", 
            "Moderate risk", 
            "Low risk", 
            " ", 
            "High Confidence Positive", 
            "Low Confidence Positive", 
            "Negative", 
            "  ", 
            "HGD, IMC, OAC", 
            "LGD/Crypt dysplasia", 
            "IND", 
            "BE"
        ),
        labels = c(
            "High risk", 
            "Moderate risk", 
            "Low risk", 
            " ", 
            "Atypia and aberrant p53", 
            "Other positive", 
            "Negative", 
            "  ", 
            "HGD, IMC, OAC", 
            "LGD, Crypt dysplasia", 
            "Indefinite", 
            "NDBO"
        ),
        values = progression_colours)+
    scale_x_discrete(
        labels = c(
            "Clinical Risk Group", 
            "Biomarker Status",
            "Capsule Sponge  & Clinical Risk Group",
            "Histopathology at Post Sponge Endoscopy") %>% str_wrap(width = 15), 
        position = "top")+
    theme(
        legend.position="right", 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 16),
        axis.text=element_text(size = 16), 
        legend.text = element_text(size=19), 
        legend.title = element_text(size=20)) + 
    ylab("Number of patients")
print(lode_plot)
ggsave("Plots/Lode_plot_alt_col.pdf", lode_plot, dpi=600, width=12.7, height=9.6)
#================= Patient Timelines ===========================================
tl_subset <- df %>% filter(Any_dysplasia == "T")
#tl_subset <- df %>% filter(Final_histopathology %in% c("HGD", "IMC", "LGD", "OAC", "IND", "Crypt dysplasia"))
BlowUp <- function(
    ID,
    CYT_ID,
    Date_sponge, 
    Date_endo1, 
    Histopathology_endo1, 
    Date_endo2, 
    Histopathology_endo2, 
    Date_endo3, 
    Histopathology_endo3, 
    Date_endo4, 
    Histopathology_endo4, 
    Sponge_risk_strat, 
    ...){
        tibb <- tibble()

        baseline <- as.Date(Date_sponge, format="%d/%m/%Y")
        endo_dates <- c(
            Date_endo1, 
            Date_endo2, 
            Date_endo3, 
            Date_endo4
        )
        endo_paths <- c(
            Histopathology_endo1, 
            Histopathology_endo2, 
            Histopathology_endo3, 
            Histopathology_endo4
        )
        final_path <- NA
        for(index in seq_along(endo_paths) %>% rev){
            path <- endo_paths[index]
            if(is.na(path)){
                next
            }else{
                if(path %in% c("HGD", "IMC")){
                    path <- "HGD, IMC"
                }
                if(is.na(final_path)){
                    final_path <- path
                }
                date <- endo_dates[index] %>% as.Date(format="%d/%m/%Y")
                time <- (date - baseline) %>% as.numeric
                time <- time * 12/365.25
                row <- tibble(
                    ID = ID, 
                    Time_from_sponge = time, 
                    Histopathology = path, 
                    Final_histopathology = final_path,
                    Sponge_risk_strat = Sponge_risk_strat
                )
                tibb <- rbind(tibb, row)
            }
        }
        row <- tibble(
            ID = ID, 
            Time_from_sponge = 0, 
            Histopathology = NA, 
            Final_histopathology = final_path,
            Sponge_risk_strat = Sponge_risk_strat
        )
        tibb <- rbind(tibb, row)

        return(tibb)
}
tl_df <- pmap_dfr(tl_subset, BlowUp)
tl_df$Sponge_risk_strat <- tl_df$Sponge_risk_strat %>% factor(levels = c(
    "High risk", 
    "Moderate risk", 
    "Low risk"
) %>% rev)
tl_df <- tl_df %>% 
    group_by(ID) %>% 
    mutate(time_to_end = max(Time_from_sponge)) %>% 
    ungroup(ID)
tl_df$ID <- factor(
    tl_df$ID, 
    levels = tl_df$ID[order(tl_df$time_to_end)] %>% unique
)
tl_df$Mask = "F"
tl_df$Histopathology <- sapply(tl_df$Histopathology, function(x) ifelse(x == "NDBE", "BE", x))
tl_df$Final_histopathology <- sapply(tl_df$Final_histopathology, function(x) ifelse(x == "NDBE", "BE", x))

n_hr <- nrow(tl_df %>% filter(Sponge_risk_strat == 'High risk') %>% count(ID))
n_mr <- nrow(tl_df %>% filter(Sponge_risk_strat == "Moderate risk")%>% count(ID)) 
n_lr <- nrow(tl_df %>% filter(Sponge_risk_strat == "Low risk")%>% count(ID)) 
lr_new_rows <- data.frame(
    ID = sapply(1:(n_hr - n_lr), function(x){
        randID <- paste0(sample(LETTERS, 10), collapse="")
        paste0("*", randID)
    }), 
    Time_from_sponge = 0, 
    Histopathology = NA, 
    Final_histopathology = NA, 
    Sponge_risk_strat = "Low risk",
    time_to_end = 10000, 
    Mask = "T"
)
mr_new_rows <- data.frame(
    ID = sapply(1:(n_hr - n_mr), function(x){
    randID <- paste0(sample(LETTERS, 10), collapse="")
    paste0("*", randID)
    }),  
    Time_from_sponge = 0, 
    Histopathology = NA, 
    Final_histopathology = NA, 
    Sponge_risk_strat = "Moderate risk",
    time_to_end = 10000, 
    Mask = "T"
)
tl_df <- rbind(tl_df, rbind(lr_new_rows, mr_new_rows))
# plot_df <- data.frame()
# for(risk in unique(tl_df$Sponge_risk_strat)){
#     new_rows <- tl_df %>% 
#         filter(Sponge_risk_strat == risk)%>% 
#         arrange(time_to_end, ID) %>%  
#         group_by(time_to_end, DELTA_ID) %>%
#         mutate(ID = cur_group_id()) %>% 
#         ungroup() 
#     if(risk == "Low risk"){
#         r <- "LR"
#     }else if(risk == "Moderate risk"){
#         r <- "MR"
#     }else{
#         r <- "HR"
#     }
#     new_rows$ID <- sapply(new_rows$ID, function(x) paste0(r, x))

#     plot_df <- rbind(plot_df, new_rows)
# }
# plot_df$ID <- plot_df$ID %>% factor(levels = 1:max(plot_df$ID)%>% as.character())
plot_df <- tl_df
plot_df$Histopathology <- plot_df$Histopathology %>% factor(levels = c(
    "OAC", 
    "HGD, IMC", 
    "LGD", 
    "IND", 
    "Crypt dysplasia", 
    "BE"
))
lr_names <- plot_df %>% filter(Sponge_risk_strat == "Low risk") %>% arrange(time_to_end) %>% pull(ID) %>% as.character() %>% unique()
lr_names <- sapply(lr_names, function(x) ifelse(grepl("\\*", x), "", x)) %>% unname
mr_names <- plot_df %>% filter(Sponge_risk_strat == "Moderate risk") %>% arrange(time_to_end) %>% pull(ID) %>% as.character() %>% unique()
mr_names <- sapply(mr_names, function(x) ifelse(grepl("\\*", x), "", x)) %>% unname
hr_names <- plot_df %>% filter(Sponge_risk_strat == "High risk") %>% arrange(time_to_end) %>% pull(ID) %>% as.character() %>% unique()
hr_names <- sapply(hr_names, function(x) ifelse(grepl("\\*", x), "", x)) %>% unname

lr_plot <- plot_df %>% filter(Sponge_risk_strat == "Low risk") %>%  
    ggplot(
        aes(
            x = Time_from_sponge, 
            y = ID,
            group = ID)) +
        theme_minimal() +
        geom_line(aes(col = Final_histopathology), lwd=2) +
        scale_colour_manual(name = "Final Endoscopy Status", values = progression_colours, guide= 'none') +
        ggnewscale::new_scale_colour() + 
        geom_point(aes(col=Histopathology), 
                    pch=19, cex=4, alpha=0.9) +
        scale_colour_manual(name = "Endoscopy Status", values = progression_colours, na.translate=FALSE, guide='none') +
        ggnewscale::new_scale_colour() + 
        geom_point(aes(col = Mask), alpha = 0) +
        scale_colour_manual(values = c("white", "white"), guide = "none") + 
        guides(fill=guide_legend(override.aes=list(shape=22))) + 
        xlab("Months Since Sponge") + 
        ylab("Patient ID") + 
        geom_vline(xintercept = 12, lty = 2, alpha = 0.6) + 
        scale_x_continuous(breaks = c(0, 12, 24, 36)) + 
        scale_y_discrete(labels = lr_names %>% rev, limits = rev) +
        ggtitle("Low Risk") +
        theme(
            legend.position="right", 
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12), 
            strip.text = element_text(size = 12),
            axis.text.y = element_text(size = 8)
        )
mr_plot <- plot_df %>% filter(Sponge_risk_strat == "Moderate risk") %>%  
    ggplot(
        aes(
            x = Time_from_sponge, 
            y = ID,
            group = ID)) +
        theme_minimal() +
        geom_line(aes(col = Final_histopathology), lwd=2) +
        scale_colour_manual(name = "Final Endoscopy Status", values = progression_colours, guide= 'none') +
        ggnewscale::new_scale_colour() + 
        geom_point(aes(col=Histopathology), 
                    pch=19, cex=4, alpha=0.9) +
        scale_colour_manual(name = "Endoscopy Status", values = progression_colours, na.translate=FALSE, guide='none') +
        ggnewscale::new_scale_colour() + 
        geom_point(aes(col = Mask), alpha = 0) +
        scale_colour_manual(values = c("white", "white"), guide = "none") + 
        guides(fill=guide_legend(override.aes=list(shape=22))) + 
        xlab("Months since capsule sponge") + 
        ylab("Patient ID") + 
        geom_vline(xintercept = 12, lty = 2, alpha = 0.6) + 
        scale_x_continuous(breaks = c(0, 12, 24, 36)) + 
        scale_y_discrete(labels = mr_names %>% rev, limits = rev) +
        ggtitle("Moderate Risk") +
        theme(
            legend.position="right", 
            axis.title.x = element_text(vjust=-2), 
            axis.title.y = element_blank() ,
            strip.text = element_text(size = 12),
            axis.text.y = element_text(size = 8)
        )
hr_plot <- plot_df %>% filter(Sponge_risk_strat == "High risk") %>%  
    ggplot(
        aes(
            x = Time_from_sponge, 
            y = ID,
            group = ID)) +
        theme_minimal(base_size=10) +
        geom_line(aes(col = Final_histopathology), lwd=2) +
        scale_colour_manual(name = "Final Endoscopy Status", values = progression_colours, guide= 'none') +
        ggnewscale::new_scale_colour() + 
        geom_point(aes(col=Histopathology), 
                    pch=19, cex=4, alpha=0.9, show.legend=TRUE) +
        scale_colour_manual(
            name = "Endoscopy Status", 
            values = progression_colours, 
            labels = c(
                "OAC", 
                "HGD, IMC",
                "LGD", 
                "Indefinite", 
                "Crypt dysplasia", 
                "NDBO"
            ),
            na.translate=FALSE, 
            drop=FALSE, 
            limits = unique(plot_df$Histopathology[!is.na(plot_df$Histopathology)]) %>% sort) +
        ggnewscale::new_scale_colour() + 
        geom_point(aes(col = Mask), alpha = 0) +
        scale_colour_manual(values = c("white", "white"), guide = "none") + 
        xlab("") + 
        ylab("Patient ID") + 
        geom_vline(xintercept = 12, lty = 2, alpha = 0.6) + 
        scale_x_continuous(breaks = c(0, 12, 24, 36)) + 
        scale_y_discrete(labels = hr_names %>% rev, limits = rev) +
        ggtitle("High Risk") +
        theme(
            legend.position="right", 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            strip.text = element_text(size = 12),
            axis.text.y = element_text(size = 8)
        )
timeline_plot <- lr_plot + mr_plot + hr_plot
plot(timeline_plot)
ggsave("Plots/dysplasia_timelines_any_dys.jpg", timeline_plot, dpi = 600, width=12.7, height=9.6)

#================= Survival Analysis ===========================================
kf$Age <- as.double(kf$Age)
surv_df <- left_join(df, kf, by = join_by(
    Age == Age, 
    Sex == Gender, 
    DELTA_ID == Study_Number,
    Date_endo0 == Date_of_Endo0
))
KM <- survfit(Surv(time_to, event) ~ SurvStrat, data = df)
cox <- coxph(Surv(time_to, event) ~ SurvStrat, data = df) # full cox regression 

hx_df <- inner_join(kf2, df,
    by = join_by(
        Study_number == DELTA_ID, 
        Age == Age, 
        Gender == Sex, 
        `Date_of_follow-up_endoscopy` == Date_endo1,
        C == Prague_C0, 
        M == Prague_M0, 
        Date_of_Cytosponge == Date_sponge
    )
)

KM_clin <- survfit(Surv(time_to, event)~Clin_risk_strat, data = df)
cox_clin <- coxph(Surv(time_to, event) ~ Clin_risk_strat, data = df)
cox_hx <- coxph(Surv(time_to, event) ~ SurvStrat + Hx_dysplasia + SurvStrat:Hx_dysplasia, data = surv_df) # subset cox post-hoc subgroup
# Check proportional hazards assumptions 
df_prop_check <- surv_df %>% 
    mutate(
        SurvStratMod = model.matrix(cox)[,"SurvStratModerate risk"], 
        SurvStratLowConf = model.matrix(cox)[,"SurvStratLow conf positive"],
        SurvStratHiConf = model.matrix(cox)[,"SurvStratHigh conf positive"], 
        time_to_SpE1 = time_to_E0E1 - time_to
    )
cox_prop_check <- coxph(Surv(time_to, event) ~ SurvStratMod + SurvStratLowConf + SurvStratHiConf, data = df_prop_check)
check <- cox.zph(cox_prop_check) # failed, constant betas not met 

KM_to_e1 <- survfit(Surv(time_to_spongeE1, rep(1, 910)) ~ SurvStrat, data = df)
df$endo_event <- 1
KM_to_e1<- survfit(Surv(time_to_spongeE1, endo_event)~SurvStrat, data =df)
# generate KM plot
progression_colours <- c(
    "Low risk" = "#525252", 
    "Moderate risk" = "#57d35b", 
    "High risk" = "#18a96d", 
    "High Confidence Positive" = "#a5149a",#8365c8"
    "Low Confidence Positive" = "#e693dc" #c453ba", 
)
survp <- ggsurvplot(
    KM,
        data = df,#[df$SurvStrat !="NY",],
        conf.int = TRUE,
        ggtheme = theme_bw(),
         legend.labs=c(
            "Low risk",
            "Moderate risk",
            "High risk tier 2", 
            "High risk tier 1"
        ),
        palette = c( progression_colours[["Low risk"]], progression_colours[["Moderate risk"]], progression_colours[["Low Confidence Positive"]], progression_colours[["High Confidence Positive"]]), 
        risk.table = TRUE, 
        pval=TRUE,
        pval.coord = c(-1,-0.05), 
        pval.method.coord = c(5,0),
        log.rank.weights = "n",
        ylab = "Proportion with endoscopy", 
        xlab = "Months", 
        xlim = c(-1,36), 
        ylim = c(-0.05, 1),
        break.x.by= 12, 
        fontsize=20,
        tables.y.text.col =F,
        risk.table.fontsize=8
) 
survp$plot <- survp$plot+
    #scale_y_continuous(labels = c("1", "0.75", "0.5", "0.25", "0")) + # modifies y-axis to be [1,0] rather than [0,1]
    theme(
        legend.title=element_blank(), 
        legend.text=element_text(size =18), 
        axis.text =element_text(size=20), 
        axis.title = element_text(size = 20))
survp$table <- survp$table + 
    scale_y_discrete(
        labels = c(
            "Low risk", 
            "Moderate risk",
            "High risk tier 2", 
            "High risk tier 1"
        ) %>%rev,
        position="right"
    ) + 
    theme(
        axis.title.x=element_text(size=18), 
        axis.title.y =element_blank(), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=18),
        plot.title = element_blank()
    )
layout <- c(
  area(t = 0, l = 0, b = 20, r = 15), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 22, l = 0, b = 26, r = 15) # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3
)
surv_plot <- survp$plot/survp$table + plot_layout(design=layout)
print(surv_plot)
ggsave("Plots/survplot_reviwed_col2.png", surv_plot, dpi=600, height = 12.7, width = 15)
# generate combined forest plot 
dysplasia <- c("Crypt dysplasia", "HGD", "LGD", "OAC", "IMC")
dysplasia <- c("HGD", "IMC", "OAC")
risk_table <- table(
    df$SurvStrat, 
    df$General_histopathology %in% dysplasia )
forest_df <- riskratio(risk_table)$measure[-1,] %>% data.frame() %>% round(digits = 1)
forest_df$pval <- riskratio(risk_table)$p.val[-1,2] # use Fischer's exact test 
forest_df$adj_pval <- ifelse(
    forest_df$pval < 1*10^(-16), 
    paste0("p < 2.2 e-16"), 
    paste0("p = ", scientific(forest_df$pval, digits = 2))
)
forest_df$Group <- c(
    "Moderate risk" ,
    "High risk tier 2" , 
    "High risk tier 1"
    # "History of dysplasia",
    # "Moderate risk : History of dysplasia",
    # "High risk tier 2 : History of dysplasia",
    # "High risk tier 1 : History of dysplasia"
)
colnames(forest_df) <- c("RR", "ci_low", "ci_high", "pval", "adj_pval", "Group")
forest_df <- forest_df %>% mutate(
                            annotation = paste0(round(RR, 2), " ", "(", round(ci_low, 2), ", ", round(ci_high, 2), ")")) 

# plot portion 
forest_df$Group <- forest_df$Group %>% factor(
    levels = c(
        "Risk Stratification",
        "Moderate risk" , 
        "High risk tier 2" , 
        "High risk tier 1", 
        "History of dysplasia",
        "Moderate risk : History of dysplasia",
        "High risk tier 2 : History of dysplasia",
        "High risk tier 1 : History of dysplasia"
    ) %>% rev
)
forest_df$annotation <- forest_df$annotation %>% factor(
    levels = forest_df$annotation %>% rev
)
forest_df <- forest_df[-4,]
p <- ggplot(aes(y = Group), data = forest_df) + 
        theme_minimal() + 
        geom_point(aes(x = log(RR) ), shape = 15, size = 8)+
        geom_linerange(aes(xmin = log(ci_low), xmax = log(ci_high)), linewidth = 2.5) + 
        geom_vline(xintercept = log(1), linetype = "dashed") + 
        labs(x = "Log Relative Risk") + 
        theme(
            axis.line.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            axis.title.x = element_text(size = 20), 
            axis.text.x = element_text(size = 18)
            )
forest_df <- rbind(
    data.frame(
        RR = NA, 
        ci_low = NA, 
        ci_high = NA, 
        Group = "Risk Stratification", 
        pval = "p-value", 
        adj_pval = "Adjusted p-value",
        annotation = "Relative Risk of Any Dysplasia (95% CI)" %>% str_wrap(width = 22)), forest_df
) 
forest_df$Group <- forest_df$Group %>% factor(
    levels = c(
        "Risk Stratification",
        "Moderate risk" , 
        "High risk tier 2" , 
        "High risk tier 1", 
        "History of dysplasia",
        "Moderate risk : History of dysplasia",
        "High risk tier 2 : History of dysplasia",
        "High risk tier 1 : History of dysplasia"
    ) %>% rev
)
p_left <- ggplot(aes(y = Group), data = forest_df)+ 
            geom_text(
                aes(
                    x = 0, 
                    label = Group, 
                    hjust = 0, 
                    fontface = "bold"), 
                    size = 7) + 
            geom_text(
                aes(
                    x = 1.5, 
                    label = annotation, 
                    fontface = ifelse(grepl("95%", annotation), "bold", "plain")), 
                    hjust = 0, 
                    size = 7) + 
            theme_void() + 
            coord_cartesian(xlim = c(0,4.5)) 
layout <- c(
  area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 8, l = 3, b = 30, r = 4) # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3
)
forest_plot <- p_left + p + plot_layout(design = layout)
print(forest_plot)
ggsave("Plots/forest_plot_RR.pdf", forest_plot, dpi = 600, width = 13, height = 4) 

names = c(
    paste0("Low risk \n (n = ", sum(df$SurvStrat == "Low risk"), ")"), 
    paste0("Moderate risk \n (n = ", sum(df$SurvStrat == "Moderate risk"), ")"), 
    paste0("Low confidence \n biomarker positive \n (n = ", sum(df$SurvStrat == "Low conf positive"), ")"), 
    paste0("High confidence \n biomarker positive \n (n = ", sum(df$SurvStrat == "High conf positive"), ")")
) 
forest_df <- data.frame(
    names = names , 
    estimate = c(0, cox$coefficients),
    se = c(0, sqrt(diag(cox$var))), 
    pval = c(NA, summary(cox)$coefficients[,5])
)
forest_plot <- forestplot(
    forest_df, 
    name = names, 
    estimate = estimate, 
    se = se, 
    pvalue = pval, 
    logodds=TRUE, 
    xlab="Hazard Ratio",
) + 
    theme(
        axis.text = element_text(size=13), 
        plot.margin=unit(c(0.1,1,0.1,0.1), "cm")
)
print(forest_plot)
ggsave("Plots/forest_plot.png", forest_plot, dpi=600)
surv_plot + forest_plot + plot_layout(ncol=2)
KM_raw <- survfit(Surv(time_to, event) ~ SurvStratAnyAtypia, data = df)
cox_raw<- coxph(Surv(time_to, event) ~ SurvStratAnyAtypia, data = df)

names = c(
    paste0("Low risk (n = ", sum(df$SurvStratAnyAtypia == "Low risk"), ")"), 
    paste0("Moderate risk (n = ", sum(df$SurvStratAnyAtypia == "Moderate risk"), ")"), 
    paste0("High risk: any atypia (no AUS, n = ", sum(df$SurvStratAnyAtypia == "High Atypia"), ")"), 
    paste0("High risk: no atypia (n = ", sum(df$SurvStratAnyAtypia == "High other"), ")")
)

forest_df <- data.frame(
    names = names, 
    estimate = c(0, cox_raw$coefficients),
    se = c(0, sqrt(diag(cox_raw$var))), 
    pval = c(NA, summary(cox_raw)$coefficients[,5])
)
forestplot(
    forest_df, 
    name = names, 
    estimate = estimate, 
    se = se, 
    pvalue = pval, 
    logodds=TRUE, 
    xlab="Hazard Ratio"
)
#================= Other Analyses ==============================================
reg_df <- df 
reg_df$Primary_outcome <- sapply(reg_df$Histopathology_endo1, function(x){
    ifelse(x %in% c("HGD", "OAC", "IMC"), 1, 0)
})
reg_df$Secondary_outcome <- sapply(reg_df$Histopathology_endo1, function(x){
    ifelse(x %in% c("HGD", "OAC", "IMC", "LGD"), 1, 0)
})
reg_df$Any_dysplasia <- reg_df$Any_dysplasia %>% factor(
    levels = c('F', 'T')
)
reg_df$P53 <- sapply(reg_df$P53, function(x) ifelse(toupper(x) == "E", "E", x)) %>% factor(
    levels = c(
        "N", 
        "E", 
        "Y"
    )
)
reg_df$Atypia <- reg_df$Atypia %>% factor(
    levels= c(
        "N", 
        "AUS", 
        "Y"
    )
)
reg_df$Sex <- reg_df$Sex %>% factor(
    levels = c(
        "Female", "Male"
    )
)
FitUnivariate <- function(predictor) {
    eq <- paste0("Any_dysplasia ~ ", predictor) %>% as.formula()
    model <- glm(eq, family = binomial(link = "logit"), data = reg_df)
    summary <- summary(model)

    n_coefs <- nrow(summary$coefficients)
    coefs <- rownames(summary$coefficients)

    p_values <- c()

    for (coef in coefs) {
        p_values <- c(
            p_values,
            summary$coefficients[coef, "Estimate"],
            summary$coefficients[coef, "Pr(>|z|)"]
        )
    }
    zero_vec <- rep(NA, 6)
    zero_vec[1:(2 * n_coefs)] <- p_values
    return(zero_vec)
}
univariate_df <- sapply(c(
    "Age",
    "Sex",
    "Prague_C0",
    "Prague_M0",
    "Atypia",
    "P53"
), FitUnivariate) %>%
    t() %>%
    data.frame() %>%
    setnames(c(
        "Est1",
        "Pval1",
        "Est2",
        "Pval2",
        "Est3",
        "Pval3"
    ))

Clinical_model <-  glm(
    Secondary_outcome ~ Age + Sex + Prague_C0 + Prague_M0 +
        Age:Prague_C0 + Age:Prague_M0 +
        Prague_C0:Prague_M0 + Age:Sex + Sex:Prague_C0 + Sex:Prague_M0,
    family = binomial("logit"), data = reg_df
)
Clinical_model_2 <- glm(
    Secondary_outcome ~ Age + Sex + Prague_C0 + Prague_M0,
    family = binomial("logit"), data = reg_df
)
lrtest(Clinical_model, Clinical_model_2) # assess performance of two models (no diff)

Biomarker_model <- glm(
    Secondary_outcome ~ Atypia + P53 + Atypia:P53, 
    family = binomial("logit"), data = reg_df
)
Biomarker_model_2 <- glm(
    Secondary_outcome ~ Atypia + P53, 
    family = binomial("logit"), data = reg_df
)
lrtest(Biomarker_model_2, Biomarker_model) # assess performance of two models (no diff)

All_model <- glm(
    Secondary_outcome ~ Atypia + P53 + Age + Sex + Prague_C0 + Prague_M0, 
    family = binomial("logit"), data = reg_df
)
# given no diff between the pairs of models, we choose the simpler ones 

clin_out <- reg_df$Secondary_outcome[
    !is.na(reg_df$Age) & !is.na(reg_df$Sex)]

clinical_roc <- roc(clin_out, predict(Clinical_model_2, type="response"), ci=TRUE)

biomark_roc <- roc(reg_df$Secondary_outcome, predict(Biomarker_model_2, type = "response"), ci=TRUE)

all_roc <- roc(clin_out, predict(All_model, type="response"), ci=TRUE)

roc_curvs <- list(Clinical = clinical_roc, Biomarker = biomark_roc, All = all_roc)
roc_plot <- ggroc(roc_curvs) + 
    theme_minimal() + 
    geom_line(lwd = 1.5, alpha = 0.8)+
    scale_colour_manual(
        values = c(
            Clinical = "#a9d9ab", 
            Biomarker = "#a5149a",
            All = "#ffd43b"
        ), 
        labels = c(
            paste0("Clinical Markers\n AUC = ", clinical_roc$auc %>% as.numeric %>% round(digits = 2), "\n"),
            paste0("Biomarkers \n AUC = ", biomark_roc$auc %>% as.numeric %>% round(digits = 2), "\n"), 
            paste0("Both \n AUC = ", all_roc$auc %>% as.numeric %>% round(digits = 2), "\n")
        )
    ) + 
    guides(col=guide_legend(title = "Risk Stratification \n Methods"))
print(roc_plot)
ggsave("Plots/ROC_curve.jpg", roc_plot, dpi=600 )