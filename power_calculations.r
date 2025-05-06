library(dplyr) 
library(data.table) 
library(rstatix)
library(clinfun)

DELTA_full <- fread("/home/somers01/ECI/OptimalScreening/Data/DELTA_all_agg.csv", na.strings=c(""))
cf_df <- fread("/home/somers01/ECI/OptimalScreening/Data/DELTA_counterfactual.csv",  na.strings=c(""))

filter_num <- function(Study_num){
    out <- gsub('[[:digit:]]+', '', Study_num)
    out <- gsub(" ", "", out)
    return(out)
}

sapply(DELTA_full$Study_Number, filter_num) %>% unique

district <- c(
    "ENH", # lister
    "Bed", # bedford
    "NUH",  # nottingham 
    "CCG", # lister, 
    "PAH", # harlow 
    "WSH", # west suffolk 
    "HAR", #Harrogate 
    "ANT", #Antrom 
    "ALT" #craigavon
)

is_district <- function(Study_num){
    (sapply(district, function(x) grepl(x, Study_num)) %>% sum) == 1
}
DELTA_full$is_district <- sapply(DELTA_full$Study_Number, is_district)
DELTA_district <- DELTA_full %>% filter(is_district)

CreateFinalDiagnosis <- function(row){
    diagnoses <- c(
        row[["Endo3_histopathology"]], 
        row[["Endo2_histopathology"]], 
        row[["Endo1_histopathology"]]
    )
    diagnosis <- NA
    for(i in seq_along(diagnoses)){
        if(!is.na(diagnoses[i])){
            diagnosis <- diagnoses[i]
            break
        }
    }
    return(diagnosis)
}
ConvertShorthand <- function(diagnosis){
    mapping <- c(
        "non-dysplastic barrett's oesophagus" = "NDBE", 
        "ndbe" = "NDBE",
        "mastric metaplasia (no im)" = "GM (no IM)", 
        "indefinite for dysplasia" = "IND", 
        "ind" = "IND", 
        "low grade dysplasia" = "LGD", 
        "crypt dysplasia" = "Crypt Dysplasia", 
        "high grade dysplasia" = "HGD", 
        "lgd" = "LGD", 
        "hgd" = "HGD", 
        "goj im" = "GOJ IM", 
        "adenocarcinoma" = "EAC", 
        "adenoca" = "EAC", 
        "intramucosal carcinoma" = "IMC", 
        "gastric metaplasia (no im)" = "NDBE"
    )
    if( tolower(diagnosis) %in% names(mapping)){
        shorthand <- mapping[tolower(diagnosis)] %>% unname 
    }else{
        if( grepl("^IND*", diagnosis)){
            shorthand <- "IND"
        }else if( grepl("^NDBE*", diagnosis)){
            shorthand <- "NDBE"
        }else if( grepl("^LGD*", diagnosis)){
            shorthand <- "LGD"
        }else{ 
            shorthand <- diagnosis
        }
    }  
    return(shorthand)
}

DELTA_district$diagnosis <- apply(DELTA_district, 1, CreateFinalDiagnosis)
DELTA_district$diagnosis <- sapply(DELTA_district$diagnosis, ConvertShorthand)
cf_df$diagnosis <- apply(cf_df, 1, CreateFinalDiagnosis)
cf_df$diagnosis <- sapply(cf_df$diagnosis, ConvertShorthand)


library(pwr)

pwr.p.test(
    h=ES.h(p1=0.035, p2=0.06), 
    sig.level = 0.05, 
    n=478,
    alternative="less"
)

pwr.p.test( # direct calculation of sample size for dysplasia rate
    h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05, 
    power = 0.9,
    alternative = "less"
)
pwr.p.test( # inferring power for HGD+ given above sample size
    h=ES.h(p1=0.0065, p2=0.02),
    sig.level=0.05, 
    n = 609,
    alternative="less"
)
pwr.p.test( # inferring power for dysplasia given protocol's sample size
    h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05,
    n=589,
    alternative="less"
)
pwr.p.test( # inferring power for HGD+ given protocol's sample size
    h=ES.h(p1=0.0065, p2=0.02),
    sig.level=c(0.05, 0.1),
    n=589,
    alternative="less"
)
pwr.p.test( # direct calculation of sample size for HGD+
    h=ES.h(p1=0.0065, p2=0.03),
    sig.level=0.05,
    power=0.9,
    alternative="less"
)
pwr.p.test( 
    h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05,
    power=0.8,
    alternative="less"
)
pwr.p.test( 
    h=ES.h(p1=0.0065, p2=0.03),
    sig.level=0.05,
    n=440,
    alternative="less"
)

power_levels <- c(0.8, 0.825, 0.85, 0.875, 0.9)
sample_size <- c()
HGD_power <- c()
for(power in power_levels){
    a <- pwr.p.test(
        h=ES.h(p1=0.035, p2=0.06),
        sig.level=0.05,
        power=power,
        alternative="less"
    )
    ss <- a$n
    sample_size[length(sample_size)+1] <- ss
    b <- pwr.p.test(
        h=ES.h(p1=0.0065, p2=0.03),
        sig.level=0.05,
        n=ss,
        alternative="less"
    )
}

sample_size <- seq(460, 550)
power_HGD <- c() 
power_dys <- c()
for(sample in sample_size){
    HGD <- pwr.p.test(
        h=ES.h(p1=0.0065, p2=0.03),
        sig.level=0.05,
        n=sample,
        alternative="less"
    )
    dys <- pwr.p.test(
      h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05,
    n =sample,
    alternative="less"
    )
    power_HGD[length(power_HGD)+1] <- HGD$power
    power_dys[length(power_dys)+1] <- dys$power
}
plot(sample_size, power_dys, lty=1, ylim = c(0.8,1))
lines(sample_size, power_HGD, lty=1, col ="blue")
abline(h=0.85, lty=3)



pwr.p.test( # direct calculation of sample size for dysplasia rate
    h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05, 
    n = 484,
    alternative = "less"
)
pwr.p.test( # direct calculation of sample size for dysplasia rate
    h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05, 
    power=0.85,
    alternative = "less"
)
pwr.p.test( # direct calculation of sample size for dysplasia rate
    h=ES.h(p1=0.0065, p2=0.03),
    sig.level=0.05, 
    n = 484,
    alternative = "less"
)

power_table <- data.frame(
    power1 = numeric() ,
    sample_size = numeric(),
    power2 = numeric(), 
    power3 = numeric()
)
power1 <- seq(0.8, 0.95, 0.05)
for(power in power1){
    calc1 <- pwr.p.test(
        h=ES.h(p1=0.035, p2=0.06),
        sig.level=0.05,
        power=power,
        alternative="less"
    )
    sample_size <- calc1$n

    calc2 <- pwr.p.test(
        h=ES.h(p1=0.0065, p2=0.03),
        sig.level=0.05,
        n = sample_size, 
        alternative="less"
    )
    power2 <- calc2$power 
    calc3 <- pwr.p.test(
        h=ES.h(p1=0.0065, p2=0.025),
        sig.level=0.05,
        n = sample_size, 
        alternative="less"
    )
    power3 <- calc3$power
    row <- c(power, sample_size, power2, power3) 
    power_table <- rbind(power_table, row)
}


pwr.p.test(
    h=ES.h(p1=0.035, p2=0.06),
    sig.level=0.05,
    n = 448, 
    alternative="less"
)
