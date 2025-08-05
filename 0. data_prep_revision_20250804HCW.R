
#######################################################################################
# Code layout ---------------------------------------------------------
# Health effects of wildfire occurred on 2023-04-11 for 8 hours in Gangneung city
# Wrote on 2025-07-17 during the manuscript revision by Changwoo Han
# 0. Dataset preparation for GSCM
# 1. GSCM function generation
# 2. Run final model
#######################################################################################

#######################################################################################
# 0. Dataset preparation for GSCM ------------------------------------------------------------
# This section is for preparing a dataset for main analysis and simple statistical analysis for Table S2.  

rm(list=ls())
# load the packages 
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'/code/packages/packages_to_load_20231205HCW.R'))

## Data loading ------------------------------------------------------------
# Loading exported NHIS data (already in weekly format) for main analysis 
# Disease prevalence data includes number of hospital visits from each district 

dir(here::here("data/raw_data"))                                       # check the data in the folder
name_1   <- paste0("data/raw_data/week_emd_prev_grwf_final.csv")       # raw data is already weekly and emd data 
a_1      <-read.csv(here::here(name_1)) %>%                            # import the data
  mutate(REGION_EMD= as.factor(REGION_EMD))                            # change the region code to the factor variable for analysis    

head(a_1)
tail(a_1)                    # 2019-12-31 to 2023-09-05 period is available for each district of Gangneung City 
a_1[is.na(a_1)] <-0          # change NA to 0

a_1 %>% distinct(REGION_EMD) # Total of 21 district confirmed in Gangreung-si


## Merge to the district name data ---------------------------------------
# load the sgg coding file (for labelling)
dir(here("data/other_data"))
sgg <-read_excel(here::here("data/other_data/hangan_sgg_code_20230522HCW.xlsx"), sheet=1)

# extract the code and name of each Gangneung-si district
gr_emd <- sgg %>% mutate(REGION_SIDO = substr(ADDR_CD,1,2),
                         REGION_SGG  = substr(ADDR_CD,1,5),
                         REGION_EMD  = substr(ADDR_CD,1,8),
                         KOR_name    = as.factor(SGG_NM)) %>% 
  select(REGION_SGG, REGION_SIDO, REGION_SGG, KOR_name, REGION_EMD, EMD_NM) %>% 
  distinct(REGION_EMD, REGION_SIDO, REGION_SGG, EMD_NM) %>% 
  filter(REGION_SGG %in% "42150")    # limit to Gangneung City

a_2 <- a_1 %>% 
  group_by(REGION_EMD) %>% 
  arrange(DT, .by_group = TRUE) %>% 
  mutate(serial =row_number()) %>%         # add serial variable
  left_join(gr_emd, by="REGION_EMD") %>%   # join the region label 
  select( -REGION_SGG, -REGION_SIDO)


## Add variables for GSC analysis ------------------------------------

# check the event (wildfire) serial 
a_2 %>% 
  filter(DT == as.Date("2023-04-11")) %>% 
  select(serial)  # its 172

a_3 <- a_2 %>% 
  mutate(eq =ifelse(DT<as.Date("2023-04-11"), 0, 1),           # define the pre and post wildfire period
         intt=ifelse(REGION_EMD %in% c(42150665)==1 , 1, 0),   # only 42150665 Gyeongpo district is the intervention region 
         int =ifelse(intt==1 & eq==1, 1, 0)) %>%               # make 'int' variable for GSCM 
  mutate(time=serial-172) %>%           # make time variable 
  filter(DT>as.Date("2022-04-11")) %>%  # data period 2022-04-12 to 2023-09-26
  filter(time<25)                       # week 24 is the longest post wildfire period

table(a_3$time)                         # pre period  : -52 to -1 week before the wildfire
                                        # post period : 1 to 24 week after the wildfire

# change logical variable to numeric: if there are many 0, the R think its logical variable 
cols         <- sapply(a_3, is.logical)
a_3[,cols]   <- lapply(a_3[,cols], as.numeric)

# generate the _W (women) and _Y (young) variables automatically: Currently, total (_T), men (_M), and elderly (_E) counts are exported from NHIS server
SPE_T <-  colnames(a_3)[str_detect(colnames(a_3), c("_T"))]
SPE_T <-  gsub("_T", "", SPE_T)

for (var_name in SPE_T) {
  new_var_name <- paste0(var_name, "_W")
  new_var_name2 <- paste0(var_name, "_Y")
  a_3 <- a_3 %>%
    mutate(!!sym(new_var_name) := .data[[paste0(var_name, "_T")]] - .data[[paste0(var_name, "_M")]],
           !!sym(new_var_name2) := .data[[paste0(var_name, "_T")]] - .data[[paste0(var_name, "_A")]] - .data[[paste0(var_name, "_E")]])
}

a_3[is.na(a_3)] <-0   # change NA to 0

# data out for main analysis 

dir.create(here::here("data/use_data/revision"))
write.csv(a_3, here::here("data/use_data/revision/grwf_main.csv"),
          fileEncoding='cp949', row.names = F)                         # main dataset ready
#######################################################################################


#######################################################################################
# Data preparation for falsification test  --------------------------------
# 2 years before the wildfire
b_3 <- a_2 %>% 
  mutate(eq =ifelse(DT<as.Date("2021-04-13"), 0, 1),
         intt=ifelse(REGION_EMD %in% c(42150665)==1 , 1, 0),
         int =ifelse(intt==1 & eq==1, 1, 0)) %>% 
  mutate(time=serial-67) %>% 
  filter(DT>as.Date("2020-04-06")) %>%  # 2022-04-12 to 2023-09-26
  filter(time<25)

# change logical variable to numeric
cols         <- sapply(b_3, is.logical)
b_3[,cols]   <- lapply(b_3[,cols], as.numeric)

# generate the _W and _Y variable automatically
SPE_T <-  colnames(b_3)[str_detect(colnames(a_3), c("_T"))]
SPE_T <-  gsub("_T", "", SPE_T)

for (var_name in SPE_T) {
  new_var_name  <- paste0(var_name, "_W")
  new_var_name2 <- paste0(var_name, "_Y")
  b_3 <- b_3 %>%
    mutate(!!sym(new_var_name) := .data[[paste0(var_name, "_T")]] - .data[[paste0(var_name, "_M")]],
           !!sym(new_var_name2) := .data[[paste0(var_name, "_T")]] - .data[[paste0(var_name, "_A")]] - .data[[paste0(var_name, "_E")]])
}

b_3[is.na(b_3)] <-0   # change NA to 0

# data out for falsification analysis 
dir.create(here::here("data/use_data/revision"))
write.csv(b_3, here::here("data/use_data/revision/grwf_falsi.csv"),
          fileEncoding='cp949', row.names = F)


#######################################################################################
# Basic analysis for Table 1 ----------------------------------------------
# Descriptive analysis for 8 weeks before and after the wildfire

a_4 <- a_3 %>% 
  filter(time<8 & time>-9) %>% 
  mutate(t1 = ifelse (time<0, "Before", "After"),
         t2 = ifelse (REGION_EMD %in% c(42150665), "Exposed", "Control")) 
  
SPE_T   <- c("ALL_J_T", "J_1_T", "J_2_T", "J_3_T", "J_4_T", "J_5_T", "J_6_T", "J_7_T", "J_8_T", "J_9_T",
             "ALL_I_T", "I_3_T", "I_4_T","I_5_T", "I_6_T", "I_7_T", "I_8_T", "I_10_T")   # Target disease categories

# disease name loading
name_1   <- paste0("data/other_data/disease_name_20231231HCW.xlsx")
d_name   <- read_excel(here::here(name_1), sheet=1)
d_name$z <- gsub(";", "\n", d_name$z )

in.list <- list()

for (i in seq_along(SPE_T)) {
  go <- SPE_T[i]
  tab_tit <- d_name[d_name$xy == SPE_T[i], ]$y
  
  a_5 <- a_4 %>%
    group_by(REGION_EMD, t1) %>%
    summarise(mean_var = mean(!!sym(SPE_T[i])), .groups = "drop")
  
  a_6 <- a_5 %>%
    pivot_wider(names_from = t1, values_from = mean_var) %>%
    mutate(Diff = After - Before) %>%
    mutate(int = ifelse(REGION_EMD == "42150665", "Exposed", "Control"))
  
  a_7 <- a_6 %>%
    group_by(int) %>%
    summarise(
      m_Before = mean(Before),
      m_After = mean(After),
      m_Diff = mean(Diff),
      sd_Before = if_else(first(int) == "Control", sd(Before), NA_real_),
      sd_After = if_else(first(int) == "Control", sd(After), NA_real_),
      sd_Diff = if_else(first(int) == "Control", sd(Diff), NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      percent_change = round(( (m_After-m_Before) / m_Before) * 100, 2)
    ) %>%
    mutate_if(is.numeric, round, digits = 1) %>%
    mutate(
      before = if_else(int == "Exposed",
                       as.character(m_Before),
                       paste0(m_Before, " (", sd_Before, ")")),
      after = if_else(int == "Exposed",
                      as.character(m_After),
                      paste0(m_After, " (", sd_After, ")")),
      diff = if_else(int == "Exposed",
                     as.character(m_Diff),
                     paste0(m_Diff, " (", sd_Diff, ")"))
    ) %>%
    select(int, before, after, diff, percent_change) %>%
    pivot_wider(
      names_from = int,
      values_from = c(before, after, diff, percent_change)
    ) %>%
    mutate(disease = tab_tit) %>%
    select(
      disease,
      before_Exposed, after_Exposed, diff_Exposed, percent_change_Exposed,
      before_Control, after_Control, diff_Control, percent_change_Control
    )
  
  in.list[[i]] <- a_7
}

descriptive_tableS2 <- bind_rows(in.list)

dir.create(here::here("output/revision"))
write.csv(descriptive_tableS2, here::here("output/revision/tableS2_20250717.csv"))

#######################################################################################