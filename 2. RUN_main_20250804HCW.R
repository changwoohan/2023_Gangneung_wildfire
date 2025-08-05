#######################################################################################
# 2. Run main analysis ------------------------------------
# This section is for running main GSCM  

# Package install and function load
rm(list=ls())

# load the needed packages and functions
project.folder = paste0(print(here::here()),'/')
source(paste0(project.folder,'/code/packages/packages_to_load_20231205HCW.R'))
source(paste0(project.folder,'/code/upload/1. GSCM_functions_20250804HCW.R'))

# Data load ---------------------------------------------------------------
dir(here::here("data/use_data/revision"))

# main data load
name_1    <- paste0("data/use_data/revision/grwf_main.csv")
prev      <- read.csv(here::here(name_1), fileEncoding = 'cp949') 

# disease name loading
name_1    <- paste0("data/other_data/disease_name_20231231HCW.xlsx")
d_name    <- read_excel(here::here(name_1), sheet=1)
d_name$z  <- gsub(";", "\n", d_name$z )
d_name$xy <- gsub("(.{2})$", "", d_name$xy)                    # make the universal key for merging

# falsification data load : 2 year before the event 
name_2    <- paste0("data/use_data/revision/grwf_falsi.csv")
falsi     <- read.csv(here::here(name_2), fileEncoding = 'cp949') 

# Analysis example for single disease outcome: ALL_J_T  -------------------------------------

basic   <- ALL_J_T ~ int                                  # basic model setting 

nint <- prev %>% 
  filter(int==1) %>% 
  distinct(REGION_EMD) %>% 
  nrow()                                                  # to get the number of intervention region, in this study, it is 1

nam <- prev %>% 
  filter(int==1) %>% 
  distinct(REGION_EMD)                                    # to get the intervention region code REGION_EMD = 42150665

# GSCM analysis 
val.name  <- "ALL_J_T"
c.name    <- gsub("(.{2})$", "",   val.name)              # get the disease code 

fig_tit   <- d_name[d_name$xy ==  c.name,]$z              # get the figure label - "Diseases of the respiratory system \n(J00-J99)"
tab_tit   <- d_name[d_name$xy ==  c.name,]$y              # get the table label  - "Diseases of the respiratory system (J00-J99)"
uf        <- update.formula(basic, as.formula(paste(val.name, "~ int")))  # for future analysis (will keep change val.name using function)
out       <- gsynth(uf, data = prev, index = c("REGION_EMD", "time"), EM= F,
                    CV = TRUE, r = c(0, 5), force = "two-way",
                    nboots = 1000, inference = "parametric", se = TRUE, parallel = TRUE, seed=1234) # apply GSCM using gsynth package

bef <- prev %>% 
  ungroup() %>%
  filter(REGION_EMD %in% nam$REGION_EMD & eq==0) %>%    # to pick the period before for intervention region
  filter(time>-9 & time<0) %>%                          # limit the pre-treatment period to match post wildfire period (8 weeks)  
  select(val.name)                                      # select the target disease outcome 

bef_sum <- round(colSums(bef)/nint, 1)                    # average visit 1-8 weeks before the wildfire in intervention region 

# Figure generation based on GSCM results, Extracting observed and synthetic control values from the GSCM results
# There are also other ways to plot using basic plot function within gsynth package
k1 <-as.data.frame(out$est.att) %>% 
  mutate(time=as.numeric(row.names(out$est.att)))
k2 <-as.data.frame(out$Y.bar)  %>% 
  mutate(time=as.numeric(row.names(out$Y.bar)) +1)

k3 <-k2 %>% left_join(k1, by=c("time")) %>%
  mutate(Exposed   = Y.ct.bar + ATT,
         Control   = Y.ct.bar,
         exp_up    = Y.ct.bar + CI.upper - ATT,
         exp_low   = Y.ct.bar + CI.lower - ATT)

k4 <- k3 %>% 
  pivot_longer(
    cols =c("Exposed", "Control"),
    names_to="value",
    values_to="rate") %>%  
  mutate(value=as.factor(value)) %>% 
  mutate(value=relevel(value,"Exposed"))

ggplot(k4, aes_string(x="time", y="rate", group=as.factor(k4$value))) + 
  geom_line(aes(lty = as.factor(value)), size=1.0) +
  labs(title=fig_tit,
       x="Weeks from wildfire", 
       y="Number of hospital visits (n)") +
  coord_cartesian(ylim=c(0, max(k4$rate)*2), xlim = c(-32, 24)) +   # changed
  geom_vline(xintercept =0, linetype="longdash", size=1.0) +
  geom_vline(xintercept =4, linetype="longdash", size=1.0) +
  geom_vline(xintercept =8, linetype="longdash", size=1.0) +
  scale_x_continuous(breaks = c(-28,-24,-20,-16,-12,-8,-4 ,0, 4, 8, 12, 16, 20)) +
  geom_ribbon(aes(x = time, ymin = exp_low, ymax = exp_up),
              fill = "grey", alpha = 0.4) +      
  theme(
    plot.title = element_text(hjust=0.5, size=20, face="bold"),  
    panel.background = element_rect(fill = "white"),         # Set plot background to white
    legend.key  = element_rect(fill = "white"),              # Set legend item backgrounds to white
    axis.line.x = element_line(colour = "black", size = 1),  # Add line to x axis
    axis.line.y = element_line(colour = "black", size = 1),   # Add line to y axis
    legend.position = c(.95, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(2, 2, 2, 2),
    legend.title = element_blank(),
    legend.text=element_text(size=20),
    legend.key.height = unit(0.4,"cm"),
    legend.key.width = unit(3,"cm"),
    axis.title=element_text(size=15),
    axis.title.x=element_text(size=15),
    axis.title.y=element_text(size=15),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15)) 

# To calculate cumulative number of excess visits during 0-3 and 0-7 weeks using cumuEff function  
in.list=NULL
out.list=NULL
pe      <- c(0, 3, 0, 7)                                  # cumulative period setting 4 weeks (0-3), 8 weeks (0-7) after the wildfire 
for (k in (1:2)) {
  
  cumu <- cumuEff(out, cumu = TRUE, id = NULL, period = c(pe[2*k -1],pe[2*k]))
  cumu <- as.data.frame(cumu$est.catt)  
  cumu <- cumu %>% mutate (dd   =out$Y,
                           CATT =round(CATT,1),
                           low  =round(CI.lower,1),
                           high =round(CI.upper,1), 
                           NO   =row.names(cumu)) %>% 
    select(c("dd","CATT","low","high","p.value","NO")) %>% 
    filter(row_number()==n())
  
  cumu_1<- cbind(cumu, bef_sum) %>% 
    mutate(pe   = round(CATT/bef_sum *100,1),
           l_pe = round(low/bef_sum *100,1),
           h_pe = round(high/bef_sum *100,1),
           DIS  = tab_tit) %>% 
    mutate(pre_event = bef_sum) %>% 
    mutate(diff_out = paste0(CATT," (",low,", ",high,")")) %>% 
    mutate(diff_per = paste0(pe," (",l_pe,", ",h_pe,")")) %>% 
    mutate(DIS2 = val.name) %>% 
    mutate(pvalue =p.value) %>% 
    select(DIS2, DIS, pre_event, diff_out, diff_per)  
  
  in.list[[k]]  <- cumu_1  
} 

print(in.list[[1]])



# Now do it with function: han_gsynth_95() --------------------------------------------------

# Get variable target disease outcomes
MAIN       <- c("ALL_J_T", "J_1_T", "J_2_T", "J_3_T", "J_4_T", "J_5_T", "J_6_T", "J_7_T", "J_8_T", "J_9_T",
                "ALL_I_T", "I_3_T", "I_4_T","I_5_T", "I_6_T", "I_7_T", "I_8_T", "I_10_T")     # For main analysis  

FIG_2      <- c("ALL_J_T", "ALL_I_T", "J_1_T", "J_2_T", "J_3_T", "J_5_T", "I_3_T", "I_4_T")   # For Figure 2 plotting

FIG_S10    <- c("ALL_J_T", "ALL_I_T", "J_1_T", "J_2_T", "J_3_T", "J_4_T", "J_5_T", "J_7_T", "J_8_T", "J_9_T",
                "I_3_T", "I_4_T","I_5_T", "I_6_T", "I_7_T", "I_8_T", "I_10_T")                # For falsification analysis

# Run the analysis using function -----------------------------------------
han_gsynth_95(dat.use=prev,      val.name=MAIN,     folder.name="20250804", table.name="TABLE2",  figure.name="TABLE2")
han_gsynth_95(dat.use=prev,      val.name=FIG_2,    folder.name="20250804", table.name="FIG_2",   figure.name="FIG_2")
han_gsynth_95(dat.use=falsi,     val.name=FIG_S10,  folder.name="20250804", table.name="FIG_S10", figure.name="FIG_S10")


# Draw figure -------------------------------------------------------------

# Figure 2
fig.list        <- readRDS(file.path(here::here("figures/20250804/FIG_2.rds")))
folder.name     <-"20250804"
path            <- here::here(paste0("figures/", folder.name)) # to create the figure output folder

png(file=paste0(path,"/Figure2.tiff"), 
    res=50, width=25, height=30, units="in",  pointsize=70)

grid.arrange(fig.list[[1]],  fig.list[[2]],
             fig.list[[3]],  fig.list[[4]],
             fig.list[[5]],  fig.list[[6]],
             fig.list[[7]],  fig.list[[8]], nrow=4, ncol=2)
dev.off()


#FIG S10 falsification test
fig.list        <- readRDS(file.path(here::here("figures/20250804/FIG_S10.rds")))
folder.name     <-"20250804"
path            <- here::here(paste0("figures/", folder.name)) # to create the figure output folder

png(file=paste0(path,"/FigureS10_1.tiff"), 
    res=50, width=25, height=35, units="in",  pointsize=70)

grid.arrange(fig.list[[1]],  fig.list[[2]],
             fig.list[[3]],  fig.list[[4]],
             fig.list[[5]],  fig.list[[6]],
             fig.list[[7]],  fig.list[[8]],
             fig.list[[9]],  fig.list[[10]], nrow=5, ncol=2)
dev.off()

png(file=paste0(path,"/FigureS10_2.tiff"), 
    res=50, width=25, height=35, units="in",  pointsize=70)

grid.arrange(fig.list[[11]],  fig.list[[12]],
             fig.list[[13]],  fig.list[[14]],
             fig.list[[15]],  fig.list[[16]],   nrow=5, ncol=2)
dev.off()

