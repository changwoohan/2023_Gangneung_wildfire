
#######################################################################################
# 1. Making functions based on GSCM gsynth package ------------------------------------
# This section is for making a function for the main analysis 
# Store the figure in RDS file format 

han_gsynth_95 = function(dat.use, val.name, folder.name, table.name, figure.name) {
  
  fig.list=NULL
  in.list=NULL
  out.list=NULL

  pe     <- c(0, 3, 0, 7)                               # cumulative period setting 4 weeks, 8 weeks after the wildfire 
  path   <- here::here(paste0("figures/", folder.name)) # to create the figure output folder
  dir.create(path)                                      # to create the figure output folder 
  basic <- COUNT_A_1_d0 ~ int                           # basic model setting
  
  nint <- dat.use %>% 
    filter(int==1) %>% 
    distinct(REGION_EMD) %>% 
    nrow()                       # to get the number of intervention region, in this study, it is 1
  
  nam <- dat.use %>% 
    filter(int==1) %>% 
    distinct(REGION_EMD)         # to get the intervention region code

  for (i in (1:length(val.name))){ 
    val.name[i]                                              # target for this round
    c.name    <- gsub("(.{2})$", "", val.name)

    fig_tit   <- d_name[d_name$xy ==  c.name[i],]$z          # get the figure label
    tab_tit   <- d_name[d_name$xy ==  c.name[i],]$y          # get the table label
    uf        <- update.formula(basic, as.formula(paste(val.name[i], "~ int")))
    out       <- gsynth(uf, data = dat.use, index = c("REGION_EMD", "time"), EM= F,
                  CV = TRUE, r = c(0, 5), force = "two-way",
                  nboots = 1000, inference = "parametric", se = TRUE, parallel = TRUE, seed=1234)
    
    bef <- dat.use %>% 
      ungroup() %>%
      filter(REGION_EMD %in% nam$REGION_EMD & eq==0) %>%  # pick the period before for intervention region
      filter(time>-9 & time<0) %>%                        # limit the pre-treatment period to match post wildfire period (8 weeks)  
      select(val.name[i])                                 # select the target 

    bef_sum <- round(colSums(bef)/nint, 1)                # average visit before the wildfire in intervention region 
    
    # figure generation
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
    
    fig.list[[i]] <- ggplot(k4, aes_string(x="time", y="rate", group=as.factor(k4$value))) + 
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
        mutate(DIS2 = val.name[i]) %>% 
        mutate(pvalue =p.value) %>% 
        select(DIS2, DIS, pre_event, diff_out, diff_per)  
      
      in.list[[k]]  <- cumu_1  
      out.list[[i]] <- as.data.frame(do.call(cbind, in.list))
      
    } 
    
  }
  
  final_out   <- as.data.frame(do.call(rbind, out.list)) 
  saveRDS(fig.list, paste0(path,"/",figure.name,".rds")) 
  
  output_ads  <- here::here(paste0("output/",folder.name))
  dir.create(output_ads)
  write.csv(final_out, paste0(output_ads,"/", table.name,"grwffire.csv"))
  
}

