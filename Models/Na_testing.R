#Model notes------
#Based on Poh et al (2012)


#Check parval values and values in the rate equations for markov map
#ensure equations are correct - testing the P total
# correct the end equations as they will need to be adjusted for Na
#compare predictions, then implement all into the full model :D

#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
Na_I <- read.csv("data/Na_IV.csv")
Na_I$x <- round(Na_I$x,0)

#normalisation
norm_stat_val <- function(b){
  sum_b <- sum(b)
  if(sum_b == 0) return(b)
  b / sum_b
}

# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)
      
      # Rate equations
      # Nernst Potentials----
      E_Ca <- ((R * Temp) / (2 * Faraday)) * log(Ca_o / Ca_i_total)
      #Single slow wave-----
      if (times >= 0 & times < t_peak_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (times/f2)
      } 
      if (times >= t_peak_ICC & times < t_plateau_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (1 + exp(-f1/(2 * t_slope))) * (1/(1 + exp((times - f2 - 0.5 * f1)/t_slope)))
      } 
      
      if (times < clamp_start | times > clamp_end ){
        Vm_eq <- Hv
      } else{
        Vm_eq <- Vm
      } 
      # Voltage-dependent Sodium Current (INa)-----
      
      
      KOI1 <- parval[23] * exp(parval[1] + parval[2] * Vm_eq) * T_correction_Na  # (S-80)
      KI1I2 <- parval[24] * exp(parval[3] + parval[4] * Vm_eq) * T_correction_Na # (S-81)
      KC3C2 <- parval[25] * exp(parval[5] + parval[6] * Vm_eq) * T_correction_Na # (S-82)
      KC2C1 <- parval[26] * exp(parval[7]+ parval[8] * Vm_eq) * T_correction_Na
      KC1O <- parval[27] * exp(parval[9] + parval[10] * Vm_eq)* T_correction_Na
      KI2I1 <- parval[28] * exp(parval[11] + parval[12] * Vm_eq)* T_correction_Na
      KC2C3 <- parval[29] * exp(parval[13] + parval[14] * Vm_eq)* T_correction_Na
      KC1C2 <- parval[30] * exp(parval[15] + parval[16] * Vm_eq)* T_correction_Na
      KOC1 <- parval[31] * exp(parval[17] + parval[18] * Vm_eq)* T_correction_Na
      KI1C1 <- parval[32] * exp(parval[19] + parval[20] * Vm_eq)* T_correction_Na
      KC1I1 <- parval[33] * exp(parval[21] + parval[22] * Vm_eq)* T_correction_Na
      KI1O <- parval[34] * exp(parval[35] + parval[36] * Vm_eq)* T_correction_Na
      #ODE-----
      # Temp-type Calcium Channel ODEs-----
      d[1] <- KI2I1 * y["I2_Na"] + KOI1 * y["O_Na"] + KC1I1 * y["C1_Na"] - (KI1I2 + KI1O + KI1C1) * y["I1_Na"]  # Inactivated State 1
      d[2] <- KI1I2 * y["I1_Na"] - KI2I1 * y["I2_Na"] # Inactivated 2
      d[3] <- KI1O * y["I1_Na"] + KC1O * y["C1_Na"] - (KOI1 + KOC1) * y["O_Na"]  # Open 
      d[4] <- KOC1 * y["O_Na"] + KI1C1 * y["I1_Na"] + KC2C1 * y["C2_Na"] - (KC1O + KC1C2 + KC1I1) * y["C1_Na"] # Closed state 1
      d[5] <- KC1C2 * y["C1_Na"] + KC3C2 * y["C3_Na"] - (KC2C1 + KC2C3) * y["C2_Na"] # Closed state 2
      d[6] <- KC2C3 * y["C2_Na"] - KC3C2 * y["C3_Na"] # Closed state 3
      return(list(d))
    })
  }
  
  #Fixed parameters
  fixed <- list()
  R <- 8.314  # J/(molK) - Ideal gas constant
  Faraday <- 96.48534  # C/mmol - Faraday’s constant
  Temp <- 310  # K - Temperature
  Texp <- 297.0 #experimental temperatures
  Cm <- 50  # pF - Cell membrane Capacitance
  Vcell <- 3.5e-12  # l - Cell volume
  Ca_o <- 2  # mM - Extracellular Calcium concentration
  K_o <- 5.4  # mM - Extracellular potassium concentration
  Na_o <- 140  # mM - Extracellular sodium concentration
  Ca_i_total <- 0.004914  # mM - Initial total intracellular Calcium concentration
  Ca_i_free <- 1.26e-4  # mM - Initial free intracellular Calcium concentration
  K_i <- 150  # mM - Intracellular potassium concentration
  Na_i <- 10.5  # mM - Intracellular sodium concentration
  CaQ10 <- 2.1  # - Q10 for Calcium channels
  KQ10 <- 3.1  # - Q10 for potassium channels
  T_correction_Ca <- 1.0 * CaQ10^((Temp - Texp) / 10.0) # experimental correction
  NaQ10 <- 2.45  # - Q10 for sodium channels
  T_correction_Na <- 1.0 * NaQ10^((Temp-Texp)  /  10.0)
  Gcouple <- 2.6  # nS - Coupling conductance between ICC and SMC
  V_rest_ICC <- -57  # mV - Resting membrane potential of ICC
  V_peak_ICC <- -23.5  # mV - Peak membrane potential of ICC
  V_amp_ICC <- 33.5  # mV - Amplitude of ICC membrane potential
  t_period <- 10000  # ms - Period of single ICC slow wave
  t_peak_ICC <- 300  # ms - Duration of ICC slow wave upstroke
  t_plateau_ICC <- 9700  # ms - Duration of ICC slow wave plateau phase
  t_slope <- 600  # ms - Slope factor in ICC membrane equation
  f1 <- 12000  # ms - ICC conditioning factor 1
  f2 <- 300  # ms - ICC conditioning factor 2
  CRT_total <- 0.034  # mM - Total Calreticulin concentration
  CRT_n <- 1  # - Hill coefficient for Calreticulin
  CRT_KD <- 0.0009  # mM - Dissociation constant for Calreticulin
  CaM_total <- 0.012  # mM - Total Calmodulin concentration
  CaM_n <- 4  # - Hill coefficient for Calmodulin
  CaM_KD <- 0.0001  # mM - Dissociation constant for Calmodulin
  G_CaL <- 1.44  # nS - Maximum conductance of CaL
  G_CaT <- 0.0425  # nS - Maximum conductance of CaT
  G_Kv <- 1.0217  # nS - Maximum conductance of Kv
  tau_x_Kv <- 4.7803  # ms - Time constant for Kvx
  tau_y_Kv <- 763.7564  # ms - Time constant for Kvy
  G_BK <- 80  # nS - Maximum conductance of BK
  G_Na <- 25.1  # nS - Maximum conductance of Na
  P_NCX <- 39.8437 # pA/pF - Maximum NCX
  K_mCa <- 1.38  # mM - Ca_i half saturation constant of NCX
  K_mNa <- 87.5  # mM - Na_i half saturation constant of NCX
  k_sat <- 0.1 # - Saturation factor for NCX
  gamma <- 0.35  # - Voltage dependence parameter of NCX
  P_NaK <- 0.1852  # pA/pF - Maximum NaK
  K_mK <- 1  # mM - K_o half saturation constant of NaK
  K_mNa <- 40  # mM - Na_i half saturation constant of NaK
  G_NSNa <- 0.022488  # nS - Maximum conductance of non-selective Na current
  G_NSK <- 0.017512  # nS - Maximum conductance of non-selective K current
  tau_d_CaT <- 1.9058 / T_correction_Ca #for ICaL
  sigma_LCaL <- 0# 0.01 #FOR ICaL set to 0 to replicate EGTA
  k_on <- 40633 #BK on rate
  k_c_off <- 11 #BK closed off rate
  k_o_off <- 1.1 #BK open off rate
  SMC_resting <- -60 #
  #list of values in C code.
  parval <- c(0.307629865993612,0.00605347675235982,0.0514901233349612,-0.0468650997035773,-0.0691017278134239,0.00319450166548765,
              -0.156598571141704,0.0583527292113407,0.00931925729714207,0.0410752869450889,2.67926519075303,0.00614677585983095,
              -0.0990741811465619,0.0364411841506491,0.363520774026992,0.0771931771113635,-13.3349821299039,-0.252888485042226,
              -2.48395798104711,0.0204061599230419,-0.0634382282818863,0.00466834101392052,1.61639957785747,0.0277351044118711,
              0.000525484598535321,1.4496348286393,1.53287364507221,0.00392388356048677,0.554315057682815,3.15656890165025,
              2.39152928164775,1.90461121472293,0.00021688112860783,0.120515381956888,-9.60280306718205,0.0830250234337321 )
  # Initial values for A, B, and C
  initial <- c(
    I1_Na = 0.000258498,
    I2_Na = 0.0252312,
    O_Na = 1.25534e-05,
    C1_Na = 0.0127485,
    C2_Na = 0.0343829,
    C3_Na = 0.136166
  )
  
  # Time sequence for the simulation
  times <- seq(0, 1000, by = 1)
  
  #history
  
  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$P_total <- rowSums(output_df[2:7])
  output_df$E_Na <- ((R * Temp) / (2 * Faraday)) * log(Na_o / Na_i)
  output_df$I_Na <- G_Na * output_df$O_Na *  
    (ifelse(output_df$time < parms[["clamp_start"]] | output_df$time > parms[["clamp_end"]], parms[["Hv"]], parms[["Vm"]]) - output_df$E_Na)
  output_df$Vm <- ifelse(
    output_df$time < parms[["clamp_start"]] | output_df$time > parms[["clamp_end"]], 
    parms[["Hv"]],  # Use Hv before clamp_start or after clamp_end
    parms[["Vm"]]   # Use Vm between clamp_start and clamp_end
  )
  
  return(output_df)
}

#Parameters to fit -----
mv_tests <- seq(-80, 35, by = 5)
sim_v <- list()
peak_store <- numeric(length(mv_tests))
i <- 1
for(i in seq_along(mv_tests)) {  # Correct looping over mv_tests
  parms <- list(Vm = mv_tests[i], Hv = -100, clamp_start = 475, clamp_end = 525)  # Set voltage parameter
  sim_temp <- Model(parms = parms)  # Run simulation
  sim_temp$Vm_identity <- parms[["Vm"]]
  sim_v[[i]] <- as.data.frame(sim_temp)
  sim_temp_peak <- subset(sim_temp, Vm == mv_tests[i])
  peak_store[i] <- min(sim_temp_peak$I_Na)
}

x <- 1
for(x in 1:length(sim_v)) {  # Correct looping over mv_tests
  sim_v[[x]]$I_Na_norm <- sim_v[[x]]$I_Na / min(peak_store)
  #sim_v[[x]]$I_Na_norm <- norm_stat_val(sim_v[[x]]$I_Na)
}


sim_v_df <- bind_rows(sim_v)  # Convert list to dataframe

# Plot all conditions with color mapped to Vm
S6.plot <- ggplot(sim_v_df, aes(x = time, y = I_Na_norm, color = factor(Vm_identity), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "Time vs Normalized Current",
    x = "Time (msec)",
    y = "Normalized Current",
    color = "Voltage (mV)"  # Legend label
  ) +
  scale_y_reverse()+
  xlim(470,520) +
  theme_minimal()  

# Display the plot
print(S6.plot)



# Create dataframe with results
IV.df <- data.frame(mV = mv_tests, I = peak_store)
IV.df$I_norm <- IV.df$I / min(IV.df$I)

IV.plot <- ggplot() +
  geom_point(data = Na_I, aes(x = x, y = y), color = "red", size = 3) +  # Scatter points for experimental data
  geom_line(data = IV.df, aes(x = mV, y = I_norm), color = "blue", linewidth = 1) +  # Model I-V curve
  labs(
    title = "Normalized I-V Curve for I_Na",
    x = "Membrane Voltage (mV)",
    y = "Normalized Peak I_Na"
  ) +
  scale_y_reverse() +
  xlim(-90,50) +
  theme_minimal()  # Clean theme

# Display the plot
print(IV.plot)