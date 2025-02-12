#Model notes------
#Based on Poh et al (2012)
#This one investigates probability of open over 3 intracellular calcium
#So this model investigates that as well


#addd in calcium buffering as it gives Ca2+ ree ? may be important?
#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
BK_0.0001<- read.csv("data/BK_0.0001.csv")
BK_0.0003<- read.csv("data/BK_0.0003.csv")
BK_0.001<- read.csv("data/BK_0.001.csv")
BK_0.0001$x <- round(BK_0.0001$x,0)
BK_0.0003$x <- round(BK_0.0003$x,0)
BK_0.001$x <- round(BK_0.001$x,0)
# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)
      
      # Rate equations
      # Nernst Potentials----
      
      # Calcium & Voltage-activated Potassium Current (IBK)-----
      a_BK <- exp((0.73 * Faraday * Vm)/(Temp * R))
      b_BK <- exp((-0.67 * Faraday * Vm)/(Temp * R))
      
      #voltage dependent transitions
      KC0O0_BK <- a_BK * 21.6200042030946 * 10^-3
      KC1O1_BK <- a_BK * 0.868781002082468 * 10^-3
      KC2O2_BK <- a_BK * 0.028063354826539 * 10^-3
      KC3O3_BK <- a_BK * 0.781455969444612 * 10^-3
      KC4O4_BK <- a_BK * 44.323856772419 * 10^-3
      
      KO0C0_BK <- b_BK * 318108.401217897 * 10^-3
      KO1C1_BK <- b_BK * 144173.645394299 * 10^-3
      KO2C2_BK <- b_BK * 32659.4033438573 * 10^-3
      KO3C3_BK <- b_BK * 95.3119836990264 * 10^-3
      KO4C4_BK <- b_BK * 0.105563199075352 * 10^-3
      
      # Define rate equations for calcium-dependent transitions
      KC0C1_BK <- (Ca_i_free * 4 * k_on)
      KC1C2_BK <- (Ca_i_free * 3 * k_on)
      KC2C3_BK <- (Ca_i_free * 2 * k_on)
      KC3C4_BK <- (Ca_i_free * 1 * k_on)
      
      KC4C3_BK <- (4 * k_c_off)
      KC3C2_BK <- (3 * k_c_off)
      KC2C1_BK <- (2 * k_c_off)
      KC1C0_BK <- (1 * k_c_off)
      
      KO0O1_BK <- (Ca_i_free * 4 * k_on)
      KO1O2_BK <- (Ca_i_free * 3 * k_on)
      KO2O3_BK <- (Ca_i_free * 2 * k_on)
      KO3O4_BK <- (Ca_i_free * 1 * k_on)
      
      KO4O3_BK <- (4 * k_o_off)
      KO3O2_BK <- (3 * k_o_off)
      KO2O1_BK <- (2 * k_o_off)
      KO1O0_BK <- (1 * k_o_off)
      
      #ODE-----
      # Temp-type Calcium Channel ODEs-----
      #Closed
      d[1] <- KC1C0_BK * y["C1_BK"] + KO0C0_BK * y["O0_BK"] - (KC0C1_BK + KC0O0_BK) * y["C0_BK"]  # closed state of BK channel 0
      d[2] <- KC0C1_BK * y["C0_BK"] + KC2C1_BK * y["C2_BK"] + KO1C1_BK * y["O1_BK"] - (KC1C0_BK + KC1O1_BK + KC1C2_BK) * y["C1_BK"]  # closed state of BK channel 1
      d[3] <- KC1C2_BK * y["C1_BK"] + KC3C2_BK * y["C3_BK"] + KO2C2_BK * y["O2_BK"] - (KC2C1_BK + KC2O2_BK + KC2C3_BK) * y["C2_BK"]  # closed state of BK channel 2
      d[4] <- KC2C3_BK * y["C2_BK"] + KC4C3_BK * y["C4_BK"] + KO3C3_BK * y["O3_BK"] - (KC3C2_BK + KC3O3_BK + KC3C4_BK) * y["C3_BK"] # closed state of BK channel 3
      d[5] <- KC3C4_BK * y["C3_BK"] + KO4C4_BK * y["O4_BK"] - (KC4C3_BK + KC4O4_BK) * y["C4_BK"]  # closed state of BK channel 4
      #Open
      d[6] <- KO1O0_BK * y["O1_BK"] + KC0O0_BK * y["C0_BK"] - (KO0O1_BK + KO0C0_BK) * y["O0_BK"]  # Open state of BK channel 0
      d[7] <- KO0O1_BK * y["O0_BK"] + KO2O1_BK * y["O2_BK"] + KC1O1_BK * y["C1_BK"] - (KO1O0_BK + KO1C1_BK + KO1O2_BK) * y["O1_BK"]  # Open state of BK channel 1
      d[8] <- KO1O2_BK * y["O1_BK"] + KO3O2_BK * y["O3_BK"] + KC2O2_BK * y["C2_BK"] - (KO2O1_BK + KO2C2_BK + KO2O3_BK) * y["O2_BK"] # Open state of BK channel 2
      d[9] <- KO2O3_BK * y["O2_BK"] + KO4O3_BK * y["O4_BK"] + KC3O3_BK * y["C3_BK"] - (KO3O2_BK + KO3C3_BK + KO3O4_BK) * y["O3_BK"]  # Open state of BK channel 3
      d[10] <- KO3O4_BK * y["O3_BK"] + KC4O4_BK * y["C4_BK"] - (KO4O3_BK + KO4C4_BK) * y["O4_BK"]  # Open state of BK channel 4
      return(list(d))
    })
  }
  
  #Fixed parameters
  fixed <- list()
  R <- 5.189*1e19  # J/(molK) - Ideal gas constant
  Faraday <- 6.022*1e23  # C/mmol - Faraday’s constant
  Temp <- 310  # K - Temperature
  Cm <- 50  # pF - Cell membrane Capacitance
  Vcell <- 3.5e-12  # l - Cell volume
  Ca_o <- 2  # mM - Extracellular Calcium concentration
  K_o <- 5.4  # mM - Extracellular potassium concentration
  Na_o <- 140  # mM - Extracellular sodium concentration
  Ca_i_total <- 0.004914  # mM - Initial total intracellular Calcium concentration
  #Ca_i_free <- 1.26e-4  # mM - Initial free intracellular Calcium concentration
  K_i <- 150  # mM - Intracellular potassium concentration
  Na_i <- 10.5  # mM - Intracellular sodium concentration
  CaQ10 <- 2.1  # - Q10 for Calcium channels
  KQ10 <- 3.1  # - Q10 for potassium channels
  NaQ10 <- 2.45  # - Q10 for sodium channels
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
  tau_d_CaT <- 1.9058 #for ICaL
  sigma_LCaL <- 0# 0.01 #FOR ICaL set to 0 to replicate EGTA
  k_on <- 40.6326276428765 * 10^9 * 10^-3#BK on rate
  k_c_off <- 1.0e9*11.0*1.0e-6 * 10^-3 #BK closed off rate
  k_o_off <- 1.0e9*1.1*1.0e-6 * 10^-3 #BK open off rate
  SMC_resting <- -60 #
  
  
  
  # Initial values for A, B, and C
  
  initial <- c(
    C0_BK = 0.1,   # Closed, 0 Ca²⁺
    C1_BK = 0.1,  # Closed, 1 Ca²⁺
    C2_BK = 0.1,  # Closed, 2 Ca²⁺
    C3_BK = 0.1, # Closed, 3 Ca²⁺
    C4_BK = 0.1, # Closed, 4 Ca²⁺
    O0_BK = 0.1,  # Open, 0 Ca²⁺
    O1_BK = 0.1,  # Open, 1 Ca²⁺
    O2_BK = 0.1,  # Open, 2 Ca²⁺
    O3_BK = 0.1,  # Open, 3 Ca²⁺
    O4_BK = 0   # Open, 4 Ca²⁺ (conducting state)
  )
  
  # Time sequence for the simulation
  times <- seq(0, 200, by = 1)
  
  #history
  
  
  
  # Solve the ODE system
  output <- lsoda(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$P_total <- rowSums(output_df[2:11])
  output_df$P_open <-  output_df$O4_BK 
  output_df$E_K <- (R * Temp / Faraday) * log(K_o / K_i)
  output_df$I_BK <- G_BK * output_df$O4_BK * (parms[["Vm"]] - output_df$E_K)
  output_df$Vm <-  parms[["Vm"]] /  10^-3   # Use Vm between clamp_start and clamp_end
  output_df$Ca_i_free_ID <- parms[["Ca_i_free"]]
  return(output_df)
}

#Parameters to fit -----
BK_test.func <- function(Ca_i_free,BK){
  i <- 1
  mv_tests <- BK$x#seq(-90, 35, by = 5)
  sim_v <- list()
  for(i in seq_along(mv_tests)) {  # Correct looping over mv_tests
    parms <- list(Vm = -mv_tests[i] * 10^-3, Ca_i_free = Ca_i_free * 10^-3)  # Set voltage parameter
    sim_temp <- Model(parms = parms)  # Run simulation
    sim_temp$Vm_identity <- parms[["Vm"]]
    sim_v[[i]] <- as.data.frame(sim_temp)
  }  
  return(sim_v)
}

sim_0.0001 <- BK_test.func(Ca_i_free = 0.0001, BK = BK_0.0001)
sim_0.0003 <- BK_test.func(Ca_i_free = 0.0003, BK = BK_0.0003)
sim_0.001 <- BK_test.func(Ca_i_free = 0.001, BK = BK_0.001)

sim_0.0001_df <- bind_rows(sim_0.0001)  # Convert list to dataframe
sim_0.0003_df <- bind_rows(sim_0.0003)
sim_0.001_df <- bind_rows(sim_0.001)
# Plot all conditions with color mapped to Vm
S7.0.0001 <- ggplot(sim_0.0001_df, aes(x = time, y = P_open, color = factor(Vm_identity), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "100nm",
    x = "Time (msec)",
    y = "Open Probability",
    color = "Voltage (mV)"  # Legend label
  ) +
  xlim(0,200) +
  theme_minimal() 

S7.0.0003 <- ggplot(sim_0.0003_df, aes(x = time, y = P_open, color = factor(Vm_identity), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "300nm",
    x = "Time (msec)",
    y = "Open Probability",
    color = "Voltage (mV)"  # Legend label
  ) +
  xlim(0,200) +
  theme_minimal()

S7.0.001 <- ggplot(sim_0.001_df, aes(x = time, y = P_open, color = factor(Vm_identity), group = Vm)) +
  geom_line(linewidth = 1) +  
  labs(
    title = "1000nm",
    x = "Time (msec)",
    y = "Open Probability",
    color = "Voltage (mV)"  # Legend label
  ) +
  xlim(0,200) +
  theme_minimal()

# Display the plot
grid.arrange(grobs = list(S7.0.001,S7.0.0003,S7.0.0001), ncol = 1, nrow = 3)



# Create dataframe with results
peak.func <- function(df) {
  unique_Vm <- unique(df$Vm)  # Get unique Vm values
  peak_store <- numeric(length(unique_Vm))  # Initialize storage vector
  
  for (i in seq_along(unique_Vm)) {
    sim_temp_peak <- subset(df, Vm == unique_Vm[i])  # Filter by Vm
    peak_store[i] <- mean(sim_temp_peak$O4_BK, na.rm = TRUE)  # Store max O4_BK
  }
  
  # Return results as a data frame (or list if needed)
  return(data.frame(Vm = unique_Vm, Max_O4_BK = peak_store))
}

OV_0.0001 <- peak.func(sim_0.0001_df)

OV_0.0003 <- peak.func(sim_0.0003_df)

OV_0.001 <- peak.func(sim_0.001_df)

combined_plot <- ggplot() +
  # Scatter plots (experimental data)
  geom_point(data = BK_0.0001, aes(x = x, y = y), color = "red", size = 3) +
  geom_point(data = BK_0.0003, aes(x = x, y = y), color = "blue", size = 3) +
  geom_point(data = BK_0.001, aes(x = x, y = y), color = "green", size = 3) +
  
  # Line plots (model I-V curves)
  geom_line(data = OV_0.0001, aes(x = Vm, y = Max_O4_BK), color = "red", linewidth = 1) +
  geom_line(data = OV_0.0003, aes(x = Vm, y = Max_O4_BK), color = "blue", linewidth = 1) +
  geom_line(data = OV_0.001, aes(x = Vm, y = Max_O4_BK), color = "green", linewidth = 1) +
  
  # Labels and theme
  labs(
    title = "Open Probability vs Membrane Voltage for Different Ca_I Concentrations",
    x = "Membrane Voltage (mV)",
    y = "Open Probability"
  ) +
  theme_minimal()  # Clean theme

# Display the plot
print(combined_plot)
