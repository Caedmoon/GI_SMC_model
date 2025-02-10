#Model notes------
#Based on Poh et al (2012)
#Testing summaritive performance then investigate the individual ion current predictions using Fig 1 2 3 etc - 
#Create additional scripts to solve for each voltage supplied 
#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")
#Import observed data----
# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)

      # Rate equations
      
      #Single slow wave-----
      if (times >= 0 & times < t_peak_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (times/f2)
      } 
      if (times >= t_peak_ICC & times < t_plateau_ICC){
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (1 + exp(-f1/(2 * t_slope))) * (1/(1 + exp((times - f2 - 0.5 * f1)/t_slope)))
      } 
      
      # Nernst Potentials----
      E_Ca <- (R * Temp / (2 * Faraday)) * log(y["Ca_o"] / y["Ca_i_free"])
      E_K <- (R * Temp / Faraday) * log(y["K_o"] / y["K_i"])
      E_Na <- (R * Temp / Faraday) * log(y["Na_o"] / y["Na_i"])
      
      # L-type Calcium Current (ICaL)----
      # Voltage-dependent opening/closing rates
      a_CaL <- 0.7310 * exp(y["Vm"] / 30)
      b_CaL <- 0.2149 * exp(-y["Vm"] / 40)
      a0_LCaL <- 4 * a_CaL  # S-15
      a1_LCaL <- 3 * a_CaL  # S-16
      a2_LCaL <- 2 * a_CaL  # S-17
      a3_LCaL <- a_CaL  # S-18
      
      b0_LCaL <- b_CaL  # S-19
      b1_LCaL <- 2 * b_CaL  # S-20
      b2_LCaL <- 3 * b_CaL  # S-21
      b3_LCaL <- 4 * b_CaL  # S-22
      # Fast and slow inactivation rates (S-23, S-24)
      phi_f_LCaL <- 0.4742 * exp(y["Vm"] / 10)   # Fast inactivation rate (S-23)
      phi_s_LCaL <- 0.05956 * exp(-y["Vm"] / 40)  # Slow inactivation rate (S-24)
      
      # Additional inactivation rate equations (S-25, S-26, S-27, S-28)
      xi_f_LCaL <- 0.01407 * exp(-y["Vm"] / 300)  # Fast inactivation recovery rate (S-25)
      xi_s_LCaL <- 0.01213 * exp(y["Vm"] / 500)   # Slow inactivation recovery rate (S-26)
      psi_f_LCaL <- 0.02197 * exp(y["Vm"] / 500)  # Fast inactivation development rate (S-27)
      psi_s_LCaL <- 0.00232 * exp(-y["Vm"] / 280) # Slow inactivation development rate (S-28)
      
      # Transition probabilities (S-29, S-30)
      omega_f_LCaL <- (b3_LCaL * xi_f_LCaL * phi_f_LCaL) / (a3_LCaL * psi_f_LCaL)  # Fast inactivation transition probability (S-29)
      omega_s_LCaL <- (b3_LCaL * xi_s_LCaL * phi_s_LCaL) / (a3_LCaL * psi_s_LCaL)  # Slow inactivation transition probability (S-30)
      
      # Transition between fast and slow inactivation (S-31, S-32)
      omega_sf_LCaL <- xi_s_LCaL * psi_f_LCaL / xi_f_LCaL  # Transition from fast to slow inactivation (S-31)
      omega_fs_LCaL <- psi_s_LCaL  # Transition from slow to fast inactivation (S-32)
      
      # Calcium-dependent inactivation rate
      theta_LCaL <- 4/(1 + (1 / Ca_i_free))
      
      # Temp-type Calcium Current (ICaT)----
      I_CaT <- G_CaT * y["d_CaT"] * y["f_CaT"] * (y["Vm"] - E_Ca)
      
      d_CaT_inf <- 1/ (1 + exp(-((y["Vm"] + 60.5)/5.3)))
      f_CaT_inf <- 1/ (1 + exp(((y["Vm"] + 75.5)/4.0)))
      
      tau_f_CaT <- 0.38117 * (8.6 + 14.7 * exp(-((y["Vm"] +50)^2/900)))
      # Voltage-dependent Potassium Current (IKv)----
      I_Kv <- G_Kv * y["x_Kv"] * y["y_Kv"] * (y["Vm"] - E_K)
      
      x_Kv_inf <- 1/ (1 + exp(-((y["Vm"] + 43)/17.36)))
      y_Kv_inf <- 1/ (1 + exp(((y["Vm"] + 44.9)/12.0096)))
      # Calcium & Voltage-activated Potassium Current (IBK)-----
      I_BK <- G_BK * y["O4_BK"] * (y["Vm"] - E_K)
      a_BK <- exp((8.47188 * y["Vm"])/Temp)
      b_BK <- exp((-7.77556 * y["Vm"])/Temp)
      # Voltage-dependent Sodium Current (INa)-----
      I_Na <- G_Na * y["O_Na"] * (y["Vm"] - E_Na)
      
      KOI1 <- 1.6164 * exp(0.30763 + 0.0060535 * y["Vm"])  # (S-80)
      KI1I2 <- 0.027735 * exp(0.051490 - 0.046865 * y["Vm"])  # (S-81)
      KC3C2 <- 0.00052548 * exp(-0.069102 + 0.0031945 * y["Vm"])  # (S-82)
      KC2C1 <- 1.4496 * exp(-0.15660 + 0.058353 * y["Vm"]) 
      KC1O <- 1.5329 * exp(0.0093193 + 0.041075 * y["Vm"])
      KI2I1 <- 0.0039239 * exp(2.6793 + 0.0061468 * y["Vm"])
      KC2C3 <- 0.55432 * exp(-0.099074 + 0.036441 * y["Vm"])
      KC1C2 <- 3.1566 * exp(0.36352 + 0.077193 * y["Vm"])
      KOC1 <- 2.3915 * exp(-13.335 - 0.25289 * y["Vm"])
      KI1C1 <- 1.9046 * exp(-2.4840 + 0.020406 * y["Vm"])
      KC1I1 <- 0.00021688 * exp(-0.063438 + 0.0046683 * y["Vm"])
      KI1O <- 0.12052 * exp(-9.6028 + 0.083025 * y["Vm"])
      # Sodium-Calcium Exchanger (NCX)------
      I_NCX <- P_NCX * ((exp(gamma * y["Vm"] * Faraday / (R * Temp)) * (y["Na_i"]^3) * y["Ca_o"]) -
                          (2.5 * exp(((gamma - 1) * y["Vm"] * Faraday) / (R * Temp)) * (y["Na_o"]^3) * y["Ca_i_free"])) /
        ((1 + k_sat * exp((gamma - 1 * y["Vm"] * Faraday) / (R * Temp))) *
           ((K_mNa^3 + y["Na_o"]^3) * (K_mCa + y["Ca_o"])))
      
      # Sodium-Potassium Pump (NaK)-----
      I_NaK <- P_NaK * ((y["K_o"] * y["Na_i"]) / ((K_mK + y["K_o"]) * (K_mNa + y["Na_i"])* (1 + 0.1245 * exp(-((0.1 * y["Vm"] * Faraday)/(R * Temp)))) + 0.0353 * exp(-((y["Vm"] * Faraday)/(R * Temp)))))
      
      # Non-Selective Leak Current (INS)-----
      I_NaNS <- G_NSNa * (y["Vm"] - E_Na)
      I_KNS <- G_NSK * (y["Vm"] - E_K)
      I_NS <- I_NaNS + I_KNS
      
      # Ionic Current-----
      I_ion <- I_CaL + I_CaT + I_Kv + I_BK + I_Na + I_NCX + I_NaK + I_NS
      
      # Ionic Stimulation-----
      I_stim <- Gcouple * (y["Vm"] - Vm_ICC)
      
      #ODE-----
      #Single hJSMC electrophysiology -dVm
      d[1] <- (I_ion + I_stim) / Cm
      #Ion concentration tracked in mM
      # Intracellular Calcium Concentration - d[Ca_i_total]
      d[2] <- (- (I_CaL + I_CaT - I_NCX) / (2 * Faraday * Vcell))
      # Intracellular Sodium Concentration - d[Na_i]
      d[3] <- (- (I_Na - I_NCX - I_NaK - I_NaNS) / (Faraday * Vcell))
      # Intracellular Potassium Concentration - d[K_i]
      d[4] <- (- (I_Kv + I_BK + I_stim - I_NaK - I_KNS) / (Faraday * Vcell))
      
      # Calcium Buffering ODEs -----
      d[5] <- y["Ca_i_total"] / (1 + ((CRT_n * CRT_total * CRT_KD * y["Ca_i_free"]^CRT_n)/(y["Ca_i_free"]^CRT_n + CRT_KD)^2) + ((CaM_n * CaM_total * CaM_KD * y["Ca_i_free"]^CaM_n)/(y["Ca_i_free"]^CaM_n + CaM_KD)^2)) # d[Ca_i_free] - Free intracellular Calcium
      
      # Closed states without Calcium binding
      d[6] <- b0_LCaL * y["C1_I_LCaL"] + sigma_LCaL * y["C0Ca_I_LCaL"] - (theta_LCaL + a0_LCaL) * y["C0_I_LCaL"]
      d[7] <- a0_LCaL * y["C0_I_LCaL"] + b1_LCaL * y["C2_I_LCaL"] + sigma_LCaL * y["C1Ca_I_LCaL"] - (b0_LCaL + a1_LCaL + theta_LCaL) * y["C1_I_LCaL"]
      d[8] <- a1_LCaL * y["C1_I_LCaL"] + b2_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["C2Ca_I_LCaL"] - (b1_LCaL + a2_LCaL + theta_LCaL ) * y["C2_I_LCaL"]
      d[9] <- a2_LCaL * y["C2_I_LCaL"] + b3_LCaL * y["O_I_LCaL"] + sigma_LCaL * y["C3Ca_I_LCaL"] + omega_s_LCaL * y["IVS_I_LCaL"] + omega_f_LCaL * y["IVF_I_LCaL"] - (b2_LCaL + a3_LCaL + phi_s_LCaL + phi_f_LCaL + theta_LCaL) * y["C3_I_LCaL"]
      
      # Closed states with Calcium binding
      d[10] <- b0_LCaL * y["C1Ca_I_LCaL"] + theta_LCaL * y["C0_I_LCaL"] - (a0_LCaL + sigma_LCaL) * y["C0Ca_I_LCaL"]
      d[11] <- a0_LCaL * y["C0Ca_I_LCaL"] + b1_LCaL * y["C2Ca_I_LCaL"] + theta_LCaL * y["C1_I_LCaL"] - (b0_LCaL + a1_LCaL + sigma_LCaL) * y["C1Ca_I_LCaL"]
      d[12] <- a1_LCaL * y["C1Ca_I_LCaL"] + b2_LCaL * y["C3Ca_I_LCaL"] + theta_LCaL * y["C2_I_LCaL"] - (b1_LCaL + a2_LCaL + sigma_LCaL) * y["C2Ca_I_LCaL"]
      d[13] <- a2_LCaL * y["C2Ca_I_LCaL"] + omega_f_LCaL * y["IVF_Ca_I_LCaL"]  + theta_LCaL * y["C3_I_LCaL"] + omega_s_LCaL * y["IVS_Ca_I_LCaL"] + b3_LCaL * y["ICa_I_LCaL"] - (b2_LCaL + a3_LCaL + phi_f_LCaL + phi_s_LCaL + sigma_LCaL) * y["C3Ca_I_LCaL"]
      
      # Open state
      d[14] <- a3_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["ICa_I_LCaL"]  + xi_s_LCaL * y["IVS_I_LCaL"] + xi_f_LCaL * y["IVF_I_LCaL"] - (b3_LCaL + theta_LCaL + psi_s_LCaL + psi_f_LCaL) * y["O_I_LCaL"]
      
      # Inactivation states
      d[15] <- psi_f_LCaL * y["O_I_LCaL"] + phi_f_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["IVF_Ca_I_LCaL"] + omega_sf_LCaL * y["IVS_I_LCaL"] - (xi_f_LCaL + theta_LCaL + omega_fs_LCaL + omega_f_LCaL) * y["IVF_I_LCaL"] # Voltage dependent - fast 
      d[16] <- phi_s_LCaL * y["C3_I_LCaL"] + sigma_LCaL * y["IVS_Ca_I_LCaL"] + omega_fs_LCaL * y["IVF_I_LCaL"] + psi_s_LCaL * y["O_I_LCaL"] - (omega_s_LCaL + xi_s_LCaL + omega_sf_LCaL + theta_LCaL) * y["IVS_I_LCaL"] # Voltage-dependent - slow
      d[17] <- a3_LCaL * y["C3Ca_I_LCaL"] + theta_LCaL * y["O_I_LCaL"] + xi_f_LCaL * y["IVF_Ca_I_LCaL"] + xi_s_LCaL * y["IVS_Ca_I_LCaL"] - (b3_LCaL + sigma_LCaL + psi_f_LCaL + psi_s_LCaL) * y["ICa_I_LCaL"] # Calcium dependent
      d[18] <- phi_f_LCaL * y["C3Ca_I_LCaL"] + psi_f_LCaL * y["ICa_I_LCaL"] + theta_LCaL * y["IVF_I_LCaL"] + omega_sf_LCaL * y["IVS_Ca_I_LCaL"] - (sigma_LCaL + xi_f_LCaL + omega_f_LCaL + omega_fs_LCaL) * y["IVF_Ca_I_LCaL"] # fast voltage-dependent Calcium-dependent
      d[19] <- phi_s_LCaL * y["C3Ca_I_LCaL"] + psi_s_LCaL * y["ICa_I_LCaL"] + theta_LCaL * y["IVS_I_LCaL"] + omega_fs_LCaL * y["IVF_Ca_I_LCaL"] - (sigma_LCaL + omega_s_LCaL + xi_s_LCaL + omega_sf_LCaL) * y["IVS_Ca_I_LCaL"] # slow voltage dependent Calcium dependent
      # Temp-type Calcium Channel ODEs-----
      d[20] <- (d_CaT_inf - y["d_CaT"]) / tau_d_CaT
      d[21] <- (f_CaT_inf - y["f_CaT"]) / tau_f_CaT
      
      # Voltage-dependent Potassium Channel ODEs
      d[22] <- (x_Kv_inf - y["x_Kv"]) / tau_x_Kv
      d[23] <- (y_Kv_inf - y["y_Kv"]) / tau_y_Kv
      
      # BK Channel
      #Closed
      d[24] <- k_c_off * y["Ca_i_free"] * y["C1_BK"] + 318.1084 * b_BK * y["O0_BK"] - (4 * k_on * y["Ca_i_free"] + 0.02162 * a_BK) * y["C0_BK"]  # closed state of BK channel 0
      d[25] <- 4 * k_on * y["Ca_i_free"] * y["C0_BK"] + 2 * k_c_off * y["Ca_i_free"] * y["C2_BK"] + 144.1736 * b_BK * y["O1_BK"] - (k_c_off * y["Ca_i_free"] + 0.000869 * a_BK + 3 * k_on * y["Ca_i_free"]) * y["C1_BK"]  # closed state of BK channel 1
      d[26] <- 3 * k_on * y["Ca_i_free"] * y["C1_BK"] + 3 * k_c_off * y["Ca_i_free"] * y["C3_BK"] + 32.6594 * b_BK * y["O2_BK"] - (2 * k_c_off * y["Ca_i_free"] + 0.0000281 * a_BK + 2 * k_on * y["Ca_i_free"]) * y["C2_BK"]  # closed state of BK channel 2
      d[27] <- 2 * k_on * y["Ca_i_free"] * y["C2_BK"] + 4 * k_c_off * y["Ca_i_free"] * y["C4_BK"] + 0.095312 * b_BK * y["O3_BK"] - (3 * k_c_off * y["Ca_i_free"] + 0.000781 * a_BK + k_on * y["Ca_i_free"]) * y["C3_BK"] # closed state of BK channel 3
      d[28] <- k_on * y["Ca_i_free"] * y["C3_BK"] + 0.000106 * b_BK * y["O4_BK"] - (4 * k_c_off * y["Ca_i_free"] + 0.044324 * a_BK) * y["C4_BK"]  # closed state of BK channel 4
      #Open
      d[29] <- k_o_off * y["Ca_i_free"] * y["O1_BK"] + 0.02162 * a_BK * y["C0_BK"] - (4 * k_on * y["Ca_i_free"] + 318.1084 * b_BK) * y["O0_BK"]  # Open state of BK channel 0
      d[30] <- 4 * k_on * y["Ca_i_free"] * y["O0_BK"] + 2 * k_o_off * y["Ca_i_free"] * y["O2_BK"] + 0.000869 * a_BK * y["C1_BK"] - (k_o_off * y["Ca_i_free"] + 144.1736 * b_BK + 3 * k_on * y["Ca_i_free"]) * y["O1_BK"]  # Open state of BK channel 1
      d[31] <- 3 * k_on * y["Ca_i_free"] * y["O1_BK"] + 3 * k_o_off * y["Ca_i_free"] * y["O3_BK"] + 0.0000281 * a_BK * y["C2_BK"] - (2 * k_o_off * y["Ca_i_free"] + 32.6594 * b_BK + 2 * k_on * y["Ca_i_free"]) * y["O2_BK"] # Open state of BK channel 2
      d[32] <- 2 * k_on * y["Ca_i_free"] * y["O2_BK"] + 4 * k_o_off * y["Ca_i_free"] * y["O4_BK"] + 0.000781 * a_BK * y["C3_BK"] - (3 * k_o_off * y["Ca_i_free"] + 0.095312 * b_BK + k_on * y["Ca_i_free"]) * y["O3_BK"]  # Open state of BK channel 3
      d[33] <- k_on * y["Ca_i_free"] * y["O3_BK"] + 0.044324 * a_BK * y["C4_BK"] - (4 * k_o_off * y["Ca_i_free"] + 0.000106 * b_BK) * y["O4_BK"]  # Open state of BK channel 4
      # Sodium Channel ODE
      d[34] <- KI2I1 * y["I2_Na"] + KOI1 * y["O_Na"] + KC1I1 * y["C1_Na"] - (KI1I2 + KI1O + KI1C1) * y["I1_Na"]  # Inactivated State 1
      d[35] <- KI1I2 * y["I1_Na"] - KI2I1 * y["I2_Na"] # Inactivated 2
      d[36] <- KI1O * y["O_Na"] + KC1O * y["C1_Na"] - (KOI1 + KOC1) * y["O_Na"]  # Open 
      d[37] <- KOC1 * y["O_Na"] + KI1C1 * y["I1_Na"] + KC2C1 * y["C2_Na"] - (KC1O + KC2C1 + KI1C1) * y["C1_Na"] # Closed state 1
      d[38] <- KC1C2 * y["C1_Na"] + KC3C2 * y["C3_Na"] - (KC2C1 + KC2C3) * y["C2_Na"] # Closed state 2
      d[39] <- KC2C3 * y["C2_Na"] - KC3C2 * y["C3_Na"] # Closed state 3
      # Return the rates of change
      return(list(d))
    })
  }
  
  #Fixed parameters
  #Fixed parameters
  fixed <- list()
  R <- 8.314  # J/(molK) - Ideal gas constant
  Faraday <- 96.48534  # C/mmol - Faraday’s constant
  Temp <- 310  # K - Temperature
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
  sigma_LCaL <- 0.01 #FOR ICaL
  k_on <- 40633 #BK on rate
  k_c_off <- 11 #BK closed off rate
  k_o_off <- 1.1 #BK open off rate
  SMC_resting <- -60 #
    # Initial values for A, B, and C
    initial <- c(
      Vm = SMC_resting,  # Initial membrane potential
      Ca_i_total = Ca_i_total,  # mM - Initial total intracellular Calcium concentration
      Na_i = 10.5,  # mM - Initial intracellular sodium concentration
      K_i = 150,  # mM - Initial intracellular potassium concentration
      Ca_i_free = Ca_i_free,  # mM - Initial free intracellular Calcium concentration
      C0_I_LCaL = 0.923,
      C1_I_LCaL = 0.049,
      C2_I_LCaL = 0,
      C3_I_LCaL = 0,
      C0Ca_I_LCaL = 0.025,
      C1Ca_I_LCaL = 0.003,
      C2Ca_I_LCaL = 0,
      C3Ca_I_LCaL = 0,
      O_I_LCaL = 0,
      IVF_I_LCaL = 0,
      IVS_I_LCaL = 0,
      ICa_I_LCaL = 0,
      IVF_Ca_I_LCaL = 0,
      IVS_Ca_I_LCaL = 0,
      d_CaT = 1/ (1 + exp(-((SMC_resting + 60.5)/5.3))), #maybe should be -70mV
      f_CaT = 1/ (1 + exp(((SMC_resting + 75.5)/4.0))),
      x_Kv = 1/ (1 + exp(-((SMC_resting + 43)/17.36))),
      y_Kv = 1/ (1 + exp(((SMC_resting + 44.9)/12.0096))),
      C0_BK = 1,
      C1_BK = 0,
      C2_BK = 0,
      C3_BK = 0,
      C4_BK = 0,
      O0_BK = 0,
      O1_BK = 0,
      O2_BK = 0,
      O3_BK = 0,
      O4_BK = 0,
      I1_Na = 0,
      I2_Na = 0,
      O_Na = 0,
      C1_Na = 1,
      C2_Na = 0,
      C3_Na = 0
    )
  
  # Time sequence for the simulation
  times <- seq(0, 1800, by = 1)
  
  #history

  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$E_Ca <- ((R * Temp) / (2 * Faraday)) * log(Ca_o / Ca_i_total)
  output_df$I_CaL <- G_CaL * output_df$O_I_LCaL * (Vm - output_df$E_Ca) # check equations again error is here? - maybe in P_O
  return(output_df)
}

#Parameters to fit -----
parms <- list(
)



#initial estimate -----
Initial_out <- Model(parms = parms)


Voltage_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time, y = Vm, colour = "Vm")) +
  sCale_colour_manual(values = c("Vm" = "black")) +
  xlim(1760,1800) +
  xlab("time (s)") +
  ylab("Voltage (mV)") +
  ggtitle("hJMC voltage simulation")

Ca_i_free_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time, y = Ca_i_free, colour = "Ca_i_free")) +
  sCale_colour_manual(values = c("Ca_i_free" = "black")) +
  xlim(1760,1800) +
  xlab("time (s)") +
  ylab("Free Intracellular\n calcium (nM)") +
  ggtitle("hJMC Free Intracellular\n calcium (nM) simulation")

print(Voltage_plot)
print(Ca_i_free_plot)
