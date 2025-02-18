#Model notes------
#Go through and print each changing variable
#Na_i goes negative at a certain point - investigate variables involved
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
      if (times >= t_period * 0.0 && times <= t_period){
        stim_start <- t_period * 0.0
      }else if(times >= t_period * 1 && times <= t_period * 2){
        stim_start <- t_period * 1.0 
      }else if(times >= t_period * 2 && times <= t_period * 3){
        stim_start <- t_period * 2.0 
      }else if(times >= t_period * 3 && times <= t_period * 4){
        stim_start <- t_period * 3.0 
      }else{
        stim_start <- 0.0
      }
      local_time <- times - stim_start
      
      if (local_time >= 0 && local_time < t_peak_ICC) {
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (local_time / f2)
      } else if (local_time >= t_peak_ICC && local_time < t_plateau_ICC) {
        Vm_ICC <- V_rest_ICC + V_peak_ICC * (1 + exp(-f1 / (2 * t_slope))) * (1 / (1 + exp((local_time - f2 - 0.5 * f1) / t_slope)))
      }
      # Nernst Potentials----
      E_Ca <- (R * Temp / (2 * Faraday)) * log(Ca_o / y["Ca_i_free"])
      E_K <- (R * Temp / Faraday) * log(K_o / y["K_i"])
      E_Na <- (R * Temp / Faraday) * log(Na_o / y["Na_i"])
      # L-type Calcium Current (ICaL)----
      I_CaL <- G_CaL * y["O_LCaL"] * (y["Vm"] - E_Ca)
      
      # Voltage-dependent opening/closing rates
      a_CaL <- 0.7310 * exp(y["Vm"] / 30) * T_correction_Ca_1
      b_CaL <- 0.2149 * exp(-y["Vm"] / 40) * T_correction_Ca_1
      a0_LCaL <- 4 * a_CaL  # S-15
      a1_LCaL <- 3 * a_CaL  # S-16
      a2_LCaL <- 2 * a_CaL  # S-17
      a3_LCaL <- a_CaL  # S-18
      
      b0_LCaL <- b_CaL  # S-19
      b1_LCaL <- 2 * b_CaL  # S-20
      b2_LCaL <- 3 * b_CaL  # S-21
      b3_LCaL <- 4 * b_CaL  # S-22
      # Fast and slow inactivation rates (S-23, S-24)
      phi_f_LCaL <- 0.4742 * exp(y["Vm"] / 10) * T_correction_Ca_1   # Fast inactivation rate (S-23)
      phi_s_LCaL <- 0.005956 * exp(-y["Vm"] / 40) * T_correction_Ca_1  # Slow inactivation rate (S-24)
      
      # Additional inactivation rate equations (S-25, S-26, S-27, S-28)
      xi_f_LCaL <- 0.01407 * exp(-y["Vm"] / 300) * T_correction_Ca_1  # Fast inactivation recovery rate (S-25)
      xi_s_LCaL <- 0.01213 * exp(y["Vm"] / 500)  * T_correction_Ca_1  # Slow inactivation recovery rate (S-26)
      psi_f_LCaL <- 0.02197 * exp(y["Vm"] / 500) * T_correction_Ca_1  # Fast inactivation development rate (S-27)
      psi_s_LCaL <- 0.00232 * exp(-y["Vm"] / 280) * T_correction_Ca_1 # Slow inactivation development rate (S-28)
      
      # Transition probabilities (S-29, S-30)
      omega_f_LCaL <- (b3_LCaL * xi_f_LCaL * phi_f_LCaL) / (a3_LCaL * psi_f_LCaL)  # Fast inactivation transition probability (S-29)
      omega_s_LCaL <- (b3_LCaL * xi_s_LCaL * phi_s_LCaL) / (a3_LCaL * psi_s_LCaL)  # Slow inactivation transition probability (S-30)
      
      # Transition between fast and slow inactivation (S-31, S-32)
      omega_sf_LCaL <- xi_s_LCaL * psi_f_LCaL / xi_f_LCaL  # Transition from fast to slow inactivation (S-31)
      omega_fs_LCaL <- psi_s_LCaL  # Transition from slow to fast inactivation (S-32)
      
      # Calcium-dependent inactivation rate
      if(CalciumDependence == 1){
        theta_LCaL <- 4/(1 + (1 / y["Ca_i_free"]))  * T_correction_Ca_1
      }else{
        theta_LCaL <- 0 #- this is set to 0 to replicate EGTA present
      }
      
      
      # Temp-type Calcium Current (ICaT)----
      I_CaT <- G_CaT * y["d_CaT"] * y["f_CaT"] * (y["Vm"] - E_Ca)
      
      d_CaT_inf <- d_CaT_inf(y["Vm"])
      f_CaT_inf <- f_CaT_inf(y["Vm"])
      
      tau_f_CaT <- 0.38117 * (8.6 + 14.7 * exp(-((y["Vm"] +50)^2/900))) / T_correction_Ca_2
      # Voltage-dependent Potassium Current (IKv)----
      I_Kv <- G_Kv * y["x_Kv"] * y["y_Kv"] * (y["Vm"] - E_K)
      
      x_Kv_inf <- x_Kv_inf(y["Vm"])
      y_Kv_inf <- y_Kv_inf(y["Vm"])
      # Calcium & Voltage-activated Potassium Current (IBK)-----
      I_BK <- G_BK * y["O4_BK"] * (y["Vm"] - E_K)
      a_BK <- exp((0.73 * Faraday * y["Vm"])/(Temp * R))
      b_BK <- exp((-0.67 * Faraday * y["Vm"])/(Temp * R))
      
      #voltage dependent transitions (converted to ms)
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
      KC0C1_BK <- (y["Ca_i_free"] * 4 * k_on)
      KC1C2_BK <- (y["Ca_i_free"] * 3 * k_on)
      KC2C3_BK <- (y["Ca_i_free"] * 2 * k_on)
      KC3C4_BK <- (y["Ca_i_free"] * 1 * k_on)
      
      KC4C3_BK <- (4 * k_c_off)
      KC3C2_BK <- (3 * k_c_off)
      KC2C1_BK <- (2 * k_c_off)
      KC1C0_BK <- (1 * k_c_off)
      
      KO0O1_BK <- (y["Ca_i_free"] * 4 * k_on)
      KO1O2_BK <- (y["Ca_i_free"] * 3 * k_on)
      KO2O3_BK <- (y["Ca_i_free"] * 2 * k_on)
      KO3O4_BK <- (y["Ca_i_free"] * 1 * k_on)
      
      KO4O3_BK <- (4 * k_o_off)
      KO3O2_BK <- (3 * k_o_off)
      KO2O1_BK <- (2 * k_o_off)
      KO1O0_BK <- (1 * k_o_off)
      # Voltage-dependent Sodium Current (INa)-----
      I_Na <- G_Na * y["O_Na"] * (y["Vm"] - E_Na)
      
      KOI1 <- parval[23] * exp(parval[1] + parval[2] * y["Vm"]) * T_correction_Na  # (S-80)
      KI1I2 <- parval[24] * exp(parval[3] + parval[4] * y["Vm"]) * T_correction_Na   # (S-81)
      KC3C2 <- parval[25] * exp(parval[5] + parval[6] * y["Vm"]) * T_correction_Na   # (S-82)
      KC2C1 <- parval[26] * exp(parval[7]+ parval[8] * y["Vm"]) * T_correction_Na  
      KC1O <- parval[27] * exp(parval[9] + parval[10] * y["Vm"]) * T_correction_Na 
      KI2I1 <- parval[28] * exp(parval[11] + parval[12] * y["Vm"]) * T_correction_Na 
      KC2C3 <- parval[29] * exp(parval[13] + parval[14] * y["Vm"]) * T_correction_Na 
      KC1C2 <- parval[30] * exp(parval[15] + parval[16] * y["Vm"]) * T_correction_Na 
      KOC1 <- parval[31] * exp(parval[17] + parval[18] * y["Vm"]) * T_correction_Na 
      KI1C1 <- parval[32] * exp(parval[19] + parval[20] * y["Vm"]) * T_correction_Na 
      KC1I1 <- parval[33] * exp(parval[21] + parval[22] * y["Vm"]) * T_correction_Na 
      KI1O <- parval[34] * exp(parval[35] + parval[36] * y["Vm"]) * T_correction_Na 
      # Sodium-Calcium Exchanger (NCX)------
      I_NCX <- P_NCX * ((exp(gamma * y["Vm"] * Faraday / (R * Temp)) * (y["Na_i"]^3) * Ca_o) -
                          (2.5 * exp(((gamma - 1) * y["Vm"] * Faraday) / (R * Temp)) * (Na_o^3) * y["Ca_i_free"])) /
        ((1 + k_sat * exp((gamma - 1 * y["Vm"] * Faraday) / (R * Temp))) *
           ((K_mNa^3 + Na_o^3) * (K_mCa + Ca_o)))
      
      # Sodium-Potassium Pump (NaK)-----
      I_NaK <- P_NaK * ((K_o * y["Na_i"]) / ((K_mK + K_o) * (K_mNa + y["Na_i"])* (1 + 0.1245 * exp(-((0.1 * y["Vm"] * Faraday)/(R * Temp)))) + 0.0353 * exp(-((y["Vm"] * Faraday)/(R * Temp)))))
      
      # Non-Selective Leak Current (INS)-----
      I_NaNS <- G_NSNa * (y["Vm"] - E_Na)
      I_KNS <- G_NSK * (y["Vm"] - E_K)
      I_NS <- I_NaNS + I_KNS
      
      # Ionic Current-----
      I_ion <- I_CaL + I_CaT + I_Kv + I_BK + I_Na + I_NCX + I_NaK + I_NS
      
      # Ionic Stimulation-----
      I_stim <- Gcouple * (y["Vm"] - Vm_ICC)
      
      #ODE-----
      # BK Channel
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
      
      #L - type channels
      d[11] <- b0_LCaL * y["C1_LCaL"] + sigma_LCaL * y["C0Ca_LCaL"] - (theta_LCaL + a0_LCaL) * y["C0_LCaL"]
      d[12] <- b0_LCaL * y["C1Ca_LCaL"] + theta_LCaL * y["C0_LCaL"] - (a0_LCaL + sigma_LCaL) * y["C0Ca_LCaL"]
      d[13] <- a0_LCaL * y["C0_LCaL"] + b1_LCaL * y["C2_LCaL"] + sigma_LCaL * y["C1Ca_LCaL"] - (b0_LCaL + a1_LCaL + theta_LCaL) * y["C1_LCaL"]
      d[14] <- a0_LCaL * y["C0Ca_LCaL"] + b1_LCaL * y["C2Ca_LCaL"] + theta_LCaL * y["C1_LCaL"] - (b0_LCaL + a1_LCaL + sigma_LCaL) * y["C1Ca_LCaL"]
      d[15] <- a1_LCaL * y["C1_LCaL"] + b2_LCaL * y["C3_LCaL"] + sigma_LCaL * y["C2Ca_LCaL"] - (b1_LCaL + a2_LCaL + theta_LCaL ) * y["C2_LCaL"]
      d[16] <- a1_LCaL * y["C1Ca_LCaL"] + b2_LCaL * y["C3Ca_LCaL"] + theta_LCaL * y["C2_LCaL"] - (b1_LCaL + a2_LCaL + sigma_LCaL) * y["C2Ca_LCaL"]
      d[17] <- a2_LCaL * y["C2_LCaL"] + b3_LCaL * y["O_LCaL"] + sigma_LCaL * y["C3Ca_LCaL"] + omega_s_LCaL * y["IVS_LCaL"] + omega_f_LCaL * y["IVF_LCaL"] - (b2_LCaL + a3_LCaL + phi_s_LCaL + phi_f_LCaL + theta_LCaL) * y["C3_LCaL"]
      d[18] <- a2_LCaL * y["C2Ca_LCaL"] + omega_f_LCaL * y["IVF_Ca_LCaL"]  + theta_LCaL * y["C3_LCaL"] + omega_s_LCaL * y["IVS_Ca_LCaL"] + b3_LCaL * y["ICa_LCaL"] - (b2_LCaL + a3_LCaL + phi_f_LCaL + phi_s_LCaL + sigma_LCaL) * y["C3Ca_LCaL"]
      d[19] <- a3_LCaL * y["C3Ca_LCaL"] + theta_LCaL * y["O_LCaL"] + xi_f_LCaL * y["IVF_Ca_LCaL"] + xi_s_LCaL * y["IVS_Ca_LCaL"] - (b3_LCaL + sigma_LCaL + psi_f_LCaL + psi_s_LCaL) * y["ICa_LCaL"] # Calcium dependent
      d[20] <- psi_f_LCaL * y["O_LCaL"] + phi_f_LCaL * y["C3_LCaL"] + sigma_LCaL * y["IVF_Ca_LCaL"] + omega_sf_LCaL * y["IVS_LCaL"] - (xi_f_LCaL + theta_LCaL + omega_fs_LCaL + omega_f_LCaL) * y["IVF_LCaL"] # Voltage dependent - fast 
      d[21] <- phi_f_LCaL * y["C3Ca_LCaL"] + psi_f_LCaL * y["ICa_LCaL"] + theta_LCaL * y["IVF_LCaL"] + omega_sf_LCaL * y["IVS_Ca_LCaL"] - (sigma_LCaL + xi_f_LCaL + omega_f_LCaL + omega_fs_LCaL) * y["IVF_Ca_LCaL"] # fast voltage-dependent Calcium-dependent
      d[22] <- phi_s_LCaL * y["C3_LCaL"] + sigma_LCaL * y["IVS_Ca_LCaL"] + omega_fs_LCaL * y["IVF_LCaL"] + psi_s_LCaL * y["O_LCaL"] - (omega_s_LCaL + xi_s_LCaL + omega_sf_LCaL + theta_LCaL) * y["IVS_LCaL"] # Voltage-dependent - slow
      d[23] <- phi_s_LCaL * y["C3Ca_LCaL"] + psi_s_LCaL * y["ICa_LCaL"] + theta_LCaL * y["IVS_LCaL"] + omega_fs_LCaL * y["IVF_Ca_LCaL"] - (sigma_LCaL + omega_s_LCaL + xi_s_LCaL + omega_sf_LCaL) * y["IVS_Ca_LCaL"] # slow voltage dependent Calcium dependent
      d[24] <- a3_LCaL * y["C3_LCaL"] + sigma_LCaL * y["ICa_LCaL"]  + xi_s_LCaL * y["IVS_LCaL"] + xi_f_LCaL * y["IVF_LCaL"] - (b3_LCaL + theta_LCaL + psi_s_LCaL + psi_f_LCaL) * y["O_LCaL"]
      
      # Sodium Channel ODE
      d[25] <- KOC1 * y["O_Na"] + KI1C1 * y["I1_Na"] + KC2C1 * y["C2_Na"] - (KC1O + KC1C2 + KC1I1) * y["C1_Na"] # Closed state 1
      d[26] <- KC1C2 * y["C1_Na"] + KC3C2 * y["C3_Na"] - (KC2C1 + KC2C3) * y["C2_Na"] # Closed state 2
      d[27] <- KC2C3 * y["C2_Na"] - KC3C2 * y["C3_Na"] # Closed state 3
      d[28] <- KI2I1 * y["I2_Na"] + KOI1 * y["O_Na"] + KC1I1 * y["C1_Na"] - (KI1I2 + KI1O + KI1C1) * y["I1_Na"]  # Inactivated State 1
      d[29] <- KI1I2 * y["I1_Na"] - KI2I1 * y["I2_Na"] # Inactivated 2
      d[30] <- KI1O * y["I1_Na"] + KC1O * y["C1_Na"] - (KOI1 + KOC1) * y["O_Na"]  # Open 
      
      # Temp-type Calcium Channel ODEs
      d[31] <- (d_CaT_inf - y["d_CaT"]) / tau_d_CaT
      d[32] <- (f_CaT_inf - y["f_CaT"]) / tau_f_CaT
      
      #Ion concentration tracked in mM
      # Intracellular Calcium Concentration - d[Ca_i_total]
      d[33] <- (- (I_CaL + I_CaT - I_NCX) / (2 * Faraday * Vcell))
      # Intracellular Potassium Concentration - d[K_i]
      d[34] <- (- (I_Kv + I_BK + I_stim - I_NaK - I_KNS) / (Faraday * Vcell))
      # Intracellular Sodium Concentration - d[Na_i]
      d[35] <- (- (I_Na - I_NCX - I_NaK - I_NaNS) / (Faraday * Vcell))

      #Single hJSMC electrophysiology -dVm
      d[36] <- (I_ion + I_stim) / Cm
      # Voltage-dependent Potassium Channel ODEs
      d[37] <- (x_Kv_inf - y["x_Kv"]) / tau_x_Kv
      d[38] <- (y_Kv_inf - y["y_Kv"]) / tau_y_Kv
      
      
      
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
  Texp <- 297
  Cm <- 50  # pF - Cell membrane Capacitance
  Vcell <- 3.5e-12  # l - Cell volume
  Ca_o <- 2  # mM - Extracellular Calcium concentration
  K_o <- 5.4  # mM - Extracellular potassium concentration
  Na_o <- 140  # mM - Extracellular sodium concentration
  Ca_i_total <- 0.004914  # mM - Initial total intracellular Calcium concentration
  Ca_i_free <- 1.26e-4  # mM - Initial free intracellular Calcium concentration
  K_i <- 150  # mM - Intracellular potassium concentration
  Na_i <- 10.5  # mM - Intracellular sodium concentration
  Q10Ca <- 2.1  # - Q10 for L type channels
  Q10K <- 3.1  # - Q10 for potassium channels
  Q10Na <- 2.45  # - Q10 for sodium channels
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
  k_on <- 40633 #BK on rate
  k_c_off <- 11 #BK closed off rate
  k_o_off <- 1.1 #BK open off rate
  SMC_resting <- -60 #
  T_correction_Ca_1 <- 1.0 * Q10Ca^((Temp - 310) / 10.0) # experimental correction
  T_correction_Ca_2 <- 1.0 * Q10Ca^((Temp - Texp) / 10.0) # experimental correction
  T_correction_Na <- 1.0 * Q10Na^((Temp - Texp) / 10.0) # experimental correction
  T_correction_K <- 1.0 * Q10K^((Temp - Texp) / 10.0) # experimental correction
  
  if(parms[["CalciumDependence"]] == 1){
    sigma_LCaL <- 0.01 * T_correction_Ca_1 #FOR ICaL set to 0 to replicate EGTA 
  }else{
    sigma_LCaL <- 0 #- this is set to 0 to replicate EGTA present
  }
  tau_d_CaT <- 1.9058 / T_correction_Ca_2 #for ICaL
  tau_x_Kv <- 4.7803 / T_correction_K  # ms - Time constant for Kvx
  tau_y_Kv <- 763.7564 / T_correction_K  # ms - Time constant for Kvy
  #Function for CaT
  d_CaT_inf <- function(Vm) {
    return(1 / (1 + exp(-((Vm + 60.5) / 5.3))))
  }
  
  f_CaT_inf <- function(Vm) {
    return(1 / (1 + exp(((Vm + 75.5) / 4.0))))
  }
  #Function for IKV
  x_Kv_inf <- function(Vm) {
    return(1/ (1 + exp(-((Vm + 43)/17.36))))
  }
  
  y_Kv_inf <- function(Vm) {
    return(1/ (1 + exp(((Vm + 44.9)/12.0096))))
  }
  #Na channels parameters from code
  parval <- c(0.307629865993612,0.00605347675235982,0.0514901233349612,-0.0468650997035773,-0.0691017278134239,0.00319450166548765,
              -0.156598571141704,0.0583527292113407,0.00931925729714207,0.0410752869450889,2.67926519075303,0.00614677585983095,
              -0.0990741811465619,0.0364411841506491,0.363520774026992,0.0771931771113635,-13.3349821299039,-0.252888485042226,
              -2.48395798104711,0.0204061599230419,-0.0634382282818863,0.00466834101392052,1.61639957785747,0.0277351044118711,
              0.000525484598535321,1.4496348286393,1.53287364507221,0.00392388356048677,0.554315057682815,3.15656890165025,
              2.39152928164775,1.90461121472293,0.00021688112860783,0.120515381956888,-9.60280306718205,0.0830250234337321 )
  
  #BK has special units?
  
  # Initial values for A, B, and C
  initial <- c(
    C0_BK = 0.48379087935899,   # Closed, 0 Ca²⁺
    C1_BK = 0.385183559520031,  # Closed, 1 Ca²⁺
    C2_BK = 0.115002824567753,  # Closed, 2 Ca²⁺
    C3_BK = 0.0152602714149774, # Closed, 3 Ca²⁺
    C4_BK = 0.000759264410974374, # Closed, 4 Ca²⁺
    O0_BK = 6.94960798375172e-7,  # Open, 0 Ca²⁺
    O1_BK = 5.55636826398253e-8,  # Open, 1 Ca²⁺
    O2_BK = 2.85143702125325e-8,  # Open, 2 Ca²⁺
    O3_BK = 1.59832806123435e-6,  # Open, 3 Ca²⁺
    O4_BK = 1.82113764497095e-6,  # Open, 4 Ca²⁺ (conducting state)
    C0_LCaL = 0.815464741971086,
    C0Ca_LCaL = 0.0175888495282545,
    C1_LCaL = 0.152399266235657,
    C1Ca_LCaL = 0.00328711668724504,
    C2_LCaL = 0.0106805060777161,
    C2Ca_LCaL = 0.000230369020877669,
    C3_LCaL = 0.000332673548872087,
    C3Ca_LCaL = 7.1754726923539e-6,
    ICa_LCaL = 8.38123983500905e-8,
    IVF_LCaL = 4.0998751301597e-6,
    IVF_Ca_LCaL = 8.84306615061238e-8,
    IVS_LCaL = 1.1193313274705e-6,
    IVS_Ca_LCaL = 2.41429816075123e-8,
    O_LCaL = 3.88576045134351e-6,
    C1_Na = 0.0119443135223679,
    C2_Na = 0.0109545368437155,
    C3_Na = 0.973782548650071,
    I1_Na = 0.000126882921013389,
    I2_Na = 0.00318975045717667,
    O_Na = 1.96760342050475e-6,
    d_CaT = 0.0791635737410974, 
    f_CaT = 0.377831534375835,
    Ca_i_free = Ca_i_free,  # mM - Initial free intracellular Calcium concentration
    K_i = K_i,  # mM - Initial intracellular potassium concentration
    Na_i = Na_i,  # mM - Initial intracellular sodium concentration
    Vm = -73.5049651455872,  # Membrane potential - millivolts
    x_Kv = 0.14714161078933,
    y_Kv = 0.99994773314105
  )
  
  
  # Time sequence for the simulation ms
  times <- seq(0, t_period, by = 1)
  
  #history
  
  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$Ca_i_free_nm <- output_df$Ca_i_free * 10^6
  return(output_df)
}

#Parameters to fit -----
parms <- list(
  CalciumDependence = 1
)



#initial estimate -----
Initial_out <- Model(parms = parms)


Voltage_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time / 1000, y = Vm, colour = "Vm")) +
  scale_colour_manual(values = c("Vm" = "black")) +
  xlim(1760,1800) +
  xlab("time (s)") +
  ylab("Voltage (mV)") +
  ggtitle("hJMC voltage simulation")

Ca_i_free_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time / 1000, y = Ca_i_free_nm, colour = "Ca_i_free")) +
  scale_colour_manual(values = c("Ca_i_free" = "black")) +
  xlim(1760,1800) +
  xlab("time (s)") +
  ylab("Free Intracellular\n calcium (nM)") +
  ggtitle("hJMC Free Intracellular\n calcium (nM) simulation")

print(Voltage_plot)
print(Ca_i_free_plot)
grid.arrange(grobs = list(Voltage_plot, Ca_i_free_plot), ncol = 1, nrow = 2)