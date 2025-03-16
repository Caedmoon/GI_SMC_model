#Model notes------
#Something wrong with Calcium channels. keep going over C version
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
      stim_cycle <- times %/% period
      stim_start <- stim_cycle * period
      local_time <- times - stim_start
      #print(local_time)
      if (local_time >= 0 && local_time < t_ICCpeak) {
        Vm_ICC <- V_ICCrest + V_ICCamp * local_time / f_2
      } else if (local_time >= t_ICCpeak && local_time < t_ICCplateau) {
        Vm_ICC <- V_ICCrest + V_ICCamp * (1 + exp(-f_1 / (2 * t_slope))) * 1 / (1 + exp((local_time - f_2 - 0.5 * f_1) / t_slope))
      }
      else{
        Vm_ICC <- V_ICCrest
      }
      I_stim <- Gcouple * (y["Vm"] - Vm_ICC)
      # Nernst Potentials----
      E_Ca <- (R * Temp / (2 * Faraday)) * log(Cao / y["Ca_i_free"])
      E_K <- (R * Temp / Faraday) * log(Ko / y["K_i"])
      E_Na <- (R * Temp / Faraday) * log(Nao / y["Na_i"])
      # L-type Calcium Current (ICaL)----
      I_CaL <- g_CaL * y["O_LCaL"] * (y["Vm"] - E_Ca)
      
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
      phi_f <- 0.4742 * exp(y["Vm"] / 10) * T_correction_Ca_1   # Fast inactivation rate (S-23)
      phi_s <- 0.005956 * exp(-y["Vm"] / 40) * T_correction_Ca_1  # Slow inactivation rate (S-24)
      
      # Additional inactivation rate equations (S-25, S-26, S-27, S-28)
      xi_f <- 0.01407 * exp(-y["Vm"] / 300) * T_correction_Ca_1  # Fast inactivation recovery rate (S-25)
      xi_s <- 0.01213 * exp(y["Vm"] / 500)  * T_correction_Ca_1  # Slow inactivation recovery rate (S-26)
      psi_f <- 0.02197 * exp(y["Vm"] / 500) * T_correction_Ca_1  # Fast inactivation development rate (S-27)
      psi_s <- 0.00232 * exp(-y["Vm"] / 280) * T_correction_Ca_1 # Slow inactivation development rate (S-28)
      
      # Transition probabilities (S-29, S-30)
      omega_f <- (b3_LCaL * xi_f * phi_f) / (a3_LCaL * psi_f)  # Fast inactivation transition probability (S-29)
      omega_s <- (b3_LCaL * xi_s * phi_s) / (a3_LCaL * psi_s)  # Slow inactivation transition probability (S-30)
      
      # Transition between fast and slow inactivation (S-31, S-32)
      omega_sf <- xi_s * psi_f / xi_f  # Transition from fast to slow inactivation (S-31)
      omega_fs <- psi_s  # Transition from slow to fast inactivation (S-32)
      
      # Calcium-dependent inactivation rate
      if(CalciumDependence == 1){
        theta_LCaL <- T_correction_Ca_1 * 4/(1 + (1 / y["Ca_i_free"]))
      }else{
        theta_LCaL <- 0 #- this is set to 0 to replicate EGTA present
      }
      
      norm2 <- y["C3_LCaL"]+y["C2_LCaL"]+y["C1_LCaL"]+y["C0_LCaL"]+y["C3Ca_LCaL"]+y["C2Ca_LCaL"]+y["C1Ca_LCaL"]+y["C0Ca_LCaL"]+y["O_LCaL"]+y["ICa_LCaL"]+y["IVS_LCaL"]+y["IVF_LCaL"]+y["IVS_Ca_LCaL"]+y["IVF_Ca_LCaL"]
      # Temp-type Calcium Current (ICaT)----
      I_CaT <- g_CaT * y["d_CaT"] * y["f_CaT"] * (y["Vm"] - E_Ca)
      
      d_CaT_inf <- d_CaT_inf(y["Vm"])
      f_CaT_inf <- f_CaT_inf(y["Vm"])
      
      tau_f_CaT <- 0.38117*(8.6+14.7*exp(-(y["Vm"]+50.0)*(y["Vm"]+50.0)/900.0)) / T_correction_Ca_2
      # Voltage-dependent Potassium Current (IKv)----
      I_Kv <- g_Kv * y["x_Kv"] * y["y_Kv"] * (y["Vm"] - E_K)
      I_NsK <- g_Kv * (y["Vm"] - E_K)
      x_Kv_inf <- x_Kv_inf(y["Vm"])
      y_Kv_inf <- y_Kv_inf(y["Vm"])
      # Calcium & Voltage-activated Potassium Current (IBK)-----
      I_BK <- g_BK * y["O4_BK"] * (y["Vm"] - E_K)
      a_BK <- 1.0*exp((8.47188*y["Vm"])/(1.0*Temp))
      b_BK <- 1.0*exp(-7.77556*y["Vm"]/(1.0*Temp))
      
      #voltage dependent transitions (ms)
      KC0O0_BK <- a_BK * 0.02162
      KC1O1_BK <- a_BK * 0.000869
      KC2O2_BK <- a_BK * 0.0000281
      KC3O3_BK <- a_BK * 0.000781
      KC4O4_BK <- a_BK * 0.044324
      
      KO0C0_BK <- b_BK * 318.1084
      KO1C1_BK <- b_BK * 144.1736
      KO2C2_BK <- b_BK * 32.6594
      KO3C3_BK <- b_BK * 0.095312
      KO4C4_BK <- b_BK * 0.000106
      
      # Define rate equations for calcium-dependent transitions
      KC0C1_BK <- (y["Ca_i_free"] * 4 * k_on)
      KC1C2_BK <- (y["Ca_i_free"] * 3 * k_on)
      KC2C3_BK <- (y["Ca_i_free"] * 2 * k_on)
      KC3C4_BK <- (y["Ca_i_free"] * 1 * k_on)
      
      KC4C3_BK <- (4 * k_off_C)
      KC3C2_BK <- (3 * k_off_C)
      KC2C1_BK <- (2 * k_off_C)
      KC1C0_BK <- (1 * k_off_C)
      
      KO0O1_BK <- (y["Ca_i_free"] * 4 * k_on)
      KO1O2_BK <- (y["Ca_i_free"] * 3 * k_on)
      KO2O3_BK <- (y["Ca_i_free"] * 2 * k_on)
      KO3O4_BK <- (y["Ca_i_free"] * 1 * k_on)
      
      KO4O3_BK <- (4 * k_off_O)
      KO3O2_BK <- (3 * k_off_O)
      KO2O1_BK <- (2 * k_off_O)
      KO1O0_BK <- (1 * k_off_O)
      norm1 <- y["C0_BK"] + y["C1_BK"] + y["C2_BK"] + y["C3_BK"] + y["C4_BK"] + y["O0_BK"] + y["O1_BK"] + y["O2_BK"] + y["O3_BK"] + y["O4_BK"]
      # Voltage-dependent Sodium Current (INa)-----
      I_Na <- g_Na * y["O_Na"] * (y["Vm"] - E_Na)
      
      k_I2I1_Na <- T_correction_Na*0.0039239*exp(2.6793+0.0061468*y["Vm"])
      k_I1O_Na <- T_correction_Na*0.12052*exp(-9.6028+0.083025*y["Vm"])
      k_OC1_Na <- T_correction_Na*2.391*exp(-13.335-0.25289*y["Vm"])
      k_C1C2_Na <- T_correction_Na*3.1566*exp(0.36352+0.077193*y["Vm"])
      k_C2C3_Na <- T_correction_Na*0.55432*exp(-0.099074+0.036441*y["Vm"])
      k_C3C2_Na <- T_correction_Na*0.00052548*exp(-0.069102+0.0031945*y["Vm"])
      k_C2C1_Na <- T_correction_Na*1.4496*exp(-0.1566+0.058353*y["Vm"])
      k_C1O_Na <- T_correction_Na*1.5329*exp(0.0093193+0.041075*y["Vm"])
      k_OI1_Na <- T_correction_Na*1.6164*exp(0.30763+0.0060535*y["Vm"])
      k_I1I2_Na <- T_correction_Na*0.027735*exp(0.05149-0.046865*y["Vm"])
      k_I1C1_Na <- T_correction_Na*1.9046*exp(-2.484+0.020406*y["Vm"])
      k_C1I1_Na <- T_correction_Na*0.00021688*exp(-0.063438+0.0046683*y["Vm"])
      norm3 <- y["O_Na"] + y["C1_Na"] + y["C2_Na"] + y["C3_Na"] + y["I1_Na"] + y["I2_Na"]
      # Sodium-Calcium Exchanger (NCX)------
      I_NCX <- P_NCX*(exp(gamma*y["Vm"]*Faraday/(R*Temp))*y["Na_i"]^3.0*Cao-2.5*exp((gamma-1.0)*y["Vm"]*Faraday/(R*Temp))*Nao^3.0*y["Ca_i_free"])/((K_mNa2^3.0+Nao^3.0)*(K_mCa+Cao)*(1.0+k_sat*exp((gamma-1.0)*y["Vm"]*Faraday/(R*Temp))))
      # Sodium-Potassium Pump (NaK)-----
      I_NaK <- P_NaK * ((Ko * y["Na_i"]) / ((K_mK + Ko) * (K_mNa + y["Na_i"])* (1 + 0.1245 * exp(((-0.1 * y["Vm"] * Faraday)/(R * Temp)))) + 0.0353 * exp(((-y["Vm"] * Faraday)/(R * Temp)))))
      
      # Non-Selective Leak Current (INS)-----
      I_NsNa <- g_NsNa * (y["Vm"] - E_Na)
      I_KNS <- g_NsK * (y["Vm"] - E_K)
      I_NS <- I_NsNa + I_KNS
      
      # Ionic Stimulation-----
      
      #Calcium buffering----
      theta_1_cabuff = (n_CRT*CRT_total*(y["Ca_i_free"]^(n_CRT - 1))*
                          (K_D_CRT))/((((y["Ca_i_free"]^n_CRT) + K_D_CRT)^2))
      theta_2_cabuff = (n_CaM*CaM_total*((y["Ca_i_free"]^(n_CaM - 1)))*
                          (K_D_CaM))/((((y["Ca_i_free"]^n_CaM) + K_D_CaM)^2))
      
      theta_cabuff = 1.0 + theta_1_cabuff + theta_2_cabuff
      
      #ODE-----
      # BK Channel
      #Closed
      d[1] <- -(KC0C1_BK+KC0O0_BK)*y["C0_BK"]/norm1+KC1C0_BK*y["C1_BK"]/norm1+KO0C0_BK*y["O0_BK"]/norm1 # closed state of BK channel 0
      d[2] <- -(KC1C0_BK+KC1O1_BK+KC1C2_BK)*y["C1_BK"]/norm1+KC0C1_BK*y["C0_BK"]/norm1+KO1C1_BK*y["O1_BK"]/norm1+KC2C1_BK*y["C2_BK"]/norm1  # closed state of BK channel 1
      d[3] <- -(KC2C1_BK+KC2O2_BK+KC2C3_BK)*y["C2_BK"]/norm1+KC1C2_BK*y["C1_BK"]/norm1+KO2C2_BK*y["O2_BK"]/norm1+KC3C2_BK*y["C3_BK"]/norm1  # closed state of BK channel 2
      d[4] <- -(KC3C2_BK+KC3O3_BK+KC3C4_BK)*y["C3_BK"]/norm1+KC2C3_BK*y["C2_BK"]/norm1+KO3C3_BK*y["O3_BK"]/norm1+KC4C3_BK*y["C4_BK"]/norm1 # closed state of BK channel 3
      d[5] <- -(KC4C3_BK+KC4O4_BK)*y["C4_BK"]/norm1+KC3C4_BK*y["C3_BK"]/norm1+KO4C4_BK*y["O4_BK"]/norm1  # closed state of BK channel 4
      #Open
      d[6] <- -(KO0O1_BK+KO0C0_BK)*y["O0_BK"]/norm1+KO1O0_BK*y["O1_BK"]/norm1+KC0O0_BK*y["C0_BK"]/norm1  # Open state of BK channel 0
      d[7] <- -(KO1O0_BK+KO1C1_BK+KO1O2_BK)*y["O1_BK"]/norm1+KO0O1_BK*y["O0_BK"]/norm1+KC1O1_BK*y["C1_BK"]/norm1+KO2O1_BK*y["O2_BK"]/norm1  # Open state of BK channel 1
      d[8] <- -(KO2O1_BK+KO2C2_BK+KO2O3_BK)*y["O2_BK"]/norm1+KO1O2_BK*y["O1_BK"]/norm1+KC2O2_BK*y["C2_BK"]/norm1+KO3O2_BK*y["O3_BK"]/norm1 # Open state of BK channel 2
      d[9] <- -(KO3O2_BK+KO3C3_BK+KO3O4_BK)*y["O3_BK"]/norm1+KO2O3_BK*y["O2_BK"]/norm1+KC3O3_BK*y["C3_BK"]/norm1+KO4O3_BK*y["O4_BK"]/norm1  # Open state of BK channel 3
      d[10] <- -(KO4O3_BK+KO4C4_BK)*y["O4_BK"]/norm1+KO3O4_BK*y["O3_BK"]/norm1+KC4O4_BK*y["C4_BK"]/norm1  # Open state of BK channel 4
      
      #L - type channels
      d[11] <- -(a0_LCaL+theta_LCaL)*y["C0_LCaL"]/norm2+b0_LCaL*y["C1_LCaL"]/norm2+delta_LCaL*y["C0Ca_LCaL"]/norm2
      d[12] <- theta_LCaL*y["C0_LCaL"]/norm2-(a0_LCaL+delta_LCaL)*y["C0Ca_LCaL"]/norm2+b0_LCaL*y["C1Ca_LCaL"]/norm2
      d[13] <- a0_LCaL*y["C0_LCaL"]/norm2-(a1_LCaL+b0_LCaL+theta_LCaL)*y["C1_LCaL"]/norm2+b1_LCaL*y["C2_LCaL"]/norm2+delta_LCaL*y["C1Ca_LCaL"]/norm2
      d[14] <- theta_LCaL*y["C1_LCaL"]/norm2+a0_LCaL*y["C0Ca_LCaL"]/norm2-(a1_LCaL+b0_LCaL+delta_LCaL)*y["C1Ca_LCaL"]/norm2+b1_LCaL*y["C2Ca_LCaL"]/norm2
      d[15] <- a1_LCaL*y["C1_LCaL"]/norm2-(a2_LCaL+b1_LCaL+theta_LCaL)*y["C2_LCaL"]/norm2+b2_LCaL*y["C3_LCaL"]/norm2+delta_LCaL*y["C2Ca_LCaL"]/norm2
      d[16] <- theta_LCaL*y["C2_LCaL"]/norm2+a1_LCaL*y["C1Ca_LCaL"]/norm2-(a2_LCaL+b1_LCaL+delta_LCaL)*y["C2Ca_LCaL"]/norm2+b2_LCaL*y["C3Ca_LCaL"]/norm2
      d[17] <- a2_LCaL*y["C2_LCaL"]/norm2-(a3_LCaL+b2_LCaL+phi_s+phi_f+theta_LCaL)*y["C3_LCaL"]/norm2+b3_LCaL*y["O_LCaL"]/norm2+omega_f*y["IVF_LCaL"]/norm2+omega_s*y["IVS_LCaL"]/norm2+delta_LCaL*y["C3Ca_LCaL"]/norm2
      d[18] <- theta_LCaL*y["C3_LCaL"]/norm2+a2_LCaL*y["C2Ca_LCaL"]/norm2-(a3_LCaL+b2_LCaL+phi_f+phi_s+delta_LCaL)*y["C3Ca_LCaL"]/norm2+b3_LCaL*y["ICa_LCaL"]/norm2+omega_f*y["IVF_Ca_LCaL"]/norm2+omega_s*y["IVS_Ca_LCaL"]/norm2
      d[19] <- theta_LCaL*y["O_LCaL"]/norm2+a3_LCaL*y["C3Ca_LCaL"]/norm2-(b3_LCaL+psi_f+psi_s+delta_LCaL)*y["ICa_LCaL"]/norm2+xi_f*y["IVF_Ca_LCaL"]/norm2+xi_s*y["IVS_Ca_LCaL"]/norm2 # Calcium dependent
      d[20] <- phi_f*y["C3_LCaL"]/norm2+psi_f*y["O_LCaL"]/norm2-(omega_f+xi_f+omega_fs+theta_LCaL)*y["IVF_LCaL"]/norm2+omega_sf*y["IVS_LCaL"]/norm2+delta_LCaL*y["IVF_Ca_LCaL"]/norm2 # Voltage dependent - fast 
      d[21] <- theta_LCaL*y["IVF_LCaL"]/norm2+phi_f*y["C3Ca_LCaL"]/norm2+psi_f*y["ICa_LCaL"]/norm2-(omega_f+xi_f+omega_fs+delta_LCaL)*y["IVF_Ca_LCaL"]/norm2+omega_sf*y["IVS_Ca_LCaL"]/norm2 # fast voltage-dependent Calcium-dependent
      d[22] <- phi_s*y["C3_LCaL"]/norm2+psi_s*y["O_LCaL"]/norm2+omega_fs*y["IVF_LCaL"]/norm2-(omega_s+xi_s+omega_sf+theta_LCaL)*y["IVS_LCaL"]/norm2+delta_LCaL*y["IVS_Ca_LCaL"]/norm2 # Voltage-dependent - slow
      d[23] <- theta_LCaL*y["IVS_LCaL"]/norm2+phi_s*y["C3Ca_LCaL"]/norm2+psi_s*y["ICa_LCaL"]/norm2+omega_fs*y["IVF_Ca_LCaL"]/norm2-(omega_s+xi_s+omega_sf+delta_LCaL)*y["IVS_Ca_LCaL"]/norm2 # slow voltage dependent Calcium dependent
      d[24] <- a3_LCaL*y["C3_LCaL"]/norm2-(b3_LCaL+psi_f+psi_s+theta_LCaL)*y["O_LCaL"]/norm2+xi_f*y["IVF_LCaL"]/norm2+xi_s*y["IVS_LCaL"]/norm2+delta_LCaL*y["ICa_LCaL"]/norm2
      
      # Sodium Channel ODE
      d[25] <- -(k_C1C2_Na+k_C1O_Na+k_C1I1_Na)*y["C1_Na"]/norm3+k_OC1_Na*y["O_Na"]/norm3+k_C2C1_Na*y["C2_Na"]/norm3+k_I1C1_Na*y["I1_Na"]/norm3 # Closed state 1
      d[26] <- -(k_C2C1_Na+k_C2C3_Na)*y["C2_Na"]/norm3+k_C1C2_Na*y["C1_Na"]/norm3+k_C3C2_Na*y["C3_Na"]/norm3 # Closed state 2
      d[27] <- -k_C3C2_Na*y["C3_Na"]/norm3+k_C2C3_Na*y["C2_Na"]/norm3 # Closed state 3
      d[28] <- -(k_I1O_Na+k_I1I2_Na+k_I1C1_Na)*y["I1_Na"]/norm3+k_I2I1_Na*y["I2_Na"]/norm3+k_C1I1_Na*y["C1_Na"]/norm3+k_OI1_Na*y["O_Na"]/norm3  # Inactivated State 1
      d[29] <- -k_I2I1_Na*y["I2_Na"]/norm3+k_I1I2_Na*y["I1_Na"]/norm3 # Inactivated 2
      d[30] <- -(k_OC1_Na+k_OI1_Na)*y["O_Na"]/norm3+k_C1O_Na*y["C1_Na"]/norm3+k_I1O_Na*y["I1_Na"]/norm3  # Open 
      
      # Temp-type Calcium Channel ODEs
      d[31] <- (d_CaT_inf - y["d_CaT"]) / tau_d_CaT
      d[32] <- (f_CaT_inf - y["f_CaT"]) / tau_f_CaT
      
      #Ion concentration tracked in mM
      # Intracellular Calcium Concentration - d[Ca_i_free]
      d[33] <- -1.0e-15*((I_CaL+I_CaT-I_NCX*2.0)/(2.0*V_myo*Faraday)/theta_cabuff)
      # Intracellular Potassium Concentration - d[K_i]
      d[34] <- -1e-15 * (I_Kv + I_BK + I_NsK + I_stim - I_NaK * 2)/(V_myo * Faraday)
      # Intracellular Sodium Concentration - d[Na_i]
      d[35] <- -1/1e15 * (I_Na+I_NsNa+I_NCX*3.0+I_NaK*3.0)/(V_myo * Faraday)
      #Single hJSMC electrophysiology -dVm
      d[36] <- -1.0/Cm * (I_Na+I_CaL+I_CaT+I_Kv+I_BK+I_NCX+I_NaK+I_NsK+I_NsNa+I_stim)
      # Voltage-dependent Potassium Channel ODEs
      d[37] <- (x_Kv_inf - y["x_Kv"]) / tau_x_Kv
      d[38] <- (y_Kv_inf - y["y_Kv"]) / tau_y_Kv
      print(times)
      # Return the rates of change
      return(list(d))
    })
  }
  #Fixed parameters
  #Fixed parameters
  fixed <- list()
  g_BK <- 80.0   # nanoS (in Calcium_voltage_activated_potassium_channel)
  Gcouple <- 2.6   # nanoS (in I_stim)
  V_ICCamp <- 33.5   # millivolt (in I_stim)
  V_ICCrest <- -57.0   # millivolt (in I_stim)
  f_1 <- 12000.0   # millisecond (in I_stim)
  f_2 <- 300.0   # millisecond (in I_stim)
  period <- 10000    # millisecond (in I_stim)
  t_ICCpeak <- 300.0   # millisecond (in I_stim)
  t_slope <- 600.0   # millisecond (in I_stim)
  t_ICCplateau <- 9700 #millisecond (in I_stim)
  Q10Ca1 <- 2.1   # dimensionless (Q10Ca1 in L_type_Ca_channel_states)
  g_CaL <- 1.44   # nanoS (in L_type_Ca_channel)
  K_mCa <- 1.38   # millimolar (in Na_Ca_exchanger)
  K_mNai <- 87.5   # millimolar (in Na_Ca_exchanger)
  P_NCX <- 1992.335   # picoA (in Na_Ca_exchanger)
  gamma <- 0.35   # dimensionless (in Na_Ca_exchanger)
  k_sat <- 0.1   # dimensionless (in Na_Ca_exchanger)
  Q10Na <- 2.45   # dimensionless (in Na_channel_states)
  Q10Ca2 <- 2.1   # dimensionless (Q10Ca1 in T_type_Ca_channel_d_gate)
  Q10Ca3 <- 2.1   # dimensionless (Q10Ca1 in T_type_Ca_channel_f_gate)
  g_CaT <- 0.0425   # nanoS (in T_type_Ca_channel)
  CRT_total <- 0.034   # millimolar (in ionic_concentrations)
  CaM_total <- 0.012   # millimolar (in ionic_concentrations)
  Cao <- 2.5   # millimolar (in ionic_concentrations)
  K_D_CRT <- 0.0009   # millimolar (in ionic_concentrations)
  K_D_CaM <- 0.0001   # millimolar4 (in ionic_concentrations)
  Ko <- 5.4   # millimolar (in ionic_concentrations)
  Nao <- 140.0   # millimolar (in ionic_concentrations)
  n_CRT <- 1.0   # dimensionless (in ionic_concentrations)
  n_CaM <- 4.0   # dimensionless (in ionic_concentrations)
  Cm <- 50.0   # picoF (in membrane)
  Faraday <- 96.48534   # coulomb_per_millimole (in membrane)
  R <- 8.314   # joule_per_mole_kelvin (in membrane)
  Temp <- 310.0   # kelvin (in membrane)
  Texp <- 297
  g_NsK <- 0.017512   # nanoS (in non_specific_K_current)
  g_NsNa <- 0.022488   # nanoS (in non_specific_Na_current)
  g_Na <- 25.1   # nanoS (in sodium_current)
  K_mK <- 1.0   # millimolar (in sodium_potassium_pump)
  K_mNa <- 40.0   # millimolar (in sodium_potassium_pump)
  K_mNa2 <- 87.5
  P_NaK <- 9.26   # picoA (in sodium_potassium_pump)
  Q10K1 <- 3.1   # dimensionless (Q10K1 in voltage_dep_K_channel_x_gate)
  Q10K2 <- 3.1   # dimensionless (Q10K1 in voltage_dep_K_channel_y_gate)
  g_Kv <- 1.0217   # nanoS (in voltage_dep_K_channel)
  k_on <- 40633.0
  k_off_C <- 11.0
  k_off_O <- 1.1
  SMC_resting <- -60 #
  V_myo <- 3.5e-12
  T_correction_Ca_1 <- 1.0 * Q10Ca1^((Temp - 310) / 10.0) # experimental correction
  T_correction_Ca_2 <- 1.0 * Q10Ca1^((Temp - Texp) / 10.0) # experimental correction
  T_correction_Na <- 1.0 * Q10Na^((Temp - Texp) / 10.0) # experimental correction
  T_correction_K <- 1.0 * Q10K1^((Temp - Texp) / 10.0) # experimental correction
  
  if(parms[["CalciumDependence"]] == 1){
    delta_LCaL <- 0.01 * T_correction_Ca_1 #FOR ICaL set to 0 to replicate EGTA 
  }else{
    delta_LCaL <- 0 #- this is set to 0 to replicate EGTA present
  }
  tau_d_CaT <- 1.9058 / T_correction_Ca_2 #for ICaL
  tau_x_Kv <- 4.7803 / T_correction_K  # ms - Time constant for Kvx
  tau_y_Kv <- 763.7564 / T_correction_K  # ms - Time constant for Kvy
  #Function for CaT
  d_CaT_inf <- function(Vm) {
    return(1 / (1 + exp((-(Vm + 60.5) / 5.3))))
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
  
  
  # Initial values from MATLAB
  Y = c(0.48379087935899, 0.385183559520031, 0.115002824567753, 0.0152602714149774, 0.000759264410974374, 6.94960798375172e-7, 5.55636826398253e-8, 2.85143702125325e-8, 1.59832806123435e-6, 1.82113764497095e-6, 0.815464741971086, 0.0175888495282545, 0.152399266235657, 0.00328711668724504, 0.0106805060777161, 0.000230369020877669, 0.000332673548872087, 7.1754726923539e-6, 8.38123983500905e-8, 4.0998751301597e-6, 8.84306615061238e-8, 1.1193313274705e-6, 2.41429816075123e-8, 3.88576045134351e-6, 0.0119443135223679, 0.0109545368437155, 0.973782548650071, 0.000126882921013389, 0.00318975045717667, 1.96760342050475e-6, 0.0791635737410974, 0.377831534375835, 9.6e-5, 153.604280337996, 10.5731241425458, -60, 0.14714161078933, 0.99994773314105)
  initial <- c(
    C0_BK = Y[1],   # Closed, 0 Ca²⁺
    C1_BK = Y[2],  # Closed, 1 Ca²⁺
    C2_BK = Y[3],  # Closed, 2 Ca²⁺
    C3_BK = Y[4],  # Closed, 3 Ca²⁺
    C4_BK = Y[5],  # Closed, 4 Ca²⁺
    O0_BK = Y[6],  # Open, 0 Ca²⁺
    O1_BK = Y[7],  # Open, 1 Ca²⁺
    O2_BK = Y[8],  # Open, 2 Ca²⁺
    O3_BK = Y[9],  # Open, 3 Ca²⁺
    O4_BK = Y[10], # Open, 4 Ca²⁺ (conducting state)
    C0_LCaL = Y[11],
    C0Ca_LCaL = Y[12],
    C1_LCaL = Y[13],
    C1Ca_LCaL = Y[14],
    C2_LCaL = Y[15],
    C2Ca_LCaL = Y[16],
    C3_LCaL = Y[17],
    C3Ca_LCaL = Y[18],
    ICa_LCaL = Y[19],
    IVF_LCaL = Y[20],
    IVF_Ca_LCaL = Y[21],
    IVS_LCaL = Y[22],
    IVS_Ca_LCaL = Y[23],
    O_LCaL = Y[24],
    C1_Na = Y[25],
    C2_Na = Y[26],
    C3_Na = Y[27],
    I1_Na = Y[28],
    I2_Na = Y[29],
    O_Na = Y[30],
    d_CaT = Y[31], 
    f_CaT = Y[32],
    Ca_i_free = Y[33],  # mM - Initial free intracellular Calcium concentration
    K_i = Y[34],  # mM - Initial intracellular potassium concentration
    Na_i = Y[35],  # mM - Initial intracellular sodium concentration
    Vm = Y[36],  # Membrane potential - millivolts
    x_Kv = Y[37],
    y_Kv = Y[38]
  )
  
  
  # Time sequence for the simulation ms
  times <- seq(1760000, 1800000, by = 1)
  #history
  
  
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  output_df$stim_cycle <- output_df$time %/% period
  output_df$stim_start <- output_df$stim_cycle * period
  output_df$local_time <- output_df$time - output_df$stim_start  # Compute local time
  output_df$Vm_ICC <- ifelse(output_df$local_time >= 0 & output_df$local_time < t_ICCpeak,
                             V_ICCrest + V_ICCamp * output_df$local_time / f_2,
                             ifelse(output_df$local_time >= t_ICCpeak & output_df$local_time < t_ICCplateau,
                                    V_ICCrest + V_ICCamp * (1 + exp(-f_1 / (2 * t_slope))) * 1 / (1 + exp((output_df$local_time - f_2 - 0.5 * f_1) / t_slope)),
                                    V_ICCrest))
  return(output_df)
}
parms <- c(CalciumDependence = 1)
Initial_out <- Model(parms = parms)

ICC_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time / 1000, y = Vm_ICC, colour = "Vm")) +
  scale_colour_manual(values = c("Vm" = "black")) +
  geom_hline(yintercept = c(-57, -23.5), linetype = "dashed", color = "red") +  # Horizontal lines
  geom_vline(xintercept = c(300 / 1000, 9700 / 1000), linetype = "dashed", color = "blue") +  # Vertical lines (convert ms to s)
  xlab("time (s)") +
  ylab("Voltage (mV)") +
  ggtitle("ICC voltage sim")


Voltage_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time /1000, y = Vm, colour = "Vm")) +
  scale_colour_manual(values = c("Vm" = "black")) +
  #xlim(1760,1800) +
  #ylim(-65,-30) +
  xlab("time (s)") +
  ylab("Voltage (mV)") +
  ggtitle("hJMC voltage simulation")

Ca_i_free_plot <- ggplot() +
  geom_line(data = Initial_out, aes(x = time /1000, y = Ca_i_free * 10^6, colour = "Ca_i_free")) +
  scale_colour_manual(values = c("Ca_i_free" = "black")) +
  #xlim(1760,1800) +
  #ylim(0,300) +
  xlab("time (s)") +
  ylab("Free Intracellular\n calcium (nM)") +
  ggtitle("hJMC Free Intracellular\n calcium (nM) simulation")

grid.arrange(grobs = list(Voltage_plot, Ca_i_free_plot), ncol = 1, nrow = 2)