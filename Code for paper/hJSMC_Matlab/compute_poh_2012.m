%===============================================================================
%   A Quantitative Model of Human Jejunal Smooth Muscle Cell Electrophysiology
%   Paper details:                                                              
%   Authors: Poh YC, Corrias A, Cheng N and Buist ML                           
%   PLoS ONE 7(8): e42385. doi:10.1371/journal.pone.0042385
%-------------------------------------------------------------------------------
%   Conversion from CellML 1.0 to MATLAB code (init) was done using COR (0.9.31.1409)
%   Copyright 2002-2012 Dr Alan Garny
%   http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%
%   Main function compute_poh_2012 was written by YC Poh 
%   (email address:  poh.yongcheng@gmail.com)
%-------------------------------------------------------------------------------
%   Feedback & comments are welcomed
% 
%   The CellML & C versions of this model are available for download from  
%   the Computational Bioengineering Laboratory (NUS) website:                   
%   http://www.bioeng.nus.edu.sg/compbiolab/
%===============================================================================

function compute_poh_2012

%clear command window, and close all pre-existing figures
clc
close all

%initial conditions
Y = [0.48379087935899, 0.385183559520031, 0.115002824567753, 0.0152602714149774, 0.000759264410974374, 6.94960798375172e-7, 5.55636826398253e-8, 2.85143702125325e-8, 1.59832806123435e-6, 1.82113764497095e-6, 0.815464741971086, 0.0175888495282545, 0.152399266235657, 0.00328711668724504, 0.0106805060777161, 0.000230369020877669, 0.000332673548872087, 7.1754726923539e-6, 8.38123983500905e-8, 4.0998751301597e-6, 8.84306615061238e-8, 1.1193313274705e-6, 2.41429816075123e-8, 3.88576045134351e-6, 0.0119443135223679, 0.0109545368437155, 0.973782548650071, 0.000126882921013389, 0.00318975045717667, 1.96760342050475e-6, 0.0791635737410974, 0.377831534375835, 5.38843941249284e-5, 153.604280337996, 10.5731241425458, -73.5049651455872, 0.14714161078933, 0.99994773314105];

%solve dY using ode fn
[T,Y]=ode15s(@poh_2012,[0 40000],Y);

%plot Vm and Cai; refer to State variables in poh_2012 for the list of
%indices
subplot(2,1,1)
plot(T,Y(:,36)); ylabel('Membrane voltage (mV)');

subplot(2,1,2)
plot(T,Y(:,33)); ylabel('Free intracellular [Ca^{2+}] (mM)'); xlabel('Time (ms)')

return

function dY = poh_2012(time, Y)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.48379087935899, 0.385183559520031, 0.115002824567753, 0.0152602714149774, 0.000759264410974374, 6.94960798375172e-7, 5.55636826398253e-8, 2.85143702125325e-8, 1.59832806123435e-6, 1.82113764497095e-6, 0.815464741971086, 0.0175888495282545, 0.152399266235657, 0.00328711668724504, 0.0106805060777161, 0.000230369020877669, 0.000332673548872087, 7.1754726923539e-6, 8.38123983500905e-8, 4.0998751301597e-6, 8.84306615061238e-8, 1.1193313274705e-6, 2.41429816075123e-8, 3.88576045134351e-6, 0.0119443135223679, 0.0109545368437155, 0.973782548650071, 0.000126882921013389, 0.00318975045717667, 1.96760342050475e-6, 0.0791635737410974, 0.377831534375835, 5.38843941249284e-5, 153.604280337996, 10.5731241425458, -73.5049651455872, 0.14714161078933, 0.99994773314105];

% YNames = {'C01', 'C11', 'C21', 'C31', 'C4', 'O0', 'O1', 'O2', 'O3', 'O4', 'C02', 'C0Ca', 'C12', 'C1Ca', 'C22', 'C2Ca', 'C32', 'C3Ca', 'ICa', 'IVf', 'IVfCa', 'IVs', 'IVsCa', 'O_CaL', 'C13', 'C23', 'C33', 'I1', 'I2', 'O_Na', 'd_CaT', 'f_CaT', 'Cai', 'Ki', 'Nai', 'V', 'x_Kv', 'y_Kv'};
% YUnits = {'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'millivolt', 'dimensionless', 'dimensionless'};
% YComponents = {'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'Calcium_voltage_activated_potassium_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'L_type_Ca_channel_states', 'Na_channel_states', 'Na_channel_states', 'Na_channel_states', 'Na_channel_states', 'Na_channel_states', 'Na_channel_states', 'T_type_Ca_channel_d_gate', 'T_type_Ca_channel_f_gate', 'ionic_concentrations', 'ionic_concentrations', 'ionic_concentrations', 'membrane', 'voltage_dep_K_channel_x_gate', 'voltage_dep_K_channel_y_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: C01 (dimensionless) (C0 in Calcium_voltage_activated_potassium_channel_states)
% 2: C11 (dimensionless) (C1 in Calcium_voltage_activated_potassium_channel_states)
% 3: C21 (dimensionless) (C2 in Calcium_voltage_activated_potassium_channel_states)
% 4: C31 (dimensionless) (C3 in Calcium_voltage_activated_potassium_channel_states)
% 5: C4 (dimensionless) (in Calcium_voltage_activated_potassium_channel_states)
% 6: O0 (dimensionless) (in Calcium_voltage_activated_potassium_channel_states)
% 7: O1 (dimensionless) (in Calcium_voltage_activated_potassium_channel_states)
% 8: O2 (dimensionless) (in Calcium_voltage_activated_potassium_channel_states)
% 9: O3 (dimensionless) (in Calcium_voltage_activated_potassium_channel_states)
% 10: O4 (dimensionless) (in Calcium_voltage_activated_potassium_channel_states)
% 11: C02 (dimensionless) (C0 in L_type_Ca_channel_states)
% 12: C0Ca (dimensionless) (in L_type_Ca_channel_states)
% 13: C12 (dimensionless) (C1 in L_type_Ca_channel_states)
% 14: C1Ca (dimensionless) (in L_type_Ca_channel_states)
% 15: C22 (dimensionless) (C2 in L_type_Ca_channel_states)
% 16: C2Ca (dimensionless) (in L_type_Ca_channel_states)
% 17: C32 (dimensionless) (C3 in L_type_Ca_channel_states)
% 18: C3Ca (dimensionless) (in L_type_Ca_channel_states)
% 19: ICa (dimensionless) (in L_type_Ca_channel_states)
% 20: IVf (dimensionless) (in L_type_Ca_channel_states)
% 21: IVfCa (dimensionless) (in L_type_Ca_channel_states)
% 22: IVs (dimensionless) (in L_type_Ca_channel_states)
% 23: IVsCa (dimensionless) (in L_type_Ca_channel_states)
% 24: O_CaL (dimensionless) (in L_type_Ca_channel_states)
% 25: C13 (dimensionless) (C1 in Na_channel_states)
% 26: C23 (dimensionless) (C2 in Na_channel_states)
% 27: C33 (dimensionless) (C3 in Na_channel_states)
% 28: I1 (dimensionless) (in Na_channel_states)
% 29: I2 (dimensionless) (in Na_channel_states)
% 30: O_Na (dimensionless) (in Na_channel_states)
% 31: d_CaT (dimensionless) (in T_type_Ca_channel_d_gate)
% 32: f_CaT (dimensionless) (in T_type_Ca_channel_f_gate)
% 33: Cai (millimolar) (in ionic_concentrations)
% 34: Ki (millimolar) (in ionic_concentrations)
% 35: Nai (millimolar) (in ionic_concentrations)
% 36: V (millivolt) (in membrane)
% 37: x_Kv (dimensionless) (in voltage_dep_K_channel_x_gate)
% 38: y_Kv (dimensionless) (in voltage_dep_K_channel_y_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

g_BK = 80.0;   % nanoS (in Calcium_voltage_activated_potassium_channel)
Gcouple = 2.6;   % nanoS (in I_stim)
V_ICCamp = 33.5;   % millivolt (in I_stim)
V_ICCrest = -57.0;   % millivolt (in I_stim)
f_1 = 12000.0;   % millisecond (in I_stim)
f_2 = 300.0;   % millisecond (in I_stim)
period = 10000.0;   % millisecond (in I_stim)
t_ICCpeak = 300.0;   % millisecond (in I_stim)
t_slope = 600.0;   % millisecond (in I_stim)
Q10Ca1 = 2.1;   % dimensionless (Q10Ca in L_type_Ca_channel_states)
g_CaL = 1.44;   % nanoS (in L_type_Ca_channel)
K_mCa = 1.38;   % millimolar (in Na_Ca_exchanger)
K_mNai = 87.5;   % millimolar (in Na_Ca_exchanger)
P_NCX = 1992.335;   % picoA (in Na_Ca_exchanger)
gamma = 0.35;   % dimensionless (in Na_Ca_exchanger)
k_sat = 0.1;   % dimensionless (in Na_Ca_exchanger)
Q10Na = 2.45;   % dimensionless (in Na_channel_states)
Q10Ca2 = 2.1;   % dimensionless (Q10Ca in T_type_Ca_channel_d_gate)
Q10Ca3 = 2.1;   % dimensionless (Q10Ca in T_type_Ca_channel_f_gate)
g_CaT = 0.0425;   % nanoS (in T_type_Ca_channel)
CRT_total = 0.034;   % millimolar (in ionic_concentrations)
CaM_total = 0.012;   % millimolar (in ionic_concentrations)
Cao = 2.0;   % millimolar (in ionic_concentrations)
K_D_CRT = 0.0009;   % millimolar (in ionic_concentrations)
K_D_CaM = 0.0001;   % millimolar4 (in ionic_concentrations)
Ko = 5.4;   % millimolar (in ionic_concentrations)
Nao = 140.0;   % millimolar (in ionic_concentrations)
n_CRT = 1.0;   % dimensionless (in ionic_concentrations)
n_CaM = 4.0;   % dimensionless (in ionic_concentrations)
Cm = 50.0;   % picoF (in membrane)
F = 96.48534;   % coulomb_per_millimole (in membrane)
R = 8.314;   % joule_per_mole_kelvin (in membrane)
T = 310.0;   % kelvin (in membrane)
g_NsK = 0.017512;   % nanoS (in non_specific_K_current)
g_NsNa = 0.022488;   % nanoS (in non_specific_Na_current)
g_Na = 25.1;   % nanoS (in sodium_current)
K_mK = 1.0;   % millimolar (in sodium_potassium_pump)
K_mNa = 40.0;   % millimolar (in sodium_potassium_pump)
P_NaK = 9.26;   % picoA (in sodium_potassium_pump)
Q10K1 = 3.1;   % dimensionless (Q10K in voltage_dep_K_channel_x_gate)
Q10K2 = 3.1;   % dimensionless (Q10K in voltage_dep_K_channel_y_gate)
g_Kv = 1.0217;   % nanoS (in voltage_dep_K_channel)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% alpha1 (per_millisecond) (alpha in Calcium_voltage_activated_potassium_channel_states)
% beta1 (per_millisecond) (beta in Calcium_voltage_activated_potassium_channel_states)
% k_C0C1 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C0O0 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C1C0 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C1C2_1 (per_millisecond) (k_C1C2 in Calcium_voltage_activated_potassium_channel_states)
% k_C1O1 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C2C1_1 (per_millisecond) (k_C2C1 in Calcium_voltage_activated_potassium_channel_states)
% k_C2C3_1 (per_millisecond) (k_C2C3 in Calcium_voltage_activated_potassium_channel_states)
% k_C2O2 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C3C2_1 (per_millisecond) (k_C3C2 in Calcium_voltage_activated_potassium_channel_states)
% k_C3C4 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C3O3 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C4C3 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_C4O4 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O0C0 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O0O1 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O1C1 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O1O0 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O1O2 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O2C2 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O2O1 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O2O3 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O3C3 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O3O2 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O3O4 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O4C4 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_O4O3 (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_off_C (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_off_O (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% k_on (per_millisecond) (in Calcium_voltage_activated_potassium_channel_states)
% norm1 (dimensionless) (norm in Calcium_voltage_activated_potassium_channel_states)
% E_K_1 (millivolt) (E_K in Calcium_voltage_activated_potassium_channel)
% i_BK (picoA) (in Calcium_voltage_activated_potassium_channel)
% I_stim (picoA) (in I_stim)
% local_time (millisecond) (in I_stim)
% stim_start (millisecond) (in I_stim)
% T_correction_Ca_1 (dimensionless) (T_correction_Ca in L_type_Ca_channel_states)
% alpha2 (per_millisecond) (alpha in L_type_Ca_channel_states)
% alpha_0 (per_millisecond) (in L_type_Ca_channel_states)
% alpha_1 (per_millisecond) (in L_type_Ca_channel_states)
% alpha_2 (per_millisecond) (in L_type_Ca_channel_states)
% alpha_3 (per_millisecond) (in L_type_Ca_channel_states)
% beta2 (per_millisecond) (beta in L_type_Ca_channel_states)
% beta_0 (per_millisecond) (in L_type_Ca_channel_states)
% beta_1 (per_millisecond) (in L_type_Ca_channel_states)
% beta_2 (per_millisecond) (in L_type_Ca_channel_states)
% beta_3 (per_millisecond) (in L_type_Ca_channel_states)
% delta (per_millisecond) (in L_type_Ca_channel_states)
% norm2 (dimensionless) (norm in L_type_Ca_channel_states)
% omega_f (per_millisecond) (in L_type_Ca_channel_states)
% omega_fs (per_millisecond) (in L_type_Ca_channel_states)
% omega_s (per_millisecond) (in L_type_Ca_channel_states)
% omega_sf (per_millisecond) (in L_type_Ca_channel_states)
% phi_f (per_millisecond) (in L_type_Ca_channel_states)
% phi_s (per_millisecond) (in L_type_Ca_channel_states)
% psi_f (per_millisecond) (in L_type_Ca_channel_states)
% psi_s (per_millisecond) (in L_type_Ca_channel_states)
% theta (per_millisecond) (in L_type_Ca_channel_states)
% xi_f (per_millisecond) (in L_type_Ca_channel_states)
% xi_s (per_millisecond) (in L_type_Ca_channel_states)
% E_Ca_1 (millivolt) (E_Ca in L_type_Ca_channel)
% i_Ca_L (picoA) (in L_type_Ca_channel)
% i_NCX (picoA) (in Na_Ca_exchanger)
% T_correction_Na (dimensionless) (in Na_channel_states)
% k_C1C2_2 (per_millisecond) (k_C1C2 in Na_channel_states)
% k_C1I1 (per_millisecond) (in Na_channel_states)
% k_C1O (per_millisecond) (in Na_channel_states)
% k_C2C1_2 (per_millisecond) (k_C2C1 in Na_channel_states)
% k_C2C3_2 (per_millisecond) (k_C2C3 in Na_channel_states)
% k_C3C2_2 (per_millisecond) (k_C3C2 in Na_channel_states)
% k_I1C1 (per_millisecond) (in Na_channel_states)
% k_I1I2 (per_millisecond) (in Na_channel_states)
% k_I1O (per_millisecond) (in Na_channel_states)
% k_I2I1 (per_millisecond) (in Na_channel_states)
% k_OC1 (per_millisecond) (in Na_channel_states)
% k_OI1 (per_millisecond) (in Na_channel_states)
% norm3 (dimensionless) (norm in Na_channel_states)
% T_correction_Ca_2 (dimensionless) (T_correction_Ca in T_type_Ca_channel_d_gate)
% d_CaT_inf (dimensionless) (in T_type_Ca_channel_d_gate)
% tau_d_CaT (millisecond) (in T_type_Ca_channel_d_gate)
% T_correction_Ca_3 (dimensionless) (T_correction_Ca in T_type_Ca_channel_f_gate)
% f_CaT_inf (dimensionless) (in T_type_Ca_channel_f_gate)
% tau_f_CaT (millisecond) (in T_type_Ca_channel_f_gate)
% E_Ca_2 (millivolt) (E_Ca in T_type_Ca_channel)
% i_Ca_T (picoA) (in T_type_Ca_channel)
% V_myo (litre) (in ionic_concentrations)
% E_K_2 (millivolt) (E_K in non_specific_K_current)
% i_NsK (picoA) (in non_specific_K_current)
% E_Na_1 (millivolt) (E_Na in non_specific_Na_current)
% i_NsNa (picoA) (in non_specific_Na_current)
% E_Na_2 (millivolt) (E_Na in sodium_current)
% i_Na (picoA) (in sodium_current)
% i_NaK (picoA) (in sodium_potassium_pump)
% T_correction_K_1 (dimensionless) (T_correction_K in voltage_dep_K_channel_x_gate)
% tau_x_Kv (millisecond) (in voltage_dep_K_channel_x_gate)
% x_Kv_inf (dimensionless) (in voltage_dep_K_channel_x_gate)
% T_correction_K_2 (dimensionless) (T_correction_K in voltage_dep_K_channel_y_gate)
% tau_y_Kv (millisecond) (in voltage_dep_K_channel_y_gate)
% y_Kv_inf (dimensionless) (in voltage_dep_K_channel_y_gate)
% E_K_3 (millivolt) (E_K in voltage_dep_K_channel)
% i_Kv (picoA) (in voltage_dep_K_channel)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (millisecond)

E_K_1 = R*T/F*log(Ko/Y(34));
i_BK = g_BK*Y(10)*(Y(36)-E_K_1);
alpha1 = 1.0*exp(8.47188*Y(36)/(1.0*T));
beta1 = 1.0*exp(-7.77556*Y(36)/(1.0*T));
k_on = 40633.0;
k_off_C = 11.0;
k_off_O = 1.1;
k_C0O0 = 0.02162*alpha1;
k_C1O1 = 0.000869*alpha1;
k_C2O2 = 0.0000281*alpha1;
k_C3O3 = 0.000781*alpha1;
k_C4O4 = 0.044324*alpha1;
k_O0C0 = 318.1084*beta1;
k_O1C1 = 144.1736*beta1;
k_O2C2 = 32.6594*beta1;
k_O3C3 = 0.095312*beta1;
k_O4C4 = 0.000106*beta1;
k_C0C1 = 4.0*k_on*Y(33);
k_C1C2_1 = 3.0*k_on*Y(33);
k_C2C3_1 = 2.0*k_on*Y(33);
k_C3C4 = 1.0*k_on*Y(33);
k_C4C3 = 4.0*k_off_C;
k_C3C2_1 = 3.0*k_off_C;
k_C2C1_1 = 2.0*k_off_C;
k_C1C0 = 1.0*k_off_C;
k_O0O1 = 4.0*k_on*Y(33);
k_O1O2 = 3.0*k_on*Y(33);
k_O2O3 = 2.0*k_on*Y(33);
k_O3O4 = 1.0*k_on*Y(33);
k_O4O3 = 4.0*k_off_O;
k_O3O2 = 3.0*k_off_O;
k_O2O1 = 2.0*k_off_O;
k_O1O0 = 1.0*k_off_O;
norm1 = Y(1)+Y(2)+Y(3)+Y(4)+Y(5)+Y(6)+Y(7)+Y(8)+Y(9)+Y(10);
dY(5, 1) = -(k_C4C3+k_C4O4)*Y(5)/norm1+k_C3C4*Y(4)/norm1+k_O4C4*Y(10)/norm1;
dY(4, 1) = -(k_C3C2_1+k_C3O3+k_C3C4)*Y(4)/norm1+k_C2C3_1*Y(3)/norm1+k_O3C3*Y(9)/norm1+k_C4C3*Y(5)/norm1;
dY(3, 1) = -(k_C2C1_1+k_C2O2+k_C2C3_1)*Y(3)/norm1+k_C1C2_1*Y(2)/norm1+k_O2C2*Y(8)/norm1+k_C3C2_1*Y(4)/norm1;
dY(2, 1) = -(k_C1C0+k_C1O1+k_C1C2_1)*Y(2)/norm1+k_C0C1*Y(1)/norm1+k_O1C1*Y(7)/norm1+k_C2C1_1*Y(3)/norm1;
dY(1, 1) = -(k_C0C1+k_C0O0)*Y(1)/norm1+k_C1C0*Y(2)/norm1+k_O0C0*Y(6)/norm1;
dY(10, 1) = -(k_O4O3+k_O4C4)*Y(10)/norm1+k_O3O4*Y(9)/norm1+k_C4O4*Y(5)/norm1;
dY(9, 1) = -(k_O3O2+k_O3C3+k_O3O4)*Y(9)/norm1+k_O2O3*Y(8)/norm1+k_C3O3*Y(4)/norm1+k_O4O3*Y(10)/norm1;
dY(8, 1) = -(k_O2O1+k_O2C2+k_O2O3)*Y(8)/norm1+k_O1O2*Y(7)/norm1+k_C2O2*Y(3)/norm1+k_O3O2*Y(9)/norm1;
dY(7, 1) = -(k_O1O0+k_O1C1+k_O1O2)*Y(7)/norm1+k_O0O1*Y(6)/norm1+k_C1O1*Y(2)/norm1+k_O2O1*Y(8)/norm1;
dY(6, 1) = -(k_O0O1+k_O0C0)*Y(6)/norm1+k_O1O0*Y(7)/norm1+k_C0O0*Y(1)/norm1;

if ((time >= period*0.0) && (time <= period*1.0))
   stim_start = period*0.0;
elseif ((time >= period*1.0) && (time <= period*2.0))
   stim_start = period*1.0;
elseif ((time >= period*2.0) && (time < period*3.0))
   stim_start = period*2.0;
elseif ((time >= period*3.0) && (time < period*4.0))
   stim_start = period*3.0;
else
   stim_start = 0.0;
end;

local_time = time-stim_start;

if (local_time < t_ICCpeak)
   I_stim = Gcouple*(Y(36)-(V_ICCrest+V_ICCamp*local_time/f_2));
else
   I_stim = Gcouple*(Y(36)-(V_ICCrest+V_ICCamp*(1.0+exp(-f_1/(2.0*t_slope)))*1.0/(1.0+exp((local_time-f_2-0.5*f_1)/t_slope))));
end;

E_Ca_1 = R*T/(2.0*F)*log(Cao/Y(33));
i_Ca_L = g_CaL*Y(24)*(Y(36)-E_Ca_1);
T_correction_Ca_1 = Q10Ca1^((T-310.0)/10.0);
alpha2 = T_correction_Ca_1*0.731*exp(Y(36)/30.0);
beta2 = T_correction_Ca_1*0.2149*exp(-Y(36)/40.0);
alpha_0 = 4.0*alpha2;
alpha_1 = 3.0*alpha2;
alpha_2 = 2.0*alpha2;
alpha_3 = 1.0*alpha2;
beta_0 = 1.0*beta2;
beta_1 = 2.0*beta2;
beta_2 = 3.0*beta2;
beta_3 = 4.0*beta2;
phi_f = T_correction_Ca_1*0.4742*exp(Y(36)/10.0);
phi_s = T_correction_Ca_1*0.05956*exp(-Y(36)/40.0);
xi_f = T_correction_Ca_1*0.01407*exp(-Y(36)/300.0);
xi_s = T_correction_Ca_1*0.01213*exp(Y(36)/500.0);
psi_f = T_correction_Ca_1*0.02197*exp(Y(36)/500.0);
psi_s = T_correction_Ca_1*0.00232*exp(-Y(36)/280.0);
omega_f = beta_3*xi_f*phi_f/(alpha_3*psi_f);
omega_s = beta_3*xi_s*phi_s/(alpha_3*psi_s);
omega_sf = xi_s*psi_f/xi_f;
omega_fs = psi_s;
delta = T_correction_Ca_1*0.01;
theta = T_correction_Ca_1*4.0/(1.0+1.0/Y(33));
norm2 = Y(17)+Y(15)+Y(13)+Y(11)+Y(18)+Y(16)+Y(14)+Y(12)+Y(24)+Y(19)+Y(22)+Y(20)+Y(23)+Y(21);
dY(17, 1) = alpha_2*Y(15)/norm2-(alpha_3+beta_2+phi_s+phi_f+theta)*Y(17)/norm2+beta_3*Y(24)/norm2+omega_f*Y(20)/norm2+omega_s*Y(22)/norm2+delta*Y(18)/norm2;
dY(15, 1) = alpha_1*Y(13)/norm2-(alpha_2+beta_1+theta)*Y(15)/norm2+beta_2*Y(17)/norm2+delta*Y(16)/norm2;
dY(13, 1) = alpha_0*Y(11)/norm2-(alpha_1+beta_0+theta)*Y(13)/norm2+beta_1*Y(15)/norm2+delta*Y(14)/norm2;
dY(11, 1) = -(alpha_0+theta)*Y(11)/norm2+beta_0*Y(13)/norm2+delta*Y(12)/norm2;
dY(18, 1) = theta*Y(17)/norm2+alpha_2*Y(16)/norm2-(alpha_3+beta_2+phi_f+phi_s+delta)*Y(18)/norm2+beta_3*Y(19)/norm2+omega_f*Y(21)/norm2+omega_s*Y(23)/norm2;
dY(16, 1) = theta*Y(15)/norm2+alpha_1*Y(14)/norm2-(alpha_2+beta_1+delta)*Y(16)/norm2+beta_2*Y(18)/norm2;
dY(14, 1) = theta*Y(13)/norm2+alpha_0*Y(12)/norm2-(alpha_1+beta_0+delta)*Y(14)/norm2+beta_1*Y(16)/norm2;
dY(12, 1) = theta*Y(11)/norm2-(alpha_0+delta)*Y(12)/norm2+beta_0*Y(14)/norm2;
dY(24, 1) = alpha_3*Y(17)/norm2-(beta_3+psi_f+psi_s+theta)*Y(24)/norm2+xi_f*Y(20)/norm2+xi_s*Y(22)/norm2+delta*Y(19)/norm2;
dY(19, 1) = theta*Y(24)/norm2+alpha_3*Y(18)/norm2-(beta_3+psi_f+psi_s+delta)*Y(19)/norm2+xi_f*Y(21)/norm2+xi_s*Y(23)/norm2;
dY(22, 1) = phi_s*Y(17)/norm2+psi_s*Y(24)/norm2+omega_fs*Y(20)/norm2-(omega_s+xi_s+omega_sf+theta)*Y(22)/norm2+delta*Y(23)/norm2;
dY(20, 1) = phi_f*Y(17)/norm2+psi_f*Y(24)/norm2-(omega_f+xi_f+omega_fs+theta)*Y(20)/norm2+omega_sf*Y(22)/norm2+delta*Y(21)/norm2;
dY(23, 1) = theta*Y(22)/norm2+phi_s*Y(18)/norm2+psi_s*Y(19)/norm2+omega_fs*Y(21)/norm2-(omega_s+xi_s+omega_sf+delta)*Y(23)/norm2;
dY(21, 1) = theta*Y(20)/norm2+phi_f*Y(18)/norm2+psi_f*Y(19)/norm2-(omega_f+xi_f+omega_fs+delta)*Y(21)/norm2+omega_sf*Y(23)/norm2;
i_NCX = P_NCX*(exp(gamma*Y(36)*F/(R*T))*Y(35)^3.0*Cao-2.5*exp((gamma-1.0)*Y(36)*F/(R*T))*Nao^3.0*Y(33))/((K_mNai^3.0+Nao^3.0)*(K_mCa+Cao)*(1.0+k_sat*exp((gamma-1.0)*Y(36)*F/(R*T))));
T_correction_Na = 1.0*Q10Na^((T-297.0)/10.0);
k_I2I1 = T_correction_Na*0.0039239*exp(2.6793+0.0061468*Y(36));
k_I1O = T_correction_Na*0.12052*exp(-9.6028+0.083025*Y(36));
k_OC1 = T_correction_Na*2.391*exp(-13.335-0.25289*Y(36));
k_C1C2_2 = T_correction_Na*3.1566*exp(0.36352+0.077193*Y(36));
k_C2C3_2 = T_correction_Na*0.55432*exp(-0.099074+0.036441*Y(36));
k_C3C2_2 = T_correction_Na*0.00052548*exp(-0.069102+0.0031945*Y(36));
k_C2C1_2 = T_correction_Na*1.4496*exp(-0.1566+0.058353*Y(36));
k_C1O = T_correction_Na*1.5329*exp(0.0093193+0.041075*Y(36));
k_OI1 = T_correction_Na*1.6164*exp(0.30763+0.0060535*Y(36));
k_I1I2 = T_correction_Na*0.027735*exp(0.05149-0.046865*Y(36));
k_I1C1 = T_correction_Na*1.9046*exp(-2.484+0.020406*Y(36));
k_C1I1 = T_correction_Na*0.00021688*exp(-0.063438+0.0046683*Y(36));
norm3 = Y(30)+Y(25)+Y(26)+Y(27)+Y(28)+Y(29);
dY(27, 1) = -k_C3C2_2*Y(27)/norm3+k_C2C3_2*Y(26)/norm3;
dY(26, 1) = -(k_C2C1_2+k_C2C3_2)*Y(26)/norm3+k_C1C2_2*Y(25)/norm3+k_C3C2_2*Y(27)/norm3;
dY(25, 1) = -(k_C1C2_2+k_C1O+k_C1I1)*Y(25)/norm3+k_OC1*Y(30)/norm3+k_C2C1_2*Y(26)/norm3+k_I1C1*Y(28)/norm3;
dY(30, 1) = -(k_OC1+k_OI1)*Y(30)/norm3+k_C1O*Y(25)/norm3+k_I1O*Y(28)/norm3;
dY(29, 1) = -k_I2I1*Y(29)/norm3+k_I1I2*Y(28)/norm3;
dY(28, 1) = -(k_I1O+k_I1I2+k_I1C1)*Y(28)/norm3+k_I2I1*Y(29)/norm3+k_C1I1*Y(25)/norm3+k_OI1*Y(30)/norm3;
E_Ca_2 = R*T/(2.0*F)*log(Cao/Y(33));
i_Ca_T = g_CaT*Y(31)*Y(32)*(Y(36)-E_Ca_2);
T_correction_Ca_2 = Q10Ca2^((T-297.0)/10.0);
d_CaT_inf = 1.0/(1.0+exp(-(Y(36)+60.5)/5.3));
tau_d_CaT = 1.9058/T_correction_Ca_2;
dY(31, 1) = (d_CaT_inf-Y(31))/tau_d_CaT;
T_correction_Ca_3 = Q10Ca3^((T-297.0)/10.0);
f_CaT_inf = 1.0/(1.0+exp((Y(36)+75.5)/4.0));
tau_f_CaT = 0.38117*(8.6+14.7*exp(-(Y(36)+50.0)*(Y(36)+50.0)/900.0))/T_correction_Ca_3;
dY(32, 1) = (f_CaT_inf-Y(32))/tau_f_CaT;
V_myo = 3.5e-12;
E_Na_2 = R*T/F*log(Nao/Y(35));
i_Na = g_Na*Y(30)*(Y(36)-E_Na_2);
E_Na_1 = R*T/F*log(Nao/Y(35));
i_NsNa = g_NsNa*(Y(36)-E_Na_1);
i_NaK = P_NaK*Ko*Y(35)/((Ko+K_mK)*(K_mNa+Y(35))*(1.0+0.1245*exp(-0.1*Y(36)*F/(R*T))+0.0353*exp(-Y(36)*F/(R*T))));
dY(35, 1) = -1.0/1.0e15*(i_Na+i_NsNa+i_NCX*3.0+i_NaK*3.0)/(V_myo*F);
E_K_3 = R*T/F*log(Ko/Y(34));
i_Kv = g_Kv*Y(37)*Y(38)*(Y(36)-E_K_3);
E_K_2 = R*T/F*log(Ko/Y(34));
i_NsK = g_NsK*(Y(36)-E_K_2);
dY(34, 1) = -1.0e-15*(i_Kv+i_BK+i_NsK+I_stim-i_NaK*2.0)/(V_myo*F);
dY(33, 1) = -1.0e-15*(i_Ca_L+i_Ca_T-i_NCX*2.0)/(2.0*V_myo*F)/(1.0+n_CRT*CRT_total*K_D_CRT*Y(33)^(n_CRT-1.0)/(Y(33)^n_CRT+K_D_CRT)^2.0+n_CaM*CaM_total*K_D_CaM*Y(33)^(n_CaM-1.0)/(Y(33)^n_CaM+K_D_CaM)^2.0);
dY(36, 1) = -1.0/Cm*(i_Na+i_Ca_L+i_Ca_T+i_Kv+i_BK+i_NCX+i_NaK+i_NsK+i_NsNa+I_stim);
T_correction_K_1 = Q10K1^((T-297.0)/10.0);
x_Kv_inf = 1.0/(1.0+exp(-(Y(36)+43.0)/17.36));
tau_x_Kv = 4.7803/T_correction_K_1;
dY(37, 1) = (x_Kv_inf-Y(37))/tau_x_Kv;
T_correction_K_2 = Q10K2^((T-297.0)/10.0);
y_Kv_inf = 1.0/(1.0+exp((Y(36)-44.9)/12.0096));
tau_y_Kv = 763.7564/T_correction_K_2;
dY(38, 1) = (y_Kv_inf-Y(38))/tau_y_Kv;

return

%===============================================================================
% End of file
%===============================================================================
