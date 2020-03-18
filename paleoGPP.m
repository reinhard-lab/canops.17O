%
function GPP_paleo = paleoGPP(D17O_sulfate,f_O_sulfate,PALO2,PALCO2,X_h,X_l,rho_0,strat_co2_m,strat_co2_b,gamma,theta)
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% +++ mass balance calculation of biospheric productivity based on 17O ++ %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% NOTE: see Cao + Bao [2013] and Crockford et al. [2018]
%
%% set some input parameters
%
pO2_ref          = 209500; % modern atmospheric pO2 [ppmv]
pCO2_ref         = 280;    % modern atmospheric pCO2 [ppmv]
X_h_ref          = 146;    % d18O of CO2_strat relative to O2_strat at high O2 [permil]
X_l_ref          = 64;     % d18O of CO2_strat relative to O2_strat at low O2 [permil]
rho_0_ref        = 1.23;   % atmospheric O2/CO2 ratio
strat_co2_m_ref  = 0.5167; % slope of 17O relationship between O2_strat and CO2_strat
strat_co2_b_ref  = -8.025; % intercept of 17O relationship between O2_strat and CO2_strat
gamma_ref        = 0.0426; % stratosphere-troposphere exchange rate for O2 [yr-1]
theta_ref        = 0.1156; % active stratospheric O2 fraction [dimensionless]
tau_ref          = 1244.0; % O2 residence time in preindustrial atmosphere [yr]
D17_O2_modern    = -0.546; % calcuated from data in Pack et al (2017)
GPP_mod          = 220.0;  % modern gross primary production [GtC yr-1]
%
%% estimate biospheric productivity based on assumed pCO2, pO2, and D17O_sulfate
%
rho_ref    = pO2_ref / pCO2_ref ;
%
rho        = rho_ref * PALO2 ./ PALCO2 ; 
%
D17O_O2    = D17O_sulfate ./ f_O_sulfate ;
%
d18O_CO2   = (X_l + X_h .* (rho ./ rho_0)) ./ (1. + (rho ./ rho_0));
D17O_CO2   = strat_co2_m .* d18O_CO2 + strat_co2_b;
%
tau        = D17O_O2 .* (1. + rho) ./ (-1. * (D17O_CO2 + D17O_O2) .* gamma .* theta) ;
%
GPP_scaled = PALO2 ./ (tau / tau_ref);
%
GPP_paleo  = GPP_scaled.*GPP_mod;
%
end

