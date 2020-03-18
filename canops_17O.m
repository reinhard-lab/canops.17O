%
function out = canops_17O(num_run,f_O_sulfate_min,f_O_sulfate_max,logPALO2_min,logPALO2_max,CO2_min,CO2_max,GPP_land,f_ex,p_fit)
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% +++ filtered inversion for estimating atmospheric pO2 based on 17O ++++ %
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% num_run         --> number of samples in Monte Carlo routine
% f_O_sulfate_min --> O atom incorporation into sulfate, minimum
% f_O_sulfate_max --> O atom incorporation into sulfate, maximum
% logPALO2_min    --> log atmospheric pO2, minimum [PAL]
% logPALO2_max    --> log atmospheric pO2, maximum [PAL]
% CO2_min         --> atmospheric pCO2, minimum [PAL]
% CO2_max         --> atmospheric pCO2, maximum [PAL]
% GPP_mod         --> modern global gross primary productivity [GtC y-1]
% GPP_land        --> assumed terrestrial gross primary productivity [GtC y-1]
% f_ex            --> export ratio relative to marine GPP
% p_fit           --> prediction interval for CANOPS filter
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% define default parameters if none are given
if exist('num_run','var')==0,         num_run = 100000; end
if exist('f_O_sulfate_min','var')==0, f_O_sulfate_min = 0.08; end
if exist('f_O_sulfate_max','var')==0, f_O_sulfate_max = 0.15; end
if exist('logPALO2_min','var')==0,    logPALO2_min = -4.0; end
if exist('logPALO2_max','var')==0,    logPALO2_max = 0.0; end
if exist('CO2_min','var')==0,         CO2_min = 2.0; end
if exist('CO2_max','var')==0,         CO2_max = 500.0; end
if exist('GPP_land','var')==0,        GPP_land = 1.74; end
if exist('f_ex','var')==0,            f_ex = 0.1; end
if exist('p_fit','var')==0,           p_fit = 0.9; end
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% initial isotope mass balance inversion
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% preallocate arrays
f_O_sulfate  = zeros(num_run,1);
PALO2        = zeros(num_run,1);
PALCO2       = zeros(num_run,1);
X_h          = zeros(num_run,1);
X_l          = zeros(num_run,1);
rho_0        = zeros(num_run,1);
strat_co2_m  = zeros(num_run,1);
strat_co2_b  = zeros(num_run,1);
gamma        = zeros(num_run,1);
theta        = zeros(num_run,1);
GPP_scaled   = zeros(num_run,1);
%
% import 17O data, if any 
D17O_sulfate_measured = importdata('_input/D17O_sulfate_measured.txt');
%
% sulfate data distribution from measured values
D17O_sulfate = min(D17O_sulfate_measured) * ones(num_run,1);
%
% fraction of O2 incorporated into sulfate
f_O_sulfate  = unifrnd(f_O_sulfate_min,f_O_sulfate_max,num_run,1);
%
% atmospheric pO2
logPALO2     = unifrnd(logPALO2_min,logPALO2_max,num_run,1);
PALO2        = 10.^(logPALO2);
%
% atmospheric pCO2
PALCO2       = unifrnd(CO2_min,CO2_max,num_run,1);
%
% THESE DISTRIBUTIONS DON'T CHANGE [see Crockford et al. 2018]
X_h          = 2.0    * randn(num_run,1) + 146.0;
X_l          = 5.5    * randn(num_run,1) + 64.0;
rho_0        = 0.325  * randn(num_run,1) + 1.23;
strat_co2_m  = 0.0336 * randn(num_run,1) + 0.5167;
strat_co2_b  = 3.673  * randn(num_run,1) - 8.052;
gamma        = 0.01   * randn(num_run,1) + 0.0426;
%
% active stratospheric O2 fraction assumed constant given an O3 layer
theta        = 0.1156 * ones(num_run,1);
%
% estimate GPP from isotope mass balance
GPP_paleo    = paleoGPP(D17O_sulfate,f_O_sulfate,PALO2,PALCO2,X_h,X_l,rho_0,strat_co2_m,strat_co2_b,gamma,theta);
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% filtered inversion results using CANOPS output
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% import CANOPS results [see Ozaki et al. 2019]
% atmospheric pO2 [PAL]
CANOPS_pO2    = importdata('CANOPS_out_pO2_PAL.txt');
% oceanic export production [GtC y-1]
CANOPS_EX     = importdata('CANOPS_out_EX_GtC.txt');
%
% calculate net primary productivity of marine biosphere based on export
GPP_sea       = CANOPS_EX./f_ex;
%
% calculate total biospheric NPP based on assumed terrestrial productivity
CANOPS_GPP    = GPP_sea+GPP_land;
%
% convert global GPP to PAL based on modern marine NPP
%GPP_PAL       = GPP_paleo./GPP_mod;
%
% construct filter based on assumed prediction interval (default is 90%)
fit_result    = fit(CANOPS_pO2,CANOPS_GPP,'poly1');
pO2_in        = logspace(-5,0,1000)';
filter        = predint(fit_result,pO2_in,p_fit);
%
% slope and intercept values for filter
m_filter      = fit_result.p1;
b_filter_low  = filter(1,1);
b_filter_high = filter(1,2);
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% MC filter based on CANOPS results
%
% preallocate filtered arrays
GPP_paleo_filtered = zeros(num_run,1);
PALO2_filtered     = zeros(num_run,1);
PALCO2_filtered    = zeros(num_run,1);
%
% initialize 
j=0;
%
% filter initial inversion output
for i=1:num_run
%
    GPP_low = m_filter*PALO2(i)+b_filter_low;
    GPP_up  = m_filter*PALO2(i)+b_filter_high;
%
    if (GPP_paleo(i)>=GPP_low) && (GPP_paleo(i)<=GPP_up)
        j = j+1;
        temp(j,:)=[GPP_paleo(i),PALO2(i),PALCO2(i)];        
    end
%    
end
%
% compile filtered arrays 
GPP_paleo_filtered  = temp(:,1);   
PALO2_filtered      = temp(:,2);  
PALCO2_filtered     = temp(:,3);  
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% do some downstream calculations
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
% central tendency and dispersion of filtered pO2 values
median_pO2 = median(PALO2_filtered);               % median filtered pO2
q90_pO2    = quantile(PALO2_filtered,[0.10 0.90]); % credible interval [90%]
%
% explicitly calculate filter bounds for plotting
pO2_in     = logspace(-4,0,1000);
GPP_low    = m_filter * pO2_in + b_filter_low;
GPP_up     = m_filter * pO2_in + b_filter_high;
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% save output
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
cd ~/Documents/MATLAB/canops.17O/
if ~exist('_output','dir')
    mkdir('_output')
end
%
cd _output
%
filename = strcat('17O_inversion_',datestr(now,'yyyymmdd'),'.mat');
save(filename);
cd ..
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
% make (and save) some plots
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ %
%
cd ~/Documents/MATLAB/canops.17O/
if ~exist('_figures','dir')
    mkdir('_figures')
end
%
cd _figures
%
figure(1);
%
edges=linspace(-1,0,10);
h=histogram(D17O_sulfate_measured,edges);
h.FaceColor='k';
h.EdgeColor='w';
xlim([-1 0]);
ylim([0 30]);
xlabel('D^{17}O [o/oo]');
ylabel('n');
grid on;
%
print(gcf,'-dpdf','Planavsky_Astrobiology_Fig_1.pdf');
%
figure(2);
%
subplot(1,3,[1,2]);
scatter(CANOPS_pO2,CANOPS_GPP,25,[0.5 0.5 0.5],'+');
set(gca,'xscale','log','yscale','log');
xlim([1.e-4 2.e-1]);
ylim([1.e0 1.e2]);
xlabel('pO_2 [PAL]');
ylabel('GPP [GtC y^{-1}]');
hold on;
plot(pO2_in,GPP_low,'r','LineWidth',1.5);
plot(pO2_in,GPP_up,'r','LineWidth',1.5);
grid on;
%
subplot(1,3,3);
edges=logspace(0,2,25);
h=histogram(CANOPS_GPP,edges);
set(gca,'yscale','log');
h.Normalization='probability';
h.FaceColor=[0.5 0.5 0.5];
h.EdgeColor='w';
h.Orientation='horizontal';
%xlim([1.e-4 2.e-1]);
ylim([1.e0 1.e2]);
xlabel('f(x)');
ylabel('GPP [GtC y^{-1}]');
grid on;
hold on;
%
print(gcf,'-dpdf','Planavsky_Astrobiology_Fig_4.pdf');
%
figure(3);
edges=logspace(-4,0,50);
h=histogram(PALO2_filtered,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
set(gca,'xscale','log');
xlabel('pO_2 [PAL]');
ylabel('f(x)');
%
hold on;
%
xline(median_pO2,'k','LineWidth',1.5);
xline(q90_pO2(:,1),'r','LineWidth',1.5);
xline(q90_pO2(:,2),'r','LineWidth',1.5);
grid on;
%
print(gcf,'-dpdf','Planavsky_Astrobiology_Fig_5.pdf');
%
figure (4);
%
subplot (3,2,1);
edges=logspace(-4,0,50);
h=histogram(PALO2,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
set(gca,'xscale','log');
xlabel('pO_2 [PAL]');
ylabel('f(x)');
title('unfiltered');
grid on;
hold on;
%
subplot(3,2,2);
edges=logspace(-4,0,50);
h=histogram(PALO2_filtered,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
set(gca,'xscale','log');
xlabel('pO_2 [PAL]');
ylabel('f(x)');
title('filtered');
%
hold on;
%
xline(median_pO2,'k','LineWidth',1.5);
xline(q90_pO2(:,1),'r','LineWidth',1.5);
xline(q90_pO2(:,2),'r','LineWidth',1.5);
grid on;
%
subplot(3,2,3);
edges=logspace(-2,4,50);
h=histogram(GPP_paleo,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
set(gca,'xscale','log');
xlabel('GPP [GtC y^{-1}]');
ylabel('f(x)');
title('unfiltered');
grid on;
hold on;
%
subplot(3,2,4);
edges=logspace(-2,4,50);
h=histogram(GPP_paleo_filtered,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
set(gca,'xscale','log');
xlabel('GPP [GtC y^{-1}]');
ylabel('f(x)');
title('filtered');
grid on;
hold on;
%
subplot(3,2,5);
edges=linspace(0,500,50);
h=histogram(PALCO2,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
xlabel('pCO_2 [PAL]');
ylabel('f(x)');
title('filtered');
grid on;
hold on;
%
subplot(3,2,6);
edges=linspace(0,500,50);
h=histogram(PALCO2_filtered,edges);
h.Normalization='probability';
h.FaceColor='k';
h.EdgeColor='w';
xlabel('pCO_2 [PAL]');
ylabel('f(x)');
title('filtered');
grid on;
hold on;
%
print(gcf,'-dpdf','filter_comparison.pdf');
%
end