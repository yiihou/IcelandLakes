% ODE MonteCarlo for HVT
clear
BiogenicSi = readtable ('IcelandBiogenicSi.xlsx');
tspan = [0 100];

Q_0 = 88.9; % m3/s
Q_0 = Q_0 * 60 * 60 * 24 * 365; % m3/yr


Jb_0 = mean(BiogenicSi.bSiFlux(9:15)); % g/cm2/yr
Jb_0 = Jb_0 * 10000; % g/m2/yr

lake_area = 28.9; % km^2
Jb_0 = Jb_0 * lake_area*1000*1000; % g/yr for whole lake
Jb_0 = Jb_0./67; % mol/yr, assuming molar mass of bSi is 67g/mol, 10% water content 

Volume = 6.31E8; % m3, calculated from lake bathymetry

C_lake = 0.2; % mmol/L
C_lake = C_lake ./ 1000; % mol/L
C_lake = C_lake .* 1000; % mol/m3

M_Si_0 = Volume * C_lake; % mol

y0 = [M_Si_0 2]; % M_Si_0 (mol), d30Si_0 (unitless)

delta_in = 1.5; % Si isotope value of river inputs
epsilon = -1.5; % diatom fractionation factor
f_Q = 1 * 2*pi; % period = 1 year
f_bSi = 1 * 2*pi; % period = 1 year

N = 1E4; % number of simulation
rng("default")
% r = a + (b-a).*rand(N,1).
A_bSi = 0.1 + (0.9-0.1).* rand(N,1);
A_Q = 0.1 + (0.9-0.1).* rand(N,1);
tau = 0 + (0.15-0).*rand(N,1);
a = 0.9 + (1.6-0.9).*rand(N,1);
b = -0.15 + (-0.05-(-0.15)).*rand(N,1);
F_true = nan(N,1);
F_calculated = nan(N,1);

t = linspace(4,5,10000);
for i = 1:N
    % f_bSi(i)= (mod(i,2)+1)*0.5*2*pi; % half 2pi, half 1pi, 1 or 2 blooms per year
    sol = lakemodel(tspan, y0, Q_0, A_Q(i),A_bSi(i),Volume,...
        Jb_0,delta_in,epsilon,tau(i),f_Q,f_bSi, a(i),b(i));
    if max(sol.x) == 100
        Q_t = Q_0 *( 1 + A_Q(i) *sin(f_Q * t)); % discharge time series 
        Q_Si_t = Q_t.* a(i) .*Q_t.^ b(i); % Si river input time series
    
        J_t = Jb_0 * (1+A_bSi(i)*sin(f_bSi*(t-tau(i)))); % bSi production time series

        d30Si_lake = deval(sol,t,2); % lake Si isotope ratio time series
        d30Si_bSi = d30Si_lake + epsilon; % diatom Si isotope ratio time series

        True_prod = mean(J_t); 
        F_true (i) = 1 - True_prod./mean(Q_Si_t); % True fraction nutrient remaining

        d30Si_bSi_ave = sum(J_t.* d30Si_bSi)./sum(J_t); % weighted average bSi Si isotope ratio
        F_calculated (i) = (d30Si_bSi_ave -delta_in)./epsilon; % calculated fraction of nutrient remaining

    end
end

%%
g = [185,146,196]./256; % glacial
ng = [134,189,112]./256; % non-glacial
%%
F_true(or(F_true<0, F_true>1)) = NaN;
F_calculated(or(F_calculated<0, F_calculated>1)) = NaN;
N_fail = sum(isnan(F_true));
diff = abs(F_true-F_calculated);
fraction_not_open = length(find(diff>0.05))./(N-N_fail);
%%
figure ('Renderer', 'painters','Position', [10 10 500 400])
histogram(abs(F_true-F_calculated),'Normalization','pdf',...
    'FaceColor',g,'FaceAlpha',1,'LineWidth',1.5)
hold on
set(gca,'Xlim',[0 0.15])
set(gca,'Ylim',[0 45])
xticks([0  0.05 0.1 0.15])
yticks([0 20 40])
makepretty_axes('abs(F_{true} - F_{calculated})','pdf')

%%

figure ('Renderer', 'painters','Position', [10 10 500 400])
histogram(tau(find(diff>0.05)),'Normalization','pdf',...
    'FaceColor',g,'FaceAlpha',1,'LineWidth',1.5)
xticks([0 0.0769 0.153])
xticklabels({'0','4','8'})
set(gca,'Xlim',[0 0.153])
makepretty_axes('\tau (week)','pdf')


