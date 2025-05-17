%% Calculate weathering flux from bSi data
%% CLOSED Rayleigh

% Read data
clear
addpath('/Users/yihou/Documents/IcelandLakePaperRepo/york_curve_fit/york_curve_fit_0_01')
% Read in chemistry data
Dissolved = readtable ('IcelandRiverLakeDissolved.xlsx');
BiogenicSi = readtable ('IcelandBiogenicSi.xlsx');

%% Regress River Array with York fit
d30_Si_err = max(Dissolved.d30Si_err(Dissolved.Type == "river"),0.11);
% regress only river samples
[a, b, sigma_a, sigma_b,b_save] = york_fit(Dissolved.Ge_Si_pM_uM_(Dissolved.Type =="river")',...
    Dissolved.d30Si(Dissolved.Type == "river")',...
    Dissolved.Ge_Si_err(Dissolved.Type == "river")',...
    d30_Si_err',0);

%% Inversion (HAK)
% For HAK, no Ge/Si fractionation for diatoms
N = 1e5; % number of iteration
% initialize arrays
HAK_GeSi_dist = nan(8,N); 
HAK_d30Si_dist = nan(8,N);
HAK_d30Si_infer = nan(8,N);
HAK_fractionation_dist = nan(8,N);
HAK_frac_lim = nan(1,8);
HAK_F_remain_real = nan(8,N);

for i = 1: 8  
     HAK_GeSi_dist (i, :) = normrnd(BiogenicSi.Ge_Si_pM_uM_(i),...
         BiogenicSi.Ge_Si_err(i)./2,...
         [1 N]); 
     HAK_d30Si_dist (i,:)= normrnd(BiogenicSi.d30Si(i),...
        max(BiogenicSi.d30Si_err(i),0.11)./2,...
        [1 N]);    
    
    a_dist = normrnd(a,sigma_a,[1 N]);
    b_dist = normrnd(b,sigma_b,[1 N]);    
    HAK_d30Si_infer(i, :) = HAK_GeSi_dist(i,:) .* b_dist + a_dist;  
    % -1.1+/-0.41 fractionation factor De La Rocha 1997
    HAK_fractionation_dist(i,:) = normrnd(-1.1,0.41,[1 N]);
    
    %----------------- truncate fractionation factor dist -----------------
    %--- HAK_frac_lim(i) = max (HAK_d30Si_dist (i,:) -  HAK_d30Si_infer(i, :));
    %--- HAK_fractionation_dist(i,:) = -1.51 + (HAK_frac_lim(i) + 1.51).*rand(1,N);  
    %------------------(don't use, 2023/8/9)------------------------------- 
    
    % HAK_F_remain_real (i,:) = (HAK_d30Si_dist (i,:) -  HAK_d30Si_infer(i, :))./...
    %     HAK_fractionation_dist(i,:);    
end
options = optimset('Display','off');
for i = 1:8
    for j = 1:N
        HAK_F_remain_real (i,j) = fsolve(@(x)HAK_d30Si_infer(i,j)-HAK_d30Si_dist (i,j)-...
            HAK_fractionation_dist(i,j).*(x.*log(x)/(1-x)),0.5, options);
    end
end


HAK_F_remain = HAK_F_remain_real;
HAK_F_remain_real(HAK_F_remain_real>1)=NaN;
HAK_F_remain_real(HAK_F_remain_real<0)=NaN;

HAK_F_ave = mean(HAK_F_remain_real,2, "omitnan");
HAK_F_std = std(HAK_F_remain_real,0,2, "omitnan");

HAK_lakearea = 3.2848;
HAK_catcharea = 172; % km^2
conversion = (28/60) ./ 28.0855 * HAK_lakearea * (100000^2); 
HAK_W_flux = BiogenicSi.bSiFlux(BiogenicSi.Location=="Efri Haukadalsá") ./(1-HAK_F_ave) .* conversion./HAK_catcharea; % in mol Si yr^1 km-2
HAK_W_flux_low = BiogenicSi.bSiFlux(BiogenicSi.Location=="Efri Haukadalsá") ./(1-(HAK_F_ave-HAK_F_std)).* conversion./HAK_catcharea;
HAK_W_flux_high = BiogenicSi.bSiFlux(BiogenicSi.Location=="Efri Haukadalsá") ./(1-(HAK_F_ave+HAK_F_std)).* conversion./HAK_catcharea;


%% Inverstion (HVT)
% for HVT, Ge/Si fractionation for the youngest 3 samples
HVT_GeSi_dist = nan(7,N); HVT_d30Si_dist = nan(7,N); HVT_d30Si_infer = nan(7,N);
HVT_fractionation_dist = nan(7,N); HVT_frac_lim = nan(1,8); HVT_F_remain_real = nan(7,N);

% older sample, no Ge/Si fractionation
for i = 12: 15
    j =i-8;
    HVT_GeSi_dist (j, :) = normrnd(BiogenicSi.Ge_Si_pM_uM_(i),...
        BiogenicSi.Ge_Si_err(i)./2,...
        [1 N]);    
    HVT_d30Si_dist (j,:)= normrnd(BiogenicSi.d30Si(i),...
        max(BiogenicSi.d30Si_err(i),0.11)./2,...
        [1 N]);
    a_dist = normrnd(a,sigma_a,[1 N]); b_dist = normrnd(b,sigma_b,[1 N]);
    
    HVT_d30Si_infer(j, :) = HVT_GeSi_dist(j,:) .* b_dist + a_dist;
    
    % -1.1+/-0.41 fractionation factor De La Rocha 1997
    HVT_fractionation_dist(j,:) = normrnd(-1.1,0.41,[1 N]);
    
    %--------------------------------------------------------------------
    %-- HVT_frac_lim(j) = max (HVT_d30Si_dist (j,:) -  HVT_d30Si_infer(j, :));
    %-- HVT_fractionation_dist(j,:) = -1.51 + (HVT_frac_lim(j) +1.51).*rand(1,N);
end    %--------------------------------------------------------------------
    
% %     HVT_F_remain_real (j,:) = (HVT_d30Si_dist (j,:) -  HVT_d30Si_infer(j, :))./...
% %         HVT_fractionation_dist(j,:);    

% Younger samples (1-3) where Ge/Si fractionation is large
offset = nan(3,N);
for i = 9:11
    j = i-8;
    offset(j,:) = abs(normrnd(BiogenicSi.Ge_Si_pM_uM_(9),...
        BiogenicSi.Ge_Si_err(9)./2,...
        [1 N]) - normrnd(Dissolved.Ge_Si_pM_uM_(30),...
        Dissolved.Ge_Si_err(30)./2,...
        [1 N]));   
    
    HVT_GeSi_dist (j, :) = normrnd(BiogenicSi.Ge_Si_pM_uM_(i),...
        BiogenicSi.Ge_Si_err(i)./2,...
        [1 N]);    
    HVT_d30Si_dist (j,:)= normrnd(BiogenicSi.d30Si(i),...
        max(BiogenicSi.d30Si_err(i),0.11)./2,...
        [1 N]);
    a_dist = normrnd(a,sigma_a,[1 N]); b_dist = normrnd(b,sigma_b,[1 N]);
    
    HVT_d30Si_infer(j, :) = (HVT_GeSi_dist(j,:)+offset(j,:)) .* b_dist...
        + a_dist;
    
    % -1.1+/-0.41 fractionation factor De La Rocha 1997
    HVT_fractionation_dist(j,:) = normrnd(-1.1,0.41,[1 N]);
    
    %----------------------------------------------------------------------
    %--- HVT_frac_lim(j) = max (HVT_d30Si_dist (j,:) -  HVT_d30Si_infer(j, :));
    %--- HVT_fractionation_dist(j,:) = -1.51 + (HVT_frac_lim(j) +1.51).*rand(1,N);
    %----------------------------------------------------------------------
end    
% %     HVT_F_remain_real (j,:) = (HVT_d30Si_dist (j,:) -  HVT_d30Si_infer(j, :))./...
% %         HVT_fractionation_dist(j,:);  

options = optimset('Display','off');
for i = 1:7
    for j = 1:N
        HVT_F_remain_real (i,j) = fsolve(@(x)HVT_d30Si_infer(i,j)-HVT_d30Si_dist (i,j)-...
            HVT_fractionation_dist(i,j).*(x.*log(x)/(1-x)),0.5, options);
    end
end

HVT_F_remain = HVT_F_remain_real;
HVT_F_remain_real(HVT_F_remain_real>1)=NaN;
HVT_F_remain_real(HVT_F_remain_real<0)=NaN;

HVT_F_ave = mean(HVT_F_remain_real,2, "omitnan");
HVT_F_std = std(HVT_F_remain_real,0,2, "omitnan");

HVT_lakearea = 29.8206;
HVT_catcharea = 820;
conversion = (28/60) ./ 28.0855 * HVT_lakearea * (100000^2); 
HVT_W_flux = BiogenicSi.bSiFlux(BiogenicSi.Location=="Hvíta Catchment") ./(1-HVT_F_ave) .* conversion ./HVT_catcharea; % in mol Si yr^1 km^-2
HVT_W_flux_low = BiogenicSi.bSiFlux(BiogenicSi.Location=="Hvíta Catchment") ./(1-(HVT_F_ave-HVT_F_std)).* conversion ./HVT_catcharea;
HVT_W_flux_high = BiogenicSi.bSiFlux(BiogenicSi.Location=="Hvíta Catchment") ./(1-(HVT_F_ave+HVT_F_std)).* conversion ./HVT_catcharea;

%% Save workspace to .mat file and write to excel
WFlux_closed = table (BiogenicSi.Sample_ID,...
    BiogenicSi.Age_BP2k_,...
    BiogenicSi.bSiFlux,...
    cat(1,HAK_F_ave, HVT_F_ave),...
    cat(1,HAK_F_std, HVT_F_std),...
    cat(1,HAK_W_flux,HVT_W_flux),...
    cat(1,HAK_W_flux_low,HVT_W_flux_low),...
    cat(1,HAK_W_flux_high,HVT_W_flux_high));
    
WFlux_closed.Properties.VariableNames = ["Sample_ID","Age_BP2k_","bSiFlux",...
    "f_remain","f_remain_std",...
    "W_ave","W_low","W_high"];

writetable(WFlux_closed,'Iceland_WFlux_closed.xlsx')
