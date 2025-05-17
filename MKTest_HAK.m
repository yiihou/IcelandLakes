%% Conduct Mann-Kendall Test for weathering fluxes at HAK
% 
clear
addpath('/Users/yihou/Documents/IcelandLakePaperRepo/Mann_Kendall')

WFlux = readtable ('Iceland_WFlux.xlsx');

N = 1E5;

%% resample weathering flux at HAK from f_remain
HAK_lakearea = 3.2848; % km^2
HAK_catcharea = 172; % km^2
conversion = (28/60) ./ 28.0855 * HAK_lakearea * (100000^2); 
Resample_f = nan(8,N);
Re_W = nan(8,N);
for i = 1:8
    Resample_f (i,:) = normrnd(WFlux.f_remain(i),WFlux.f_remain_std(i)./2,[1,N]);
    Resample_f(Resample_f>=1)=nan;
    Re_W (i,:) = WFlux.bSiFlux(i)./(1-Resample_f(i,:)).*conversion./HAK_catcharea;
end
H = nan(1,N);
alpha = 0.05;
for i = 1:N
    if sum(isnan(Re_W(:,i)))==0
        [H(i),p_value(i)]=Mann_Kendall(Re_W(:,i),alpha);
    end
end

% calculate the percentage of tests that have significant value < alpha
M = p_value <alpha;
Percent = sum(M)/N;

disp('Percent of tests that have monotonic trend: ')
disp(Percent)