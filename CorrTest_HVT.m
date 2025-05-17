%% Conduct Statistical Test for weathering fluxes at HVT
% 
clear

WFlux = readtable ('Iceland_WFlux.xlsx');
BiogenicSi = readtable ('IcelandBiogenicSi.xlsx');
A = readmatrix ('HVT_IceArea.txt');
IceArea.age = A(:,1); % in kyr
IceArea.area = A(:,2);
B = readmatrix('HVT_Q.txt');
Discharge.age = B(:,1); % in kyr
Discharge.q = B(:,2);
clearvars A B

N = 1E5;

%% find ice area and discharge corresponding to sample ages
index = nan(1,7);
index2 = nan(1,7);

for i = 9:15
    j = i-8; 
    [~, index(j)] = min(abs(BiogenicSi.Age_BP2k_(i) - IceArea.age .*1000));
    [~, index2(j)] = min(abs(BiogenicSi.Age_BP2k_(i) - Discharge.age .*1000));
end

%% resample weathering flux at HVT from f_remain

HVT_lakearea = 29.8206; % km^2
HVT_catcharea = 820; % km^2

conversion = (28/60) ./ 28.0855 * HVT_lakearea * (100000^2); 
Resample_f = nan(7,N);
Re_W = nan(7,N);

for j = 9:15
    i= j-8;
    Resample_f (i,:) = normrnd(WFlux.f_remain(j),WFlux.f_remain_std(j)./2,[1,N]);
    Resample_f(Resample_f>=1)=nan;
    Re_W (i,:) = WFlux.bSiFlux(j)./(1-Resample_f(i,:)).*conversion./HVT_catcharea;
end

%% Pearson Correlation
r= nan(1,N);
p = r;

alpha = 0.05; % significance level

% Ice area
for i = 1:N
    if sum(isnan(Re_W(:,i)))==0
        [r(i),p(i)] = corr(IceArea.area(index),Re_W(:,i),'Type','Pearson');        
    end
end
H1 = p < alpha;
Percent_IceArea = sum(H1)/N;

r= nan(1,N);
p = r;
% Discharge
for i = 1:N
    if sum(isnan(Re_W(:,i)))==0
       [r(i),p(i)] = corr(Discharge.q(index2),Re_W(:,i),'Type','Pearson');  
    end
end

H2 = p < alpha;
Percent_Discharge = sum(H2)/N;

% clearvars -except Percent*

disp('Percent of tests with significant correlation with Ice Area: ')
disp(Percent_IceArea)
disp('Percent of tests with significant correlation with Glacial Discharge: ')
disp(Percent_Discharge)
