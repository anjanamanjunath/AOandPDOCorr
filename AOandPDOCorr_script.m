%% load data
load("AO_data.txt") %Arctic Oscillation
load("meiv2.data.txt") %ENSO
load("PDO.txt") %Pacific Decadal Oscillation
load("Sc_county_precipitation.txt") %SC precipitation

%% reformat oscillation index data and put all the indices within the same time frame:
% latest earliest start date 1979, and the end date: 2/2021

AODataClean = [datenum(AO_data(:,1), AO_data(:,2), 1), AO_data(:,3)];
AODataClean = AODataClean(349:end-8,:);

ENSODataClean = [];
cntr = 0;
for i = 1:length(meiv2_data)
    for j = 1:12
        cntr = cntr+1;
        ENSODataClean(cntr,1) = datenum(meiv2_data(i,1), j, 1);
        ENSODataClean(cntr,2) = meiv2_data(i,j+1);
    end
end
ENSODataClean = ENSODataClean(1:end-10, :);

PDODataClean = [];
cntr = 0;
for i = 1:length(PDO)
    for j = 1:12
        cntr = cntr+1;
        PDODataClean(cntr,1) = datenum(PDO(i,1), j, 1);
        PDODataClean(cntr,2) = PDO(i,j+1);
    end
end
PDODataClean = PDODataClean(373:end-10,:);


%creating one matrix with all the data
projectData = [Sc_county_precipitation AODataClean(:,2) ENSODataClean(:,2) PDODataClean(:,2)];

%going from monthly to annual signals
precipAnnual = [];
AOAnnual = [];
ENSOAnnual = [];
PDOAnnual = [];
year = []; 

cntr = 0;
for i = 9:12:492
    cntr = cntr + 1;
    year(cntr, 1) = cntr + 1979 - 1; 
    precipAnnual(cntr,1) = sum(projectData(i:i+12,2))/12;
    AOAnnual(cntr,1) = max(projectData(i:i+12,3),[],'ComparisonMethod','abs'); 
    ENSOAnnual(cntr,1) = max(projectData(i:i+12,4),[],'ComparisonMethod','abs'); 
    PDOAnnual(cntr,1) = max(projectData(i:i+12,5),[],'ComparisonMethod','abs'); 
end 

projectDataAnnual = [year precipAnnual AOAnnual ENSOAnnual PDOAnnual];

%% plotting all the oscillation index datasets as a function of time and then the Santa Cruz County Precipitation
figure(1)
subplot(3, 1, 1);
hold on;
plot(year, AOAnnual, 'b -o', year, precipAnnual, 'm -*')
xlabel("Time (year)")
title("AO index and Santa Cruz County Precipitation")

subplot(3, 1, 2);
hold on;
plot(year, ENSOAnnual, 'b -o', year, precipAnnual, 'm -*')
xlabel("Time (year)")
title("ENSO index and Santa Cruz County Precipitation")

subplot(3, 1, 3);
hold on;
plot(year, PDOAnnual, 'b -o', year, precipAnnual, 'm -*')
xlabel("Time (year)")
title("PDO index and Santa Cruz County Precipitation")
%% plotting the precipitation data and then applying a fourier transform on it

SCPrecip = fft(detrend(Sc_county_precipitation(:,2)));

N = length(SCPrecip);
dt = 1/12; %units: years
f = [0:N-1]/(N*dt); % units: cycles/year

figure(2)
plot(f(1:end/2), abs(SCPrecip(1:end/2)));
xlabel("Frequency (cycles per year)")

%Taking a look at the frequencies of the Oscillation Indicies themselves,
%we repeat the same thing: 

AO = fft(AODataClean(:,2)); 
PDOFourier = fft(PDODataClean(:,2)); 
ENSO = fft(ENSODataClean(:,2)); 

figure(3)
subplot(3,1,1)
plot(f(1:end/2), abs(AO(1:end/2)), 'b');
xlabel("Frequency (cycles per year)")

subplot(3, 1, 2)
plot(f(1:end/2), abs(PDOFourier(1:end/2)), 'b');
xlabel("Frequency (cycles per year)")

subplot(3, 1, 3)
plot(f(1:end/2), abs(ENSO(1:end/2)), 'b');
xlabel("Frequency (cycles per year)")
%% Determining Correlation between the Oscillation Indicies
% first demean the data
precipAnnual = precipAnnual - mean(precipAnnual);
AOAnnual = AOAnnual - mean(AOAnnual); 
ENSOAnnual = ENSOAnnual - mean(ENSOAnnual); 
PDOAnnual = PDOAnnual - mean(PDOAnnual); 

% comparing the Arctic Oscillation with SC Precipitation
AOPosPhaseLog = AOAnnual > 0;
AOPosPhase = AOPosPhaseLog .* AOAnnual;
AONegPhaseLog = AOAnnual < 0;
AONegPhase = AONegPhaseLog .* AOAnnual;

[Phi, lags] = xcorr(precipAnnual, AOPosPhase, 'unbiased');
Phi0x = sum(AOPosPhase.^2);
Phi0y = sum(precipAnnual.^2);
PhiNorm = Phi./(sqrt(Phi0x + Phi0y));
max(PhiNorm) % 0.2477

[Phi2, lags2] = xcorr(precipAnnual, AONegPhase, 'unbiased');
Phi0x2 = sum(AONegPhase.^2);
PhiNorm2 = Phi2./(sqrt(Phi0x2 + Phi0y));
max(PhiNorm2) % 0.1540

[Phi3, lags3] = xcorr(precipAnnual, AOAnnual, 'unbiased');
Phi0x3 = sum(AOAnnual.^2);
PhiNorm3 = Phi3./(sqrt(Phi0x3 + Phi0y));
max(PhiNorm3) % 0.1724

%plotting these correlations as a function of lag
figure(4)
subplot(3, 1, 1);
hold on;
plot(lags, PhiNorm, 'b')
title('lags v. normalized correlation for AO Positive Phases')

subplot(3, 1, 2);
hold on;
plot(lags2, PhiNorm2, 'b')
title('lags v. normalized correlation for AO Negative Phases')

subplot(3, 1, 3);
hold on;
plot(lags3, PhiNorm3, 'b')
title('lags v. normalized correlation for AO Original Index')

% comparing the ENSO Index with SC Precipitation - repeating the same procedure
ENSOPosPhaseLog = ENSOAnnual > 0;
ENSOPosPhase = ENSOPosPhaseLog .* ENSOAnnual;
ENSONegPhaseLog = ENSOAnnual < 0;
ENSONegPhase = ENSONegPhaseLog .* ENSOAnnual;

[Phi, lags] = xcorr(precipAnnual, ENSOPosPhase, 'unbiased');
Phi0x = sum(ENSOPosPhase.^2);
Phi0y = sum(precipAnnual.^2);
PhiNorm = Phi./(sqrt(Phi0x + Phi0y));
max(PhiNorm) % 0.1310

[Phi2, lags2] = xcorr(precipAnnual, ENSONegPhase, 'unbiased');
Phi0x2 = sum(ENSONegPhase.^2);
PhiNorm2 = Phi2./(sqrt(Phi0x2 + Phi0y));
max(PhiNorm2) % 0.0587

[Phi3, lags3] = xcorr(precipAnnual, ENSOAnnual, 'unbiased');
Phi0x3 = sum(ENSOAnnual.^2);
PhiNorm3 = Phi3./(sqrt(Phi0x3 + Phi0y));
max(PhiNorm3) % 0.1309

figure(5)
subplot(3, 1, 1);
hold on;
plot(lags, PhiNorm, 'b')
title('lags v. normalized correlation for ENSO Positive Phases')

subplot(3, 1, 2);
hold on;
plot(lags2, PhiNorm2, 'b')
title('lags v. normalized correlation for ENSO Negative Phases')

subplot(3, 1, 3);
hold on;
plot(lags3, PhiNorm3, 'b')
title('lags v. normalized correlation for ENSO Original Index')

% for the PDO values

PDOPosPhaseLog = PDOAnnual > 0;
PDOPosPhase = PDOPosPhaseLog .* PDOAnnual;
PDONegPhaseLog = PDOAnnual < 0;
PDONegPhase = PDONegPhaseLog .* PDOAnnual;

[Phi, lags] = xcorr(precipAnnual, PDOPosPhase, 'unbiased');
Phi0x = sum(PDOPosPhase.^2);
Phi0y = sum(precipAnnual.^2);
PhiNorm = Phi./(sqrt(Phi0x + Phi0y));
max(PhiNorm) % 0.0904


[Phi2, lags2] = xcorr(precipAnnual, PDONegPhase, 'unbiased');
Phi0x2 = sum(PDONegPhase.^2);
PhiNorm2 = Phi2./(sqrt(Phi0x2 + Phi0y));
max(PhiNorm2) % 0.0570

[Phi3, lags3] = xcorr(precipAnnual, PDOAnnual, 'unbiased');
Phi0x3 = sum(PDOAnnual.^2);
PhiNorm3 = Phi3./(sqrt(Phi0x3 + Phi0y));
max(PhiNorm3) % 0.0729

figure(6)
subplot(3, 1, 1);
hold on;
plot(lags, PhiNorm, 'b')
title('lags v. normalized correlation for PDO Positive Phases')

subplot(3, 1, 2);
hold on;
plot(lags2, PhiNorm2 ,'b')
title('lags v. normalized correlation for PDO Negative Phases')

subplot(3, 1, 3);
hold on;
plot(lags3, PhiNorm3, 'b')
title('lags v. normalized correlation for PDO Original Index')
