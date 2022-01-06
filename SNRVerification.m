clear all
close all
import HypWright.*
import HypWright.Models.*
clc
%% Variable Setup
N = 60;
TR = 2;
TRList = (0:(N-1))*TR;
timeAx = (0:(N-1)).*TR;
T1a = 43;
T1b = 33;
Kpl = 0.1;
alpha = 2.5;
beta = 4.5;
noise_fact_vect = [1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2];
for j = 1:numel(noise_fact_vect)
noiseFact = noise_fact_vect(j);
nAverages = 10000;
M0 = [0,0];
kve = 0.02;
ve = 0.95;
VIFScaleFact = [1;0];
opts = optimset('lsqcurvefit');
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';
params = struct('t0',[0;0],'gammaPdfA',[alpha;1],'gammaPdfB',[beta;1],...
    'scaleFactor',VIFScaleFact,'T1s',[T1a,T1b],'ExchangeTerms',[0,Kpl;0,0],...
    'TRList',TRList,'PerfusionTerms',[kve,0],'volumeFractions',ve,...
    'fitOptions', opts);
A = MultiPoolTofftsGammaVIF();
%% Choose Excitation Angle
flipAngles = zeros(2,length(TRList))+20*pi/180;
params.FaList = flipAngles;
%% Fitting
[~,Mxy,~] = A.compile(M0.',params);
for k = 1:nAverages
    noise = noiseFact.*(rand(size(Mxy))-0.5);
    SNR(j,k) = max(max(Mxy))/mean(std(noise));
end
end
errorbar(noise_fact_vect,mean(SNR,2),std(SNR,[],2))
set(gca, 'XScale','log', 'YScale','log')
grid on
xlabel('noise factor (AU)')
ylabel('SNR (AU)')
