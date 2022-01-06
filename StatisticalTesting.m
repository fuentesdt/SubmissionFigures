%% Run statistical tests on the VFA paper data
% Date: 8/23/2018
% Author: Chris Walker

clear all
close all
clc

%% Run t-test on the means
significance = 0.01;
% Run on Kpl when perfusion is known
load('tmpPerffusedPub2')
[~,p(1,:),~,~] = ttest(squeeze(fits(:,1,:)).',squeeze(fits(:,2,:)).');
[~,p(2,:),~,~] = ttest(squeeze(fits(:,2,:)).',squeeze(fits(:,3,:)).');
[~,p(3,:),~,~] = ttest(squeeze(fits(:,1,:)).',squeeze(fits(:,3,:)).');
subplot(2,2,1)
semilogy(N_vect,p,N_vect,zeros(size(N_vect))+significance,'g--')
ylim([1e-4,2]),xlim([5,60])
title('k_P_L difference t-test p value')
legend('CEA-bb Vs VEA-const','VEA-const Vs VEA-max','CEA-bb Vs VEA-max',...
    sprintf('p = %f',significance))
xlabel(' Number of Excitations '),ylabel('p value')
grid on
sprintf('Kpl Only Mean Kpl average: %f, Min: %f',mean(p(3,:)),min(p(3,:)))

%% Run 
% Run on Kpl when perfusion is known
load('tmpPerffusedPub2')
[~,p(1,:),~,~] = vartest2(squeeze(fits(:,1,:)).',squeeze(fits(:,2,:)).');
[~,p(2,:),~,~] = vartest2(squeeze(fits(:,2,:)).',squeeze(fits(:,3,:)).');
[~,p(3,:),~,~] = vartest2(squeeze(fits(:,1,:)).',squeeze(fits(:,3,:)).');
subplot(2,2,2)
semilogy(N_vect,p,N_vect,zeros(size(N_vect))+significance,'g--')
ylim([1e-4,2]),xlim([5,60])
title('k_P_L varianece F-test p value')
legend('CEA-bb Vs VEA-const','VEA-const Vs VEA-max','CEA-bb Vs VEA-max',...
    sprintf('p = %f',significance))
xlabel(' Number of Excitations '),ylabel('p value')
grid on
sprintf('Kpl Only Variance in Kpl average: %f, Min: %f',mean(p(3,:)),min(p(3,:)))

% Run on Kpl when perfusion is also fit
load('tmpPerffusedPubPerf2')
[~,p(1,:),~,~] = ttest(squeeze(fits(:,1,:,1)).',squeeze(fits(:,2,:,1)).');
[~,p(2,:),~,~] = ttest(squeeze(fits(:,2,:,1)).',squeeze(fits(:,3,:,1)).');
[~,p(3,:),~,~] = ttest(squeeze(fits(:,1,:,1)).',squeeze(fits(:,3,:,1)).');
subplot(2,2,3)
semilogy(N_vect,p,N_vect,zeros(size(N_vect))+significance,'g--')
ylim([1e-4,2]),xlim([5,60])
title('k_P_L difference t-test p value')
legend('CEA-bb Vs VEA-const','VEA-const Vs VEA-max','CEA-bb Vs VEA-max',...
    sprintf('p = %f',significance))
xlabel(' Number of Excitations '),ylabel('p value')
grid on
sprintf('Kpl and Perffusion Mean Kpl average: %f, Min: %f',mean(p(3,:)),min(p(3,:)))

%% Run 
% Run on Kpl when perfusion is known
load('tmpPerffusedPubPerf2')
[~,p(1,:),~,~] = vartest2(squeeze(fits(:,1,:,1)).',squeeze(fits(:,2,:,1)).');
[~,p(2,:),~,~] = vartest2(squeeze(fits(:,2,:,1)).',squeeze(fits(:,3,:,1)).');
[~,p(3,:),~,~] = vartest2(squeeze(fits(:,1,:,1)).',squeeze(fits(:,3,:,1)).');
subplot(2,2,4)
semilogy(N_vect,p,N_vect,zeros(size(N_vect))+significance,'g--')
ylim([1e-4,2]),xlim([5,60])
title('k_P_L varianece F-Test p value')
legend('CEA-bb Vs VEA-const','VEA-const Vs VEA-max','CEA-bb Vs VEA-max',...
    sprintf('p = %f',significance))
xlabel(' Number of Excitations '),ylabel('p value')
grid on
sprintf('Kpl and Perffusion Variance in Kpl average: %f, Min: %f',mean(p(3,:)),min(p(3,:)))