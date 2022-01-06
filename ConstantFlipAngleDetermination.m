% SNRLimit determins the SNR threshold that the improved lactate SNR from
% VEA-max results in more accurate Kpl
% Author: Chris Walker
% Date: 8/3/2018

clear all
close all
import HypWright.*
import HypWright.Models.*
clc

M0 = [0,0];
kve = 0.02;
ve = 0.95;
VIF_scale_fact = [1;0];
excitation_angle_vect = [5,10,15,20,25,30,35,40]; %Vector to store all the noise factors tested, noise increase with noise factor
saveID = 'CEAOptTmp';
verbose = false;
n_averages = 30;

for m = 1:length(excitation_angle_vect)
%% Variable Setup
N = 40;
TR = 2;
TR_list = (0:(N-1))*TR;
T1a = 43;
T1b = 33;
Kpl = 0.1;
alpha = 2.5;
beta = 4.5;
noise_fact = 5e-4;
bb_flip_angle = excitation_angle_vect(m);
opts = optimset('lsqcurvefit');
opts.TolFun = 1e-09;
opts.TolX = 1e-09;
opts.Display = 'off';
params = struct('t0',[0;0],'gammaPdfA',[alpha;1],'gammaPdfB',[beta;1],...
    'scaleFactor',VIF_scale_fact,'T1s',[T1a,T1b],'ExchangeTerms',[0,Kpl;0,0],...
    'TRList',TR_list,'PerfusionTerms',[kve,0],'volumeFractions',ve,...
    'fitOptions', opts);
model = MultiPoolTofftsGammaVIF();
%% Get true Mz
%% Choose Excitation Angle
FAType = {'BroadBand'};
clear flips
for i = 1:numel(FAType)
    switch (FAType{i})
        case('BroadBand') % 20 degree lactate 20 degree pyruvate
            flipAngles = zeros(2,length(TR_list))+bb_flip_angle*pi/180;
            params.FaList = flipAngles;
        case('Lac90') % 90 degree lactate 10 degree pyruvate
            flipAngles = zeros(2,length(TR_list))+[10;90]*pi/180;
            params.FaList = flipAngles;
        case('ConstantLacMB') % Constant Signal multi-band
            tic
            A = [-Kpl-1/T1a-kve/ve,Kpl;0,-1/T1b];
            b = @(t)VIF_scale_fact.*[gampdf(t,alpha,beta);0];
            calMz = @(i,Mz)Mz*expm(A*TR)+kve/ve*integral(...
                @(t2)expm(A.'.*(i*TR-t2))*(b(t2)),...
                ((i-1)*TR),i*TR,'ArrayValued',true).';
            tmp =  constSignal(calMz,M0,N,ve,b,TR,1,1);
            flipAngles = tmp;
            tmp = constSignalMultiband(calMz,M0,N,ve,b,TR,2,flipAngles(:,2));
            flipAngles(:,2) = tmp(:,2);
            flipAngles = flipAngles.';
            params.FaList = flipAngles;
        case('ConstantTotalBB') % Broad band pulse VFA constant Signal
            tic
            A = [-Kpl-1/T1a-kve/ve,Kpl;0,-1/T1b];
            b = @(t)VIF_scale_fact.*[gampdf(t,alpha,beta);0];
            calMz = @(i,Mz)Mz*expm(A*TR)+kve/ve*integral(...
                @(t2)expm(A.'.*(i*TR-t2))*(b(t2)),...
                ((i-1)*TR),i*TR,'ArrayValued',true).';
            [flipAngles] = constTotalSignal(calMz,M0,N,ve,b,TR);
            flipAngles = flipAngles.';
            params.FaList = flipAngles;
        case('UCSFMethod') % Nagashima for lactate const VFA for pyuvate
            tic
            E1(1) = exp(-TR*(1/T1a+Kpl));
            E1(2) = exp(-TR/T1b);
            for n = 1:N
                flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
            end
            flips(1,:) = vfa_const_amp(N, pi/2, E1(1));
            params.FaList = flips;
        case('UCSFMethodPerfused') % Nagashima for lactate const VFA for pyuvate but accout for perfusion
            tic
            E1(1) = exp(-TR*(1/T1a+Kpl));
            E1(2) = exp(-TR/T1b);
            for n = 1:N
                flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
            end
            A = [-Kpl-1/T1a-kve/ve,Kpl;0,-1/T1b];
            b = @(t)VIF_scale_fact.*[gampdf(t,alpha,beta);0];
            calMz = @(i,Mz)Mz*expm(A*TR)+kve/ve*integral(...
                @(t2)expm(A.'.*(i*TR-t2))*(b(t2)),...
                ((i-1)*TR),i*TR,'ArrayValued',true).';
            tmp =  constSignal(calMz,M0,N,ve,b,TR,1,1);
            flips(1,:) = tmp(:,1).';
            params.FaList = flips;
        case('UCSFMethodConst') % Nagashima for lactate const 10 pyruvate
            tic
            E1(1) = exp(-TR*(1/T1a+Kpl));
            E1(2) = exp(-TR/T1b);
            for n = 1:N
                flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
                flips(1,n) = 10*pi/180;
            end
            params.FaList = flips;
        case('Nagashima') % Nagashima for Both Pyruvate and Lactate
            tic
            E1(1) = exp(-TR*(1/T1a));
            E1(2) = exp(-TR/T1b);
            for n = 1:N
                flips(1,n) = acos(sqrt((E1(1)^2-E1(1)^(2*(N-n+1)))/(1-E1(1)^(2*(N-n+1)))));
                flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
            end
            params.FaList = flips;
        case('T1Effective') % T1 Effective for Both
            tic
            E1(1) = exp(-TR*(1/T1a+Kpl));
            E1(2) = exp(-TR*(1/T1b-Kpl));
            for n = 1:N
                flips(1,n) = acos(sqrt((E1(1)^2-E1(1)^(2*(N-n+1)))/(1-E1(1)^(2*(N-n+1)))));
                flips(2,n) = acos(sqrt((E1(2)^2-E1(2)^(2*(N-n+1)))/(1-E1(2)^(2*(N-n+1)))));
            end
            params.FaList = flips;
    end
    tic
    %% Fitting
    [~,Mxy,~] = model.compile(M0.',params);
    for n = 1:n_averages
        noise = noise_fact.*(rand(size(Mxy))-0.5);
        tmpMxy = abs(Mxy+noise);
        guessParams = struct('ExchangeTerms',[NaN,0.0;NaN,NaN]);
        [x,~,allParams,resnorm,residual,exitflag,~,~,~]=model.fitData(...
            params,guessParams,TR_list,tmpMxy,'lb',0,'ub',10,'Y0',M0.');
        fits(m,i,n) = x;
        signals(m,i,n) = sum(sum(Mxy));
        lacSig(m,i,n) = sum(Mxy(2,:));
        pyrSig(m,i,n) = sum(Mxy(1,:));
    end
    if verbose
%         figure
%         plot(flipAngles.'*180/pi)
%         legend('lactate','pryvate')
%         xlabel('Time (sec)'),ylabel('Excitation Angle')
        model.DataCompare(model,allParams,TR_list,tmpMxy)
        title(sprintf('%s Nois:%f Kpl %1.2f',FAType{i},noise_fact,x(1)))
        drawnow
    end
        toc       
end
size(fits)
save(saveID)
end
figure
subplot(2,3,1)
plot(excitation_angle_vect,abs(squeeze(mean(fits,3))-0.1)/0.1)
xlabel('Flip Angle'),ylabel('k_p_l Error'),title('k_p_l Error')
grid on
subplot(2,3,2)
plot(excitation_angle_vect,abs(squeeze(std(fits,[],3)))/0.1)
xlabel('Flip Angle'),ylabel('k_p_l Coeff. Varr.'),title('k_p_l Coeff. Varr.')
grid on
subplot(2,3,4)
plot(excitation_angle_vect,squeeze(mean(signals,3)))
xlabel('Flip Angle'),ylabel('Total Signal'),title('Total Signal')
grid on
subplot(2,3,5)
plot(excitation_angle_vect,squeeze(mean(pyrSig,3)))
xlabel('Flip Angle'),ylabel('Pyruvate Signal'),title('Pyruvate Signal')
grid on
subplot(2,3,6)
plot(excitation_angle_vect,squeeze(mean(lacSig,3)))
xlabel('Flip Angle'),ylabel('Lactate Signal'),title('Lactate Signal')
grid on
