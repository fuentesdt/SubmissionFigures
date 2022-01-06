% ShowMxyPub Shows the flip angle schedules and the resulting magnetization
% evolutions for the paper
% Author: Chris Walker
% Date: 8/6/2018

clear all
close all
import HypWright.*
import HypWright.Models.*
clc
N_vect = round(linspace(5,60,15));
%% Variable Setup
N = 23;
TR = 2;
TR_list = (0:(N-1))*TR;
T1a = 43;
T1b = 33;
Kpl = 0.1;
alpha = 2.5;
beta = 4.5;
M0 = [0,0];
kve = 0.02;
ve = 0.95;
VIF_scale_fact = [1;0];
bb_flip_angle = 20;
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
FAType = {'BroadBand','ConstantLacMB','UCSFMethodConst'};
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
                flips(1,n) = 20*pi/180;
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
    [t_axis,Mxy,Mz] = model.compile(M0.',params);
    toc
    save_Mxy{i} = Mxy;
    save_Mz{i} = Mz;
    save_t_axis{i} = t_axis;
    save_flip_angles{i} = params.FaList;
end
save('tmpShowMxyPub')

