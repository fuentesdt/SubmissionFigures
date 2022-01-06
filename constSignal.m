function [ flipAngles ] = constSignal(calMz,M0,N,ve,b2,TR,dim,scaleFact)
%CONSTSIGNAL Summary of this function goes here
%   Detailed explanation goes here
if(sum(M0) <= 1)
    tmp = 1;
else
    tmp = sum(M0);
end
if dim == 1
    altDim = 2;
else
    altDim = 1;
end
sRange = [0 tmp];  % bound product signal
% Iterative solve for maximum constant signal
flipAngles = zeros(N,size(M0,2));
Mz = zeros(N,size(M0,2));
MzTot = Mz;
while diff(sRange) > 1e-12
    tmpS = sum(sRange)/2;
    % Calculat Mz
    Mz(1,:) = M0; % init Mz
    tmp = b2(0*TR);
    MzTot(1,:) = (ve*Mz(1,:)+(1-ve)*tmp.');
    flipAngles(1,dim)= asin(tmpS/MzTot(1,dim));
    flipAngles(1,altDim) = flipAngles(1,dim);
    Mxy(1,:) = sin(real(flipAngles(1,:))).*MzTot(1,:);
    Mz(1,:) = cos(real(flipAngles(1,:))).*Mz(1,:);
    for k = 2:N
        tmp = b2((k-1)*TR);
        Mz(k,:) = calMz(k-1,Mz(k-1,:));
        MzTot(k,:) = (ve*Mz(k,:)+(1-ve)*tmp.');
        flipAngles(k,dim)= asin(tmpS/MzTot(k,dim));
        if ~isreal(flipAngles(k,dim))&&(k<N)
            flipAngles(k,dim) = 10*pi/180;
        end
        flipAngles(k,altDim) = flipAngles(k,dim);
        Mxy(k,:) = sin(real(flipAngles(k,:))).*MzTot(k,:);
        Mz(k,:) = cos(real(flipAngles(k,:))).*Mz(k,:);
    end
    % Check Conditions
    if isreal(flipAngles(N,dim)) && flipAngles(N,dim) < pi/2
        % tmpS is low set it to the new lower bound
        sRange(1) = tmpS;
    else
        % tmpS is too high, set it to the upper bound
        sRange(2) = tmpS;
    end
end
%% Change S
tmpS = scaleFact*tmpS;
Mz(1,:) = M0; % init Mz
tmp = b2(0*TR);
MzTot(1,:) = (ve*Mz(1,:)+(1-ve)*tmp.');
flipAngles(1,dim)= asin(tmpS/MzTot(1,dim));
flipAngles(1,altDim) = flipAngles(1,dim);
Mxy(1,:) = sin(real(flipAngles(1,:))).*MzTot(1,:);
Mz(1,:) = cos(real(flipAngles(1,:))).*Mz(1,:);
for k = 2:N
    tmp = b2(k*TR);
    Mz(k,:) = calMz(k,Mz(k-1,:));
    MzTot(k,:) = (ve*Mz(k,:)+(1-ve)*tmp.');
    flipAngles(k,dim)= asin(tmpS/MzTot(k,dim));
%     if ~isreal(flipAngles(k,dim))&&(k<N)
%         flipAngles(k,dim) = 10*pi/180;
%     end
    flipAngles(k,altDim) = flipAngles(k,dim);
    Mxy(k,:) = sin(real(flipAngles(k,:))).*MzTot(k,:);
    Mz(k,:) = cos(real(flipAngles(k,:))).*Mz(k,:);
end
flipAngles(N,:) = pi/2; % set last flip angle to 90
flipAngles = real(flipAngles);
end

