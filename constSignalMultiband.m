function [ flipAngles ] = constSignalMultiband(calMz,M0,N,ve,b2,TR,dim,FlipAngleOther)
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
while diff(sRange) > 1e-14
    tmpS = sum(sRange)/2;
    % Calculat Mz
    Mz(1,:) = M0; % init Mz
        tmp = b2(0*TR);
    flipAngles(1,dim) = asin(tmpS/(ve*Mz(1,dim)+(1-ve)*tmp(dim)));
    flipAngles(1,altDim) = FlipAngleOther(1);
    Mz(1,:) = Mz(1,:).*cos(real(flipAngles(1,:)));
    for k = 2:N
        Mz(k,:)= calMz(k-1,Mz(k-1,:));
        tmp = b2((k-1)*TR);
        flipAngles(k,dim)= asin(tmpS/(ve*Mz(k,dim)+(1-ve)*tmp(dim)));
        flipAngles(k,altDim) = FlipAngleOther(k);
%         if ~isreal(flipAngles(k,dim))&&(k<N)
%             flipAngles(k,dim) = 10*pi/180;
%         end
        Mz(k,:) = Mz(k,:).*cos(real(flipAngles(k,:))); % loss due to RF
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
% tmpS = sRange(1);
% Mz(1,:) = M0; % init Mz
% flipAngles(1,1) = asin(tmpS/Mz(1,1));
% flipAngles(1,2) = flipAngles(1,1);
% Mz(1,:) = Mz(1,:).*cos(real(flipAngles(1,:)));
% for k = 2:N
%     Mz(k,:)= calMz(k,Mz(k-1,:));
%     flipAngles(k,1)= asin(tmpS/Mz(k,1));
%     flipAngles(k,2) = flipAngles(k,1);
%     Mz(k,:) = Mz(k,:).*cos(real(flipAngles(k,:))); % loss due to RF
% end
flipAngles(N,:) = pi/2; % set last flip angle to 90
flipAngles = real(flipAngles);
end

