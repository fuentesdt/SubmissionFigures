function val = fwhm(x,y,varargin)
%% Computes the FWHM of the largest peak in the passed in data

p = inputParser();
p.addOptional('verbose',[],@isstr)
p.parse(varargin{:})
if (strcmp(p.Results.verbose,'verbose')), 
    verbose = 1;
else
    verbose = 0; 
end
%Find peak max
[maxHeight, maxHeightI] = max(y);
% find half height value
halfHeight = maxHeight/2;
% find left half height
[~,tmpI] = find(halfHeight>y(1:maxHeightI),1,'last');
% return error if not found
if(isempty(tmpI)||tmpI==maxHeightI)
%     fprintf(['Half Height not found to left of peak! HalfHeight %d'...
%         ' minimum value left of peak %d\n'],halfHeight,min(y(1:maxHeightI)));
    val = 0;
    return
end
% interpolate to exact point of left half height
leftHalfHeightI = interp1(y(tmpI:tmpI+1),x(tmpI:tmpI+1),halfHeight);
% store left side of integration waveform
integrationY = [halfHeight,y(tmpI+1:maxHeightI)];
integrationX = [leftHalfHeightI,x(tmpI+1:maxHeightI)];
% find right half height
[~,tmpI] = find(halfHeight>y(maxHeightI:end),1,'first');
% return error if not found
if(isempty(tmpI)||tmpI==0)
%     fprintf(['Half Height not found to right of peak! HalfHeight %d'...
%         ' minimum value right of peak %d\n'],halfHeight,min(y(maxHeightI:end)));
    val = 0;
    return
end
% adjust search results
tmpI = tmpI+maxHeightI-1;
% interpolate to exact point of right half height
rightHalfHeightI = interp1(y(tmpI-1:tmpI),x(tmpI-1:tmpI),halfHeight);
% store right side of integration waveform
integrationY = [integrationY,y(maxHeightI+1:tmpI-1),halfHeight,];
integrationX = [integrationX,x(maxHeightI+1:tmpI-1),rightHalfHeightI,];
% Assemble half height waveform
val = trapz(integrationX,integrationY);
if(verbose)
    figure
    plot(x,y,'k',integrationX,integrationY,'ro')
    title('FWHM results')
    xlabel('x'),ylabel('y'),legend('input','FWHM')
    text(rightHalfHeightI*1.02,halfHeight,num2str(val));
end
