classdef TwoSiteEmpirical < HypWright.Models.Base
    %TWOSITEEMPIRICAL model for a two site excange system some defined pool a
    % signal
    %   This is the basic two site exchange model. This model has a linear flip
    %   angle correction and will not work for non-linear sampling.
    %   The parameters for this model follow
    %   Kab - exchange reate from pool a to b: default 0
    %   Kba - exchange reate from pool b to a: default 0
    %   T1b - T1 decay constant for pool b: default 31
    %   flipAngle - excitation angle in radians: default 0
    %   TR - repetition time (again this is a linear TR model): default 1
    %   a - the measured pool a signal: default 0
    %   aTimes - the time points for a: default 0
    properties (Access = private)
    end
    methods (Static)
        function [Y,T,fun] = evaluate(params,tSpan,Y0)
            % EVALUATE: solves this model over some time span (tSpan), with an
            % initial Y (Y0) and some parameters (params).
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % a - the measured pool a signal: default 0
            % aTimes - the time points for a: default 0
            params = HypWright.Models.TwoSiteEmpirical.parseParams(params);
            Kab = params.Kab;
            A = -1/params.T1b-(1-cos(params.flipAngle))/params.TR-params.Kba;
            a = params.a;
            aTime = params.aTime;
            fun = @(t,y)exp(A*t)*Y0+Kab.*trapz(...
                linspace(0,t,1000),exp(-A.*((linspace(0,t,1000))-t)).*...
                interp1(aTime,a(1,:),linspace(0,t,1000),'linear',0));
            T = tSpan;
            Y = zeros(size(T));
            for i = 1:length(T)
                Y(i) = fun(T(i));
            end
        end
        function [x,resultParams,resnorm,residual] = fitData(params,guess,...
                xdata,ydata)
            % FITDATA: fits some data set (xdata, ydata) with some constant
            % parameters (params) and variable parameters (guess) using the two
            % site exchange model. Which ever parameters are in the guess struct
            % will be fit, any parameters in the params struct will be held
            % constant, any parameters not specfied will be set to their
            % defaults and held constant
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % a - the measured pool a signal: default 0
            % aTimes - the time points for a: default 0
            xNames = fieldnames(guess);
            x0 = zeros(numel(xNames),1);
            for i = 1:numel(xNames)
                params.(xNames{i}) = guess.(xNames{i});
                x0(i) = guess.(xNames{i});
            end
            params = HypWright.Models.TwoSiteEmpirical.parseParams(params);
            Y0 = ydata(:,1);
            fun = @(x,xdata)HypWright.Models.TwoSiteEmpirical.fitFunction(...
                params,x,xNames,xdata,Y0);
            opts = optimset('lsqcurvefit');
            opts = optimset(opts,'display','off','tolFun',1e-15);
            [x,resnorm,residual] = ...
                lsqcurvefit(fun,x0,xdata,ydata,[],[],opts);
            resultParams = params;
            for i = 1:numel(xNames)
                resultParams.(xNames{i}) = x(i);
            end
        end
        function dataCompare(params,xdata,ydata,varargin)
            % DATACOMPARE: displays the model with the parameters parameters
            % (params) against the data (xdata, ydata). optiona 4th argument for
            % a figure axis in which to draw the plot
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % a - the measured pool a signal: default 0
            % aTimes - the time points for a: default 0
            p = inputParser();
            p.addOptional('axis',[])
            p.parse(varargin{:})
            Y0 = ydata(:,1);
            Y = HypWright.Models.TwoSiteEmpirical.evaluate(params,xdata,Y0);
            resNorm = sum(sum((Y-ydata).^2));
            if(isempty(p.Results.axis))
                figure;
                curAxis = gca;
            else
                curAxis = p.Results.axis;
            end
            plot(curAxis,xdata,Y','b',xdata,ydata,'bo',params.aTime,params.a,'g')
            xlabel('Time')
            ylabel('Signal Intensity')
            legend('Model Pool b','Data Pool b')
            title('Comparison of data with two site exchage model')
            fprintf('The norm of the residual is: %d\n',resNorm)
        end
    end
    methods (Static, Access = protected)
        function paramsOut = parseParams(paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            if isfield(paramsIn,'Kab')
                paramsOut.Kab = paramsIn.Kab;
            else
                paramsOut.Kab = 0;
            end
            if isfield(paramsIn,'Kba')
                paramsOut.Kba = paramsIn.Kba;
            else
                paramsOut.Kba = 0;
            end
            if isfield(paramsIn,'T1b')
                paramsOut.T1b = paramsIn.T1b;
            else
                paramsOut.T1b = 31;
            end
            if isfield(paramsIn,'flipAngle')
                paramsOut.flipAngle = paramsIn.flipAngle;
            else
                paramsOut.flipAngle = 0;
            end
            if isfield(paramsIn,'TR')
                paramsOut.TR = paramsIn.TR;
            else
                paramsOut.TR = 1;
            end
            if isfield(paramsIn,'a')
                paramsOut.a = paramsIn.a;
            else
                paramsOut.a = 0;
            end
            if isfield(paramsIn,'aTime')
                paramsOut.aTime = paramsIn.aTime;
            else
                paramsOut.a = 0:length(paramsOut.a)-1;
            end
        end
        function Y = fitFunction(params,x,xNames,tSpan,Y0)
            % fitFunction packs the parameter in params and x up and evaluates
            % using the evaluate funnction over some time (tSpan) with some
            % initial value (Y0)
            for i = 1:numel(xNames)
                params.(xNames{i}) = x(i);
            end
            Y = HypWright.Models.TwoSiteEmpirical.evaluate(params,tSpan,Y0);
        end
    end
end