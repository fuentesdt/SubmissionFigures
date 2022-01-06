classdef TwoSiteExchange < HypWright.Models.Base
    %TWOSITEEXCHANGE model for a two site excange system with a time dependant
    % input
    %   This is the basic two site exchange model. This model has a linear flip
    %   angle correction and will not work for non-linear sampling.
    %   The parameters for this model follow
    %   Kab - exchange reate from pool a to b: default 0
    %   Kba - exchange reate from pool b to a: default 0
    %   T1a - T1 decay constant for pool a: default 56
    %   T1b - T1 decay constant for pool b: default 31
    %   flipAngle - excitation angle in radians: default 0
    %   TR - repetition time (again this is a linear TR model): default 1
    %   b - input function, some function of time that returns a change in pool
    %   a and b must return a 2 row vector: default [0;0]
    properties (Access = private)
    end
    methods (Static)
        function [Y,T,sol] = evaluate(params,tSpan,Y0)
            % EVALUATE: solves this model over some time span (tSpan), with an
            % initial Y (Y0) and some parameters (params).
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1a - T1 decay constant for pool a: default 56
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % b - input function, some function of time that returns a change in
            % pool a and b must return a 2 row vector: default [0;0]
            params = HypWright.Models.TwoSiteExchange.parseParams(params);
            A = [-(params.Kab+1/params.T1a+((1-cos(params.flipAngle))/params.TR)),...
                params.Kba; params.Kab,...
                -(params.Kba+1/params.T1b+((1-cos(params.flipAngle))/params.TR))];
            fun = @(t,y)A*y+params.b(t);
            sol = ode45(fun,tSpan,Y0);
            T = tSpan;
            if length(T) == 2, T = sol.x; end
            Y = deval(sol,T);
        end
        function [x,resultParams,resnorm,residual] = fitData(params,guess,...
                xdata,ydata,varargin)
            % FITDATA: fits some data set (xdata, ydata) with some constant
            % parameters (params) and variable parameters (guess) using the two
            % site exchange model. Which ever parameters are in the guess struct
            % will be fit, any parameters in the params struct will be held
            % constant, any parameters not specfied will be set to their
            % defaults and held constant. The last number is a scaling factor
            % applied to the fit data, to better match the magnitude of the
            % input function
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1a - T1 decay constant for pool a: default 56
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % b - input function, some function of time that returns a change in
            % pool a and b must return a 2 row vector: default [0;0]
            p = inputParser();
            p.addOptional('lb',[])
            p.addOptional('ub',[])
            p.parse(varargin{:})
            xNames = fieldnames(guess);
            x0 = zeros(numel(xNames),1);
            for i = 1:numel(xNames)
                params.(xNames{i}) = guess.(xNames{i});
                x0(i) = guess.(xNames{i});
            end
            params = HypWright.Models.TwoSiteExchange.parseParams(params);
            Y0 = ydata(:,1);
            fun = @(x,xdata)HypWright.Models.TwoSiteExchange.fitFunction(...
                params,x,xNames,xdata,Y0);
            opts = params.fitOptions;
            [x,resnorm,residual] = ...
                lsqcurvefit(fun,[x0;1],xdata,ydata,...
                [p.Results.lb,0],[p.Results.ub,1e20],opts);
            if x(end) > 1.1 || x(end) < 0.9
                fprintf(['Warning a significant scalling factor was applied'...
                    ' to the data to yield this fit\n']);
            end
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
            % T1a - T1 decay constant for pool a: default 56
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % b - input function, some function of time that returns a change in
            % pool a and b must return a 2 row vector: default [0;0]
            p = inputParser();
            p.addOptional('axis',[])
            p.parse(varargin{:})
            Y0 = ydata(:,1);
            Y = HypWright.Models.TwoSiteExchange.evaluate(params,xdata,Y0);
            resNorm = sum(sum((Y-ydata).^2));
            if(isempty(p.Results.axis))
                figure;
                curAxis = gca;
            else
                curAxis = p.Results.axis;
            end
            plot(curAxis,xdata,Y(1,:)','g',xdata,Y(2,:),'b',...
                xdata,ydata(1,:),'go',xdata,ydata(2,:),'bo')
            xlabel('Time')
            ylabel('Signal Intensity')
            legend('Model Pool A','Model Pool B')
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
            if isfield(paramsIn,'T1a')
                paramsOut.T1a = paramsIn.T1a;
            else
                paramsOut.T1a = 56;
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
            if isfield(paramsIn,'b')
                paramsOut.b = paramsIn.b;
            else
                paramsOut.b = @(t)0.*[t;t];
            end
            if isfield(paramsIn,'fitOptions')
                paramsOut.fitOptions = paramsIn.fitOptions;
            else
                paramsOut.fitOptions = optimset('lsqcurvefit');
            end
        end
        function Y = fitFunction(params,x,xNames,tSpan,Y0)
            % fitFunction packs the parameter in params and x up and evaluates
            % using the evaluate funnction over some time (tSpan) with some
            % initial value (Y0)
            for i = 1:numel(xNames)
                params.(xNames{i}) = x(i);
            end
            % Uses the last value of x as a scaling factor.
            Y = x(end)*HypWright.Models.TwoSiteExchange.evaluate(params,tSpan,Y0);
        end
    end
end