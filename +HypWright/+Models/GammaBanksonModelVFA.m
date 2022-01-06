classdef GammaBanksonModelVFA < HypWright.Models.GammaBanksonModel
    %BANKSONMODEL model for a two site excange system with perfused by a
    %vascular pool. b defines the vascular input function wich is assumed to be
    %uneffected.
    %   This is the basic two site exchange model. This model has a linear flip
    %   angle correction and will not work for non-linear sampling.
    %   The parameters for this model follow
    %   Kab - exchange reate from pool a to b: default 0
    %   Kba - exchange reate from pool b to a: default 0
    %   T1a - T1 decay constant for pool a: default 56
    %   T1b - T1 decay constant for pool b: default 31
    %   flipAngle - A vector of excitation angle in radians
    %   TR - A vector of repetition time
    %   kve - vascular extraction fraction: default 0.122
    %   ve - vascular volume fraction: defaul 0.91
    %   b - input function, some function of time that returns a change in
    %   pool a and b must return a 2 row vector: default [0;0]
    properties (Access = private)
    end
    methods (Static)
        function [Mxy,TRs,Mz] = evaluate(params,tSpan,Y0)
            % EVALUATE: solves this model over some time span (tSpan), with an
            % initial Y (Y0) and some parameters (params).
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1a - T1 decay constant for pool a: default 56
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - A vector of excitation angle in radians
            % TR - A vector of repetition time
            % kve - vascular extraction fraction: default 0.122
            % ve - vascular volume fraction: defaul 0.91
            % b - input function, some function of time that returns a change in
            % poola and b must return a 2 row vector: default [0;0]
            params = HypWright.Models.GammaBanksonModel.parseParams(params);
            flipAngles = params.flipAngle;
            TRs = params.TR;
            A = [-params.Kab-1/params.T1a-params.kve/params.ve,-params.Kba;...
                params.Kab,-params.Kba-1/params.T1b];
            b = params.b;
            M0 = Y0;
            fun = @(t,y)A*y+params.kve/params.ve*b(t);
            Mz = zeros(size(flipAngles));
            Mxy = zeros(size(flipAngles));
            Mz(:,1) = M0.*cos(flipAngles(:,1));
            Mxy(:,1) = (params.ve*M0+(1-params.ve)*params.b(TRs(1))).*sin(flipAngles(:,1));
            for i = 2:length(TRs)
                [~,Y] = ode45(fun,[TRs(i-1),TRs(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(flipAngles(:,i)).*(params.ve*Mz(:,i)+...
                    (1-params.ve)*b(TRs(i)));
                Mz(:,i) = cos(flipAngles(:,i)).*Mz(:,i);
            end
        end
        function [x,resultParams,resnorm,residual,exitflag,output,lambda,jacobian] = fitData(params,guess,...
                xdata,ydata,varargin)
            % FITDATA: fits some data set (xdata, ydata) with some constant
            % parameters (params) and variable parameters (guess) using the 
            % perfused two site exchange model. Which ever parameters are in the
            % guess struct will be fit, any parameters in the params struct will
            % be held constant, any parameters not specfied will be set to their
            % defaults and held constant. The fit will return one more argument
            % than the number of guesses. The last number is a scaling factor
            % applied to the fit data, to better match the magnitude og the
            % input
            % Params is a struct with the values
            % Kab - exchange reate from pool a to b: default 0
            % Kba - exchange reate from pool b to a: default 0
            % T1a - T1 decay constant for pool a: default 56
            % T1b - T1 decay constant for pool b: default 31
            % flipAngle - excitation angle in radians: default 0
            % TR - repetition time (again this is a linear TR model): default 1
            % kve - vascular extraction fraction: default 0.122
            % ve - vascular volume fraction: defaul 0.91
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
            params = HypWright.Models.GammaBanksonModel.parseParams(params);
            Y0 = ydata(:,1);
            fun = @(x,xdata)HypWright.Models.GammaBanksonModelVFA.fitFunction(...
                params,x,xNames,xdata,Y0);
            opts = params.fitOptions;
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(fun,x0,xdata,ydata,...
                [p.Results.lb],[p.Results.ub],opts);
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
            % kve - vascular extraction fraction: default 0.122
            % ve - vascular volume fraction: defaul 0.91
            % b - input function, some function of time that returns a change in
            % pool a and b must return a 2 row vector: default [0;0]
            % axis - axis handle to plot data
            p = inputParser();
            p.addOptional('axis',[])
            p.parse(varargin{:})
            Y0 = ydata(:,1);
            Y = HypWright.Models.GammaBanksonModel.evaluate(params,xdata,Y0);
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
            default = struct('Kab',0,'Kba',0,'T1a',56,'T1b',31,'flipAngle',0,...
                'TR',1,'kve',0.02,'ve',0.91,'t0',0,'gammaPdfA',2.8,...
                'gammaPdfB',4.5,'scaleFactor',1,'fitOptions', optimset('lsqcurvefit'));
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            paramsOut.b = @(t)params.scaleFactor*[...
                gampdf(t-params.t0,params.gammaPdfA,params.gammaPdfB);0];
        end
        function Y = fitFunction(params,x,xNames,tSpan,Y0)
            % fitFunction packs the parameter in params and x up and evaluates
            % using the evaluate funnction over some time (tSpan) with some
            % initial value (Y0)
            for i = 1:numel(xNames)
                params.(xNames{i}) = x(i);
            end
            % Uses the last value of x as a scaling factor.
            Y = HypWright.Models.GammaBanksonModelVFA.evaluate(params,tSpan,Y0);
        end
    end
end