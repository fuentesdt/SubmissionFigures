classdef MultiPool 
    %TWOPOOL A simple chemical exchange model assuming no inpu functions
    %   parameters Values
    %   * ExchangeTerms - A Matrix defining chemical Exchange. Defalt: 0
    %   * T1s - A row vector of T1 decay terms. Default: 100
    %   * FaList - A matrix of excitation angles. Default: 0
    %   * TRList - A matrix of excitation times. Default: 0
    %   * fitOptions - A matlab fit option structure. Default: optimset(''lsqcurvefit'')
    %   There is NO imput validation for the parameters passed in, for more
    %   detail on the assumed data structur of these parameters use the
    %   defaults function
    
    properties
    end
    
    methods (Static)
        function defaults()
            % DEFAULTS explains the default values for each parameter
            names = {'ExchangeTerms','T1s','FaList','TRList','fitOptions'};
            discriptions = {'A  NxN Matrix of Exchange Terms, where N is the number of chemical pools. The From pools should be along the rRows With the To pool along the Columns. Diagnal elemets will be set to zero'...
                ' A  Row vector of T1 decay times for each chemical pool.'...
                ' A  NxM of matrix of flip angles in radians, where N is the number of excitations and M is the number of chemical Pools'...
                ' A  NxM of Excitation Times in seconds, where N is the number of excitations and M is the number of chemical Pools'...
                ' Matlab FitOptions object'};
            defaultsVals = {'0','100','0','0','optimset(''lsqcurvefit'')'};
            for i = 1:numel(names)
                fprintf('''%s'': %s\n Default Vaule: %s\n',...
                    names{i},discriptions{i},defaultsVals{i});
            end
        end
        function paramsOut = parseParams(paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('ExchangeTerms',0,'T1s',100,'FaList',0,...
                'TRList',0,'fitOptions', optimset('lsqcurvefit'));
            tmpNames = fieldnames(default);
            paramsOut = paramsIn;
            for i = 1:numel(tmpNames)
                if ~isfield(paramsOut,tmpNames{i})
                    paramsOut.(tmpNames{i}) = default.(tmpNames{i});
                end
            end
            % Fill all flip angles with a value if only one flip angle is passed in
            if size(paramsOut.FaList,2)==1
                paramsOut.FaList = repmat(paramsOut.FaList(:,1),...
                    1,length(paramsOut.TRList));
            end
%             % Assuming input validation will put too much computational burdon on the fitting
%             % Hopeing user will supply valid input!
%             % Validte input
%             if (size(params.ExchangeTerms,1)~=size(params.ExchangeTerms,1))
%                 error('Exchange matrix not square.')
%             end
            N = size(paramsIn.ExchangeTerms,1);
            K = triu(paramsIn.ExchangeTerms)+tril(paramsIn.ExchangeTerms);
            T1 = paramsIn.T1s;
            A = zeros(N);
            for i = 1:N
                for j = 1:N
                    if(i == j)
                        A(i,i) = -sum(K(i,:))-1/T1(i);
                    else
                        A(i,j) = K(j,i);
                    end
                end
            end
            paramsOut.A = A;
        end
        function [TRList,Mxy,Mz] = compile(M0,params)
            % EVALUATE: runs the model based on some input parameters
            params = HypWright.Models.MultiPool.parseParams(params);
            A = params.A;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = HypWright.Models.MultiPool.evaluate(...
                TRList,FaList,M0,A);
        end
        function [x,resultParams,allParams,resnorm,residual,exitflag,output,lambda,jacobian]...
                = fitData(params,guess,xdata,ydata,varargin)
            p = inputParser();
            p.addOptional('lb',[])
            p.addOptional('ub',[])
            p.parse(varargin{:})
            xNames = fieldnames(guess);
            j = 1;
            xIndex = cell(size(xNames));
            for i = 1:numel(xNames)
                iFits = ~isnan(guess.(xNames{i}));
                xIndex{i} = find(iFits==1);
                for k = 1:numel(xIndex{i})
                    x0(j) = guess.(xNames{i})(xIndex{i}(k));
                    j = j+1;
                end
            end
            params = HypWright.Models.MultiPool.parseParams(params);
            Y0 = ydata(:,1)./sin(params.FaList(:,1));
            fun = @(x,xdata)HypWright.Models.MultiPool.fitFunction(...
                params,x,xNames,xIndex,Y0);
            opts = params.fitOptions;
            [x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(fun,x0,xdata,ydata,...
                [p.Results.lb],[p.Results.ub],opts);
            resultParams = guess;
            allParams = params;
            j = 1;
            for i = 1:numel(xNames)
                for k = 1:numel(xIndex{i})
                    resultParams.(xNames{i})(xIndex{i}(k)) = x(j);
                    allParams.(xNames{i})(xIndex{i}(k)) = x(j);
                    j = j+1;
                end
            end
        end
        function [TRList, Mxy, Mz] = evaluate(TRList,FaList,M0,A)
            % EVALUATE: runs the model based on some input parameters
            fun = @(t,y)A*y;
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (M0.*sin(FaList(:,1)));
            for i = 2:length(TRList)
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(FaList(:,i)).*Mz(:,i);
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
        function DataCompare(A,params,M0,xdata,ydata)
            M0 = M0./sin(params.FaList(:,1));
            [TRList,Mxy,~] = A.compile(M0,params);
            figure
            for i = 1:size(Mxy,1)
                tmpLine = plot(TRList,Mxy(i,:));
                hold on
                plot(xdata,ydata(i,:),'o','MarkerEdgeColor',tmpLine.Color);
            end
            hold off
            xlabel('Time (sec)')
            ylabel('Signal (arb)')
        end
    end
    
    methods (Access = private, Static)
    function Y = fitFunction(params,x,xNames,xIndex,Y0)
            % fitFunction packs the parameter in params and x up and evaluates
            % using the evaluate funnction over some time (tSpan) with some
            % initial value (Y0)
            j = 1;
            for i = 1:numel(xNames)
                for k = 1:numel(xIndex{i})
                    params.(xNames{i})(xIndex{i}(k)) = x(j);
                    % Check if fitting flip angle (there mus be a better
                    % way to do this
                    if strcmp(xNames{i}, 'FaList')
                        params.(xNames{i}) =...
                            repmat(x(j),size(params.(xNames{i})));
                    end
                    j = j+1;
                end
            end
            [~, Y, ~] = HypWright.Models.MultiPool.compile(Y0,params);
    end
    end
end

