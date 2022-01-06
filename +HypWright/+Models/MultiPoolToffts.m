classdef MultiPoolToffts 
    %TWOPOOLTOFFTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function defaults()
            % DEFAULTS explains the default values for each parameter
            names = {'ExchangeTerms','T1s','FaList','TRList',...
                'PerfusionTerms','volumeFractions','VIF','fitOptions'};
            discriptions = {'A  NxN Matrix of Exchange Terms, where N is the number of chemical pools. The From pools should be along the rRows With the To pool along the Columns. Diagnal elemets will be set to zero'...
                ' A  Row vector of T1 decay times for each chemical pool.'...
                ' A  NxM of matrix of flip angles in radians, where N is the number of excitations and M is the number of chemical Pools'...
                ' A  NxM of Excitation Times in seconds, where N is the number of excitations and M is the number of chemical Pools'...
                ' A Row Vector of perfusion Exchange Constnats for each chemical pool.'...
                ' A Row Vector of volme fraction for each chemical pool. Only one value can be use if all pools have the same volume fraction.'...
                ' A function of a time variable (t) in seconds that returns a Row vector for the VIF of each chemical pool at the time t.'...
                ' Matlab FitOptions object'};
            defaultsVals = {'0','100','0','0','0','1','@(t)0','optimset(''lsqcurvefit'')'};
            fprintf('*Note* all terms must be a vector of size 1 x N where N is the number of chemical Pools\n')
            for i = 1:numel(names)
                fprintf('''%s'': %s\n Default Vaule: %s\n',...
                    names{i},discriptions{i},defaultsVals{i});
            end
        end
        function paramsOut = parseParams(paramsIn)
            % parseParams: a function to fill default param values if they are
            % not defined
            default = struct('ExchangeTerms',0,'T1s',100,'FaList',0,...
                'TRList',0,'PerfusionTerms',0,'volumeFractions',1,'VIF',@(t)0,...
                'fitOptions', optimset('lsqcurvefit'));
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
            kve = paramsIn.PerfusionTerms;
            if(length(paramsIn.volumeFractions)==1)
                ve = zeros(N,1)+paramsIn.volumeFractions;
            else
                ve = paramsIn.volumeFractions;
            end
            A = zeros(N);
            for i = 1:N
                for j = 1:N
                    if(i == j)
                        A(i,i) = -sum(K(i,:))-1/T1(i)-kve(i)/ve(i);
                    else
                        A(i,j) = K(j,i);
                    end
                end
            end
            paramsOut.A = A;
            paramsOut.b = paramsIn.VIF;
            paramsOut.kve = paramsIn.PerfusionTerms;
            paramsOut.ve = paramsIn.volumeFractions;
        end
        function [TRList,Mxy,Mz] = compile(M0,params)
            % EVALUATE: runs the model based on some input parameters
            params = HypWright.Models.MultiPoolToffts.parseParams(params);
            A = params.A;
            b = params.b;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = HypWright.Models.MultiPoolToffts.evaluate(...
                TRList,FaList,M0,A,b,params);
        end
        function [TRList, Mxy, Mz] = evaluate(TRList,FaList,M0,A,b,params)
            % EVALUATE: runs the model based on some input parameters
            kve = params.kve;
            ve = params.ve;
            fun = @(t,y)A*y+(kve/ve).'.*b(t);
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (params.ve*M0+(1-params.ve)*params.b(TRList(1))).*sin(FaList(:,1));
            for i = 2:length(TRList)
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(FaList(:,i)).*(params.ve.*Mz(:,i)+...
                    (1-params.ve).*b(TRList(i)));
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
        function DataCompare(A,params,xdata,ydata,varargin)
            p = inputParser();
            p.addOptional('Axes',[])
            p.addOptional('M0',[])
            p.parse(varargin{:})
            if isempty(p.Results.M0)
                M0 = ydata(:,1);
                M0 = M0./sin(params.FaList(:,1));
            else
                M0 = p.Results.M0;
            end
            [TRList,Mxy,~] = A.compile(M0,params);
            if isempty(p.Results.Axes)
                figure
            else
                axes(p.Results.Axes)
            end
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
end

