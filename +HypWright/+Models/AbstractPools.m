classdef AbstractPools 
    %ABSTRACTPOOLS Summary of this class goes here
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
                'TRList',0,'PerfusionTerms',0,'volumeFractions',1,...
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
%             % Hoping user will supply valid input!
            P = paramsOut.PerfusionTerms;
            C = paramsOut.ExchangeTerms;
            T1 = paramsIn.T1s;
            nC = size(P,1);
            nP = size(C,1);
            if numel(paramsOut.volumeFractions) == nP
                paramsOut.volumeFractions = repelem(...
                    paramsOut.volumeFractions,nC);
            end
            if numel(P) == 1;
                
            end
            if numel(T1) == nC
                T1 = repmat(T1,[nP 1]);
            end
            tmpA = zeros(nP,nP,nC,nC);
            A = zeros(nP*nC,nP*nC);
            R1 = diag(1./T1);
            for i = 1:nP
                for j = 1:nP
                    if i == j
                        tmpA(i,j,:,:) = squeeze(C(i,:,:))-...
                            diag(sum(squeeze(C(i,:,:)),1)+...
                            sum(squeeze(P(:,:,i)),2).');
                    else
                        tmpA(i,j,:,:) = diag(P(:,i,j));
                    end
                    A(i*nC-nC+1:i*nC,j*nC-nC+1:j*nC) = tmpA(i,j,:,:);
                end
            end
            A = A-R1;
            paramsOut.A = A;
        end
        function [TRList,Mxy,Mz] = compile(M0,params)
            % EVALUATE: runs the model based on some input parameters
            params = HypWright.Models.AbstractPools.parseParams(params);
            A = params.A;
            FaList = params.FaList;
            TRList = params.TRList;
            [TRList, Mxy, Mz] = HypWright.Models.AbstractPools.evaluate(...
                TRList,FaList,M0,A,params);
        end
        function [TRList, Mxy, Mz] = evaluate(TRList,FaList,M0,A,params)
            % EVALUATE: runs the model based on some input parameters
            fun = @(t,y)A*y;
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = params.volumeFractions.*M0.*sin(FaList(:,1));
            for i = 2:length(TRList)
                [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                Mz(:,i) = Y(end,:).';
                Mxy(:,i) = sin(FaList(:,i)).*params.volumeFractions.*Mz(:,i);
                Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
            end
        end
        function DataCompare(A,params,xdata,ydata)
            M0 = ydata(:,1);
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
end

