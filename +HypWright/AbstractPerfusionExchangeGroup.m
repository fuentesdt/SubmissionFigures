classdef AbstractPerfusionExchangeGroup < HypWright.SpinGroup
    % TwoSiteExchangeGroup a class that represents a set of isolated spins
    %   Detailed explanation goes here
    %   Properties
    %   M - magnetization vector
    %   M0 - Equilibrium magnetization
    %   T1 - T1 decay constant
    %   T2 - T2 Decay constant
    %   gamma - gyromagnetic ratio
    %   density - relative number of spins in this group
    %   Methods
    %   IsolatedSpinGrp(M,M0,T1,T2,gamma,density) - initializes the spin
    %   group with an initial magnetization (M), equilibrium magentization
    %   (M0), T1 decay (T2), T2 Decay (T2), gyromagnetic ratio (gamma) and
    %   spin density (density)
    %   calculationFrame() -  returns the frequencty of the rotating
    %   refrence frame that dm is calculated in
    %   dM(x,y,z,t,M) - returns a dm at some position (x,y,z), some time
    %   (t), and some initial magnetization M ([Mx;My;Mz])
    
    properties (SetAccess = private)
        Minit % more user friendly form of initial Magnetization
        M % Magnetization
        M0 % equilibrium Magnetization
        T1 % vector of T1 decay constants
        T2 % vecore of T2 decay constants
        ppm % vector of chemical shifts
        gamma % gyromagnetic ratio
        density % Spin density
        kExchange % exchange constants for chemical species
        kPerfusion % Exchange terms for perfusion of chemical species
        b % vasuclular input function
        nChemicalSpecies % number of chemical species
        meanPPM % avarge ppm, this is used as the calsulation frame
    end
    methods
        function self = AbstractPerfusionExchangeGroup(varargin)
            % Constructor - initializes the spin group
            % IsolatedSpinGrp(M,M0,T1,T2,gamma,density) - initializes the spin
            % group with an initial magnetization (M), equilibrium magentization
            % (M0), T1 decay (T2), T2 Decay (T2), gyromagnetic ratio (gamma) and
            % spin density (density)
            p = inputParser();
            p.addRequired('Minit',[],@isnumeric)
            p.addRequired('M0',[],@isnumeric)
            p.addRequired('T1',[],@isnumeric)
            p.addRequired('T2',[],@isnumeric)
            p.addRequired('ppm',[],@isnumeric)
            p.addRequired('kExchange',[],@isnumeric)
            p.addRequired('kPerfusion',[],@isnumeric)
            %% TODO refine inout checking
            p.addRequired('gamma',67.262e6,@isscalar)
            p.addRequired('density',1,@isscalar)
            p.addOptional('b',[])
            p.parse(varargin{:})
            self.Minit = p.Results.Minit;
            self.M0 = p.Results.M0;
            self.T1 = p.Results.T1;
            self.T2 = p.Results.T2;
            self.ppm = p.Results.ppm;
            self.kExchange = p.Results.kExchange;
            self.kPerfusion = p.Results.kPerfusion;
            self.density = p.Results.density;
            self.kExchangeTotals = sum(self.kExchange,2);
            if isempty(p.Results.b)
                self.b = @(t)zeros(size(self.Minit,1)*3);
            else
                self.b = p.Results.b;
            end
            for i = 1:self.nChemicalSpecies
                indexI = 1+(i-1)*3;
            Beff(indexI:indexI+2) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
                0,0,1]*B-[0;0;1].*B0*(1+self.meanPPM-self.ppm(i));
            end
            self.nChemicalSpecies = size(self.Minit,1);
            self.meanPPM = mean(self.ppm);
            self.checkInputs();
        end
        function A = getA(self,x,y,z,t,PS,B0,varargin)
            if ~isempty(varargin) == 1
                B = varargin{1};
            else
                B = repmat(B0,1,length(t))+PS.B(x,y,z,t);
            end
            theta =  -self.calculationFrame(B0)*t;
            Beff = zeros(self.nChemicalSpecies*3,1);
            for i = 1:self.nChemicalSpecies
                indexI = 1+(i-1)*3;
            Beff(indexI:indexI+2) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
                0,0,1]*B-[0;0;1].*B0*(1+self.meanPPM-self.ppm(i));
            end
            A = zeros(self.nChemicalSpecies*3,self.nChemicalSpecies*3);
            for i = 1:self.nChemicalSpecies
                indexI = 1+(i-1)*3;
                for j = 1:self.nChemicalSpecies
                    indexJ = 1+(j-1)*3;
                    % Set terms along Diagnal 
                    if i == j
                        A(indexI,indexJ) = 0-1/self.T2(i)...
                            -self.kExchangeTotals(i)-self.kPerfusion(i);
                        A(indexI,indexJ+1) = -Beff(indexI+2);
                        A(indexI,indexJ+2) = Beff(indexI+1);
                        A(indexI+1,indexJ) = Beff(indexI+2);
                        A(indexI+1,indexJ+1) = 0-1/self.T2(i)...
                            -self.kExchangeTotals(i)-self.kPerfusion(i);
                        A(indexI+1,indexJ+2) = -Beff(indexI);
                        A(indexI+2,indexJ) = -Beff(indexI+1);
                        A(indexI+2,indexJ+1) = Beff(indexI);
                        A(indexI+2,indexJ+2) = 0-1/self.T1(i)...
                            -self.kExchangeTotals(i)-self.kPerfusion(i);
                    % Set Cross terms
                    else
                        A(indexI,indexJ) = 0+self.kPerfusion(i,j);
                        A(indexI,indexJ+1) = 0;
                        A(indexI,indexJ+2) = 0;
                        A(indexI+1,indexJ) = 0;
                        A(indexI+1,indexJ+1) = 0+self.kPerfusion(i,j);
                        A(indexI+1,indexJ+2) = 0;
                        A(indexI+2,indexJ) = 0;
                        A(indexI+2,indexJ+1) = 0;
                        A(indexI+2,indexJ+2) = 0+self.kPerfusion(i,j);
                    end  
                end
            end
        end
        function dm = dM(self,x,y,z,t,M,PS,B0)
            % DM: returns the delta m at some time and location and given M
            % dM(x,y,z,t,M) - returns a dm at some position (x,y,z), some time
            % (t), and some initial magnetization M ([Mx;My;Mz])
            A = self.getA(x,y,z,t,PS,B0);
            for i = 1:self.nChemicalSpecies
                indexI = 1+(i-1)*3;
            Beff(indexI:indexI+2) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
                0,0,1]*B-[0;0;1].*B0*(1+self.meanPPM-self.ppm(i));
            end
            Recovery(1:3) = 1/self.T1a*self.M0(1:3);
            Recovery(4:6) = 1/self.T1b*self.M0(4:6);
            dm = A*M+Recovery.'+self.kve*self.b(t);
        end
        function vals = analytical(self,x,y,z,t0,M,t,PS,B0,varargin)
            % ANALYTICAL: retun a function handle to the analytical soluton
            A = self.getA(x,y,z,t0+1e-9,PS,B0,varargin{:});
            if length(t) == 1
                vals = expm(A*(t-t0))*(M+integral(@(t2)expm(-A*(t2-t0))*self.kve*self.b(t2),...
                    t0,t,'ArrayValued',true));
            else
                vals = zeros(length(M),length(t));
                tmpT = t0-0.001:0.001:t(end)+0.001;
                tmpVals = cell2mat(arrayfun(@(t2)expm(-A*(t2-t0))*...
                    self.kve*self.b(t2),tmpT,'UniformOutput',false));
                tmpIntelgral = cumtrapz(tmpT,tmpVals,2);
                ForceFunIntegral = interp1(tmpT,tmpIntelgral.',t).';
                for i = 1:length(t)
                    vals(:,i) = expm(A*(t(i)-t0))*(M+ForceFunIntegral(:,i));
                end
            end
            if (any(isnan(vals(:))) || any(isinf(vals(:))))
                odefun = @(M,t)self.dM(x,y,z,M,t,PS,B0);
                tmpSol = ode45(odefun,[t0,t(end)],M);
                vals = deval(tmpSol,t);
            end
        end
        function ret = useAnalytical(self)
            %USEANALYTICAL: determins if the Analytical Soultion shouldbe used
            %for the given spin group under the given conditions
            ret = all(self.M0 == 0);
        end
    end
    methods (Access = private)
        function checkInputs(self)
            %% Checks the variables in the object to make sure they are valid
            if size(self.M0(2)) ~= 3 % make sure each chemical species has 3 components
                error('M needs to be a n b 3 matrix where n is umber of chemical species /n %d',...
                    self.M0)
            end
            if size(self.M0) ~= size(self.Minit) % make sure M vectors are the right size
                error('M and M0 vectors are not the same size/n %d /n %d',...
                    self.Minit,self.M0)
            end
            if size(self.T1) ~= self.nChemicalSpecies % make sure T1 vecor matchs M
                error('M and T1 vectors are not the same size/n %d /n %d',...
                    self.Minit,self.T1)
            end
            if size(self.T2) ~= self.nChemicalSpecies % make sure T2 vecor matchs M
                error('M and T2 vectors are not the same size/n %d /n %d',...
                    self.Minit,self.T2)
            end
            if size(self.ppm) ~= self.nChemicalSpecies % make sure ppm vecor matchs M
                error('M and ppm vectors are not the same size/n %d /n %d',...
                    self.Minit,self.ppm)
            end
            % make sure kExchange is square matrix of correct size
            if ~isequal(size(self.kExchange),[self.nChemicalSpecies,self.nChemicalSpecies]) 
                error('kExchange must be a nxn matrix where n is the number of chemical species/n %d /n %d',...
                    size(self.Minit),self.kExchange)
            end
            if size(self.kPerfusion) ~= self.nChemicalSpecies % make sure kPerfusion vecor matchs M
                error('M and kPerfusion vectors are not the same size/n %d /n %d',...
                    self.Minit,self.kPerfusion)
            end
            try
                tmpVIF = self.b(0);
            catch
                error('b function not taking a scaler input');
            end
            %% TODO check if transposed
            if size(tmpVIF) ~= self.nChemicalSpecies*3
                error('b must return a vector 3 times longer than number of chemical species/n %d /n %d',...
                    size(self.Minit),tmpVIF)
            end
        end
    end
    
end

