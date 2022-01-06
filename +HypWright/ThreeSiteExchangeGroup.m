classdef ThreeSiteExchangeGroup < HypWright.SpinGroup
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
    
    properties
        M % magnetization vector
        M0 % Equilibrium magnetization
        T1a % T1 decay constant for spin a
        T2a % T2 Decay constant for spin a
        ppma % chemical shift for spin a
        T1b % T1 decay constant for spin b
        T2b % T2 Decay constant for spin b
        ppmb % chemical shift for spin b
        T1c % T1 decay constant for spin c
        T2c % T2 Decay constant for spin c
        ppmc % chemical shift for spin c
        gamma % gyromagnetic ratio
        density % relative number of spins in this group
        kab % a to b exchange rate
        kba % b to a exchange rate
        kbc % b to c exchange rate
        kcb % c to b exchange rate
        b % input function
    end
    methods
        function self = ThreeSiteExchangeGroup(varargin)
            % Constructor - initializes the spin group
            % IsolatedSpinGrp(M,M0,T1,T2,gamma,density) - initializes the spin
            % group with an initial magnetization (M), equilibrium magentization
            % (M0), T1 decay (T2), T2 Decay (T2), gyromagnetic ratio (gamma) and
            % spin density (density)
            p = inputParser();
            p.addOptional('M',[0;0;1;0;0;1;0;0;1],@isnumeric)
            p.addOptional('M0',[0;0;0;0;0;0;0;0;0],@isnumeric)
            p.addOptional('T1a',56,@isnumeric)
            p.addOptional('T2a',0.02,@isnumeric)
            p.addOptional('ppma',171*1e-6,@isnumeric)
            p.addOptional('T1b',56,@isnumeric)
            p.addOptional('T2b',0.02,@isnumeric)
            p.addOptional('ppmb',171*1e-6,@isnumeric)
            p.addOptional('T1c',30,@isnumeric)
            p.addOptional('T2c',0.02,@isnumeric)
            p.addOptional('ppmc',185*1e-6,@isnumeric)
            p.addOptional('gamma',67.262e6,@isnumeric)
            p.addOptional('density',1,@isnumeric)
            p.addOptional('kab',0.0,@isnumeric)
            p.addOptional('kba',0.0,@isnumeric)
            p.addOptional('kbc',0.0,@isnumeric)
            p.addOptional('kcb',0.0,@isnumeric)
            p.addOptional('b',0.0)
            p.parse(varargin{:})
            if ~isempty(p.Results.M), self.M = p.Results.M;
            else self.M = [0;0;1;0;0;1;0;0;1];end
            if ~isempty(p.Results.M0),self.M0 = p.Results.M0;
            else self.M0 = [0;0;0;0;0;0;0;0;0];end
            if ~isempty(p.Results.T1a),self.T1a = p.Results.T1a;
            else self.T1a = 56;end
            if ~isempty(p.Results.T2a),self.T2a = p.Results.T2a;
            else self.T2a = 0.02;end
            if ~isempty(p.Results.ppma),self.ppma = p.Results.ppma;
            else self.ppma = 171*1e-6;end
            if ~isempty(p.Results.T1b),self.T1b = p.Results.T1b;
            else self.T1b = 56;end
            if ~isempty(p.Results.T2b),self.T2b = p.Results.T2b;
            else self.T2b = 0.02;end
            if ~isempty(p.Results.ppmb),self.ppmb = p.Results.ppmb;
            else self.ppmb = 171*1e-6;end
            if ~isempty(p.Results.T1c),self.T1c = p.Results.T1c;
            else self.T1c = 30;end
            if ~isempty(p.Results.T2c),self.T2c = p.Results.T2c;
            else self.T2c = 0.02;end
            if ~isempty(p.Results.ppmc),self.ppmc = p.Results.ppmc;
            else self.ppmc = 185*1e-6;end
            if ~isempty(p.Results.gamma),self.gamma = p.Results.gamma;
            else self.gamma = 67.262e6;end
            if ~isempty(p.Results.density),self.density = p.Results.density;
            else self.density = 1;end
            if ~isempty(p.Results.kab),self.kab = p.Results.kab;
            else self.kab = 0.0;end
            if ~isempty(p.Results.kba),self.kba = p.Results.kba;
            else self.kba = 0.0;end
            if ~isempty(p.Results.kbc),self.kbc = p.Results.kbc;
            else self.kbc = 0.0;end
            if ~isempty(p.Results.kcb),self.kcb = p.Results.kcb;
            else self.kcb = 0.0;end
            if ~isempty(p.Results.b),self.b = p.Results.b;
            else self.b = 0.0;end
        end
        function val = calculationFrame(self,B0)
            % CALCULATIONFRAME -  returns the frequencty of the rotating
            % refrence frame that dm is calculated in
            tmp = [0;0;1].*B0*self.gamma+...
                (1+mean([self.ppma,self.ppmb]));
            val = tmp(3);
        end
        function dm = dM(self,x,y,z,t,M,PS,B0)
            % DM: returns the delta m at some time and location and given M
            % dM(x,y,z,t,M) - returns a dm at some position (x,y,z), some time
            % (t), and some initial magnetization M ([Mx;My;Mz])
            Recovery(1:3) = 1/self.T1a*self.M0(1:3);
            Recovery(4:6) = 1/self.T1b*self.M0(4:6);
            Recovery(7:9) = 1/self.T1c*self.M0(7:9);
            dm = self.getA(x,y,z,t,PS,B0)*M+Recovery.';
        end
        function A = getA(self,x,y,z,t,PS,B0)
            B = repmat(B0,1,length(t))+PS.B(x,y,z,t);
            theta =  -self.calculationFrame(B0)*t;
            Beff(1:3) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
                0,0,1]*B-[0;0;1].*B0*...
                (1+mean([self.ppma,self.ppmb])-self.ppma);
            Beff(4:6) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
                0,0,1]*B-[0;0;1].*B0*...
                (1+mean([self.ppma,self.ppmb])-self.ppmb);
            Beff(7:9) = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
                0,0,1]*B-[0;0;1].*B0*...
                (1+mean([self.ppma,self.ppmb])-self.ppmc);
            A = zeros(6);
            A(1:3,1:3) = -self.gamma*...
                [0,-Beff(3),Beff(2);Beff(3),0,-Beff(1);-Beff(2),Beff(1),0]+...
                [-1/self.T2a,0,0;0,-1/self.T2a,0;0,0,-1/self.T1a];
            A(4:6,4:6) = -self.gamma*...
                [0,-Beff(6),Beff(5);Beff(6),0,-Beff(4);-Beff(5),Beff(4),0]+...
                [-1/self.T2b,0,0;0,-1/self.T2b,0;0,0,-1/self.T1b];
            A(7:9,7:9) = -self.gamma*...
                [0,-Beff(9),Beff(8);Beff(9),0,-Beff(7);-Beff(8),Beff(7),0]+...
                [-1/self.T2c,0,0;0,-1/self.T2c,0;0,0,-1/self.T1c];
            A = A+[-self.kab,0,0,self.kba,0,0,0,0,0;...
                0,-self.kab,0,0,self.kba,0,0,0,0;...
                0,0,-self.kab,0,0,self.kba,0,0,0;...
                self.kab,0,0,-self.kba-self.kbc,0,0,self.kcb,0,0;...
                0,self.kab,0,0,-self.kba-self.kbc,0,0,self.kcb,0;...
                0,0,self.kab,0,0,-self.kba-self.kbc,0,0,self.kcb;...
                0,0,0,self.kbc,0,0,-self.kcb,0,0;...
                0,0,0,0,self.kbc,0,0,-self.kcb,0;...
                0,0,0,0,0,self.kbc,0,0,-self.kcb];
        end
        function ret = useAnalytical(self)
            %USEANALYTICAL: determins if the Analytical Soultion shouldbe used
            %for the given spin group under the given conditions
            ret = all(self.M0 == 0);
        end
        function fun = analytical(self,x,y,z,tspan,M,PS,B0)
            A = getA(self,x,y,z,t,PS,B0);
            % ANALYTICAL: retun a function handle to the analytical soluton 
            if sum(integral(self.b,tspan(1),tspan(end),'ArrayValued',true))>1e-3;
            fun = @(t)expm(A*(t-tspan(1)))*M+integral(...
                @(t2)expm(A*(t-t2))*self.b(t2),tspan(1),t,'ArrayValued',true);
            else
                fun = @(t)expm(A*(t-tspan(1)))*M;
            end
        end
    end
    
end

