classdef BanksonSpinGrp < HypWright.SpinGroup

    
    properties
        M % magnetization vector
        M0 % Equilibrium magnetization
        T1 % T1 decay constant
        T2 % T2 Decay constant
        gamma % gyromagnetic ratio
        ppm %  ppm shift of the spin
        density % relative number of spins in this group
        shapeTerms % terms that define the gamma pdf
        t0 %injection delay
    end
    methods
        function self = BanksonSpinGrp(M,M0,T1,T2,gamma,ppm,density,shapeTerms,t0)
            % Constructor - initializes the spin group
            % v(M,M0,T1,T2,gamma,ppm,density) - initializes the spin
            % group with an initial magnetization (M), equilibrium magentization
            % (M0), T1 decay (T2), T2 Decay (T2), gyromagnetic ratio (gamma),
            % some chemical shift (ppm), spin density (density), and shape terms
            % for the gamma pdf (shapeTerms)
            self.M = M;
            self.M0 = M0;
            self.T1 = T1;
            self.T2 = T2;
            self.gamma = gamma;
            self.ppm = ppm;
            self.density = density;
            self.shapeTerms = shapeTerms;
            self.t0 = t0;
            % War user if the use analytical will always be true, if so
            % this spin group will not account for perfusion
            if(~useAnalytical(self))
                warning('The Analytical Solution of a Bankson Spin group will no be used. Any Perfusion in that spin group will be ignored\n')
            end
        end
        function val = calculationFrame(self,B0)
            % CALCULATIONFRAME -  returns the frequencty of the rotating
            % refrence frame that dm is calculated in;
            tmp = [0;0;1].*B0*self.gamma*(1+self.ppm);
            val = tmp(3);
        end
        function dm = dM(self,x,y,z,t,M,PS,B0)
            % DM: returns the delta m at some time and location and given M
            % dM(x,y,z,t,M) - returns a dm at some position (x,y,z), some time
            % (t), and some initial magnetization M ([Mx;My;Mz])
            dm = self.getA(x,y,z,t,PS,B0)*M+1/self.T1*self.M0;
        end
        function A = getA(self,x,y,z,t,PS,B0,varargin)
            % GETA - gets the matrix that defines dm
            if ~isempty(varargin) == 1
                B = varargin{1};
            else
                B = repmat(B0,1,length(t))+PS.B(x,y,z,t);
            end
            theta =  -self.calculationFrame(B0)*t;
            Beff = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]*B-...
                [0;0;1].*B0;
            A = self.gamma*...
                [0,-Beff(3),Beff(2);Beff(3),0,-Beff(1);-Beff(2),Beff(1),0] + ...
                [-1/self.T2,0,0;0,-1/self.T2,0;0,0,-1/self.T1];
        end    
        function ret = useAnalytical(self)
            %USEANALYTICAL: determins if the Analytical Soultion shouldbe used
            %for the given spin group under the given conditions
            ret = all(self.M0 == 0);
        end
        function vals = analytical(self,x,y,z,t0,M,t,PS,B0,varargin)
            % ANALYTICAL: retun a function handle to the analytical soluton 
            A = self.getA(x,y,z,t0+1e-9,PS,B0,varargin{:});
            vals = cell2mat(arrayfun(@(t2)[subsref(expm(A*(t2-t0))*M,...
                struct('type','()','subs',{{1:2,1}}));...
                self.shapeTerms(3)*gampdf(...
                t2-self.t0,self.shapeTerms(1),self.shapeTerms(2))],t,...
                'UniformOutput',false));
        end
    end
end

