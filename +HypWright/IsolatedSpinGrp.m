classdef IsolatedSpinGrp < HypWright.SpinGroup
    %ISOLATEDSPINGRP a class that represents a set of isolated spins
    %   Detailed explanation goes here
    %   Properties
    %   M - magnetization vector
    %   M0 - Equilibrium magnetization
    %   T1 - T1 decay constant
    %   T2 - T2 Decay constant
    %   gamma - gyromagnetic ratio
    %   ppm -  ppm shift of the spin
    %   density - relative number of spins in this group
    %   Methods
    %   IsolatedSpinGrp(M,M0,T1,T2,gamma,density) - initializes the spin
    %   group with an initial magnetization (M), equilibrium magentization
    %   (M0), T1 decay (T2), T2 Decay (T2), gyromagnetic ratio (gamma)
    %   some chemical shift (ppm) and, spin density (density)
    %   calculationFrame() -  returns the frequencty of the rotating
    %   refrence frame that dm is calculated in
    %   dM(x,y,z,t,M) - returns a dm at some position (x,y,z), some time
    %   (t), and some initial magnetization M ([Mx;My;Mz])
    
    properties
        M % magnetization vector
        M0 % Equilibrium magnetization
        T1 % T1 decay constant
        T2 % T2 Decay constant
        gamma % gyromagnetic ratio
        ppm %  ppm shift of the spin
        density % relative number of spins in this group
    end
    methods
        function self = IsolatedSpinGrp(M,M0,T1,T2,gamma,ppm,density)
            % Constructor - initializes the spin group
            % IsolatedSpinGrp(M,M0,T1,T2,gamma,ppm,density) - initializes the spin
            % group with an initial magnetization (M), equilibrium magentization
            % (M0), T1 decay (T2), T2 Decay (T2), gyromagnetic ratio (gamma),
            % some chemical shift (ppm) and, spin density (density)
            self.M = M;
            self.M0 = M0;
            self.T1 = T1;
            self.T2 = T2;
            self.gamma = gamma;
            self.ppm = ppm;
            self.density = density;
        end
        function val = calculationFrame(self,B0)
            % CALCULATIONFRAME -  returns the frequencty of the rotating
            % refrence frame that dm is calculated in
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
%             vals = cell2mat(arrayfun(@(t2)expm(A*(t2-t0))*M,...
%                 t,'UniformOutput',false));
            vals = zeros(size(M,1),size(t,2));
            for i = 1:numel(t)
                vals(:,i) = expm(A*(t(i)-t0))*M;
            end
        end
    end
end

