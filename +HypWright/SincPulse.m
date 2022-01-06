classdef SincPulse < HypWright.RFPulse
    %SINCPULS Sinc Enveloped RF Pulse
    %   Properties
    %   bandwidth: Bandwidth of the sinc pulse
    %   amplitude: Amplitude of the Sinc
    %   lobes: number of lobes in the sinc envelope default(5)
    %   Methods
    %   SincPulse(center,bandwidth,amplitude,omega,varargin) - sets the 
    %   center time (center), the pulse bandwidth (bandwidth), the pulse
    %   amplitude (amplitude) and the carrier frequencey (omega) Will set
    %   the name to a random number and the number of lobes to 5
    %   SinPulse(...,lobes) - same as above but accepts a posotive integer
    %   as the number of lobes in this pulse
    %   SinPulse(...,name) - same as above but the last argument will be
    %   set as the pulses' name, to use the default value (5) for the
    %   number of lobes just pass in [] as the 5th argument
    %   setDurration(self,durration) -  does nothing and warns the user note
    %   that for a sinc pulse the durration is a function of the bandwidth
    properties  
    end
    properties(SetAccess = private)
        bandwidth % Bandwidth of the sinc pulse
        amplitude % Amplitude of the Sinc
        lobes % number of lobes in the sinc envelope
    end
    methods
        function self = SincPulse(center,bandwidth,amplitude,omega,varargin)
            % CONSTRUCTOR - Initializes the Sinc pulse
            % SincPulse(center,bandwidth,amplitude,omega,varargin) - sets the 
            % center time (center), the pulse bandwidth (bandwidth), the pulse
            % amplitude (amplitude) and the carrier frequencey (omega) Will set
            % the name to a random number and the number of lobes to 5
            % SinPulse(...,lobes) - same as above but accepts a posotive integer
            % as the number of lobes in this pulse
            % SinPulse(...,name) - same as above but the last argument will be
            % set as the pulses' name, to use the default value (5) for the
            % number of lobes just pass in [] as the 5th argument
            function val = lobeTest(x)
                if isempty(x)
                    val = 1;
                else
                    val = (mod(x,1) == 0 && x > 0);
                end
            end
            p = inputParser();
            p.addOptional('lobes',5,@lobeTest)
            p.addOptional('name',sprintf('Pulse%d',int16(rand(1)*10000)),@isstr)
            p.parse(varargin{:})
            self = self@HypWright.RFPulse(center,0,omega,p.Results.name);
            self.amplitude = amplitude;
            self.bandwidth = bandwidth;
            if isempty(p.Results.lobes)
                self.lobes = 5;
            else
                self.lobes = p.Results.lobes;
            end
            self.durration = ((self.lobes))*2/self.bandwidth;
            self.calB();
        end
        function setDurration(self,durration)
            % SETDURRATION - overloaded for a nLobed sincPulse as it should not
            % be changeable. Durration is a function of bandwidth
            % setDurration(self,durration) -  does nothing and warns the user
            disp(['the durration of this pulse is a function of bandwidth and'...
                ' should be altered by changing thebandwidth']);
        end
        function setAmplitude(self,value)
            % SETAMPLITUDE - sets the amplitude of the pulse
            self.amplitude = value;
            self.calB();
        end
        function setBW(self,newBW)
            % SETBW - Sets the bandwidth of the pulse
            % SetBW(self, newBW) - sets thebandwidth to some ne wbandwidth
            % (newBW)
            self.functionBW = newBW;
            self.durration = ((self.lobes)-1)*2/self.bandwidth;
            self.calB();
        end
    end
    methods (Access = protected)
        function calB(self)
            % CALB - re-calculates the function that defines the envelope for
            % this sinc pulse
            self.Bfun = @(x,y,z,t)self.amplitude.*self.bandwidth.*...
                sinc(self.bandwidth.*t).*...
                interp1(-self.durration/2:self.durration/100:self.durration/2,...
                blackman(...
                length(-self.durration/2:self.durration/100:self.durration/2))...
                ,t);
        end
    end
end

