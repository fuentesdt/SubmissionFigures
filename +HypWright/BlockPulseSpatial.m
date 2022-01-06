classdef BlockPulseSpatial < HypWright.BlockPulse
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
        uniFun
    end
    methods
        function self = BlockPulseSpatial(center,durration,omega,amplitude,uniFun,varargin)
            % CONSTRUCTOR - Initializes the Block pulse
            % BlockPulse(center,bandwidth,amplitude,omega,varargin) - sets the 
            % center time (center), the pulse durration (durration), the pulse
            % amplitude (amplitude) and the carrier frequencey (omega) 
            % BlockPulse(...,name) - same as above but the last argument will be
            % set as the pulses' name
            self = self@HypWright.BlockPulse(center,durration,omega,amplitude,varargin{:});
            self.uniFun = uniFun;
            self.calB();
        end
        function setUniformity(self,uniFun)
            self.uniFun = uniFun;
        end
    end
    methods (Access = protected)
        function calB(self)
            % CALB - re-calculates the function that defines the envelope for
            % this sinc pulse
            self.Bfun = @(x,y,z,t)self.amplitude.*self.uniFun(x,y,z);
        end
    end
end

