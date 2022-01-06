classdef AbstractRFPulse < HypWright.RFPulse
    %SLRPULSE Generates a Pulse using the Shinnar Le-Roux Algorithm
    properties  
    end
    properties(SetAccess = private)
        amplitude % Amplitude of the Sinc
        rfShape
        tAxis

    end
    methods
        function self = AbstractRFPulse(center,omega,rfShape,tAxis,varargin)
            p = inputParser();
            p.addOptional('name',sprintf('Pulse%d',int16(rand(1)*10000)),@isstr)
            p.parse(varargin{:})
            durration = tAxis(end)-tAxis(1);
            self = self@HypWright.RFPulse(center,durration,omega,p.Results.name);
            self.rfShape = rfShape;
            self.tAxis = tAxis;
            self.calB();
        end
        function modulateRF(self,tAxis,modulation)
            tAxis = tAxis-tAxis(floor(length(tAxis)/2));
            if tAxis ~= self.tAxis
                modulation = resample(tAxis,modulation,self.tAxis);
            end
            self.rfShape = self.rfShape.*modulation;
            self.calB();
        end
    end
    methods (Access = protected)
        function calB(self)
            % CALB - re-calculates the function that defines the envelope for
            % this sinc pulse
            self.Bfun = @(x,y,z,t)interp1(self.tAxis,...
                self.rfShape,t,'linear',0);
        end
    end
end

