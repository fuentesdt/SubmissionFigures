classdef BrukerPulse < HypWright.RFPulse
    %SLRPULSE Generates a Pulse using the Shinnar Le-Roux Algorithm
    properties  
    end
    properties(SetAccess = private)
        bandwidth % Bandwidth of the sinc pulse
        amplitude % Amplitude of the Pulse
        rfShape
        tAxis
        nPoints
        bwFact
        intFact
        type

    end
    methods
        function self = BrukerPulse(center,type,bandwidth,amplitude,omega,varargin)
            p = inputParser();
            p.addOptional('name',sprintf('Pulse%d',int16(rand(1)*10000)),@isstr)
            p.parse(varargin{:})
            self = self@HypWright.RFPulse(center,0,omega,p.Results.name);
            self.bandwidth = bandwidth;
            self.getShape(type);
            self.type = sprintf('Bruker.%s',type);
            self.durration = self.bwFact/self.bandwidth;
            self.amplitude = amplitude/self.durration/self.intFact/max(self.rfShape);
            self.tAxis = linspace(-self.durration/2,self.durration/2,length(self.rfShape));
            self.calB();
        end
    end
    methods (Access = protected)
        function calB(self)
            % CALB - re-calculates the function that defines the envelope for
            % this sinc pulse
            self.Bfun = @(x,y,z,t)self.amplitude*interp1(self.tAxis,...
                self.rfShape,t);
        end
        function getShape(self,name)
            %% Get Package Location
            mc = metaclass(self);
            s = what(mc.ContainingPackage.Name);
            path = s.path;
            BrukerFile = sprintf('%s/BrukerPulses.mat',path);
            if exist(BrukerFile,'file') ~= 2
                error('Bruker Pulse Shape File Corrupted (BrukerPulses.mat), check HypWright library for file.')
            end
            load(BrukerFile)
            if isfield(shapes,name)
                self.rfShape = shapes.(name);
                self.nPoints = headers.(name).NPOINTS;
                self.bwFact = headers.(name).SHAPE_BWFAC;
                self.intFact = headers.(name).SHAPE_INTEGFAC;
            else
                error('input name: %s not in Bruker shape file')
            end
        end
    end
end

