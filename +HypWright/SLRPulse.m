classdef SLRPulse < HypWright.RFPulse
    %SLRPULSE Generates a Pulse using the Shinnar Le-Roux Algorithm
    properties  
    end
    properties(SetAccess = private)
        bandwidth % Bandwidth of the sinc pulse
        amplitude % Amplitude of the Sinc
        rfShape
        tAxis

    end
    methods
        function self = SLRPulse(nPoints,bandwidth,amplitude,center,durration,omega,varargin)
            % SLRPulse - Constructor
            %    ptype -- pulse type.  Options are:
            %      st  -- small tip angle         (default)
            %      ex  -- pi/2 excitation pulse
            %      se  -- pi spin-echo pulse
            %      sat -- pi/2 saturation pulse
            %      inv -- inversion pulse
            %    ftype -- filter design method.  Options are:
            %      ms  -- Hamming windowed sinc (an msinc)
            %      pm  -- Parks-McClellan equal ripple
            %      ls  -- Least Squares           (default)
            %      min -- Minimum phase (factored pm)
            %      max -- Maximum phase (reversed min)
            %    d1 -- Passband ripple        (default = 0.01)
            %    d2 -- Stopband ripple        (default = 0.01)
            %    pclsfrac -- pcls tolerance   (default = 1.5)
            p = inputParser();
            p.addOptional('name',sprintf('Pulse%d',int16(rand(1)*10000)),@isstr)
            p.addParameter('ptype','st')
            p.addParameter('ftype','ls')
            p.addParameter('d1',0.01)
            p.addParameter('d2',0.01)
            p.addParameter('pclsfrac',1.5)
            p.parse(varargin{:})
            self = self@HypWright.RFPulse(center,durration,omega,p.Results.name);
            timeBWProduct = bandwidth*durration;
            self.rfShape = dzrf(nPoints,timeBWProduct,p.Results.ptype,...
                p.Results.ftype,p.Results.d1,p.Results.d2,p.Results.pclsfrac);
            self.amplitude = amplitude*nPoints/durration;
            self.bandwidth = bandwidth;
            self.tAxis = linspace(-durration/2,durration/2,nPoints);
            self.calB();
        end
        function setDurration(self,durration)
            % SETDURRATION - overloaded for a nLobed sincPulse as it should not
            % be changeable. Durration is a function of bandwidth
            % setDurration(self,durration) -  does nothing and warns the user
            disp(['the durration of this pulse is a function of bandwidth and'...
                ' should bealtered by changing thebandwidth']);
        end
        function setAmplitude(self,value)
            % SETAMPLITUDE - sets the amplitude of the pulse
            self.amplitude = value;
            self.calB();
        end
        function setBW(self,newBW)
                        disp(['the durration of this pulse is a function of bandwidth and'...
                ' should bealtered by changing thebandwidth']);
        end
    end
    methods (Access = protected)
        function calB(self)
            % CALB - re-calculates the function that defines the envelope for
            % this sinc pulse
            self.Bfun = @(x,y,z,t)self.amplitude*interp1(self.tAxis,...
                self.rfShape,t);
        end
    end
end

