classdef RFPulse < handle
    %RFPULSE The base Class or Radio frequency Pulses
    %   Properties
    %   Bfun: function pointer defining the pulse
    %   center: the center of the pulse
    %   durration: the length of the pulse (truncates Bfun otherwise)
    %   omega: carrier frequency of the pulse
    %   Name: a (hopefully) unique name for the pulse
    %   startTime: start time of the pulse (truncate Bfun before this point)
    %   endTime: end time of the pulse (truncates Bfun after this point)
    %   Methods
    %   RFPulse(center, durration,omega, name) initializes center time (center),
    %   pulse durration (durraton), carrier frequency (omega), and name
    %   display(self) displays in a new figure 
    %   display(self,h) displays in a passed in figure h 
    %   B1(x,y,z,t) returns a 3D vector defining B1 for this pulse at
    %   the passed in postion (x,y,z) and time (t)
    properties (SetAccess = protected)
        Bfun % function pointer defining the pulse
        center % the center of the pulse
        durration % the length of the pulse (truncates Bfun otherwise)
        omega % carrier frequency of the pulse
        name % a (hopefully) unique name for the pulse
    end
    properties (Dependent)
        startTime % start time of the pulse
        endTime % end time of the pulse
    end
    properties (Constant)
        timeDependence = true;
    end
    methods (Abstract = true, Access = protected)
        calB(self) % returns a pointer to the function defining the pulse
    end
    methods
        function setCenter(self,center),self.center = center; end
        function setDurration(self,durration), self.durration = durration; end
        function setOmega(self,omega), self.omega = omega; end
        function val = get.startTime(self),val=self.center-self.durration/2;end
        function val = get.endTime(self),val=self.center+self.durration/2;end
        function self = RFPulse(center,durration,omega,name)
            % CONSTRUCTOR - Base Constructo for all RFPulse subclasses
            % RFPulse(center, durration,omega, name) initializes center time,
            % durration, carrier frequency omega, and name
            self.center = center;
            self.durration = durration;
            self.omega = omega;
            self.name = name;
        end
        function display(self,varargin)
            % DISPLAY - Displays the RF pulse evelope, in both frequency and
            % time domaines 
            % display(self) displays in a new figure 
            % display(self,h) displays in a passed in figure h
            p = inputParser();
            p.addOptional('figure',[])
            p.parse(varargin{:})
            N = 2^10;
            t = linspace(-self.durration/2,self.durration/2,N);
            B1 = zeros(length(t),1);
            x = 0; y = 0; z = 0;
            for i = 1:length(t)
                B1(i) = self.Bfun(x,y,z,t(i));
            end
            FT = fftshift(fft(fftshift(B1)));
            freqAxis = linspace(1/(t(2)-t(1))/length(t),1/(t(2)-t(1)),...
                length(t)); % Calculate frequency axis
            if(isempty(p.Results.figure))
                figure;
            else
                figure(p.Results.figure);
            end
            subplot(2,1,1),plot(t,real(B1),'k',t,imag(B1),'r',t,abs(B1),'b')
            xlabel('Time (seconds)')
            ylabel('B1 (Tesla)')
            legend('X','Y','Magnitude')
            subplot(2,1,2),plot(freqAxis,real(FT),'k',freqAxis,imag(FT),'r',...
                freqAxis,abs(FT),'b')
            xlabel('Frequency (Hz)')
            ylabel('Magnitude (arb)')
            legend('real','Imaginary','Magnitude')
        end
        function B1 = B(self,x,y,z,t)
            % B1: gives the B1 of this pulse at a time and location
            % B1(x,y,z,t) returns a 3D vector defining B1 for this pulse at
            % the passed in postion (x,y,z) and time (t)
            B1 = zeros(1,length(t));
            pulseOnTimes = find(self.startTime < t & t < self.endTime);
            if ~isempty(pulseOnTimes)
                B1(pulseOnTimes) = self.Bfun(x,y,z,t(pulseOnTimes)-self.center)...
                    .*exp(1i*self.omega*(t(pulseOnTimes)));
            end
        end
    end
end

