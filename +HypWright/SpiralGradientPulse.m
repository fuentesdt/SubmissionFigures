classdef SpiralGradientPulse < HypWright.GradientPulse
    %GRADIENTPULSE Class to represent a gradient field
    %   Detailed explanation goes here
    properties (SetAccess = private)
        slope
        fov
        npix
        arms
        ksamp
        kMap
        dcf
        nSpiralPoints
        gMax
        slew
        gamma
        readTimes
    end
    properties (Access = protected)
        bFun
    end
    properties (Constant)
        timeDependence = true;
    end
    methods
        function self = SpiralGradientPulse(...
                start,fov,npix,arms,ksamp,gMax,slew,gamma,varargin)
            % CONSTRUCTOR - initializes  a gradient pulse objects
            % GradientPulse(startTime, endTime, slope, magnitude) initializes a
            % gradient pulse with a start time, end time, magnitude, and slope,
            % will give the pulse a randome name
            % GradientPulse(...,name) - same as above but will give the puse a
            % specified name
            p = inputParser();
            p.addOptional('name',sprintf('GradPulse%d',int16(rand(1)*10000)),...
                @isstr)
            p.parse(varargin{:})
            self = self@HypWright.GradientPulse(0,0,p.Results.name);
            if ~exist('gMax','var'), gMax = []; end
            if isempty(gMax), gMax = 40; end;           % [mT/m] max gradient amplitude
            if ~exist('slew','var'), slew = []; end
            if isempty(slew), slew = 150; end;          % [T/m/s] max slew rate
            if ~exist('gamma','var'), gamma = []; end
            if isempty(gamma), gamma = 267.513e6; end;          % [rads/sec/T] Gyromagnetic ratio for Protons
            self.fov = fov;
            self.npix = npix;
            self.arms = arms;
            self.ksamp = ksamp;
            self.gMax =  gMax;
            self.slew = slew;
            self.gamma = gamma;
            self.calB(start);
        end
        function [dataCart] = reGrid(self,dataSpiral)
            %Make Data a Row vector for fidAll function
            if size(dataSpiral,1)  > size(dataSpiral,2)
                dataSpiral = dataSpiral.';
            end
            dataCart = gridding(dataSpiral,self.kMap,self.dcf,self.npix);
        end
        function display(self)
            figure
            t = linspace(self.startTime,self.endTime,2048);
            for i = 1:length(t)
            tmp = self.bFun(1,0,0,t(i));
            Gx(i) = tmp(3);
            tmp = self.bFun(0,1,0,t(i));
            Gy(i) = tmp(3);
            tmp = self.bFun(0,0,1,t(i));
            Gz(i) = tmp(3);
            end
            plot(t,Gx,t,Gy,t,Gz);
            xlabel('Time (sec)'),ylabel('Gradient (T/m)')
            legend('Gx','Gy','Gz')
        end
    end
    methods (Access =  private)
        function calB(self,start)
            % CALB: recalculates the function defining the gradiaent pulse
                        %figure
                        [k1,g1,s1,t1,r1,theta1] = vds(0.99*self.slew*1e2,0.99*self.gMax*1e-1,self.ksamp*1e-6, ...
                            self.arms,[self.fov/10*self.gamma/267.513e6,0],self.npix/(2*self.fov/10*self.gamma/267.513e6));
                        close(gcf)
                        close(gcf)
                        self.kMap = k1./max(k1)*0.5;
                        self.dcf = voronoi_area(k1*self.npix);
                        Gx = real(g1)*1e-2;
                        Gy = imag(g1)*1e-2;
                        tAxis = t1-t1(1)+start;
                        self.readTimes = tAxis;
                        self.nSpiralPoints = size(self.kMap,2);
                        self.center = tAxis(floor(length(tAxis)/2));
                        self.durration = tAxis(end)-tAxis(1);
                        self.bFun = @(x,y,z,t)[0,0,0;...
                            0,0,0;...
                            interp1(tAxis,Gx,t,[],0),interp1(tAxis,Gy,t,[],0),0]*[x;y;z];
            
%             [self.kMap,self.dcf,t,self.ind,out]=design_spiral(self.fov,self.npix,self.arms,...
%                 self.ksamp,[],self.gMax*1e3,self.slew,self.nucleus);
%             close(gcf)
%             close(gcf)
%             self.nSpiralPoints = size(self.kMap,2);
%             C =  find(self.ind==1);
%             Gx = real(out.grad(C));
%             Gy = imag(out.grad(C));
%             tAxis = t+start;
%             figure
%             self.durration = t(end);
%             self.center = start+self.durration/2;
%             self.bFun = @(x,y,z,t)[0,0,0;...
%                 0,0,0;...
%                 interp1(tAxis,Gx,t,[],0),interp1(tAxis,Gy,t,[],0),0]*[x;y;z];
            
        end
    end
end

