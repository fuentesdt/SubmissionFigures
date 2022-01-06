classdef GradientPulse < handle
    %GRADIENTPULSE Class to represent a gradient field
    %   Detailed explanation goes here
    properties
        center
        durration
        name
    end
    properties (Abstract, Access = protected)
        bFun
    end
    properties (Dependent)
        startTime
        endTime
    end
    methods
        function self = GradientPulse(center, durration, name)
            % CONSTRUCTOR - initializes  a gradient pulse objects
            % GradientPulse(startTime, endTime, slope, magnitude) initializes a
            % gradient pulse with a start time, end time, magnitude, and slope,
            % will give the pulse a randome name
            % GradientPulse(...,name) - same as above but will give the puse a
            % specified name
            self.center = center;
            self.durration = durration;
            self.name = name;
        end
        function val = get.startTime(self),val=self.center-self.durration/2;end
        function val = get.endTime(self),val=self.center+self.durration/2;end
        function bOut = B(self,x,y,z,t)
            % B: returns the gradient Bfied at some point and time
            % bOut = B(x,y,z,t) returns the B-filed(bOut) ate somepoint (x,y,z)
            % some time t
            if length(t) == 1
                if ((t > self.startTime)&&(t < self.endTime))
                    bOut = self.bFun(x,y,z,t);
                else
                    bOut = zeros(3,1);
                end
            else
            bOut = zeros(3,length(t));
            pulseOnTimes = find(self.startTime -t < 1e-9 & t < self.endTime);
            if ~isempty(pulseOnTimes)
            bOut(:,pulseOnTimes) = self.bFun(x,y,z,t(pulseOnTimes));
            end
            end
        end
        function display(self,varargin)
            p = inputParser();
            p.addOptional('axis',[])
            p.parse(varargin{:})
            if isempty(p.Results.axis)
                figure
                curAxis = gca;
            else
                curAxis = p.Results.axis;
            end
            span = -1:0.3:1;
            [X,Y,Z] = meshgrid(span,span,span);
            U = zeros(size(X));
            V = zeros(size(Y));
            W = zeros(size(Z));
            centerTime = (self.startTime+self.endTime)/2;
            for i = 1:numel(span)
                for j = 1:numel(span)
                    for k = 1:numel(span)
                        tmpVect = self.B(i,j,k,centerTime);
                        U(i,j,k) = tmpVect(1);
                        V(i,j,k) = tmpVect(2);
                        W(i,j,k) = tmpVect(3);
                    end
                end
            end
            quiver3(curAxis,X,Y,Z,U,V,W)
        end
    end
    methods (Access =  private)
        function calB(self)
            % CALB: recalculates the function defining the gradiaent pulse
            self.bFun = @(x,y,z,t)[0,0,0;0,0,0;self.slope]*[x;y;z];
        end
    end
end

