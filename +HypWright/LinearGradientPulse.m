classdef LinearGradientPulse < HypWright.GradientPulse
    %GRADIENTPULSE Class to represent a gradient field
    %   Detailed explanation goes here
    properties (SetAccess = private)
        slope
    end
    properties (Access = protected)
        bFun
    end
    properties (Constant)
        timeDependence = false;
    end
    methods
        function self = LinearGradientPulse(center, durration, slope,varargin)
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
            self = self@HypWright.GradientPulse(center,durration,p.Results.name);
            self.slope = slope;
            self.calB();
        end
        function setSlope(self,slope)
           self.slope = slope;
           self.calB();
        end
    end
    methods (Access =  private)
        function calB(self)
            % CALB: recalculates the function defining the gradiaent pulse
            self.bFun = @(x,y,z,t)[0,0,0;0,0,0;self.slope]*[x;y;z];
        end
    end
end

