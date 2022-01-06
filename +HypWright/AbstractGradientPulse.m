classdef AbstractGradientPulse < HypWright.GradientPulse
    %GRADIENTPULSE Class to represent a gradient field
    %   Detailed explanation goes here
    properties (SetAccess = private)
        Gx
        Gy
        Gz
        tAxis
    end
    properties (Access = protected)
        bFun
    end
    properties (Constant)
        timeDependence = true;
    end
    methods
        function self = AbstractGradientPulse(center,bShape,tAxis, varargin)
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
            durration = tAxis(end)-tAxis(1);
            self = self@HypWright.GradientPulse(center, durration, p.Results.name);
            self.tAxis = tAxis+self.startTime;
            self.Gx = bShape(1,:);
            self.Gy = bShape(2,:);
            self.Gz = bShape(3,:);
            self.calB();
        end
    end
    methods (Access =  private)
        function calB(self)
                        self.bFun = @(x,y,z,t)[0,0,0;...
                            0,0,0;...
                            interp1(self.tAxis,self.Gx,t,[],0),...
                            interp1(self.tAxis,self.Gy,t,[],0),...
                            interp1(self.tAxis,self.Gz,t,[],0)]*[x;y;z]; 
        end
    end
end

