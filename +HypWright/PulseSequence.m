classdef PulseSequence < handle
    %PULSESEQUENCE Classrepresenting a MRI pusle sequence
    %   Detailed explanation goes here
    properties
        ADC = []; % A vector storing start, stop time and bandwidth of the ADC
    end
    properties (SetAccess = private)
        rfPulses = {} % Vector of rf Pulses in the sequence
        gradientPulses = {} % Vector of static Gradient Pulses in the sequence
        calcGradPulses = {} % Vector of time varying gradient pulses
        time = {}
        eventTimes % vector defining when pulses are on or off
        slewRate = 1e9 % maximal gradient slope in T/(m*sec)
        gradientVect % The gradient amplitudes
        RFVect % Grid Containg the Times for each RF Pulse
        compiled = false;
    end
    properties (Constant)
        maxSize = 1e4; % Maximum number of points used for displaying pulses
    end
    properties (Dependent)
       timeStep % Not sure whatthis does
        
    end
    methods
        function compile(self)
            if (isempty(self.gradientPulses))
                self.gradientVect = [0,1e9;0,0;0,0;0,0];
            else
                self.compileGrads();
            end
            if(isempty(self.rfPulses))
                self.RFVect = [0,1e9;0,0];
            else
                self.compileRF();
            end
            self.updateTime()
            self.compiled = true;
        end
        function val = get.timeStep(self)
            % timeStep: not quite sure whatthis does
            val = 1;
           for i = 1:length(self.eventTimes)-1
               if(self.eventTimes(i+1)-self.eventTimes(i)<val)
                   val = self.eventTimes(i+1)-self.eventTimes(i);
               end
           end
        end
        function b = B(self,x,y,z,t)
            % B: retuns the Bfield at some point defined by (x,y,z) at some time
            % t
            b = zeros(3,length(t));
            RFPulses = interp1(self.RFVect(1,:),self.RFVect(2,:),t,'nearest','extrap');
            for i = 1:length(t)
                if RFPulses(i) ~= 0
                    b(:,i) = [real(self.rfPulses{RFPulses(i)}.B(x,y,z,t(i)));...
                        imag(self.rfPulses{RFPulses(i)}.B(x,y,z,t(i)));0];
                end
            end
            % Old implementation
%             for i = 1:numel(self.rfPulses)
%                 if(t>=self.rfPulses{i}.startTime && t<=self.rfPulses{i}.endTime)
%             b = b + [real(self.rfPulses{i}.B(x,y,z,t));...
%                 imag(self.rfPulses{i}.B(x,y,z,t));zeros(1,length(t))];
%                 end
%             end
            if(size(self.gradientVect,2) >1 )
            b = b + self.BGrad(x,y,z,t);
            end
        end
        function b = BGrad(self,x,y,z,t)
            % B: retuns the Bfield of just the gradient fields at some time t
            % and some location x,y,z
            b = zeros(3,length(t));
            b(3,:) = interp1(self.gradientVect(1,:),self.gradientVect(2,:),t,...
                'nearest','extrap')*x + ...
                interp1(self.gradientVect(1,:),self.gradientVect(3,:),t,...
                'nearest','extrap')*y + ...
                interp1(self.gradientVect(1,:),self.gradientVect(4,:),t,...
                'nearest','extrap')*z;
%             b = zeros(3,length(t));
            for i = 1:numel(self.calcGradPulses)
                b = b + self.calcGradPulses{i}.B(x,y,z,t);
            end
        end
        function addPulse(self,Pulse)
            % ADDPULSE: adds a pulse to the Pulse Sequence
            RFPulseList = {'HypWright.SincPulse','HypWright.BlockPulse',...
                'HypWright.BlockPulseSpatial','HypWright.SincPulseSpatial',...
                'HypWright.SLRPulse','HypWright.BrukerPulse',...
                'HypWright.AbstractRFPulse'};
            GradientPulseList = {'HypWright.GradientPulse',...
                'HypWright.LinearGradientPulse'};
            calcGradPulseList = {'HypWright.SpiralGradientPulse',...
                'HypWright.AbstractGradientPulse'};
            if(find(ismember(class(Pulse),RFPulseList)))
                self.rfPulses{length(self.rfPulses)+1} = Pulse;
                self.compiled = false;
            else if(find(ismember(class(Pulse),GradientPulseList)))
                    self.gradientPulses{length(self.gradientPulses)+1} = Pulse;
                    self.compiled = false;
                else if(find(ismember(class(Pulse),calcGradPulseList)))
                        self.calcGradPulses{length(self.calcGradPulses)+1} = Pulse;
                        self.compiled = false;
                    else error(['Pulse passed in not a recognized pulse type.'...
                            'Consider updating the pulse lists in the PS object'])
                    end
                end
            end
        end
        function removePulse(self,n)
            % REMOVEPULSE: removes the nth pulse in the pulse sequence
            if(~(length(self.rfPulses) > n || n > 0 || isscalar(n)))
                disp('Error! Pulse index selected to delete not in the sequence')
                return
            end
            self.rfPulses(n) = [];
            self.compiled = false;
        end
        function addADC(self,startTime,bandwidth,nPoints)
            % ADDADC: adds a virtual ADC to record the MR data.
            % addADC(self,StartTime,Bandwidth,nPoints) (StartTiem), sampling rate (Bandwidth) and number of
            % points (nPoints). The lenght of time the ADC is on is defined by
            % nPoints/Bandwidth
            endTime = startTime+nPoints*(1/bandwidth);
            self.ADC(end+1,:) = [startTime,endTime,bandwidth];
        end
        function removeADC(self,n)
            % REMOVEADC - removes an ADC
            % removeADC(n) - removes the nth (n) ADC from the pulse sequence
           self.ADC(n,:) = [];
        end
        % TODO this seems like an old method and should be removed
        function [value,isterminal,direction] = events(self,t)
            % EVENTS: uses the stored pollynomilas to tell the solver when RF
            % pulses turn on and off
            value = polyval(self.eventTimes,t);
            isterminal = 1;
            direction = 0;
        end
        % TODO: it is not very intuitive to remove pulses based on some arbitray
        % number assigned when they were created. need a beeter method to
        % identify pulses
        function clearRFPulses(self)
            %CLEARPULSES: Clear all pulses from the pulse sequence
            self.rfPulses = [];
            self.updateTime();
            self.compiled = false;
        end
        % TODO: rework the ADC display feature to have ADC be accounted for in
        % the sampling time
        function display(self,varargin)
            % DISPLAY - displays all the RF,aand gradient pulses as well as when
            % ADCS are on
            % display() - displays in a new figure with thedefault time range )
            % to 100 seconds
            % display([startTime,endTime]) displays in a new figure from the
            % start time to the end time
            % display([startTime,endTime], figure) same as above but plots in
            % the passed in figure
            self.compile();
            p = inputParser();
            p.addOptional('timeRange',[0,self.eventTimes(end)])
            p.addOptional('figure',[])
            p.parse(varargin{:})
            if(isempty(p.Results.figure))
                figure('units','normalized','outerposition',[0.25 0 0.5 1]);
            else
                figure(p.Results.figure);
            end
            for i = 1:numel(self.time)
                if self.time(i)> p.Results.timeRange(1) &&...
                        self.time(i)< p.Results.timeRange(2)
                dispTime(i) = self.time(i);
                end
                if self.time(i)> p.Results.timeRange(2)
                    break
                end
            end
            if dispTime(1) > self.time(1)
                dispTime = [self.time(1),dispTime];
            end
            if dispTime(1) > p.Results.timeRange(1)
                dispTime = [p.Results.timeRange(1),dispTime];
            end
            if dispTime(end) < p.Results.timeRange(2)
                dispTime = [dispTime,p.Results.timeRange(2)];
            end
            x=0;y=0;z=0;
            GXDisp = zeros(size(dispTime));
            GYDisp = zeros(size(dispTime));
            GZDisp = zeros(size(dispTime));
            for i = 1:numel(dispTime)
                for j = 1:numel(self.gradientPulses)
                    if(dispTime(i) < self.gradientPulses{j}.endTime && ...
                            dispTime(i) > self.gradientPulses{j}.startTime)
                        tmpGVect = self.gradientPulses{j}.slope;
                        GXDisp(i) = GXDisp(i) + tmpGVect(1);
                        GYDisp(i) = GYDisp(i) + tmpGVect(2);
                        GZDisp(i) = GZDisp(i) + tmpGVect(3);
                    end
                end
                for j = 1:numel(self.calcGradPulses)
                    if(dispTime(i) < self.calcGradPulses{j}.endTime && ...
                            dispTime(i) > self.calcGradPulses{j}.startTime)
                        tmpGVect = self.calcGradPulses{j}.B(1,1,1,dispTime(i));
                        GXDisp(i) = GXDisp(i) + tmpGVect(1);
                        GYDisp(i) = GYDisp(i) + tmpGVect(2);
                        GZDisp(i) = GZDisp(i) + tmpGVect(3);
                    end
                end
            end
            subplot(5,1,1),plot(dispTime,zeros(size(dispTime)))
            hold on
            for i = 1:numel(self.rfPulses)
                plot([self.rfPulses{i}.center,self.rfPulses{i}.center],[0,1])
                text(self.rfPulses{i}.center,1,self.rfPulses{i}.name)
            end
            hold off
            xlabel('Time (seconds)'),ylabel('RF Magnitude'),title('RF Pulses')
            legend('Real','Iamginary','Magnitude')
            subplot(5,1,2),plot(dispTime,GXDisp,'k')
            xlabel('Time (seconds)'),ylabel('Gradient Slope')
            title('X Gradient'), 
            subplot(5,1,3),plot(dispTime,GYDisp,'k')
            xlabel('Time (seconds)'),ylabel('Gradient Slope')
            title('Y Gradient')
            subplot(5,1,4),plot(dispTime,GZDisp,'k')
            xlabel('Time (seconds)'),ylabel('Gradient Slope')
            title('Z Gradient')
            subplot(5,1,5), plot(self.eventTimes(:,1),self.eventTimes(:,2))
            xlabel('Time (seconds)'),ylabel('Using Analytic Solution')
            title('Solver Type')
        end
        function S = solver(self,eventTimes)
            % SOLVER - returns the time dependency in the pulse sequence
            S = false;
            for i = 1:numel(self.rfPulses)
                if (self.rfPulses{i}.startTime - eventTimes(1)) < 1e-8 &&...
                        (self.rfPulses{i}.endTime - eventTimes(2)) > -1e-8 
                    S = S || self.rfPulses{i}.timeDependence;
                end
            end
            for i = 1:numel(self.calcGradPulses)
                if (self.calcGradPulses{i}.startTime - eventTimes(1)) < 1e-8 &&...
                        (self.calcGradPulses{i}.endTime - eventTimes(2)) > -1e-8 
                    S = S || self.calcGradPulses{i}.timeDependence;
                end
            end
        end
        function eventTimes = getEventTimes(self,endTime)
            %% Determins the relevent event times for the pulse sequence
            % when evaluated over a time frame (evalTimes)
            pulseTimes = self.eventTimes(1,:); % Times when B changes in the Pulse sequence
            pulseTimes = pulseTimes(find(pulseTimes>0,1,'first'):end);
            pulseTimes = pulseTimes(pulseTimes<=endTime);
            pulseTimes = [0,pulseTimes,endTime];
            pulseTimes = sort(pulseTimes);
            pulseTimes = unique(round(pulseTimes.*1e9))./1e9;
            eventTimes = unique(pulseTimes);
        end
    end
    methods (Access =  private)
        function compileGrads(self)
            %A sub function that Converts the Gradien Pulses into a single
            % 4 by n vector that has the valuse for each gradient direction and
            % the times those valuse change
            %Initialize the gradient vector with the first gradient pulse
            % Calculate the slew time
            slewTime = max(abs(0-self.gradientPulses{1}.slope(:)))/self.slewRate;
            % Fill the time vectors for the pulse
            tmpVect(1,1) = self.gradientPulses{1}.startTime;
            tmpVect(1,2) = self.gradientPulses{1}.startTime+slewTime;
            tmpVect(1,3) = self.gradientPulses{1}.endTime-slewTime;
            tmpVect(1,4) = self.gradientPulses{1}.endTime;
            % Fill the gradient slopes for the pulse
            tmpVect(2:4,1) = 0;
            tmpVect(2:4,2) = self.gradientPulses{1}.slope(:);
            tmpVect(2:4,3) = 0;
            tmpVect(2:4,4) = -self.gradientPulses{1}.slope(:);
            % Add the rest of the gradient to the pulse
            for i=2:numel(self.gradientPulses)
                % get the gradent slop valuse at the begining and end of  the
                % pulse
                startSlope = [interp1(tmpVect(1,:),tmpVect(2,:),...
                    self.gradientPulses{1}.startTime,'nearest','extrap');...
                    interp1(tmpVect(1,:),tmpVect(3,:),...
                    self.gradientPulses{1}.startTime,'nearest','extrap');...
                    interp1(tmpVect(1,:),tmpVect(4,:),...
                    self.gradientPulses{1}.startTime,'nearest','extrap')];
                endSlope = [interp1(tmpVect(1,:),tmpVect(2,:),...
                    self.gradientPulses{1}.endTime,'nearest','extrap');...
                    interp1(tmpVect(1,:),tmpVect(3,:),...
                    self.gradientPulses{1}.endTime,'nearest','extrap');...
                    interp1(tmpVect(1,:),tmpVect(4,:),...
                    self.gradientPulses{1}.endTime,'nearest','extrap')];
                % Calculate the slew times
                slewTimeStart = max(abs(startSlope...
                    -self.gradientPulses{1}.slope(:)))/self.slewRate;
                slewTimeEnd = max(abs(endSlope...
                    -self.gradientPulses{1}.slope(:)))/self.slewRate;
                j = (i-1)*4+1; % counter for tmpVect
                % Fill the time vectors for the pulse
                tmpVect(1,j) = self.gradientPulses{i}.startTime;
                tmpVect(1,j+1) = self.gradientPulses{i}.startTime+slewTimeStart;
                tmpVect(1,j+2) = self.gradientPulses{i}.endTime-slewTimeEnd;
                tmpVect(1,j+3) = self.gradientPulses{i}.endTime;
                % Fill the gradient slopes for the pulse
                tmpVect(2:4,j) = 0;
                tmpVect(2:4,j+1) = self.gradientPulses{i}.slope(:);
                tmpVect(2:4,j+2) = 0;
                tmpVect(2:4,j+3) = -self.gradientPulses{i}.slope(:);
                [~,I]=sort(tmpVect(1,:)); % Sort Times
                tmpVect = tmpVect(:,I); % match slopes to Times
                % Combine duplicates
                tmpI = 1;
                sumI = 1;
                trashI = [];
                for k = 2:size(tmpVect,2)
                    % grab all the slop changes at a particular time
                    if (tmpVect(1,k) == tmpVect(1,tmpI))
                        sumI = [sumI,k]; % cant think of a way to pre-allocate this
                        % sum all identical slope changes then move to next time point
                    else
                        tmpVect(2:4,sumI(1)) = sum(tmpVect(2:4,sumI),2);
                        trashI = [trashI,sumI(2:end)];
                        tmpI = k;
                        sumI = k;
                    end
                end
                tmpVect(:,trashI) = [];
            end
            tmpBSlope = zeros(3,1);
            self.gradientVect = zeros(size(tmpVect));
            for i=1:size(tmpVect,2)
                self.gradientVect(1,i) = tmpVect(1,i);
                tmpBSlope = tmpBSlope+tmpVect(2:4,i);
                self.gradientVect(2:4,i) = tmpBSlope;
            end
            % Debug Plotting
%             figure('Position',[700,200,1100,800])
%             subplot(3,1,1),plot(self.gradientVect(1,:),self.gradientVect(2,:));
%             title('Slope X'),xlabel('Time (seconds)'),ylabel('Gradient Slope (T/m)')
%             subplot(3,1,2),plot(self.gradientVect(1,:),self.gradientVect(3,:));
%             title('Slope Y'),xlabel('Time (seconds)'),ylabel('Gradient Slope (T/m)')
%             subplot(3,1,3),plot(self.gradientVect(1,:),self.gradientVect(4,:));
%             title('Slope Z'),xlabel('Time (seconds)'),ylabel('Gradient Slope (T/m)')
        end
        function compileRF(self)
            self.RFVect = [-1,1e3;0,0];
            for i = 1:numel(self.rfPulses)
                % Store which pulse is on
                % note a 2 pico second buffer is added to each pulse. this will
                % result in an error is a pulse starts at teh exact time one
                % ends
                self.RFVect = [self.RFVect,[self.rfPulses{i}.startTime-1e-12;0]];
                self.RFVect = [self.RFVect,[self.rfPulses{i}.endTime+1e-12;0]];
                % Store 0 for when pulse is off
                self.RFVect = [self.RFVect,[self.rfPulses{i}.startTime;i]];
                self.RFVect = [self.RFVect,[self.rfPulses{i}.endTime;i]];
            end
            [~,I] = sort(self.RFVect(1,:));
            self.RFVect = self.RFVect(:,I);
            for i = 1:2:size(self.RFVect,2)
%                 if self.RFVect(2,i) ~= self.RFVect(2,i+1)
%                     error('RF Pulses %d and %d are overlapping with times %d, %d and %d, %d\n',...
%                         self.RFVect(1,i),self.RFVect(1,i+1),...
%                         self.rfPulses{i}.startTime,self.rfPulses{i}.endTime,...
%                         self.rfPulses{i+1}.startTime,self.rfPulses{i+1}.endTime)
%                 end
            end
        end
        function updateTime(self)
            % UPDATETIME: iterates through the pulse sequence and stores the
            % times that pulses are turned on or off for use by the solver
            allPulses = [self.rfPulses,self.gradientPulses,self.calcGradPulses];
            pointsPerPulse = self.maxSize/length(allPulses);
            self. time = 0;
            self.eventTimes = [];
            for i = 1:length(allPulses)
                self.time = [self.time,linspace(allPulses{i}.startTime,...
                    allPulses{i}.endTime,pointsPerPulse)];
                self.eventTimes(end+1) = allPulses{i}.startTime;
                self.eventTimes(end+1) = allPulses{i}.endTime;
            end
            % stores the on and off times for each pulse as a root to a
            % polynomial
            self.time = sort(self.time);
            self.eventTimes = unique(self.eventTimes);
            self.eventTimes = sort(self.eventTimes);
%             self.eventTimes = [];
%             for i = 1:size(self.RFVect,2)
%                 if(self.RFVect(2,i))
%                     self.eventTimes = [self.eventTimes,[self.RFVect(1,i);1]];
%                 else
%                     self.eventTimes = [self.eventTimes,[self.RFVect(1,i);0]];
%                 end
%             end
%             for i = 1:size(self.gradientVect,2)-1
%                 if(sum(self.gradientVect(2:4,i))~=sum(self.gradientVect(2:4,i+1)))
%                     self.eventTimes = [self.eventTimes,[self.gradientVect(1,i+1);1]];
%                 else
%                     self.eventTimes = [self.eventTimes,[self.gradientVect(1,i+1);0]];
%                 end
%             end
%             [~,I] = sort(self.eventTimes(1,:));
%             self.eventTimes = self.eventTimes(:,I);
%             I = find(diff(self.eventTimes(1,:))<=0);
%             self.eventTimes(2,I+1) = max([self.eventTimes(2,I);self.eventTimes(2,I+1)]);
%             self.eventTimes(:,I) = [];
        end
    end
end

