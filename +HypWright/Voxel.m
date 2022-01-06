classdef Voxel < handle
    %VOXEL Represents a volume of space
    %   A voxel represents a voulme of space that contains some set of spins.
    %   once the calculate method has been run the voxel stores a solution
    %   describing the evolution of the total magnetization vector up to some
    %   time. this can be evaluated with the getM method.
    %   Properties
    %   position - Vector defining the position of the Voxel
    %   Methods
    %   Voxel(position) - initializes an empty voxel at the coordinates
    %   defined by position
    %   Voxel(position, spinList)- initializes a Voxel at the coordinates
    %   defined by position and fills it with any spin groups in the
    %   optional variable spinList
    %   addSpin(spin) - adds the input spin to the voxel
    %   calculate(endTime) -  runs the solvers to the specified end time
    %   getM(t) - returns the magnetization vector for all time points in
    %   the vector t
    
    properties
        position % Vector defining the position of the Voxel
        sol = {}; % list of all the solution structurs defing the time evolution of M
        anlyticSol = {}; %functions for all the analytic solutions
        debug = false; %switches on debug mode (will save the Mz of the spin)
        solTimes % the time ranges each solution structure spans
    end
    properties (SetAccess = private)
        solution % cell that stores all the data needed to solve for this Voxel
    end
    properties (Access = private)
        spinGroups % list of all spins in the voxel
    end
    properties (Constant)
        T2Star = 2e-14; % determins B0 inhomogenaety thus T2 star in thsis voxel
        numSubSpins = 1^3; % Defines the number of spin groups in this voxel
    end
    
    methods
        function self = Voxel(position,varargin)
            % CONSTRUCTOR - initializes the voxel at some position
            % Voxel(position) - initializes an empty voxel at the coordinates
            % defined by position
            % Voxel(position, spinList)- initializes a Voxel at the coordinates
            % defined by position and fills it with any spin groups in the
            % optional variable spinList
            p = inputParser();
            p.addOptional('spinList',[])
            p.parse(varargin{:})
            self.position = position;
            for i = 1:numel(p.Results.spinList)
                self.spinGroups = {p.Results.spinList(i)};
            end
        end
        function addSpin(self,spin)
            % ADDSPIN - adds spins to the Voxel
            % addSpin(spin) - adds the input spin to the voxel
            self.spinGroups{end+1} = spin;
        end
        function retunFID = solve(self,eventTimes,evalTimes,pulseSequence,B0,ref)
            %% Evaluates the FID for all spins in a single Voxel
            retunFID = cell(size(evalTimes)); % innitialize the return variable
            % Get the voxels' position
            x = self.position(1);
            y = self.position(2);
            z = self.position(3);
            % loop over all spins in the voxel
            for j = 1:numel(self.spinGroups)
                % Set the magentization equal to the spins inital condition
                tmpM = self.spinGroups{j}.M;
                % initialize a cell to hold the spin's mangetization
                % evolution
                Mframe = cell(size(evalTimes));
                % Loop over all the pulse sequence events
                for i = 1:(length(eventTimes)-1)
                    % Fill the times to be solved, the start and end of the
                    % event as well as any evaluation points
                    tSpan = [eventTimes(i),evalTimes{i}.',eventTimes(i+1)];
                    % Check if the pulse sequence requires a numeric
                    % solution
                    ODEBool = pulseSequence.solver(tSpan([1,end]));
                    if(ODEBool || ~self.spinGroups{j}.useAnalytical())
                        % set up dm/dt for ode45
                        odefun = @(M,t)self.spinGroups{j}.dM(x,y,z,M,t,pulseSequence,B0);
                        tmpSol = ode45(odefun,[eventTimes(i),eventTimes(i+1)],tmpM);
                        tmpM = deval(tmpSol,tSpan(end));
                        if length(tSpan) > 2
                            M = deval(tmpSol,evalTimes{i}.').';
                            % If there are evaluation times, store the
                            % values in MFrame
                            if isempty(Mframe{i})
                                Mframe{i} = M.';
                            else
                                Mframe{i} = Mframe{i}+M.';
                            end
                        end
                    else
                        % if no numeric solution is required use an
                        % analytic solution
                        t = mean(tSpan([1,end])); % get the center time of the event
                        B = B0+pulseSequence.B(x,y,z,t); % get the B field. This alue should nto be changing over time when the analytic solution is used.
                        % Build the analytic expression
                        tmpAnlyticSol = @(t)self.spinGroups{j}.analytical(...
                            x,y,z,tSpan(1),tmpM,t,pulseSequence,B0,B);
                        tmpM = tmpAnlyticSol(tSpan(end));% Store the final magnetization for the next pulse sequence event
                        if length(tSpan) > 2
                            % If there are evaluation times, store the
                            % values
                            M = tmpAnlyticSol(tSpan(2:(end-1)));
                            if isempty(Mframe{i})
                                Mframe{i} = M;
                            else
                                Mframe{i} = Mframe{i}+M;
                            end
                        end
                    end
                    % Combine multiple magnetizations from coupled spins into a single signal
                    tmpSig = zeros(1,numel(tSpan(2:(end-1))));
                    tmpMz = zeros(1,numel(tSpan(2:(end-1))));
                    for k = 1:3:size(Mframe{i},1)
                        tmpSig = tmpSig+Mframe{i}(k,:)+1i*Mframe{i}(k+1,:);
                        tmpMz = tmpMz+Mframe{i}(k+2,:);
                    end
                    % Mix signal to referenec frequency
                    omegaShift = self.spinGroups{j}.calculationFrame(B0)+ref;
                    framePhase = exp(1i*omegaShift*tSpan(2:(end-1)));
                    tmpSig = tmpSig.*framePhase;
                    % Store the final FID
                    if isempty(retunFID{i})
                        retunFID{i} = self.spinGroups{j}.density*tmpSig;
                    else
                        retunFID{i} = retunFID{i}+self.spinGroups{j}.density*tmpSig;
                    end
                end
            end
        end
        function solution = calculate(self,endTime,PS,B0)
            % CALCULATE - calculates the time dependent magnetization vector for
            % this voxel
            % calculate(endTime) -  runs the solvers to the specified end time
            %% get pulse times and clean up any overlap or negatives and sort them
            pulseTimes = PS.eventTimes(1,:); % Times when B changes in the Pulse sequence
            pulseTimes = pulseTimes(find(pulseTimes>0,1,'first'):end);
            pulseTimes = pulseTimes(pulseTimes<=endTime);
            pulseTimes = [0,pulseTimes,endTime];
            pulseTimes = sort(pulseTimes);
            pulseTimes = unique(round(pulseTimes.*1e9))./1e9;
            pulseTimes = unique(pulseTimes);
            self.solTimes = pulseTimes;
            solution.solTimes = pulseTimes;
            % initialize the storage for the solutions
            self.anlyticSol = cell(length(pulseTimes)-1,numel(self.spinGroups));
            self.sol = cell(length(pulseTimes)-1,numel(self.spinGroups));
            % generate and stor the solutions for each pulse time
            for i = 1:length(pulseTimes)-1
                tSpan = [pulseTimes(i),pulseTimes(i+1)];
                x = self.position(1);
                y = self.position(2);
                z = self.position(3);
                t = mean(tSpan);
                tmpSpinGroups = self.spinGroups;
                if(i>1)
                    tmpAnlyticSol2 = tmpAnlyticSol;
                    tmpSol2 = tmpSol;
                else
                    tmpAnlyticSol2 = self.anlyticSol(i,:);
                    tmpSol2 = self.sol(i,:);
                end
                tmpAnlyticSol = self.anlyticSol(i,:);
                tmpSol = self.sol(i,:);
                ODEBool = PS.solver(tSpan);
                %tic
                for j = 1:numel(tmpSpinGroups)
                    % callculate current M
                    if i == 1,
                        tmpM = tmpSpinGroups{j}.M;
                    else
                        if isempty(tmpSol2{j})
                            tmpFun = tmpAnlyticSol2{j};
                            tmpM = tmpFun(tSpan(1));
                        else
                            tmpFun = tmpSol2{j};
                            tmpM = deval(tmpFun,tSpan(1));
                        end
                    end
                    if(ODEBool || ~tmpSpinGroups{j}.useAnalytical())
                        % used ODE solver when PS is changing
                        odefun = @(M,t)tmpSpinGroups{j}.dM(x,y,z,M,t,PS,B0);
%                       figure
%                        ode45(odefun,tSpan,tmpM)
%                        pause(waitforbuttonpress)
%                         options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Refine',6);
%                         tmpSol{j} = ode45(odefun,tSpan,tmpM,options);
                        tmpSol{j} = ode45(odefun,tSpan,tmpM);
                        solution.functions(i,:) = tmpSol;
                        solution.useAnalytical(i,:) = false;
                    else
                        B = B0+PS.B(x,y,z,t);
                        tmpAnlyticSol{j} = ...
                            @(t)tmpSpinGroups{j}.analytical(...
                            x,y,z,tSpan(1),tmpM,t,PS,B0,B);
                        tmpSol{j} = {};
                        solution.functions(i,:) = tmpAnlyticSol;
                        solution.useAnalytical(i,:) = true;
%                         odefun = @(M,t)tmpSpinGroups{j}.dM(x,y,z,M,t,PS,B0);
%                       figure
%                       ode45(odefun,tSpan,tmpM)
%                       pause(waitforbuttonpress)
                        %tmpsol{i,j} = ode45(odefun,tspan,tmpM);
                    end
                end
                for j = 1:numel(tmpSpinGroups)
                    solution.frameFreq{j} = @(B0)...
                        self.spinGroups{j}.calculationFrame(B0);
                    solution.spinDensity{j} = self.spinGroups{j}.density;
                end
                self.anlyticSol(i,:) = tmpAnlyticSol;
                self.sol(i,:) = tmpSol;
                solution.anlyticSol(i,:) = tmpAnlyticSol;
                solution.sol(i,:) = tmpSol;
                %toc
            end
        end
        function M = getM(self,t,ref,B0)
            % GETM -  retuns the magnetization vector over some time vector
            % getM(t) - returns the magnetization vector for all time points in
            % the vector t
            M = zeros(3,length(t));
            % Break up the time vector into chunks that mach the diffrent
            % solutions
            start = find(t(1)<self.solTimes,1,'first');
            devalTimes = 1;
            if self.debug
            tmp = {};
            save('tmp','tmp')
            end
            for j = start:numel(self.solTimes)
                devalTimes(end+1) = find(t<self.solTimes(j),1,'last');
                tmp = find(t>=self.solTimes(j),1,'first');
                if isempty(tmp)
                    break;
                else
                    devalTimes(end+1) = tmp;
                end
            end
            devalTimes(end+1) = length(t);
            devalTimes = sort(devalTimes);
            for n = 1:2:numel(devalTimes)-1
                tDeval = t(unique(devalTimes(n):devalTimes(n+1)));
                % find wich solution to use for this time
                j = find(self.solTimes>tDeval(1),1,'first');
                if isempty(j), j = length(self.solTimes); end
                j = j-1;
                tmpM = zeros(3,numel(tDeval));
                for i = 1:numel(self.spinGroups)
                    % calculate M0 for a spin group at the passed in time
                    if isempty(self.sol{j,i})
                        % Calculat with analytical
                        %tmp = self.anlyticSol{j,i}(tDeval(1));
                        Mframe = self.anlyticSol{j,i}(tDeval);
                        if self.debug
                        load('tmp.mat')
                        tmp{j,i} = Mframe;
                        save('tmp','tmp')
                        end
                    else
                        % calculate woth ode
                        Mframe = deval(self.sol{j,i},tDeval);
                        if self.debug
                        load('tmp.mat')
                        tmp{j,i} = Mframe;
                        save('tmp','tmp')
                        end
                    end
                    omegaShift = self.spinGroups{i}.calculationFrame(B0)+ref;
                    framePhase = exp(1i*omegaShift*tDeval);
                    tmpSig = zeros(1,numel(tDeval));
                    tmpMz = zeros(1,numel(tDeval));
                    for k = 1:3:size(Mframe,1)
                        tmpSig = tmpSig+Mframe(k,:)+1i*Mframe(k+1,:);
                        tmpMz = tmpMz+Mframe(k+2,:);
                    end
                    tmpSig = tmpSig.*framePhase;
                    tmpM = tmpM+self.spinGroups{i}.density*[real(tmpSig);imag(tmpSig);tmpMz];
%% Old Method for rotating frame correction It was slow and replaced with the aove method                    
%                     for p = 1:numel(tDeval)
%                         % rotate to refrence frame
%                         theta = (self.spinGroups{i}.calculationFrame(B0)+ref)*...
%                             tDeval(p);
%                         %loop over all the spins in a group
%                         for k = 1:3:size(Mframe,1)
%                             tmpM(:,p) = tmpM(:,p)+self.spinGroups{i}.density*...
%                                 [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;...
%                                 0,0,1]*Mframe(k:k+2,p);
%                         end
%                     end
                end
                M(:,unique(devalTimes(n):devalTimes(n+1))) =...
                    tmpM./self.numSubSpins;
            end
        end
        function [M,t] = showM(self,endTime,PS,B0)
            for j = 1:numel(self.spinGroups)
            self.spinGroups = self.spinGroups;
            x = self.position(1);
            y = self.position(2);
            z = self.position(3);
            odefun = @(M,t)self.spinGroups{j}.dM(x,y,z,M,t,PS,B0);
            options = odeset('MaxStep',1e-3);
            [t{j},M{j}] = ode45(odefun,[0,endTime],self.spinGroups{j}.M,options);
%             t{j} = tmp_t;
%             M{j} = tmp_M;
            end
            
        end
    end
end

