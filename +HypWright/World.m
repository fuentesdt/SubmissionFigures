classdef (Sealed) World < handle
    %WORLD: Hello World!
    %   storage for global system states, this is a singleton and golbal
    %   PROPERTIES
    %   B0 - The main magnetic field of the scanner
    %   pulseSequence - The MR pulse sequence
    %   Voxels - all of the active voxels
    %   init - logical for weather or not the world has been initiated
    %   calEndTime - the last timepoint for which a solution has been found
    %   METHODS
    %   B0 - The main magnetic field of the scanner, default 3T
    %   PulseSequence - The MR pulse sequence
    %   Voxels - all of the active voxels
    %   init - logical for weather or not the world has been initiated
    %   Methods
    %   setB0(B0) - sets the static B field to B0
    %   setPulseSequence(pulseSequence) - sets the sroted pulse sequence
    %   to the input PulseSequence
    %   initWorld() initializes the world and sets B0 to its default and
    %   initializes an empty pusle sequence
    %   initWorld(B0) - initializes the world with some input B0 that
    %   should be a 3x1 column vector of the form [x;y;z]
    %   addVoxel(voxelList) - adds all the voxels in voxelList to the
    %   world
    %   calculate(times) - claculates the MR signal for all the time points
    %   specifed by (time)
    %   evaluate(time) - returns a M vector for each time point in the
    %   time vector t
    %   *Note M is only defined from 0 to the end time passed in to
    %   calulate, and will only reflect the system state at the last
    %   calculate. This function will sort the time vector and remove any
    %   points outside of the range [0, calEndTime] as well as removing
    %   any redundent time points  rounding to the nearest picosecond
    properties
        parallel = true;
    end
    properties (SetAccess = private)
        B0; % The main magnetic field of the scanner
        pulseSequence;  % The MR pulse sequence
        Voxels; % all of the active voxels
        init; % logical for weather or not the world has been initiated
        calEndTime % the last timepoint for which a solution has been found
        solutions; % a cell of enums that stores the calculated solutions
    end
    properties (Constant)
        voxelSize = 1e-15; % The voxel size, I am not sure if this is used
    end
    
    methods (Access = private)
        function self = World
            % CONSTRUCTOR: Starts the Bloch Simulator. Initializes a pulse
            % sequence 
            self.init = false;
        end
    end
    methods (Static)
        function singleObj = getWorld
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = HypWright.World;
            end
            singleObj = localObj;
        end
    end
    methods
        %% Getters and settes
        function value = getB0(self),value = self.B0; end
        function setB0(self,B0),self.B0 = B0;end
        function value = getPulseSequence(self),value = self.pulseSequence;end
        function setPulseSequence(self,pulseSequence),self.pulseSequence = pulseSequence;end
        function b = getB(self,x,y,z,t)
            % Gets the combined magnetic field from all sources at a
            % position (x,y,z) and a time t.
            b = repmat(self.B0,1,length(t))+self.pulseSequence.B(x,y,z,t);
        end
        %% Other Methods
        function initWorld(self,varargin)
            % INITWORLD: initializes the a new empty world.
            % initWorld() initializes the world and sets B0 to its default and
            % initializes an empty pusle sequence and clears out any voxel
            % initWorld(B0) - initializes the world with some input B0 that
            % should be a 3x1 column vector of the form [x;y;z]
            p = inputParser;
            p.addParameter('B0',[0;0;3.0],@isnumeric);
            p.parse(varargin{:})
            self.B0 = p.Results.B0;
            self.pulseSequence = [];
            self.solutions = [];
            self.calEndTime = 0;
            self.clearVoxels();
            self.init = true;
        end
        function addVoxel(self,voxelList)
            % ADDVOXEL -  adds a voxel to the world
            % addVoxel(voxelList) - adds all the voxels in voxelList to the
            % world. They are added in the order they were passed in.
            % Currently there is no great for managing and manipulating
            % multiple voxels, this should probably be addressed if mor
            % complicated voxel geometries are to be used, probably a
            % factory object.
            for i = 1:numel(voxelList)
                self.Voxels = [self.Voxels,voxelList(i)];
            end
        end
        function clearVoxels(self)
            % CLEARVOXELS - removes all voxels from the world
            self.Voxels = [];
        end
        function FID = solve(self,times,varargin)
            %% Numerically solves the simulation world for a given set of
            % time inputs and then returns world to the base state
            
            %Parese Input
            p = inputParser;
            p.addOptional('ref',0,@isnumeric)
            p.parse(varargin{:})
            ref= p.Results.ref;
            % reshape the times matrix into a monotonically increasing vector 
            [sortedTimes, isortedTimes] = sort(times(:));
            [~, iUndoSort] = sort(isortedTimes); % get indecies to put Time and FID back in the requested order

            % Get the relevent pulse sequence event times for the time
            % frame specified by times. Assumes times is monotonically
            % increasing with times and starting at or before t=0
            eventTimes = self.pulseSequence.getEventTimes(sortedTimes(end));
            % Split the requested eval times (times) into their own
            % evalTime for each time frame
            evalTimes = cell(numel(eventTimes)-1,1);
            for i = 1:(numel(eventTimes)-2)
                evalTimes{i} = sortedTimes(((sortedTimes>=eventTimes(i))&(sortedTimes<eventTimes(i+1))));
            end
            evalTimes{end} = sortedTimes(sortedTimes>=eventTimes(end-1));
            
            %% Evaluate the FID for each Voxel
            FID = zeros(1,length(sortedTimes));
            % Solve for each Voxel in parralell if requested
            if self.parallel
                % Parallel processing without tmp variables. I think this
                % is slower
%                 parfor i = 1:numel(self.Voxels)
%                     tmpFID = self.Voxels(i).solve(eventTimes,evalTimes,self.pulseSequence,self.B0,ref);
%                     tmpFID = tmpFID(~cellfun('isempty',tmpFID)).';
%                     FID = FID+cell2mat(tmpFID);
%                 end
                %temporary variables to be passed to the workers in the
                %parfor
                tmpVoxels = self.Voxels;
                tmpPS = self.pulseSequence;
                tmpB0 = self.B0;
                parfor i = 1:numel(tmpVoxels)
                    tmpFID = tmpVoxels(i).solve(eventTimes,evalTimes,tmpPS,tmpB0,ref);
                    tmpFID = tmpFID(~cellfun('isempty',tmpFID)).';
                    FID = FID+cell2mat(tmpFID);
                end
            else
                for i = 1:numel(self.Voxels)
                    % Solve for the FID at each pulse sequence even
                    tmpFID = self.Voxels(i).solve(eventTimes,evalTimes,self.pulseSequence,self.B0,ref);
                    % strip the pulse sequence events with no concurent evaluation times
                    tmpFID = tmpFID(~cellfun('isempty',tmpFID)).'; 
                    FID = FID+cell2mat(tmpFID); % Repack the FID into a single vector
                end
            end
            % Reshape fid to match the passed in evaluation times
            FID = reshape(FID(iUndoSort),size(times));
        end
        function calculate(self,timeRange)
            % CALCULATE - calculates the MR signal from all voxels over
            % some time range, assumes a start time of zeros if only one
            % number is passed in
            
            % Compiles the pulse sequence for efficency. Need to add a
            % check to make sure the sequence has not already been
            % compiled.
            if ~self.pulseSequence.compiled
                self.pulseSequence.compile();
            end
            %             % Let each voxl calculate it's own solution,
            tmpVoxels = self.Voxels;
            tmpEndTime = timeRange(end);
            tmpPS = self.pulseSequence;
            tmpB0 = self.B0;
            if self.parallel
                parfor i = 1:numel(self.Voxels)
                    tmpVoxel = tmpVoxels(i);
                    tmpVoxel.calculate(tmpEndTime,tmpPS,tmpB0);
                    tmpVoxels(i) = tmpVoxel;
                end
                self.Voxels = tmpVoxels;
                %Sequential Code
            else
                for i = 1:numel(self.Voxels)
                    self.solutions{i} = self.Voxels(i).calculate(timeRange(end),...
                        self.pulseSequence,self.B0);
                end
            end
            %self.solutions = tmpSolutions; % store the solutions an replace any old solutions
            % stores the end time of the calculated range
            self.calEndTime = timeRange(end);
        end
        function [signal, freqAxis, timeAxis,M] = evaluate(self,times,varargin)
            % EVALUATE - returns the complext MR signal for the time points
            % passed in. The world needs to be calculated before it can be
            % evaluated.
            % [signal, freqAxis, timeAxis] = evaluate(times) -
            % evaluate(time) - returns a M vector for each time point in the
            % time vector t
            % *Note M is only defined from 0 to the end time passed in to
            % calulate, and will only reflect the system state at the last
            % calculate. This function will sort the time vector and remove any
            % points outside of the range [0, calEndTime] as well as removing
            % any redundent time points  rounding to the nearest picosecond
            %% Check to make sure time is passed in as a row vector
            if(size(times,1)~=1)
                times = times.';
            end
            times = unique(round(times,1e12)); % Remove anny "duplicate times"
            times = sort(times);
            p = inputParser;
            p.addOptional('ref',0,@isnumeric)
            p.addOptional('verbose',0,@islogical)
            p.parse(varargin{:})
            tmpM = zeros(numel(self.Voxels),3,length(times));
            for i = 1:numel(self.Voxels)
                tmpM(i,:,:) = self.Voxels(i).getM(times,p.Results.ref,self.B0);
            end
%             %% Not sure if parallel evaluation helps
%             tmpB0 = self.B0;
%             tmpRef = p.Results.ref;
%             tmpVoxels = self.Voxels;
%             if self.parallel
%                 parfor i = 1:numel(self.Voxels)
%                     tmpM(i,:,:) = tmpVoxels(i).getM(times,tmpRef,tmpB0);
%                 end
%             else
%                 for i = 1:numel(self.Voxels)
%                     tmpM(i,:,:) = tmpVoxels(i).getM(times,tmpRef,tmpB0);
%                 end
%             end
            M = squeeze(sum(tmpM,1));
            if size(M,1) == 1
                signal = M(1)+1i*M(2);
                freqAxis = [];
            else
                signal = M(1,:)+1i*M(2,:);
                BW = 1/(times(2)-times(1));
                freqAxis = linspace(-BW/2,BW/2,length(times));
            end
            signal = signal.';
            
            timeAxis = times;
            if (p.Results.verbose)
                figure
                subplot(2,1,1),plot(times,real(signal),'r',times,imag(signal),...
                    'b')
                xlabel('Time (seconds)')
                FTSig = fftshift(fft(fftshift(signal)));
                subplot(2,1,2),plot(freqAxis,real(FTSig),'r',freqAxis,...
                    imag(FTSig),'b',freqAxis,abs(FTSig),'k')
                xlabel('Frequency (HZ)')
            end
        end
        function [M,time] = showMagnetization(self,end_time)
            for i = 1:numel(self.Voxels)
                [tmp_M,tmp_time] = self.Voxels(i).showM(end_time,self.pulseSequence,self.B0);
                M{i} = tmp_M;
                time{i} = tmp_time;
            end
        end
    end
    
end