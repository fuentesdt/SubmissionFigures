classdef DataClass < handle & matlab.mixin.Copyable
    %DATACLASS This is the base class for hyperpolarized data synthetic or
    %otherwise
    %   As the core data storage class for data processing It allows for
    %   minor general manipulation, display and visualization
    %   Private Properties Accessable by Setters
    %       hypPyruvate - The hyperpolarized Pyruvate "signal" (edit with care)
    %       hypLactate -- The hyperpolarized Lactate "signal" (edit with care)
    %       tax -- The time axes associated with the hyperpolarized signal 
    %       (edit with care)
    %   Public Methods    
    %       display: Generates a plot of the data ine ether a new figure or the
    %       axis passed in
    %       getData: returns a structure containg the trimmed normalized data.
    %       normalize: normalizes the data set to some given factor 
    %       ratio: computes either normalized lactat (nLac) or the lactate to
    %       pyruvate ratio of either trimmed or un-trimmed data
    %       reset: Removes any normalization or trimming from the data
    %       trim: trims the data in some specified way 
    
    
    %% Info
    % Author: Chris Walker
    % Creation Data:6/5/2014
    % Last Updata: 6/6/2014 CMW
    %% Properties
    properties
    end
    properties (Access = protected)
        hypPyruvate;                    %Polarized Pyruvate Signal
        hypLactate;                     %Polarized Lactate Signal
        tax;                            %total time axis
        normFact = 1;                   %Normilization factor
        trimTime = 1;                   %Some cutoff time
    end
    
    %% Methods
    methods
        %% Getters and Setters
        function hypPyr = getHypPyruvate(self), hypPyr = self.hypPyruvate(); end
        function hypLac = getHypLactate(self), hypLac = self.hypPyruvate(); end
        function normFact = getNormFact(self), normFact = self.normFact; end
        function tax = getTax(self), tax = self.tax; end
        function setHypPyruvate(self,hypPyr),self.hypPyruvate = hypPyr(); end
        function setHypLactate(self,hypLac), self.hypLactate = hypLac(); end
        function setTax(self,tax), self.tax = tax; end
        %% Other Methods
        function display(self,varargin)
            % DISPLAY Displays timecourse pyruvat and lactate signal
            %   () displays data in a new figure
            %   (A) Displays the data into the passed in axis A
            p = inputParser();
            p.addOptional('axis',[])
            p.parse(varargin{:})
            self.dataCheck;
            if(isempty(p.Results.axis))
                figure;
                curAxis = gca;
            else
                curAxis = p.Results.axis;
            end
            plot(curAxis, self.tax(self.trimTime:end),...
                self.hypPyruvate(self.trimTime:end)./self.normFact,'g',...
                self.tax(self.trimTime:end),...
                self.hypLactate(self.trimTime:end)./self.normFact,'b');
            title('Polarized Pyruvate and Lactate');
            xlabel('Time (seconds)');
            ylabel('Concentration (moles)');
            legend('Pyruvate','Lactate');
        end
        function data = getData(self)
            % GETDATA returns a structure containing the normalized trimmed data.
            % The structure (result) has 3 values, result.hypPyruvate,
            % result.hypLactate and result.tax
            self.dataCheck();
            data.hypPyruvate = self.hypPyruvate(self.trimTime:end)...
                ./self.normFact;
            data.hypLactate = self.hypLactate(self.trimTime:end)./self.normFact;
            data.tax = self.tax(self.trimTime:end);
        end
        function normalize(self,normFact)
            % NORMALIZE Normalizes the data set to some given factor
            % (normFact) normalizes by a factor normFact (probably
            % should be posotive)
            %
            % see also reset
            self.dataCheck;
            self.normFact = normFact;
        end
        function ratioVal = ratio(self,varargin)
            % RATIO computse either nLac (default) or L/P
            % (specify) computes nLac if specify = 'nLac' or L/P is specify =
            % 'L/P'
            % (specify,clean) the same specify argument as above if clean =
            % 'clean' any data trimming will be removed for the ratio
            % calculation
            p = inputParser();
            p.addOptional('ratio','nLac',@isstr);
            p.addOptional('clean',' ',@isstr);
            p.parse(varargin{:});
            if strcmpi('clean',p.Results.clean)
                ratioIndexVect = 1:length(self.tax);
            else
                ratioIndexVect = self.trimTime:length(self.tax);
            end
            switch(lower(p.Results.ratio))
                case 'l/p'
                    ratioVal = sum(self.hypLactate(ratioIndexVect))/...
                        sum(self.hypPyruvate(ratioIndexVect));
                case 'nlac'
                    ratioVal = sum(self.hypLactate(ratioIndexVect))/...
                        (sum(self.hypLactate(ratioIndexVect))...
                        +sum(self.hypPyruvate(ratioIndexVect)));
                otherwise
                    ratioVal = nan;
                    disp('ERROR! inproper input argument!')
            end
        end
        function reset(self,varargin)
            % RESET removes any normalization or trimming from the data
            % (Specify) pass a string 'trimtime' or ' normfact' to only remove
            % those correction
            % 
            % see also normalize and trim
            p = inputParser();
            p.addOptional('specified',' ',@isstr)
            p.parse(varargin{:})
            specified = lower(p.Results.specified);
            switch specified
                case('trimtime')
                    self.trimTime = 1;
                case('normfact')
                    self.normFact = 1;
                otherwise
                    self.normFact = 1;
                    self.trimTime = 1;
            end
        end
        function trim(self,varargin)
            % TRIM Trims the data
            % () Defaults to prompt user for the trim time
            % (Method) Takes a string (Method) to determin the method of
            % inputting the trim time. 'max' trims to the max pyruvate signal;
            % 'start' trims to 5 data points before the make pyruvate, 'select'
            % Displays the data and propts the user to inpit the trim time and
            % 'inp' takes a secondary input see below
            % ('inp' trimTime) sets the trims the data to the highest value in
            % tax below the specified trimeTime, if trim time is not specified
            % no change is made
            %
            % see also reset
            p = inputParser();
            p.addOptional('method','select',@isstr)
            p.addOptional('trimTime',self.trimTime,@isscalar);
            p.parse(varargin{:});
            method = lower(p.Results.method); %Strip caps
            if (~strcmp(p.Results.method,'inp') && nargin > 2)
                disp('Warning! a trim time was specifed but not used!')
            end
            if (strcmp(p.Results.method,'inp') && nargin < 3)
                disp(['Warning! Input trim time was specified but not given.'...
                    ' No changes being made'])
            end
            switch(method)
                case('max')
                    [~, I] = max(self.hypPyruvate);
                    self.trimTime = I;
                case('start')
                    [~, I] = max(self.hypPyruvate);
                    self.trimTime = I-5;
                    if(self.trimTime<0)
                        self.trimTime = 0;
                    end
                case('select')
                    self.Display()
                    h = gcf;
                    selTime = input('Enter the Time you want to trim to: ');
                    self.trimTime = find(self.tax<selTime,1,'last');
                    close(h);
                case('inp')
                    % Check user input for Trim Time.
                    if ( 1 <= p.Results.trimTime &&...
                            p.Results.trimTime <= self.tax(end))
                        self.trimTime = find(self.tax<p.Results.trimTime,1,'last');
                    elseif (p.Results.trimTime < 1), self.trimTime = 1;
                    else
                        self.trimTime = length(self.tax);
                    end
            end
        end
    end
    
    methods (Access = protected)
        function checkBool = dataCheck(self)
            % DATACHECK - Checks the integrity of the data, will return and
            % error if data is corrupted (lengths not matching, trim time out of
            % bounds, etc)
            initPyr = isempty(self.hypPyruvate);
            initLac = isempty(self.hypLactate);
            initTax = isempty(self.tax);
            lengthCheck = length(self.hypPyruvate) == length(self.hypLactate)...
                && length(self.hypPyruvate) == length(self.tax);
            trimTimeCheck = 0 < self.trimTime &&...
                self.trimTime <= length(self.tax);
            checkBool = ~initPyr && ~initLac && ~initTax && lengthCheck &&...
                trimTimeCheck;
            if(~checkBool)
                error(['The data is corrupted! The trim time is outside of'...
                    'the range of the data or the vectors are not the'...
                    'same length']);
            end
        end
    end
    
end

%% Methods from older version. Kept here for refrences if similar methods need
% to be added
%         %% Compare: Calculates the error between data and
%         function Compare (self, model,params,inpParam)
%             inpParam.normFact = self.normFact;
%             model.Compare(self.GetData,inpParam,params)
%         end 
%         %% Resample: Resamples data based on given dt
%         function Resample(self,dt)
%             if(not(self.DataCheck))
%                 error('No Data defined')
%             end
%             newTax = linspace(self.tax(1),self.tax(end),...
%                 (self.tax(end)-self.tax(1))/dt);
%             self.hypPyruvate = interp1(self.tax,self.hypPyruvate,newTax);
%             self.hypLactate = interp1(self.tax,self.hypLactate,newTax);
%             self.tax = newTax;
%         end
%         %% FlipangleCorrect: Corrects for RF excitation
%         function FlipangleCorrect(self,alpha)
%             if(not(self.DataCheck()))
%                 error('not data to flip angle correct')
%             end
%             for i = 1:length(self.hypPyruvate)
%                 self.hypPyruvate(i) = self.hypPyruvate(i)/(cos(pi*alpha/180)^(i-1));
%                 self.hypLactate(i) = self.hypLactate(i)/(cos(pi*alpha/180)^(i-1));
%                 
%             end
%         end
%         %% RFLoss: applies RF losses to data
%         function RFLoss(self,alpha)
%             if(not(self.DataCheck()))
%                 error('not data to flip angle correct')
%             end
%             for i = 1:length(self.hypPyruvate)
%                 self.hypPyruvate(i) = self.hypPyruvate(i)*(cos(pi*alpha/180)^(i-1));
%                 self.hypLactate(i) = self.hypLactate(i)*(cos(pi*alpha/180)^(i-1));
%                 
%             end
%         end

