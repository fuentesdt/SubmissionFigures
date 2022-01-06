classdef ScanClass < HypWright.DataClass
    %SCANCLASS This is the data strucure for data taken off a Bruker scanner
    %   Supports the storage and relavant manipulation of data taken off a
    %   Bruker scanner. Inherits from DataClass
    %   Private Properties Accessable by Setters
    %       scanPath - That directory that the raw data is stored in
    %   Protected set acces properties
    %       pyrSignalLor - Pyruvate time corse, Lorentzian method
    %       lacSignalLor - Lactate time corse, Lorentzian method
    %       pyrSignalFWHM - Pyruvate time corse, FWHM method
    %       lacSignalFWHM - Lactate time corse, FWHM method
    %   Public Methods
    %       process - takes the data stored in the data path and process it both
    %       using the Lorentzian fit method and integrating each peak (defined
    %       by the FWHM of the summed spectrum). the FWHM method is what is
    %       stored by default
    %   see also DataClass
    
    properties
        verbose = false                 % Boolean controling verbose processing
    end
    properties (Access = protected)
        scanPath                        % The path to some data
        rawData                         % The raw data stored at at scanPath
        dataLB                          % Data after line broadening
        % HEADER - Enumeration to store header information
        % nLines - number of data lines aquired (same as npe)
        % lineLength - number of points per data line (same as nop)
        % spectralWidth - spectral width of each line (same as sw)
        % See Also readHeader, readSpflash
        header                          
        fidTax                          % Time axis of the FID
        fullTax                         % The time axis for the whole scan
        ppmAx                           % axis in ppm
        % SPECTRUMPEAKS - Enumeration to store peak information
        % names - the name of each peak ('pyr','lac')
        % center - the center index of each peak in the stored spectrum
        % params - the lorentizian parameters of each peak
        % See Also findPeak, fillPeaks
        spectrumPeaks
    end
    properties (SetAccess = protected)
        pyrSignalLor                    % Pyruvate time corse, Lorentzian method
        lacSignalLor                    % Lactate time corse, Lorentzian method
        pyrSignalFWHM                   % Pyruvate time corse, FWHM method
        lacSignalFWHM                   % Lactate time corse, FWHM method
    end
    properties (Dependent)
        spectralData                    % Fourier transformed data
        sumSpectrum                     % Sum of spectral data over time
    end
    
    methods
        %% Getters and Setters
        function ret = getScanPath(self), ret = self.scanPath; end
        function setScanPath(self, varargin)
            % SETSCANPATH: Function for changing the scanPath. This will wipe
            % any data stored in the object as it is no longer associated with
            % the stored scanPath
            % setScanPath()  - Uses GUI to set scanPath
            % setScanPath(path) - sets the scanPath to whatever sring is passed
            % in
            
            function setPath(self,inp)
                % SETPATH: helper function for setting the path
                if isempty(inp)
                    tmpPath = uigetdir('','Select Data folder');
                    if (tmpPath == 0)
                        disp(['Warning no path was set! No changed '...
                            'were made to this object.']);
                        return
                    else
                        self.clearScan()
                        self.scanPath = tmpPath;
                    end
                else
                    self.clearScan()
                    self.scanPath = inp;
                end
            end
            
            p = inputParser();
            p.addOptional('inpPath','',@isstr);
            p.parse(varargin{:});
            if ~isempty(self.header)
                inp = input(['Changing the scanPath will clear any data '...
                    'stored in this object. \nAre you sure you want to do '...
                    'this (Y/N)?\n'],'s');
                switch(lower(inp))
                    case 'y'
                        setPath(self,p.Results.inpPath)
                        disp('scanPath changed. All other data has been cleared')
                    case 'n'
                        disp('No changes were made')
                    otherwise
                        disp('Input not recognized. No Changes were madde.')
                end
            else
                disp(['Warning! This Scan is not fully initilized assuming '...
                    'its safe to go ahead and clear any data stored.']);
                setPath(self,p.Results.inpPath)
            end
        end
        function ret = getHeader(self), ret = self.header; end
        function ret = get.spectralData(self) 
            ret = fftshift(fft(self.dataLB.'),1).';
        end
        function ret = get.sumSpectrum(self)
            ret = sum(self.spectralData,1);
        end
        %% Other Methods
        function self = ScanClass(varargin)
            % CONSTRUCTOR: Creats an instance of ScanClass
            % ScanClass(scanPath) sets scanPath to path
            p = inputParser();
            p.addOptional('inpPath','',@isstr);
            p.parse(varargin{:});
            self.scanPath =  p.Results.inpPath;
            self.header = struct([]);
        end
        function process(self,varargin)
            % PROCESS: Takes the data from some scan folder in the processes it
            % to metabolite curves.
            
            %% process input
            p = inputParser();
            p.addOptional('names',[],@iscellstr)
            p.addOptional('centers',[],@isnumeric)
            p.parse(varargin{:})
            if length(p.Results.names) ~= length(p.Results.centers)
                error('Input names and centers not matching')
            end
            %% Get data path
            if isempty(self.scanPath)
                self.getScanPath()
            end
            %% Read data and header
            headerFilePath = sprintf('%s/acqp',self.scanPath);
            fidFilePath = sprintf('%s/fid',self.scanPath);
            if isempty(self.header)
                self.readHeader(headerFilePath);
                self.readSpflash(fidFilePath);
            end
            %% Filter data and compute axis
            self.fidTax = (0:self.header.lineLength-1)./...
                (self.header.spectralWidth);
            % FIXME hardCoding TR will need to pull this from the header
            self.fullTax = (0:self.header.nLines-1)*2;
            self.dataLB = HypWright.ScanClass.lineBroadening(self.rawData,...
                self.fidTax,15);
            %% Get Peaks
            if isempty(p.Results.names)
                error('Peak selector not yet functional')
                % TODO: fix the output of the peak selector to be a struct
            else
                self.spectrumPeaks =  HypWright.ScanClass.fillPeaks(...
                    p.Results.names, p.Results.centers);
            end
            %% Define ppm axis
            iPyr = HypWright.ScanClass.findPeak(self.spectrumPeaks,'pyr');
            pyruvateCenter = self.spectrumPeaks(iPyr).center;
            self.ppmAx = ((1:self.header.lineLength)-pyruvateCenter)*...
                (self.header.spectralWidth/75.53140/self.header.lineLength)+171;
            %% Display for testing
            if(self.verbose)
                figure
                subplot(2,2,1),imagesc(self.fidTax,self.fullTax,...
                    abs(self.rawData))
                title('Raw Data'),xlabel('Time(sec)'),ylabel('Time(sec)')
                subplot(2,2,2),imagesc(self.fidTax,self.fullTax,...
                    abs(self.dataLB))
                title('Linebroadened Data'),xlabel('Time(sec)')
                ylabel('Time(sec)')
                subplot(2,2,3),imagesc(self.ppmAx,self.fullTax,...
                    abs(self.spectralData))
                title('Spectral Data'),xlabel('ppm'),ylabel('Time(sec)')
                subplot(2,2,4),plot(self.ppmAx,abs(self.sumSpectrum))
                title('Sum of Specs'),xlabel('ppm'),ylabel('Signal')
            end
            %% Fit peaks
            self.fitLorentzian()
            %% Integrate peaks
            self.integrateFWHM()
            self.hypPyruvate = self.pyrSignalFWHM;
            self.hypLactate = self.lacSignalFWHM;
            self.tax = self.fullTax;
        end
    end
    methods (Access = protected)
        function clearScan(self)
            % clearScan: returns the class to a clean initilization
            self.scanPath = '';
            self.hypPyruvate = [];
            self.hypLactate = [];
            self.tax = [];
            self.header = struct([]);
            self.rawData = [];
            self.reset();
        end
        function fitLorentzian(self)
            % FITLORENTZIAN: fits the current raw data set of specral data to a
            % series of Lorentzian distributed peaks
            spectralData = self.spectralData;
            sumSpectrum = self.sumSpectrum;
            %% Find good initial guesses
            centers = zeros(length(self.spectrumPeaks),1);
            heights = zeros(length(self.spectrumPeaks),1);
            phases = zeros(length(self.spectrumPeaks),1);
            linewidths = zeros(length(self.spectrumPeaks),1);
            for i = 1:length(self.spectrumPeaks)
                centers(i) = self.spectrumPeaks(i).center;
                heights(i) = abs(sumSpectrum(centers(i)));
                phases(i) = angle(sumSpectrum(centers(i)));
                [lefti, righti] = HypWright.ScanClass.FWHM(centers(i),...
                    abs(sumSpectrum));
                linewidths(i) = righti-lefti;
                if (linewidths(i) > 40), linewidths(i) = 40; end
            end
            %% Fit the summed spectrum
            xAxis = 1:length(self.ppmAx);
            opts = optimset('lsqcurvefit');
            opts = optimset(opts,'Display','off');
            fun = @(x,xdata)HypWright.ScanClass.genSplitSpectrum(x(:,1),...
                x(:,2),x(:,3),x(:,4),xdata);
            guesses = [centers,heights,phases,linewidths];
            [fitHeights,resnorm,residual,exitflag,output] = lsqcurvefit(fun,...
                guesses, xAxis, [real(sumSpectrum) imag(sumSpectrum)]...
                ,[],[],opts);
            %% Display for testing
            if(self.verbose)
                fprintf(['Summed spectrum fit with norm: %f'...
                    ' and residual: %f as a fraction of pyruvate peak\n'],...
                    mean(resnorm)/max(heights), mean(residual)/max(heights))
                disp(exitflag)
                disp(output)
                tmpSpectrum = HypWright.ScanClass.genSpectrum(centers,...
                    heights, phases, linewidths,1:length(self.ppmAx));
                fitSpectrum = HypWright.ScanClass.genSpectrum(fitHeights(:,1),...
                    fitHeights(:,2), fitHeights(:,3), fitHeights(:,4),...
                    1:length(self.ppmAx));
                figure('units','normalized','outerposition',[0.3 0.1 .6 0.9]);
                subplot(3,1,1),plot(self.ppmAx,real(tmpSpectrum),'b',...
                    self.ppmAx,imag(tmpSpectrum),'r')
                title('input Guess'),xlabel('ppm')
                subplot(3,1,2),plot(self.ppmAx,real(sumSpectrum),'b',...
                    self.ppmAx,imag(sumSpectrum),'r')
                title('Spectrum'),xlabel('ppm')
                subplot(3,1,3),plot(self.ppmAx,real(fitSpectrum),'b',...
                    self.ppmAx,imag(fitSpectrum),'r')
                title('Fit Spectrum'),xlabel('ppm')
            end
            %% Store the fit
            for i = 1:length(self.spectrumPeaks)
                self.spectrumPeaks(i).params = fitHeights(i);
            end
            centers = fitHeights(:,1);
            phases = fitHeights(:,3);
            linewidths = fitHeights(:,4);
            guesses = zeros(length(centers),1);
            %% Fit each time point
            if(self.verbose)
                figure('units','normalized','outerposition',[0.3 0.1 .6 0.9]);
            end
            fitSpectrum = zeros(size(spectralData));
            % find the location of the pyruvate and lactate peaks
            ipyr = HypWright.ScanClass.findPeak(self.spectrumPeaks,'pyr');
            ilac = HypWright.ScanClass.findPeak(self.spectrumPeaks,'lac');
            self.pyrSignalLor = zeros(length(self.header.nLines),1);
            self.lacSignalLor = zeros(length(self.header.nLines),1);
            for i = 1:self.header.nLines
                tmpSpectrum = spectralData(i,:);
                % Guess peak height
                for j = 1:length(self.spectrumPeaks)
                    guesses(j) = abs(tmpSpectrum(self.spectrumPeaks(j).center));
                end
                fun = @(x,xdata)HypWright.ScanClass.genSplitSpectrum(centers,...
                    x,phases,linewidths,xdata);
                [fitHeights,resnorm,residual,exitflag,output] = lsqcurvefit(...
                    fun, guesses, xAxis, [real(tmpSpectrum) imag(tmpSpectrum)]...
                    ,[],[],opts);
                fitSpectrum(i,:) = HypWright.ScanClass.genSpectrum(centers,...
                    fitHeights, phases,linewidths,1:length(self.ppmAx));
                % Compute the area under each Lorentzian
                freqPerPixel = 2*pi*self.header.spectralWidth/...
                    self.header.lineLength;
                self.pyrSignalLor(i) = abs(pi*linewidths(ipyr)*...
                    fitHeights(ipyr)*freqPerPixel);
                self.lacSignalLor(i) = abs(pi*linewidths(ilac)*...
                    fitHeights(ilac)*freqPerPixel);
                %% Display for testing
                if(self.verbose)
                    subplot(2,1,1), plot(self.ppmAx,real(tmpSpectrum),'b',...
                        self.ppmAx, imag(tmpSpectrum),'r')
                    title('Data'),xlabel('ppm')
                    subplot(2,1,2), plot(self.ppmAx,real(fitSpectrum(i,:)),...
                        'b',self.ppmAx, imag(fitSpectrum(i,:)),'r')
                    title('Fit'),xlabel('ppm')
                    pause(0.1)
                end
            end
            %% Calculate the metabolit curves
            %% Display for testing
            if(self.verbose)
                figure('units','normalized','outerposition',[0.3 0.1 .6 0.9]);
                subplot(2,1,1),imagesc(self.ppmAx,self.fullTax,abs(fitSpectrum))
                title('fit spectrum'), xlabel('ppm'), ylabel('Time (sec)')
                subplot(2,1,2),imagesc(self.ppmAx,self.fullTax,abs(spectralData))
                title('Data'), xlabel('ppm'), ylabel('Time (sec)')
                figure
                plot(1:self.header.nLines,self.pyrSignalLor,'g',...
                    1:self.header.nLines,self.lacSignalLor,'b')
                title('Lorentzian integration results'), xlabel('Time (sec)')
            end
        end
        function integrateFWHM(self)
            % INTEGRATEFWHM: integrates the pyruvate and potentially lactate
            % peaks to give signal over time data
            sumSpectrum = abs(self.sumSpectrum);
            spectralData = self.spectralData;
            %% Find the peak centers
            iPyr = HypWright.ScanClass.findPeak(self.spectrumPeaks,'pyr');
            iLac = HypWright.ScanClass.findPeak(self.spectrumPeaks,'lac');
            specCenters = zeros(length(self.spectrumPeaks),1);
            specNames = cell(length(self.spectrumPeaks),1);
            for i = 1:length(self.spectrumPeaks)
                specCenters(i) = self.spectrumPeaks(i).center;
                specNames{i} = self.spectrumPeaks(i).name;
            end
            %% find the FWHM range
            if iPyr ~= 0
                pyrCenter = specCenters(iPyr);
                tmpSpectrum = real(sumSpectrum*exp(-1i*angle(sumSpectrum(...
                    pyrCenter))));
                [pyrLeftI, pyrRightI] = HypWright.ScanClass.FWHM(pyrCenter,...
                    tmpSpectrum);
            end
            if iLac ~= 0
                lacCenter = specCenters(iLac);
                tmpSpectrum = real(sumSpectrum*exp(-1i*angle(sumSpectrum(...
                    pyrCenter))));
                [lacLeftI, lacRightI] = HypWright.ScanClass.FWHM(lacCenter,...
                    tmpSpectrum);
            end
            %% Phase and sum data
            self.pyrSignalFWHM = zeros(length(self.header.nLines));
            for i = 1:self.header.nLines
                tmpSpectrum = real(spectralData*exp(-1i*angle(spectralData(...
                    i,pyrCenter))));
                self.pyrSignalFWHM(i) = HypWright.ScanClass.peakIntegration(...
                    pyrLeftI,pyrRightI,tmpSpectrum(i,:));
                tmpSpectrum = real(spectralData*exp(-1i*angle(spectralData(...
                    i,lacCenter))));
                self.lacSignalFWHM(i) = HypWright.ScanClass.peakIntegration(...
                    lacLeftI,lacRightI,tmpSpectrum(i,:));
            end
            %% Display for testing
            if(self.verbose)
                figure
                plot(1:length(sumSpectrum),sumSpectrum,pyrRightI,...
                    sumSpectrum(pyrCenter)/2,'go',...
                    pyrLeftI,sumSpectrum(pyrCenter)/2,'go',...
                    lacRightI,sumSpectrum(lacCenter)/2,'bo',...
                    lacLeftI,sumSpectrum(lacCenter)/2,'bo')
                title('FWHM Check')
                figure
                plot(self.fullTax,self.pyrSignalFWHM,'g',self.fullTax,...
                    self.lacSignalFWHM,'b')
                title('FWHM time courses'), xlabel('Time(sec)')
            end
        end
        function readHeader(self, headerFilePath)
            % READHEADER: Pulls relavant data from the Bruker header located at
            % the headerFilePath string. This method is not the best way to do
            % this and should be a priority for futher development.
            headerFileID = fopen(headerFilePath,'rt');
            self.header = struct('init',false);
            count = 1;
            while count > 0;
                [tmpStr, count] = fscanf(headerFileID,'%s=/n',1);
                tmpStrVect = strsplit(tmpStr,'=');
                switch(tmpStrVect{1})
                    case '##$ACQ_size'
                        % throw away the first 2 strings
                        fscanf(headerFileID,'%s',1);
                        fscanf(headerFileID,'%s',1);
                        % Store number of points per line and number of lines in
                        % the header
                        tmpLineLength = str2double(fscanf(headerFileID,'%s',1));
                        % Account for real and imaginary parts being storred
                        % seperatly
                        self.header.('lineLength') = tmpLineLength/2;
                        self.header.('nLines') = str2double(fscanf(...
                            headerFileID,'%s',1));
                    case '##$SW_h'
                        self.header.('spectralWidth') = str2double(...
                            tmpStrVect{2});
                end
            end
            fclose(headerFileID);
            self.header.init = true;
        end
        function readSpflash(self, fidFilePath)
            % READSPFLASH: Reads the fid file from a SPFlash scan specified by
            % the fidFilePath input string
            fidFileID = fopen(fidFilePath,'r','l'); % Open fid file, little endian
            % Read the Binary data into a array storing real and complex parts
            % seperatly, the raw data from bruker fid files in stored as a 32
            % bit int.
            m = self.header.lineLength*2;       % accounts for real and imaginary parts
            n = self.header.nLines;             % number of lines
            fileData = fread(fidFileID,[m,n],'int32');
            realData = fileData(1:2:end,:);       % pull out real part
            imaginaryData = fileData(2:2:end,:);  % pull out imaginary part
            % generate a complex array of proper size to be stored as rawData
            self.rawData = complex(realData,imaginaryData).';
            fclose(fidFileID);
        end
        
    end
    methods (Static, Access = private)
        function peaks = fillPeaks(names,centers)
            % FILLPEAKS: takes a list of names and corrosponding centers and
            % packs them into a struct
            peaks(length(names)) = struct('name',[],'center',[],'params',[]);
            for i = 1:length(names)
                peaks(i).name = names{i};
                peaks(i).center = centers(i);
            end
        end
        function iPeak = findPeak(peakStruct,peakName)
            % FINDPEAK: finds the index of the pyruvate peak in the passed in
            % peak structure
            % findPeak(peakStruct)
            switch (peakName)
                case 'pyr'
                    formats = {'pyruvate','pyr'};
                case 'lac'
                    formats = {'lactate','lac'};
            end
            iPeak = 0;
            specNames = cell(length(peakStruct),1);
            for i = 1:length(peakStruct)
                specNames{i} = peakStruct(i).name;
            end
            for i = 1:length(formats)
                if(iPeak == 0)
                    iPeak = find(strcmpi(formats{i},specNames));
                end
            end
            if (iPeak == 0)
                error('%s was not found in the peak structure!',peakName) 
            end
        end
        function [halfHeightLeft, halfHeightRight] = FWHM(center,ydata)
            % FWHM: returns the right and left indicies of the FWHM of some data
            % ydata centered at center NOTE will return fractional indicies if
            % no index is at the FWHM
            % See Also peakIntegration
            ydata = ydata./ydata(center);
            for i = center:length(ydata)
                if(ydata(i) < 0.5)
                    break
                end
                righti = i;
            end
            for i = center:-1:1
                if(ydata(i) < 0.5)
                    break
                end
                lefti = i;
            end
            halfHeightLeft = interp1([ydata(lefti-1), ydata(lefti)],...
                [lefti-1, lefti],0.5);
            halfHeightRight = interp1([ydata(righti), ydata(righti+1)],...
                [righti, righti+1],0.5);
        end
        function val = peakIntegration(leftPoint,rightPoint,ydata)
           % PEAKINTEGRATION: integrates some data from the left point to a
           % right point
           % peakIntegration(leftPoint,rightPoint,ydata) - integrates y data
           % from leftPoint to rightPoint where left and right points are the
           % index of the bounds. NOTE if non integer left and right points are
           % passed in then ydata will be interpolated at the left and right
           % edges
           % See Also FWHM
           xdata = [leftPoint,ceil(leftPoint):floor(rightPoint),rightPoint];
           leftYPoint = interp1(floor(leftPoint):ceil(leftPoint),...
               ydata(floor(leftPoint):ceil(leftPoint)),leftPoint);
           rightYPoint = interp1(floor(rightPoint):ceil(rightPoint),...
               ydata(floor(rightPoint):ceil(rightPoint)),rightPoint);
           ydata = [leftYPoint, ydata(ceil(leftPoint):floor(rightPoint)),...
               rightYPoint];
           val = trapz(xdata,ydata);
        end
        function spectrum = genSpectrum(centers,heights,phases,linewidths,xAxis)
            % GENSPECTRUM: generates a spectrum composed of Lorentzian 
            % distrobution
            % genSpectrum(center,height,phase,linewidth,xAx) Uses the shape
            % parameters centers, heights, phases and linewidths to generate a
            % a series of Lorentzian Distrobutions over xAxis
            if ~isequal(size(centers),size(heights),size(phases),...
                    size(linewidths))
               error('arg passed into genSpectrum are not of the same length') 
            end
            spectrum = zeros(length(xAxis),1);
            for i = 1:length(centers)
                spectrum = spectrum + HypWright.ScanClass.Lorentzian(...
                    centers(i),heights(i), phases(i), linewidths(i),xAxis);
            end
        end
        function lineBroadenedData = lineBroadening(rawData, sTax, bandwidth)
            % LINEBROADENING: Applies a line broadening fliter to some data
            lb=bandwidth*2*pi;     % lb: radian
            em = exp(-(sTax.*sTax)*lb*lb/2/4/log(2));
            lineBroadenedData = rawData;
            for i=1:length(rawData(:,1))
                lineBroadenedData(i,:)=rawData(i,:).*em;
            end
        end
        function LorDist = Lorentzian(center,height,phase,linewidth,xAxis)
            % LORENTZIAN: generates a Lorentzian distrobution 
            % Lorentzian(center,height,phase,linewidth,xAx) Uses the shape
            % parameters center, height, phase and linewidth to generate a
            % Lorentzian over xaxis
            LorDist = ((height*linewidth*exp(1i*phase))./...
                (linewidth+1i*(xAxis-center))).';
        end
        function spectrum = genSplitSpectrum(centers,heights,phases,linewidths,...
                xAxis)
            % GENSPLITSPECTRUM: generates a spectrum composed of Lorentzian 
            % distrobution seperated into real and imaginary parts
            % genSplitSpectrum(centers,heights,phases,linewidths,xAxis) Uses the
            % parameters, centers, heights phases and linewidths to generate
            % a series of Lorentzian Distrobutions over xAxis. NOTE: the real
            % and imaginary parts are broken up so the return will be twice as
            % long as xAxis with all real points followed by all the imaginary
            % points
            tmpSpectrum = HypWright.ScanClass.genSpectrum(centers,heights,...
                phases,linewidths,xAxis);
            spectrum = [real(tmpSpectrum); imag(tmpSpectrum)].';
        end
    end   
end

