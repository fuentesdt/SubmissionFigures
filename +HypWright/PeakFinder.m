classdef PeakFinder < handle
    %PEAKFINDER GUI to select the peaks present in a spectrum
    %   Detailed explanation goes here
    
    properties
        GUIHandle;
        data
        VisiblePeaks = [0, 0, 0, 0, 0];
        PeakNames= {'Pyruvate', 'Lactate', 'Pyruvate-Hydrate', 'Urea', 'Alanine'};
        PeakIndex = [0,0,0,0,0];
        currentPeak
        dcm_obj
        nop
        sw
        ppmAxes
    end
    
    properties(Dependent)
        PeakList
    end
    
    methods
        %% Constructor:
        function self = PeakFinder(Data,varargin)
            p = inputParser();
            p.addRequired('Data',@isnumeric)
            p.addParamValue('nop',2048,@isscalar)
            p.addParamValue('sw',4.9603e+03,@isscalar)
            p.parse(Data,varargin{:})
            self.data = p.Results.Data;
            self.nop = p.Results.nop;
            self.sw = p.Results.sw;
            self.initGUI();
            uiwait(gcf);
        end
        %% Getters and Setters
        function PeakList = get.PeakList(self)
            j = 1;
            PeakList = {};
            for i = 1:length(self.VisiblePeaks)
                if self.VisiblePeaks(i) == 1
                    PeakList(j) = self.PeakNames(i);
                    j = j+1;
                end
            end
        end
        %% Callbacks
        function chkbox(self,i)
            if self.VisiblePeaks(i) == 1
                self.VisiblePeaks(i) = 0;
            else 
                self.VisiblePeaks(i) = 1;
            end
            set(self.GUIHandle.Met_listbox,'String',self.PeakList)
            self.PeakSel();
            self.CalPeaks();
        end
        function PeakSel(self)
             val = get(self.GUIHandle.Met_listbox,'Value');
             srt = get(self.GUIHandle.Met_listbox,'String');
             self.currentPeak = srt(val);
        end
        %% CalPeaks: Calculate peak locations
        function CalPeaks(self)
            info = getCursorInfo(self.dcm_obj);
            if(isfield(info,'DataIndex'))
                switch(self.currentPeak{1})
                    case 'Pyruvate'
                        pyruvate = info.DataIndex;
                    case 'Lactate'
                        pyruvate = info.DataIndex-264;
                    case 'Pyruvate-Hydrate'
                        pyruvate = info.DataIndex-264-119;
                    case 'Urea'
                        pyruvate = info.DataIndex+208;
                    case 'Alanine'
                        pyruvate = info.DataIndex-264-119-60;
                end
                if self.VisiblePeaks(1)                                 %Pyruvate
                    [C , I] = max(self.data(pyruvate-20:pyruvate+20));
                    self.PeakIndex(1) = I+pyruvate-21;
                end
                if self.VisiblePeaks(2)                                 %Pyruvate-Hydrate
                    tmp = pyruvate+264+119;
                    [C , I] =  max(self.data(tmp-20:tmp+20));
                    self.PeakIndex(2) = I+tmp-21;
                end
                if self.VisiblePeaks(3)                                 %Lactate
                    tmp = pyruvate+264;
                    [C , I] = max(self.data(tmp-20:tmp+20));
                    self.PeakIndex(3) = I+tmp-21;
                end
                if self.VisiblePeaks(4)                                 %Urea
                    tmp = pyruvate-208;
                    [C , I] = max(self.data(tmp-20:tmp+20));
                    self.PeakIndex(4) = I+tmp-21;
                end
                if self.VisiblePeaks(5)                                 %Alanine
                    tmp = pyruvate+180;
                    [C , I] = max(self.data(tmp-20:tmp+20));
                    self.PeakIndex(5) = I+tmp-21;
                end
                self.CalPPM()
                self.ShowPeaks()
                set(self.GUIHandle.MetValsDisp,'String',self.PeakIndex)
            end
        end
        %% CalPPM: calculates and sets ppm axis
        function CalPPM(self)
            self.ppmAxes = ((1:self.nop)-self.PeakIndex(1))*...
                (self.sw/75.53140/self.nop)+171;
            plot(self.GUIHandle.axes,self.ppmAxes,self.data)
            xlabel('ppm');ylabel('Signal');title('Sum of Spectrum');
        end
        %% ShowPeaks: displays peaks
        function ShowPeaks(self)
            if(get(self.GUIHandle.ShowPeaks,'Value'))
                for i=1:length(self.VisiblePeaks)
                    if(self.VisiblePeaks(i))
                        text(self.ppmAxes(self.PeakIndex(i)),self.data(self.PeakIndex(i)),self.PeakNames{i})
                    end
                end
            else
                plot(self.GUIHandle.axes,self.ppmAxes,self.data)
                xlabel('ppm');ylabel('Signal');title('Sum of Spectrum');
            end
        end
        %% initGUI:
        function initGUI(self)
            sz=get(0,'ScreenSize');
            self.GUIHandle.PeakFig = figure('Units','pixel','Interruptible','off','Position',[10 sz(4)-640 820 580],...
                'Name','Peak Tool','Tag','Peaktool','MenuBar','none','Color',[0.83,0.81,0.78],'Resize','off');
            set(self.GUIHandle.PeakFig,'BusyAction','queue');
            self.GUIHandle.axes=axes('Position',[0.25 0.2 0.7 0.7]);
            plot(self.GUIHandle.axes,self.data)
            datacursormode on
            self.dcm_obj = datacursormode(gcf);
            xlabel('ppm');ylabel('Signal');title('Sum of Spectrum');
            self.GUIHandle.Met_panel = uipanel('Title','Metabolites','FontSize',12,...
                'BackgroundColor','white','Position',[0 0.2 0.19 0.8]);
            self.GUIHandle.Pyr_chkbox = uicontrol(self.GUIHandle.Met_panel,'Style','checkbox','String','Pyruvate','Units','pixel','Position',[10,130,120,20],...
                'Callback',@(src,event)self.chkbox(1));
            self.GUIHandle.Lac_chkbox = uicontrol(self.GUIHandle.Met_panel,'Style','checkbox','String','Lactate','Units','pixel','Position',[10,100,120,20],...
                'Callback',@(src,event)self.chkbox(2));
            self.GUIHandle.PyHy_chkbox = uicontrol(self.GUIHandle.Met_panel,'Style','checkbox','String','Pyruvate-Hydrate','Units','pixel','Position',[10,70,120,20],...
                'Callback',@(src,event)self.chkbox(3));
            self.GUIHandle.Urea_chkbox = uicontrol(self.GUIHandle.Met_panel,'Style','checkbox','String','Urea','Units','pixel','Position',[10,40,120,20],...
                'Callback',@(src,event)self.chkbox(4));
            self.GUIHandle.Alanine_chkbox = uicontrol(self.GUIHandle.Met_panel,'Style','checkbox','String','Alanine','Units','pixel','Position',[10,10,120,20],...
                'Callback',@(src,event)self.chkbox(5));
            self.GUIHandle.Met_listbox = uicontrol(self.GUIHandle.Met_panel,'Style','listbox','String',self.PeakList,'Units','pixel','Position',[10,300,120,100],...
                'Callback',@(src,event)self.PeakSel());
            self.GUIHandle.Peak = uicontrol('Style','pushbutton','String','Calculate Peaks','Position',[100,50,100,60],'Callback',@(src,event)self.CalPeaks());
            self.GUIHandle.ShowPeaks = uicontrol('Style','radiobutton','String','Show Peaks','Position',[10,10,100,20],'Callback',@(src,event)self.ShowPeaks());
            self.GUIHandle.MetValsDisp = uicontrol(self.GUIHandle.Met_panel,'Style','text','String',self.PeakIndex,'Units','pixel','Position',[90,200,50,80]);
            self.GUIHandle.MetDisp = uicontrol(self.GUIHandle.Met_panel,'Style','text','String',self.PeakNames,'Units','pixel','Position',[10,200,90,80]);
        end
        %% Set Xaxes
        %% getPeaks: returns a Peak Object
        function Specs = getPeaks(self)
           Specs = Peaks();
           for i = 1:4
              if(self.VisiblePeaks(i))
                  Specs.addPeak(self.PeakIndex(i),self.PeakNames(i))
              end
           end
        end
    end
    
end

