classdef DiffAppV02_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        DiffusonCalculatorforSIEMENSUIFigure  matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        CalculatorTab                   matlab.ui.container.Tab
        SmoothingDegreeEditField        matlab.ui.control.NumericEditField
        SmoothingDegreeEditFieldLabel   matlab.ui.control.Label
        DtUpperLimitSegEditField        matlab.ui.control.NumericEditField
        DtUpperLimitSegEditFieldLabel   matlab.ui.control.Label
        DpUpperLimitSegEditField        matlab.ui.control.NumericEditField
        DpUpperLimitSegEditFieldLabel   matlab.ui.control.Label
        mapRoiSwitch                    matlab.ui.control.Switch
        FilterDWIBetaButton             matlab.ui.control.Button
        StateField                      matlab.ui.control.EditField
        DtUpperLimitFreeEditField       matlab.ui.control.NumericEditField
        DtUpperLimitFreeEditFieldLabel  matlab.ui.control.Label
        DpUpperLimitFreeEditField       matlab.ui.control.NumericEditField
        DpUpperLimitFreeEditFieldLabel  matlab.ui.control.Label
        KernelSizeEditField             matlab.ui.control.NumericEditField
        KernelSizeEditFieldLabel        matlab.ui.control.Label
        SliceKnob                       matlab.ui.control.DiscreteKnob
        SliceKnobLabel                  matlab.ui.control.Label
        BValueKnob                      matlab.ui.control.DiscreteKnob
        BValueKnobLabel                 matlab.ui.control.Label
        GenerateMapsButton              matlab.ui.control.Button
        EditField_2                     matlab.ui.control.NumericEditField
        EditField                       matlab.ui.control.NumericEditField
        AlgorithmButtonGroup            matlab.ui.container.ButtonGroup
        BiExpSegButton                  matlab.ui.control.RadioButton
        BiexpFreeButton                 matlab.ui.control.RadioButton
        MonoExpButton                   matlab.ui.control.RadioButton
        LoadDICOMButton                 matlab.ui.control.Button
        OptionsPanel                    matlab.ui.container.Panel
        MultiRoiSwitch                  matlab.ui.control.Switch
        DrawMultipleROIButton           matlab.ui.control.Button
        DrawROIButton                   matlab.ui.control.StateButton
        ROIStyleDropDown                matlab.ui.control.DropDown
        ROIStyleDropDownLabel           matlab.ui.control.Label
        CalculateButton                 matlab.ui.control.Button
        ADCEditField                    matlab.ui.control.NumericEditField
        ADCEditFieldLabel               matlab.ui.control.Label
        pFEditField                     matlab.ui.control.NumericEditField
        pFEditFieldLabel                matlab.ui.control.Label
        DpEditField                     matlab.ui.control.NumericEditField
        DpEditFieldLabel                matlab.ui.control.Label
        DtEditField                     matlab.ui.control.NumericEditField
        DtEditFieldLabel                matlab.ui.control.Label
        ImageViewPanel                  matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
        FittingGraphicsPanel            matlab.ui.container.Panel
        RSquareEditField                matlab.ui.control.NumericEditField
        RSquareEditFieldLabel           matlab.ui.control.Label
        UIAxes2                         matlab.ui.control.UIAxes
        DicomTab                        matlab.ui.container.Tab
        DicomState                      matlab.ui.control.EditField
        PatientsPanel                   matlab.ui.container.Panel
        PatientsListBox                 matlab.ui.control.ListBox
        PatientsListBoxLabel            matlab.ui.control.Label
        SequencesPanel                  matlab.ui.container.Panel
        LoadtoBOLDPcsButton             matlab.ui.control.Button
        LoadToSegmenterButton           matlab.ui.control.Button
        FramesEditField                 matlab.ui.control.NumericEditField
        FramesEditFieldLabel            matlab.ui.control.Label
        LoadtoViewerButton              matlab.ui.control.Button
        LoadtoCalculatorButton          matlab.ui.control.Button
        SequencesListBox                matlab.ui.control.ListBox
        SequencesListBoxLabel           matlab.ui.control.Label
        ViewerTab                       matlab.ui.container.Tab
        ImageViewPanel_2                matlab.ui.container.Panel
        BritghtnessSlider               matlab.ui.control.Slider
        BritghtnessSliderLabel          matlab.ui.control.Label
        ContrastSlider                  matlab.ui.control.RangeSlider
        ContrastSliderLabel             matlab.ui.control.Label
        ColormapDropDown                matlab.ui.control.DropDown
        ColormapDropDownLabel           matlab.ui.control.Label
        InstanceNumberEditField         matlab.ui.control.NumericEditField
        InstanceNumberEditFieldLabel    matlab.ui.control.Label
        ViewSlider                      matlab.ui.control.Slider
        ViewSliderLabel                 matlab.ui.control.Label
        UIAxes8                         matlab.ui.control.UIAxes
        MapsTab                         matlab.ui.container.Tab
        DWIMAPSPanel                    matlab.ui.container.Panel
        ColorMapsDropDown               matlab.ui.control.DropDown
        ColorMapsDropDownLabel          matlab.ui.control.Label
        UIAxes3_7                       matlab.ui.control.UIAxes
        UIAxes3_6                       matlab.ui.control.UIAxes
        UIAxes3_5                       matlab.ui.control.UIAxes
        UIAxes3_4                       matlab.ui.control.UIAxes
        UIAxes3_3                       matlab.ui.control.UIAxes
        UIAxes3_2                       matlab.ui.control.UIAxes
        UIAxes3                         matlab.ui.control.UIAxes
        SegmenterBetaTab                matlab.ui.container.Tab
        CircularityEditField            matlab.ui.control.NumericEditField
        CircularityEditFieldLabel       matlab.ui.control.Label
        AreaEditField                   matlab.ui.control.NumericEditField
        AreaEditFieldLabel              matlab.ui.control.Label
        RegionStatisticsButton          matlab.ui.control.Button
        DropPointButton                 matlab.ui.control.Button
        AddPointButton                  matlab.ui.control.Button
        SegmenterStateField             matlab.ui.control.EditField
        SegmentROIButton                matlab.ui.control.Button
        SegmentButton                   matlab.ui.control.Button
        SegmenterSlider                 matlab.ui.control.Slider
        SegmenterSliderLabel            matlab.ui.control.Label
        UIAxes9_2                       matlab.ui.control.UIAxes
        UIAxes9                         matlab.ui.control.UIAxes
        BOLD_Proceesor                  matlab.ui.container.Tab
        DrawROIforFitButton             matlab.ui.control.StateButton
        BoldState                       matlab.ui.control.EditField
        GenerateR2MapButton             matlab.ui.control.Button
        BoldSliceKnob                   matlab.ui.control.DiscreteKnob
        BoldSliceLabel                  matlab.ui.control.Label
        TEKnob                          matlab.ui.control.DiscreteKnob
        TEKnobLabel                     matlab.ui.control.Label
        R2MapandFittingPanel            matlab.ui.container.Panel
        MeanR2S1EditField               matlab.ui.control.NumericEditField
        MeanR2S1EditFieldLabel          matlab.ui.control.Label
        UIAxes12                        matlab.ui.control.UIAxes
        UIAxes11                        matlab.ui.control.UIAxes
        R2ImagesPanel                   matlab.ui.container.Panel
        UIAxes10                        matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        width;
        height;
        vol ;
        instanceSize=0;
        showBval=0;
        bvals;
        sliceNum=1;
        roi;
        pos;
        posx1;
        posx2;
        posy1;
        posy2;
        M;
        N;
        Z;
        sigData = [0];
        s_s0;
        adcfit;
        ADC=0;
        segFitHigh;
        segFitAll;
        freeFit;
        sigDataHigh;
        bvalHigh;
        paramFree;
        paramSeg;
        Dt=0;
        Dp=0;
        pF=0;
        collection;
        rows;
        cols;
        rsqr;
        mapMask;
        meanImage;
        ADCmap;
        Dtmap;
        Dpmap;
        pFmap;
        fo;
        fotype;
        segfo;
        fo2;
        foadc;
        mymodel;
        DpSegmap;
        DtSegmap;
        pFSegmap;
        kernel;
        DpUpperLimitfree = 0.05;
        DtUpperLimitfree = 0.005;
        DpUpperLimitseg = 0.05;
        DtUpperLimitseg = 0.005;
        bvals_ivim;
        s_s0_ivim;
        b_map_ivim;
        patientvalue;
        sequencevalue;
        patienttable;
        viewVol;
        viewSlidervalue;
        viewM;
        viewN;
        viewZ;
        cmaps = ["gray","jet","turbo","hsv","parula","hot","cool","sky","winter"];
        combmask;
        adctype;
        frames;
        smoothDegree = 0.5;
        segmenterVol;
        segmenterM;
        segmenterN;
        segmenterZ;
        segmenterBW;
        segmenterMasked;
        % segmentROI;
        segmentIm;
        segroi;
        tempax;
        analyzeFigure;
        contrastLowLimit=0;
        contrastUpLimit=2.5;
        contrastValue;
        brightnessValue=0;
        viewImage;
        viewAdjustedImage;
        samBox;
        samObj;
        samMasks;
        samScores;
        imgEmbeddings;
        forePoints;
        backPoints;
        selectedTab;
        boldVol;
        TEvals;
        showTEval = 0
        boldInstace = 0;
        boldSlice = 1;
        R2star;
        boldmask;
        boldfittype;
        boldopts;
        bM;
        bN;
        bZ;
        OEFvol
      
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadDICOMButton
        function LoadDICOMButtonPushed(app, event)
            
            folder = uigetdir;
            try
            app.collection = dicomCollection(folder,IncludeSubfolders=true);
            app.PatientsListBox.Items = unique(app.collection.PatientName);
            % app.PatientsListBox.ItemsData = 1:length(app.PatientsListBox.Items);
            app.patientvalue = app.PatientsListBox.Value;
            app.SequencesListBox.Items = app.collection.SeriesDescription(find(app.collection.PatientName == app.patientvalue));
            app.SequencesListBox.ItemsData = 1:length(app.SequencesListBox.Items);
            app.sequencevalue = app.SequencesListBox.Items(app.SequencesListBox.Value);
            app.ColormapDropDown.Items = app.cmaps;
            app.ColorMapsDropDown.Items = app.cmaps;
            app.SmoothingDegreeEditField.Value = app.smoothDegree;
            app.StateField.Value = "DICOM Collection Installed. See Dicom Tab.";
            app.StateField.BackgroundColor = [1.0,0.48,0.51];
            catch
                app.StateField.Value = "Error Occured";
            end
            
            



            

            
             
            

            
        end

        % Value changed function: DrawROIButton
        function DrawROIButtonValueChanged(app, event)
            app.sigData = [0];
             delete(findall(app.UIAxes,'Type','images.roi'));
             % findall(app.UIAxes,'Type','images.roi.Polygon');
             % findall(app.UIAxes,'Type','images.roi.Freehand');
            
            if app.ROIStyleDropDown.Value == "Rectangle"
                app.ROIStyleDropDown.Value
                app.roi = drawrectangle('Parent',app.UIAxes);
                mask = createMask(app.roi);
                
            end
            if app.ROIStyleDropDown.Value == "Polygon"
                app.roi = drawpolygon('Parent',app.UIAxes);
                mask = createMask(app.roi);
            end
            if app.ROIStyleDropDown.Value == "Freehand"
                app.roi = drawfreehand('Parent',app.UIAxes);
                mask = createMask(app.roi);
            end
            if app.ROIStyleDropDown.Value == "Assisted"
                
                app.roi = drawassisted(findobj(app.UIAxes,'Type','Image'));
                mask = createMask(app.roi);
            end


               for i = 0:1:length(app.bvals)-1
                    tempvol = app.vol(:,:,(i+1)+((app.sliceNum-1)*length(app.bvals)));
               
                    meanData = mean(tempvol(mask));
                    app.sigData = horzcat(app.sigData,meanData);
              end
            
            app.sigData = app.sigData(2:end);
            app.s_s0 = app.sigData./app.sigData(1);
            % app.s_s0 = app.s_s0(find(app.s_s0 <= 1));
            % app.s_s0
            % app.bvals = app.bvals(find(app.s_s0 <= 1));
            % app.bvals
            plot(app.bvals,app.s_s0,'Parent',app.UIAxes2);
            title(app.UIAxes2,'Signal Data');
            xlabel(app.UIAxes2,'B-Value');
            ylabel(app.UIAxes2,'S / S0');
            
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            % try
            app.bvals_ivim = app.bvals(find(app.bvals<=1000));
            app.s_s0_ivim = app.s_s0(1:length(app.bvals_ivim));
            app.fotype = fittype('f * exp(-b * D_star) + (1 - f) * exp(-b * D)','independent', 'b', 'coefficients', {'f', 'D_star', 'D'});
            app.segfo = fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt');
            app.fo =fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[0.1,0.01,0.001],'Lower',[0,0,0],'Upper',[1,0.1,0.005]);
            app.fo2 = fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[0.01],'Lower',[0],'Upper',[0.1]);
            app.foadc =fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[0.001]);
            app.adctype =fittype( 'exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
            [app.adcfit, gadc] = fit(app.bvals_ivim.',app.s_s0_ivim.',app.adctype,app.foadc);
            [app.freeFit, gfree] = fit(app.bvals_ivim.',app.s_s0_ivim.',app.fotype,app.fo);
            app.paramFree = [app.freeFit.D app.freeFit.D_star app.freeFit.f];
            app.bvalHigh = app.bvals_ivim(find(app.bvals_ivim > 200));
            app.sigDataHigh = app.s_s0_ivim(find(app.bvals_ivim > 200));
            app.segFitHigh = fit(app.bvalHigh.',app.sigDataHigh.','exp1',app.segfo);
            dt = app.segFitHigh.b*-1;
            si = app.segFitHigh.a * app.s_s0_ivim(1);
            pf = (app.s_s0(1)-si)/app.s_s0_ivim(1);
            app.mymodel = @(dp,x) (pf*(exp(-x*dp)) + (1-pf)*exp(-x*dt));
            [app.segFitAll, gseg] = fit(app.bvals_ivim.',app.s_s0_ivim.',app.mymodel,app.fo2);
            dp = app.segFitAll.dp;
            app.paramSeg = [dp dt pf];

            app.ADC = app.adcfit.b;
            if (app.MonoExpButton.Value == 1)
                app.rsqr = gadc.rsquare;
                app.RSquareEditField.Value = app.rsqr;
                app.ADCEditField.Value = app.ADC;
                app.DpEditField.Value = 0;
                app.DtEditField.Value = 0;
                app.pFEditField.Value = 0;
                plot(app.UIAxes2,app.adcfit,app.bvals_ivim,app.s_s0_ivim);
                ylim(app.UIAxes2,[0 1]);
                title(app.UIAxes2,'Fitting Graphic');
                xlabel(app.UIAxes2,'B-Value');
                ylabel(app.UIAxes2,'S / S0');
            end
            if(app.BiexpFreeButton.Value == 1)
                app.rsqr = gfree.rsquare;
                app.RSquareEditField.Value = app.rsqr;
                app.ADCEditField.Value = 0;
                app.DpEditField.Value = app.paramFree(2);
                app.DtEditField.Value = app.paramFree(1);
                app.pFEditField.Value = app.paramFree(3);
                app.DtUpperLimitFreeEditField.Value = app.DtEditField.Value*2;
                app.DpUpperLimitFreeEditField.Value = app.DpEditField.Value*2;
                plot(app.UIAxes2,app.freeFit,app.bvals_ivim,app.s_s0_ivim);
                ylim(app.UIAxes2,[0 1]);
                title(app.UIAxes2,'Fitting Graphic');
                xlabel(app.UIAxes2,'B-Value');
                ylabel(app.UIAxes2,'S / S0');
            end
            if(app.BiExpSegButton.Value == 1)
                app.rsqr = gseg.rsquare;
                app.RSquareEditField.Value = app.rsqr;
                app.ADCEditField.Value = 0;
                app.DpEditField.Value = app.paramSeg(1);
                app.DtEditField.Value = app.paramSeg(2);
                app.pFEditField.Value = app.paramSeg(3);
                app.DtUpperLimitSegEditField.Value = app.DtEditField.Value*2;
                app.DpUpperLimitSegEditField.Value = app.DpEditField.Value*2;
                plot(app.UIAxes2,app.segFitAll,app.bvals_ivim,app.s_s0_ivim);
                ylim(app.UIAxes2,[0 1]);
                title(app.UIAxes2,'Fitting Graphic');
                xlabel(app.UIAxes2,'B-Value');
                ylabel(app.UIAxes2,'S / S0');
            end
            
            % catch
            %     app.StateField.Value = "Error Occured";
            % end
        end

        % Value changed function: BValueKnob
        function BValueKnobValueChanged(app, event)
            app.showBval = find(app.bvals == app.BValueKnob.Value);
            app.EditField.Value = app.BValueKnob.Value;
            % imshow(adapthisteq(uint8(app.vol(:,:,(app.showBval)+((app.sliceNum-1)*length(app.bvals)))),'clipLimit',0.005,'Distribution','rayleigh'),'Parent',app.UIAxes);
            imagesc(app.UIAxes,imadjust(app.vol(:,:,(app.showBval)+((app.sliceNum-1)*length(app.bvals)))));
            colormap(app.UIAxes,"gray");
            xlim(app.UIAxes,[0 ,app.N]);
            ylim(app.UIAxes,[0, app.M]);
            title(app.UIAxes,'IMAGE');
            
        end

        % Value changed function: SliceKnob
        function SliceKnobValueChanged(app, event)
            app.sliceNum = app.SliceKnob.Value;
            app.showBval = find(app.bvals == app.BValueKnob.Value);
            app.EditField_2.Value = app.SliceKnob.Value;
            % imshow(adapthisteq(uint8(app.vol(:,:,(app.showBval)+((app.sliceNum-1)*length(app.bvals)))),'clipLimit',0.005,'Distribution','rayleigh'),'Parent',app.UIAxes);
            imagesc(app.UIAxes,imadjust(app.vol(:,:,(app.showBval)+((app.sliceNum-1)*length(app.bvals)))));
            colormap(app.UIAxes,"gray");
            xlim(app.UIAxes,[0 ,app.N]);
            ylim(app.UIAxes,[0, app.M]);
            title(app.UIAxes,'IMAGE');
            
        end

        % Button pushed function: GenerateMapsButton
        function GenerateMapsButtonPushed(app, event)
                app.StateField.Value = "Processing...";
                app.StateField.BackgroundColor = [0.95,0.84,0.18];
                pause(0.5);
                try
                app.fotype = fittype('f * exp(-b * D_star) + (1 - f) * exp(-b * D)','independent', 'b', 'coefficients', {'f', 'D_star', 'D'});
                app.segfo = fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt');
                app.fo =fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[0.1,0.01,0.001],'Lower',[0,0,0],'Upper',[1,0.1,0.005]);
                app.fo2 = fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[0.01],'Lower',[0],'Upper',[0.1]);
                app.foadc =fitoptions('Method','NonLinearLeastSquares','Algorithm','Levenberg-Marquardt','StartPoint',[0.001]);
                app.adctype =fittype( 'exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
                if app.mapRoiSwitch.Value == "SingleROI"
                        app.mapMask = uint16(createMask(app.roi));
                elseif app.mapRoiSwitch.Value == "MultiROI"
                        app.mapMask = uint16(app.combmask);
                end
                [maskX, maskY] = size(app.mapMask);
                
                app.kernel = app.KernelSizeEditField.Value;
            
            for xc = round(app.kernel/2)+1:1:maskX-app.kernel
                for yc=round(app.kernel/2)+1:1:maskY-app.kernel
                    pixelSigMatrix = [0];
                    if app.mapMask(xc,yc) ~= 0
                        for i = 0:1:length(app.bvals)-1
                            calculationPixel = app.vol(xc-(round(app.kernel/2)):xc+(round(app.kernel/2)),yc-(round(app.kernel/2)):yc+(round(app.kernel/2)),(i+1)+((app.sliceNum-1)*length(app.bvals)));
                            pixelSigMatrix = horzcat(pixelSigMatrix,mean2(calculationPixel));
               
                    
                        end
                        app.b_map_ivim = app.bvals(find(app.bvals<=1000));
                        pixelSigMatrix = pixelSigMatrix(2:length(app.b_map_ivim)+1);
                        
                        pixelSigMatrix = pixelSigMatrix./pixelSigMatrix(1);
                        pixelBvalHigh = app.b_map_ivim(find(app.b_map_ivim>200));
                        pixelBvals = app.b_map_ivim;
                    
                        try
                        pixelSigMatrixHigh = pixelSigMatrix(find(app.b_map_ivim > 200));
                        pixelADCfit = fit(pixelBvals.',pixelSigMatrix.',app.adctype,app.foadc);
                        pixelFreefit = fit(pixelBvals.',pixelSigMatrix',app.fotype,app.fo);
                        pixelsegFitHigh = fit(pixelBvalHigh.',pixelSigMatrixHigh.','exp1',app.segfo);
                        DT = pixelsegFitHigh.b*-1;
                        SI = pixelsegFitHigh.a * pixelSigMatrix(1);
                        PF = (pixelSigMatrix(1)-SI)/pixelSigMatrix(1);
                        pixelmodel = @(DP,x) (PF*(exp(-x*DP)) + (1-PF)*exp(-x*DT));
                        pixelSegAllFit = fit(pixelBvals.',pixelSigMatrix.',pixelmodel,app.fo2);
                        DP = pixelSegAllFit.DP;
                        pixelParamFree = [pixelFreefit.D pixelFreefit.D_star pixelFreefit.f];
                        app.ADCmap(xc,yc) = pixelADCfit.b;
                        app.Dpmap(xc,yc) = pixelParamFree(2);
                        app.Dtmap(xc,yc) = pixelParamFree(1);
                        app.pFmap(xc,yc) = pixelParamFree(3);
                        app.DpSegmap(xc,yc) =DP;
                        app.DtSegmap(xc,yc) = DT;
                        app.pFSegmap(xc,yc) = PF;
                        catch
                        app.ADCmap(xc,yc) = 0;
                        app.Dpmap(xc,yc) = 0;
                        app.Dtmap(xc,yc) = 0;
                        app.pFmap(xc,yc) = 0;
                        app.DpSegmap(xc,yc) =0;
                        app.DtSegmap(xc,yc) = 0;
                        app.pFSegmap(xc,yc) = 0;
                        end
                    else
                        app.ADCmap(xc,yc) = 0;
                        app.Dpmap(xc,yc) = 0;
                        app.Dtmap(xc,yc) = 0;
                        app.pFmap(xc,yc) = 0;
                        app.DpSegmap(xc,yc) =0;
                        app.DtSegmap(xc,yc) = 0;
                        app.pFSegmap(xc,yc) = 0;
                    end
                    
                end
            end
                    app.DpUpperLimitfree = app.DpUpperLimitFreeEditField.Value;
                    app.DtUpperLimitfree = app.DtUpperLimitFreeEditField.Value;
                    app.DpUpperLimitseg = app.DpUpperLimitSegEditField.Value;
                    app.DtUpperLimitseg = app.DtUpperLimitSegEditField.Value
                     
                     imagesc(app.UIAxes3_7,app.ADCmap);
                     colormap(app.UIAxes3_7,"jet");
                     colorbar(app.UIAxes3_7);
                     xlim(app.UIAxes3_7,[0 ,app.N]);
                     ylim(app.UIAxes3_7,[0, app.M]);
                     
                     imagesc(app.UIAxes3,app.Dpmap,[0, app.DpUpperLimitfree]);
                     colormap(app.UIAxes3,"jet");
                     colorbar(app.UIAxes3);
                     xlim(app.UIAxes3,[0 ,app.N]);
                     ylim(app.UIAxes3,[0, app.M]);
                     
                     imagesc(app.UIAxes3_2,app.Dtmap,[0 ,app.DtUpperLimitfree]);
                     colormap(app.UIAxes3_2,"jet");
                     colorbar(app.UIAxes3_2);
                     xlim(app.UIAxes3_2,[0 ,app.N]);
                     ylim(app.UIAxes3_2,[0, app.M]);
                    
                     imagesc(app.UIAxes3_3,app.pFmap,[0 1]);
                     colormap(app.UIAxes3_3,"jet");
                     colorbar(app.UIAxes3_3);
                     xlim(app.UIAxes3_3,[0 ,app.N]);
                     ylim(app.UIAxes3_3,[0, app.M]);
                     
                     imagesc(app.UIAxes3_4,app.DpSegmap,[0, app.DpUpperLimitseg]);
                     colormap(app.UIAxes3_4,'jet');
                     colorbar(app.UIAxes3_4);
                     xlim(app.UIAxes3_4,[0 ,app.N]);
                     ylim(app.UIAxes3_4,[0, app.M]);
                   
                     imagesc(app.UIAxes3_5,app.DtSegmap,[0 ,app.DtUpperLimitseg]);
                     colormap(app.UIAxes3_5,"jet");
                     colorbar(app.UIAxes3_5);
                     xlim(app.UIAxes3_5,[0 ,app.N]);
                     ylim(app.UIAxes3_5,[0, app.M]);
                    
                     imagesc(app.UIAxes3_6,app.pFSegmap,[0 1]);
                     colormap(app.UIAxes3_6,"jet");
                     colorbar(app.UIAxes3_6);
                     xlim(app.UIAxes3_6,[0 ,app.N]);
                     ylim(app.UIAxes3_6,[0, app.M]);


                    app.StateField.Value = "Map Calculation is Done";
                    app.StateField.BackgroundColor = [0.29,0.84,0.16];
                    
                catch
                    app.StateField.Value = "Maps Could Not be Generated, Try Changing Kernel Size";
                end

            
        end

        % Button pushed function: LoadtoCalculatorButton
        function LoadtoCalculatorButtonPushed(app, event)
            app.DicomState.Value = "Loading..";
            pause(0.5);
            try
            app.frames = app.collection.Frames;
            b_val=[0];
            descriptions = app.collection.SeriesDescription;
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            I = find(app.patienttable.SeriesDescription == app.sequencevalue);
            % candidate = convertStringsToChars(descriptions(I));
            I = I(find(I == app.SequencesListBox.Value));
            % if candidate(1) == 'e'
                fileNames = app.patienttable.Filenames(I);
                app.frames = app.collection.Frames;
            % end
            
       
             app.rows = app.patienttable.Rows(I);
            app.cols = app.patienttable.Columns(I);
            app.vol = zeros(app.rows , app.cols);
            fileNames = fileNames{1,1};
            
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                info = dicominfo(f);
                seq = info.SequenceName;
                
                b_val = horzcat(b_val,str2double(extract(seq, digitsPattern)));
                svol = dicomread(f);
                app.vol = cat(3,app.vol,svol);

            
            end
            
            % app.LoadField.Value = horzcat(int2str(count),' Files Loaded');
            app.vol = app.vol(:,:,2:end);
            
            
            app.bvals = unique(b_val);
            
            [app.M, app.N, app.Z] = size(app.vol);
   
            app.instanceSize = app.Z/length(unique(b_val));
            
            app.BValueKnob.Items = string(app.bvals);
            app.BValueKnob.ItemsData = app.bvals;
            app.SliceKnob.Items = string(1:app.instanceSize);
            app.SliceKnob.ItemsData = (1:app.instanceSize);
            
            % imshow(adapthisteq(uint8(app.vol(:,:,1+(app.showBval*app.instanceSize))),'clipLimit',0.005,'Distribution','rayleigh'),'Parent',app.UIAxes);
            imagesc(app.UIAxes,imadjust(app.vol(:,:,1+(app.showBval*app.instanceSize))));
            colormap(app.UIAxes,"gray");
            title(app.UIAxes,'IMAGE');
            xlim(app.UIAxes,[0 ,app.N]);
            ylim(app.UIAxes,[0, app.M]);
            app.DtUpperLimitFreeEditField.Value = app.DtUpperLimitfree;
            app.DpUpperLimitFreeEditField.Value = app.DpUpperLimitfree;
            app.DtUpperLimitSegEditField.Value = app.DtUpperLimitseg;
            app.DpUpperLimitSegEditField.Value = app.DpUpperLimitseg;
            app.DicomState.Value = "Images Loaded";
            app.StateField.Value = "";
            catch
                app.DicomState.Value = "Error Occured";
            end
        end

        % Value changed function: SequencesListBox
        function SequencesListBoxValueChanged(app, event)
            app.sequencevalue = app.SequencesListBox.Items(app.SequencesListBox.Value);
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            I = find(app.patienttable.SeriesDescription == app.sequencevalue);
            I = I(find(I == app.SequencesListBox.Value));
            app.frames = app.patienttable.Frames;
            
            app.FramesEditField.Value = app.frames(I);
        end

        % Value changed function: PatientsListBox
        function PatientsListBoxValueChanged(app, event)
            app.patientvalue = app.PatientsListBox.Value;
            app.SequencesListBox.Items = app.collection.SeriesDescription(find(app.collection.PatientName == app.patientvalue));
            
        end

        % Button pushed function: LoadtoViewerButton
        function LoadtoViewerButtonPushed(app, event)
            app.DicomState.Value = "Loading..";
            pause(0.5);
             try
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            I = find(app.patienttable.SeriesDescription == app.sequencevalue);
            I = I(find(I == app.SequencesListBox.Value));
            
            fileNames = app.patienttable.Filenames(I);
            fileNames = fileNames{1,1};
             % app.rows = app.patienttable.Rows(I);
            % app.cols = app.patienttable.Columns(I);
            app.viewVol = [];
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                svol = dicomread(f);
                app.viewVol = cat(3,app.viewVol,svol);
            end
            [app.viewM, app.viewN, app.viewZ] = size(app.viewVol);
            
            app.viewVol = app.viewVol(:,:,1:end);
            
            imagesc(app.UIAxes8,imadjust(app.viewVol(:,:,app.ViewSlider.Limits(1))));
            colormap(app.UIAxes8,app.ColormapDropDown.Value);
            xlim(app.UIAxes8,[0 ,app.viewN]);
            ylim(app.UIAxes8,[0, app.viewM]);
            if app.viewZ > 1
            app.ViewSlider.Limits = [1,app.viewZ];
            app.ViewSlider.MajorTicks = 1:round(app.viewZ/10):app.viewZ;
            app.ViewSlider.MajorTickLabels = string(1:round(app.viewZ/10):app.viewZ);
            app.ViewSlider.MinorTicks = 1:1:app.viewZ;
            else
                app.ViewSlider.Limits = [0,1];
                app.ViewSlider.MajorTicks = 0:1;
                app.ViewSlider.MajorTickLabels = string(0:1);
                app.ViewSlider.MinorTicks = [];
            end
            app.ViewSlider.Value = 1;
            app.InstanceNumberEditField.Value = app.ViewSlider.Value;
            app.brightnessValue = 0;
            app.contrastLowLimit = 0;
            app.contrastUpLimit = 5;
            app.BritghtnessSlider.Value = 0;
            app.ContrastSlider.Value = [0 5];
            app.DicomState.Value = "Images Loaded";
              catch
                  app.DicomState.Value = "Error Occured";
              end
            
        end

        % Value changed function: ColorMapsDropDown
        function ColorMapsDropDownValueChanged(app, event)
            colormap(app.UIAxes3,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_2,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_3,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_4,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_5,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_6,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_7,app.ColorMapsDropDown.Value);
            
        end

        % Value changed function: ColormapDropDown
        function ColormapDropDownValueChanged(app, event)
            colormap(app.UIAxes8,app.ColormapDropDown.Value);
            
        end

        % Button pushed function: FilterDWIBetaButton
        function FilterDWIBetaButtonPushed(app, event)
            app.StateField.Value = "Filtering Image...";
            pause(0.5);
            app.smoothDegree = app.SmoothingDegreeEditField.Value;
            for idx = 1:1:app.Z

                patch = std2(app.vol(:,:,idx));
                app.vol(:,:,idx) = imnlmfilt(app.vol(:,:,idx),DegreeOfSmoothing=app.smoothDegree*patch);
            end
            app.StateField.Value = "Filtering Done";
        end

        % Button pushed function: DrawMultipleROIButton
        function DrawMultipleROIButtonPushed(app, event)
            app.MultiRoiSwitch.Value = "On"
            app.sigData = [0];
            delete(findall(app.UIAxes,'Type','images.roi'));
            app.combmask = zeros(app.M,app.N);
            app.combmask = double(app.combmask);
            while true
            if app.MultiRoiSwitch.Value == "On"
            if app.ROIStyleDropDown.Value == "Rectangle"
                app.roi = drawrectangle('Parent',app.UIAxes);
                mask = createMask(app.roi);
                app.combmask = imadd(app.combmask,double(mask));
                
                
                
            end
            if app.ROIStyleDropDown.Value == "Polygon"
                app.roi = drawpolygon('Parent',app.UIAxes);
                mask = createMask(app.roi);
                app.combmask = imadd(app.combmask,double(mask));
            end
            if app.ROIStyleDropDown.Value == "Freehand"
                app.roi = drawfreehand('Parent',app.UIAxes);
                mask = createMask(app.roi);
                
                app.combmask = imadd(app.combmask,double(mask));
            end
            if app.ROIStyleDropDown.Value == "Assisted"
                app.roi = drawassisted(findobj(app.UIAxes,'Type','Image'));
                mask = createMask(app.roi);

                app.combmask = imadd(app.combmask,double(mask));
            end
            else
                    pause(.0001)
                    break
            end
            end
            
            app.combmask = imbinarize(app.combmask);
            imshow(app.combmask);
            for i = 0:1:length(app.bvals)-1
                    tempvol = app.vol(:,:,(i+1)+((app.sliceNum-1)*length(app.bvals)));
               
                    meanData = mean(tempvol(app.combmask));
                    app.sigData = horzcat(app.sigData,meanData);
            end
            
            app.sigData = app.sigData(2:end);
            app.s_s0 = app.sigData./app.sigData(1);
            % app.s_s0 = app.s_s0(find(app.s_s0 <= 1));
            % app.s_s0
            % app.bvals = app.bvals(find(app.s_s0 <= 1));
            % app.bvals
            plot(app.bvals,app.s_s0,'Parent',app.UIAxes2);
            title(app.UIAxes2,'Signal Data');
            xlabel(app.UIAxes2,'B-Value');
            ylabel(app.UIAxes2,'S / S0');

        end

        % Button pushed function: LoadToSegmenterButton
        function LoadToSegmenterButtonPushed(app, event)
             app.DicomState.Value = "Loading..";
            pause(0.5);
             try
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            I = find(app.patienttable.SeriesDescription == app.sequencevalue);
            I = I(find(I == app.SequencesListBox.Value));
            
            fileNames = app.patienttable.Filenames(I);
            fileNames = fileNames{1,1};
             % app.rows = app.patienttable.Rows(I);
            % app.cols = app.patienttable.Columns(I);
            app.segmenterVol = [];
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                svol = dicomread(f);
                app.segmenterVol = cat(3,app.segmenterVol,svol);
            end
            [app.segmenterM, app.segmenterN, app.segmenterZ] = size(app.segmenterVol);
            
            app.segmenterVol = app.segmenterVol(:,:,1:end);
            
            imagesc(app.UIAxes9,app.segmenterVol(:,:,1));
            colormap(app.UIAxes9,'gray');
            xlim(app.UIAxes9,[0 ,app.segmenterN]);
            ylim(app.UIAxes9,[0, app.segmenterM]);
            if app.segmenterZ > 1
            app.SegmenterSlider.Limits = [1,app.segmenterZ];
            app.SegmenterSlider.MajorTicks = 1:round(app.segmenterZ/10):app.segmenterZ;
            app.SegmenterSlider.MajorTickLabels = string(1:round(app.segmenterZ/10):app.segmenterZ);
            app.SegmenterSlider.MinorTicks = 1:1:app.segmenterZ;
            else
                app.SegmenterSlider.Limits = [0,1];
                app.SegmenterSlider.MajorTicks = 0:1;
                app.SegmenterSlider.MajorTickLabels = string(0:1);
                app.SegmenterSlider.MinorTicks = [];
            end
            app.samObj = segmentAnythingModel;
            
            app.DicomState.Value = "Images Loaded";
              catch
                  app.DicomState.Value = "Error Occured";
              end
        end

        % Button pushed function: SegmentButton
        function SegmentButtonPushed(app, event)
            app.segmentIm = app.segmenterVol(:,:,round(app.SegmenterSlider.Value));
            app.SegmenterStateField.Value = "Embeddings are Creating";
            pause(0.5);
            app.imgEmbeddings = extractEmbeddings(app.samObj,app.segmenterVol(:,:,round(app.SegmenterSlider.Value)));
            app.SegmenterStateField.Value = "Embeddings Created";
            pause(0.5);
            app.SegmenterStateField.Value = "Creating Masks";
            pause(0.5);
            if ~isnan(app.samBox)
            [app.samMasks, app.samScores] = segmentObjectsFromEmbeddings(app.samObj, app.imgEmbeddings, size(app.segmentIm), BoundingBox = app.samBox);
            else
            [app.samMasks, app.samScores] = segmentObjectsFromEmbeddings(app.samObj, app.imgEmbeddings, size(app.segmentIm), ForegroundPoints = app.forePoints,BackgroundPoints = app.backPoints);
            end
            app.SegmenterStateField.Value = "Masks Created";
            % app.segmenterMasked(~app.segmenterBW) = 0;
            imagesc(app.UIAxes9_2,labeloverlay(imadjust(app.segmentIm),app.samMasks));
            colormap(app.UIAxes9_2,'gray');
            title(app.UIAxes9_2,["Score" + num2str(app.samScores)]);
        end

        % Button pushed function: SegmentROIButton
        function SegmentROIButtonPushed(app, event)
            delete(findall(app.UIAxes9,'Type','images.roi'));
            delete(app.tempax);
            app.segroi = drawrectangle('Parent',app.UIAxes9);
            
            
            
            app.samBox = app.segroi.Position;
            
            
        end

        % Value changing function: ContrastSlider
        function ContrastSliderValueChanging(app, event)
            app.contrastValue = event.Value;
            app.contrastLowLimit = app.contrastValue(1);
            app.contrastUpLimit = app.contrastValue(2);
            
            app.viewAdjustedImage = imadjust(app.viewVol(:,:,round(app.ViewSlider.Value)),[app.contrastLowLimit/255, app.contrastUpLimit/255],[])
            imagesc(app.UIAxes8,app.viewAdjustedImage+app.brightnessValue);
        end

        % Value changing function: ViewSlider
        function ViewSliderValueChanging(app, event)
            app.viewSlidervalue = round(event.Value);
            app.InstanceNumberEditField.Value = app.viewSlidervalue;
            app.viewImage = imadjust(app.viewVol(:,:,app.viewSlidervalue),[app.contrastLowLimit/255, app.contrastUpLimit/255],[]);
            imagesc(app.UIAxes8,app.viewImage + app.brightnessValue);
            colormap(app.UIAxes8,app.ColormapDropDown.Value);
            xlim(app.UIAxes8,[0 ,app.viewN]);
            ylim(app.UIAxes8,[0, app.viewM]);
        end

        % Value changing function: BritghtnessSlider
        function BritghtnessSliderValueChanging(app, event)
            app.brightnessValue = event.Value*(2^16);
            app.viewAdjustedImage = imadjust(app.viewVol(:,:,round(app.ViewSlider.Value)),[app.contrastLowLimit/255, app.contrastUpLimit/255],[])
            imagesc(app.UIAxes8,app.viewAdjustedImage+app.brightnessValue);
            
        end

        % Button pushed function: AddPointButton
        function AddPointButtonPushed(app, event)
            tempAddPoint = drawpoint('Parent',app.UIAxes9,Color = 'green');
            app.forePoints = vertcat(app.forePoints,tempAddPoint.Position);

        end

        % Button pushed function: DropPointButton
        function DropPointButtonPushed(app, event)
            tempDropPoint = drawpoint('Parent',app.UIAxes9,Color = 'red');
            app.backPoints = vertcat(app.backPoints,tempDropPoint.Position);
        end

        % Button pushed function: RegionStatisticsButton
        function RegionStatisticsButtonPushed(app, event)
            stats = regionprops(app.samMasks,'Area','Circularity');
            app.AreaEditField.Value = stats.Area;
          
            app.CircularityEditField.Value = stats.Circularity;
        end

        % Value changing function: SegmenterSlider
        function SegmenterSliderValueChanging(app, event)
            value = round(event.Value);
            imagesc(app.UIAxes9,app.segmenterVol(:,:,value));
            
            xlim(app.UIAxes9,[0 ,app.segmenterN]);
            ylim(app.UIAxes9,[0, app.segmenterM]);
            app.forePoints =[];
            app.backPoints =[];
            app.samBox = [];
            
        end

        % Window scroll wheel function: 
        % DiffusonCalculatorforSIEMENSUIFigure
        function DiffusonCalculatorforSIEMENSUIFigureWindowScrollWheel(app, event)
            if app.selectedTab.Title == "Viewer"
            verticalScrollCount = event.VerticalScrollCount;
            scrollvalue = round(app.ViewSlider.Value) + verticalScrollCount;
            
            
            if scrollvalue < 1
                scrollvalue = 1
            elseif scrollvalue > app.ViewSlider.Limits(2)
                    scrollvalue = app.ViewSlider.Limits(2)
            end
            app.ViewSlider.Value = scrollvalue;
            app.InstanceNumberEditField.Value = scrollvalue;
            scrollImage = imadjust(app.viewVol(:,:,scrollvalue),[app.contrastLowLimit/255, app.contrastUpLimit/255],[]);
            imagesc(app.UIAxes8,scrollImage + app.brightnessValue);
            colormap(app.UIAxes8,app.ColormapDropDown.Value);
            xlim(app.UIAxes8,[0 ,app.viewN]);
            ylim(app.UIAxes8,[0, app.viewM]);
            end
        end

        % Selection change function: TabGroup
        function TabGroupSelectionChanged(app, event)
            app.selectedTab = app.TabGroup.SelectedTab;
            
        end

        % Value changed function: mapRoiSwitch
        function mapRoiSwitchValueChanged(app, event)
            value = app.mapRoiSwitch.Value;
            
            if value == "MultiROI"
                app.GenerateMapsButton.BackgroundColor = [ 0.48,0.54,0.78];
            else
                app.GenerateMapsButton.BackgroundColor = [0.96,0.96,0.96];
            end
        end

        % Button pushed function: LoadtoBOLDPcsButton
        function LoadtoBOLDPcsButtonPushed(app, event)
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            I = find(app.patienttable.SeriesDescription == app.sequencevalue);
            % candidate = convertStringsToChars(descriptions(I));
            I = I(find(I == app.SequencesListBox.Value));
            % if candidate(1) == 'e'
                fileNames = app.patienttable.Filenames(I);
                app.frames = app.collection.Frames;
            % end
            
       
             app.rows = app.patienttable.Rows(I);
            app.cols = app.patienttable.Columns(I);
            app.boldVol = zeros(app.rows , app.cols);
            fileNames = fileNames{1,1};
            TEmatrix = [];
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                info = dicominfo(f);
                TE = info.EchoTime;
                
                TEmatrix = horzcat(TEmatrix,TE);
                svol = dicomread(f);
                app.boldVol = cat(3,app.boldVol,svol);

            
            end
            app.TEvals = unique(TEmatrix);
            app.boldVol = app.boldVol(:,:,2:end);
            [app.bM, app.bN, app.bZ] = size(app.boldVol);
            disp(app.Z)
            disp(app.TEvals)
            app.boldInstace = app.bZ/length(app.TEvals);

            app.TEKnob.Items = string(app.TEvals);
            app.TEKnob.ItemsData = app.TEvals;
            app.BoldSliceKnob.Items = string(1:app.boldInstace);
            app.BoldSliceKnob.ItemsData = (1:app.boldInstace);
            imagesc(app.UIAxes10,imadjust(app.boldVol(:,:,1+(app.showTEval*app.boldInstace))));
            colormap(app.UIAxes10,"gray");
            title(app.UIAxes,'IMAGE');
            xlim(app.UIAxes10,[0,app.bM]);
            ylim(app.UIAxes10,[0,app.bN]);


        end

        % Value changed function: TEKnob
        function TEKnobValueChanged(app, event)
            app.showTEval = find(app.TEvals == app.TEKnob.Value);
            imagesc(app.UIAxes10,imadjust(app.boldVol(:,:,(app.showTEval)+((app.boldSlice-1)*length(app.TEvals)))));
            colormap(app.UIAxes10,"gray");
            xlim(app.UIAxes10,[0,app.bM]);
            ylim(app.UIAxes10,[0,app.bN]);
            
        end

        % Value changed function: BoldSliceKnob
        function BoldSliceKnobValueChanged(app, event)
            app.boldSlice = app.BoldSliceKnob.Value;
            app.showTEval = find(app.TEvals == app.TEKnob.Value);
            imagesc(app.UIAxes10,imadjust(app.boldVol(:,:,(app.showTEval)+((app.boldSlice-1)*length(app.TEvals)))));
            colormap(app.UIAxes10,"gray");
            xlim(app.UIAxes10,[0,app.bM]);
            ylim(app.UIAxes10,[0,app.bN]);
        end

        % Button pushed function: GenerateR2MapButton
        function GenerateR2MapButtonPushed(app, event)
            
            z = app.Z;
            app.R2star = zeros(app.M,app.N);
            [Xmask,Ymask] = size(app.boldmask);
            app.BoldState.Value = "Proccesing";
            pause(0.5)
            for k = 1:1:Xmask
                for j = 1:1:Ymask
                    boldmatrix = [0];
                    if app.boldmask(k,j) ~=0 
                    for m = 0:1:length(app.TEvals)-1
                        
                            boldpixel = app.boldVol(k,j,(m+1)+((app.boldSlice-1)*length(app.TEvals)));
                            boldmatrix = horzcat(boldmatrix,boldpixel);
                        
                    end
                    try
                    boldmatrix = boldmatrix(2:end);
                    boldmatrix = boldmatrix./boldmatrix(1);
                    boldfit = fit(app.TEvals.',boldmatrix.',app.boldfittype,app.boldopts);
                    
                    app.R2star(k,j) = boldfit.b*1000;
                    catch
                        app.R2star(k,j) = 0;
                    end
                    else
                        app.R2star(k,j) = 0;
                    end
                end
            end
           
            
            
            
            app.BoldState.Value = "Done";
            imagesc(app.UIAxes11,app.R2star,[0, max(app.R2star(:))]);
            colormap(app.UIAxes11,'turbo');
            colorbar(app.UIAxes11);
            xlim(app.UIAxes11,[0,app.bM]);
            ylim(app.UIAxes11,[0,app.bN]);
            app.BoldState.Value = "Done";
        end

        % Value changed function: DrawROIforFitButton
        function DrawROIforFitButtonValueChanged(app, event)
            delete(findall(app.UIAxes10,'Type','images.roi'));
            boldroi = drawassisted(findobj(app.UIAxes10,'Type','Image'));
            app.boldmask = createMask(boldroi);
            app.boldfittype = fittype( 'exp(-b*x)', 'independent', 'x', 'dependent', 'y' );
            app.boldopts =fitoptions( 'Method', 'NonlinearLeastSquares','Algorithm','Levenberg-Marquardt');
            app.boldopts.Display = 'Off';
            app.boldopts.StartPoint = 1;
            boldFitsig = [0];
            for i=1:1:length(app.TEvals)
                fittempvol = app.boldVol(:,:,i+((app.boldSlice-1)*length(app.TEvals)));
                fitBoldMean = mean(fittempvol(app.boldmask));
                boldFitsig =horzcat(boldFitsig,fitBoldMean);
            end
            boldFitsig = boldFitsig(2:end);
            boldFitsig = boldFitsig./boldFitsig(1);
            roifit = fit(app.TEvals.',boldFitsig.',app.boldfittype,app.boldopts);
            app.MeanR2S1EditField.Value = roifit.b;
           
            plot(app.UIAxes12,roifit,app.TEvals,boldFitsig);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create DiffusonCalculatorforSIEMENSUIFigure and hide until all components are created
            app.DiffusonCalculatorforSIEMENSUIFigure = uifigure('Visible', 'off');
            app.DiffusonCalculatorforSIEMENSUIFigure.Position = [100 100 1600 900];
            app.DiffusonCalculatorforSIEMENSUIFigure.Name = 'Diffuson Calculator for SIEMENS';
            app.DiffusonCalculatorforSIEMENSUIFigure.WindowScrollWheelFcn = createCallbackFcn(app, @DiffusonCalculatorforSIEMENSUIFigureWindowScrollWheel, true);

            % Create TabGroup
            app.TabGroup = uitabgroup(app.DiffusonCalculatorforSIEMENSUIFigure);
            app.TabGroup.SelectionChangedFcn = createCallbackFcn(app, @TabGroupSelectionChanged, true);
            app.TabGroup.Position = [34 25 1546 850];

            % Create CalculatorTab
            app.CalculatorTab = uitab(app.TabGroup);
            app.CalculatorTab.Title = 'Calculator';

            % Create FittingGraphicsPanel
            app.FittingGraphicsPanel = uipanel(app.CalculatorTab);
            app.FittingGraphicsPanel.BorderType = 'none';
            app.FittingGraphicsPanel.Title = 'Fitting Graphics';
            app.FittingGraphicsPanel.Position = [1025 361 511 437];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.FittingGraphicsPanel);
            title(app.UIAxes2, 'Title')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [2 49 497 352];

            % Create RSquareEditFieldLabel
            app.RSquareEditFieldLabel = uilabel(app.FittingGraphicsPanel);
            app.RSquareEditFieldLabel.HorizontalAlignment = 'right';
            app.RSquareEditFieldLabel.Position = [41 14 56 22];
            app.RSquareEditFieldLabel.Text = 'R-Square';

            % Create RSquareEditField
            app.RSquareEditField = uieditfield(app.FittingGraphicsPanel, 'numeric');
            app.RSquareEditField.ValueDisplayFormat = '%.4f';
            app.RSquareEditField.Position = [112 13 54 24];

            % Create ImageViewPanel
            app.ImageViewPanel = uipanel(app.CalculatorTab);
            app.ImageViewPanel.BorderType = 'none';
            app.ImageViewPanel.Title = 'Image View';
            app.ImageViewPanel.Position = [76 306 940 492];

            % Create UIAxes
            app.UIAxes = uiaxes(app.ImageViewPanel);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [116 22 667 431];

            % Create OptionsPanel
            app.OptionsPanel = uipanel(app.CalculatorTab);
            app.OptionsPanel.BorderType = 'none';
            app.OptionsPanel.Title = 'Options';
            app.OptionsPanel.Position = [1015 170 518 138];

            % Create DtEditFieldLabel
            app.DtEditFieldLabel = uilabel(app.OptionsPanel);
            app.DtEditFieldLabel.HorizontalAlignment = 'right';
            app.DtEditFieldLabel.Position = [149 91 25 22];
            app.DtEditFieldLabel.Text = 'Dt';

            % Create DtEditField
            app.DtEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DtEditField.ValueDisplayFormat = '%.8f';
            app.DtEditField.Position = [185 93 59 17];

            % Create DpEditFieldLabel
            app.DpEditFieldLabel = uilabel(app.OptionsPanel);
            app.DpEditFieldLabel.HorizontalAlignment = 'right';
            app.DpEditFieldLabel.Position = [149 62 25 22];
            app.DpEditFieldLabel.Text = 'Dp';

            % Create DpEditField
            app.DpEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DpEditField.ValueDisplayFormat = '%.8f';
            app.DpEditField.Position = [185 64 59 18];

            % Create pFEditFieldLabel
            app.pFEditFieldLabel = uilabel(app.OptionsPanel);
            app.pFEditFieldLabel.HorizontalAlignment = 'right';
            app.pFEditFieldLabel.Position = [149 34 25 22];
            app.pFEditFieldLabel.Text = 'pF';

            % Create pFEditField
            app.pFEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.pFEditField.ValueDisplayFormat = '%.8f';
            app.pFEditField.Position = [185 35 59 19];

            % Create ADCEditFieldLabel
            app.ADCEditFieldLabel = uilabel(app.OptionsPanel);
            app.ADCEditFieldLabel.HorizontalAlignment = 'right';
            app.ADCEditFieldLabel.Position = [144 8 30 22];
            app.ADCEditFieldLabel.Text = 'ADC';

            % Create ADCEditField
            app.ADCEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.ADCEditField.ValueDisplayFormat = '%.8f';
            app.ADCEditField.Position = [185 10 59 18];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.OptionsPanel, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Icon = fullfile(pathToMLAPP, 'icons', 'calculate.png');
            app.CalculateButton.Position = [414 -4 98 35];
            app.CalculateButton.Text = 'Calculate';

            % Create ROIStyleDropDownLabel
            app.ROIStyleDropDownLabel = uilabel(app.OptionsPanel);
            app.ROIStyleDropDownLabel.HorizontalAlignment = 'right';
            app.ROIStyleDropDownLabel.Position = [354 83 56 22];
            app.ROIStyleDropDownLabel.Text = 'ROI Style';

            % Create ROIStyleDropDown
            app.ROIStyleDropDown = uidropdown(app.OptionsPanel);
            app.ROIStyleDropDown.Items = {'Rectangle', 'Polygon', 'Freehand', 'Assisted'};
            app.ROIStyleDropDown.Position = [425 83 84 20];
            app.ROIStyleDropDown.Value = 'Rectangle';

            % Create DrawROIButton
            app.DrawROIButton = uibutton(app.OptionsPanel, 'state');
            app.DrawROIButton.ValueChangedFcn = createCallbackFcn(app, @DrawROIButtonValueChanged, true);
            app.DrawROIButton.Icon = fullfile(pathToMLAPP, 'icons', 'draw.png');
            app.DrawROIButton.Text = 'DrawROI';
            app.DrawROIButton.Position = [420 43 94 35];

            % Create DrawMultipleROIButton
            app.DrawMultipleROIButton = uibutton(app.OptionsPanel, 'push');
            app.DrawMultipleROIButton.ButtonPushedFcn = createCallbackFcn(app, @DrawMultipleROIButtonPushed, true);
            app.DrawMultipleROIButton.Icon = fullfile(pathToMLAPP, 'icons', 'multidraw.png');
            app.DrawMultipleROIButton.Position = [256 42 158 35];
            app.DrawMultipleROIButton.Text = 'Draw Multiple ROI';

            % Create MultiRoiSwitch
            app.MultiRoiSwitch = uiswitch(app.OptionsPanel, 'slider');
            app.MultiRoiSwitch.Position = [325 14 45 20];

            % Create LoadDICOMButton
            app.LoadDICOMButton = uibutton(app.CalculatorTab, 'push');
            app.LoadDICOMButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDICOMButtonPushed, true);
            app.LoadDICOMButton.Icon = fullfile(pathToMLAPP, 'icons', 'import.png');
            app.LoadDICOMButton.Position = [1394 15 140 38];
            app.LoadDICOMButton.Text = 'Load DICOM';

            % Create AlgorithmButtonGroup
            app.AlgorithmButtonGroup = uibuttongroup(app.CalculatorTab);
            app.AlgorithmButtonGroup.BorderType = 'none';
            app.AlgorithmButtonGroup.BorderWidth = 0;
            app.AlgorithmButtonGroup.Title = 'Algorithm';
            app.AlgorithmButtonGroup.Position = [1021 169 100 110];

            % Create MonoExpButton
            app.MonoExpButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.MonoExpButton.Text = 'MonoExp';
            app.MonoExpButton.Position = [8 66 73 22];
            app.MonoExpButton.Value = true;

            % Create BiexpFreeButton
            app.BiexpFreeButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.BiexpFreeButton.Text = 'BiexpFree';
            app.BiexpFreeButton.Position = [8 45 77 22];

            % Create BiExpSegButton
            app.BiExpSegButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.BiExpSegButton.Text = 'BiExpSeg';
            app.BiExpSegButton.Position = [8 22 75 22];

            % Create EditField
            app.EditField = uieditfield(app.CalculatorTab, 'numeric');
            app.EditField.Position = [395 80 38 19];

            % Create EditField_2
            app.EditField_2 = uieditfield(app.CalculatorTab, 'numeric');
            app.EditField_2.Position = [597 80 38 19];

            % Create GenerateMapsButton
            app.GenerateMapsButton = uibutton(app.CalculatorTab, 'push');
            app.GenerateMapsButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateMapsButtonPushed, true);
            app.GenerateMapsButton.Icon = fullfile(pathToMLAPP, 'icons', 'generate.png');
            app.GenerateMapsButton.Position = [979 57 160 51];
            app.GenerateMapsButton.Text = 'Generate Maps';

            % Create BValueKnobLabel
            app.BValueKnobLabel = uilabel(app.CalculatorTab);
            app.BValueKnobLabel.HorizontalAlignment = 'center';
            app.BValueKnobLabel.FontSize = 10;
            app.BValueKnobLabel.FontWeight = 'bold';
            app.BValueKnobLabel.Position = [393 112 42 22];
            app.BValueKnobLabel.Text = 'B-Value';

            % Create BValueKnob
            app.BValueKnob = uiknob(app.CalculatorTab, 'discrete');
            app.BValueKnob.ValueChangedFcn = createCallbackFcn(app, @BValueKnobValueChanged, true);
            app.BValueKnob.FontSize = 10;
            app.BValueKnob.FontWeight = 'bold';
            app.BValueKnob.Position = [363 149 99 99];

            % Create SliceKnobLabel
            app.SliceKnobLabel = uilabel(app.CalculatorTab);
            app.SliceKnobLabel.HorizontalAlignment = 'center';
            app.SliceKnobLabel.FontSize = 10;
            app.SliceKnobLabel.FontWeight = 'bold';
            app.SliceKnobLabel.Position = [600 118 28 22];
            app.SliceKnobLabel.Text = 'Slice';

            % Create SliceKnob
            app.SliceKnob = uiknob(app.CalculatorTab, 'discrete');
            app.SliceKnob.ValueChangedFcn = createCallbackFcn(app, @SliceKnobValueChanged, true);
            app.SliceKnob.FontSize = 10;
            app.SliceKnob.FontWeight = 'bold';
            app.SliceKnob.Position = [567 152 95 95];

            % Create KernelSizeEditFieldLabel
            app.KernelSizeEditFieldLabel = uilabel(app.CalculatorTab);
            app.KernelSizeEditFieldLabel.HorizontalAlignment = 'right';
            app.KernelSizeEditFieldLabel.Position = [991 112 66 22];
            app.KernelSizeEditFieldLabel.Text = 'Kernel Size';

            % Create KernelSizeEditField
            app.KernelSizeEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.KernelSizeEditField.Position = [1099 115 38 17];
            app.KernelSizeEditField.Value = 4;

            % Create DpUpperLimitFreeEditFieldLabel
            app.DpUpperLimitFreeEditFieldLabel = uilabel(app.CalculatorTab);
            app.DpUpperLimitFreeEditFieldLabel.HorizontalAlignment = 'right';
            app.DpUpperLimitFreeEditFieldLabel.Position = [1140 147 113 22];
            app.DpUpperLimitFreeEditFieldLabel.Text = 'Dp Upper Limit Free';

            % Create DpUpperLimitFreeEditField
            app.DpUpperLimitFreeEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.DpUpperLimitFreeEditField.ValueDisplayFormat = '%.5f';
            app.DpUpperLimitFreeEditField.Position = [1268 146 64 23];

            % Create DtUpperLimitFreeEditFieldLabel
            app.DtUpperLimitFreeEditFieldLabel = uilabel(app.CalculatorTab);
            app.DtUpperLimitFreeEditFieldLabel.HorizontalAlignment = 'right';
            app.DtUpperLimitFreeEditFieldLabel.Position = [1142 119 110 22];
            app.DtUpperLimitFreeEditFieldLabel.Text = 'Dt Upper Limit Free';

            % Create DtUpperLimitFreeEditField
            app.DtUpperLimitFreeEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.DtUpperLimitFreeEditField.ValueDisplayFormat = '%.5f';
            app.DtUpperLimitFreeEditField.Position = [1267 118 64 23];

            % Create StateField
            app.StateField = uieditfield(app.CalculatorTab, 'text');
            app.StateField.FontSize = 14;
            app.StateField.BackgroundColor = [0.902 0.902 0.902];
            app.StateField.Position = [15 15 427 22];

            % Create FilterDWIBetaButton
            app.FilterDWIBetaButton = uibutton(app.CalculatorTab, 'push');
            app.FilterDWIBetaButton.ButtonPushedFcn = createCallbackFcn(app, @FilterDWIBetaButtonPushed, true);
            app.FilterDWIBetaButton.Icon = fullfile(pathToMLAPP, 'icons', 'filter.png');
            app.FilterDWIBetaButton.Position = [1394 73 140 43];
            app.FilterDWIBetaButton.Text = 'Filter DWI (Beta)';

            % Create mapRoiSwitch
            app.mapRoiSwitch = uiswitch(app.CalculatorTab, 'slider');
            app.mapRoiSwitch.Items = {'SingleROI', 'MultiROI'};
            app.mapRoiSwitch.ValueChangedFcn = createCallbackFcn(app, @mapRoiSwitchValueChanged, true);
            app.mapRoiSwitch.Position = [1041 30 45 20];
            app.mapRoiSwitch.Value = 'SingleROI';

            % Create DpUpperLimitSegEditFieldLabel
            app.DpUpperLimitSegEditFieldLabel = uilabel(app.CalculatorTab);
            app.DpUpperLimitSegEditFieldLabel.HorizontalAlignment = 'right';
            app.DpUpperLimitSegEditFieldLabel.Position = [1142 89 113 22];
            app.DpUpperLimitSegEditFieldLabel.Text = 'Dp Upper Limit Seg ';

            % Create DpUpperLimitSegEditField
            app.DpUpperLimitSegEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.DpUpperLimitSegEditField.ValueDisplayFormat = '%.5f';
            app.DpUpperLimitSegEditField.Position = [1270 89 61 20];

            % Create DtUpperLimitSegEditFieldLabel
            app.DtUpperLimitSegEditFieldLabel = uilabel(app.CalculatorTab);
            app.DtUpperLimitSegEditFieldLabel.HorizontalAlignment = 'right';
            app.DtUpperLimitSegEditFieldLabel.Position = [1145 57 106 22];
            app.DtUpperLimitSegEditFieldLabel.Text = 'DtUpper Limit Seg ';

            % Create DtUpperLimitSegEditField
            app.DtUpperLimitSegEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.DtUpperLimitSegEditField.ValueDisplayFormat = '%.5f';
            app.DtUpperLimitSegEditField.Position = [1266 57 66 20];

            % Create SmoothingDegreeEditFieldLabel
            app.SmoothingDegreeEditFieldLabel = uilabel(app.CalculatorTab);
            app.SmoothingDegreeEditFieldLabel.HorizontalAlignment = 'right';
            app.SmoothingDegreeEditFieldLabel.Position = [1346 131 105 22];
            app.SmoothingDegreeEditFieldLabel.Text = 'Smoothing Degree';

            % Create SmoothingDegreeEditField
            app.SmoothingDegreeEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.SmoothingDegreeEditField.Position = [1475 131 51 22];

            % Create DicomTab
            app.DicomTab = uitab(app.TabGroup);
            app.DicomTab.Title = 'Dicom';

            % Create SequencesPanel
            app.SequencesPanel = uipanel(app.DicomTab);
            app.SequencesPanel.Title = 'Sequences';
            app.SequencesPanel.Position = [1 66 1513 364];

            % Create SequencesListBoxLabel
            app.SequencesListBoxLabel = uilabel(app.SequencesPanel);
            app.SequencesListBoxLabel.HorizontalAlignment = 'right';
            app.SequencesListBoxLabel.Position = [2 308 65 22];
            app.SequencesListBoxLabel.Text = 'Sequences';

            % Create SequencesListBox
            app.SequencesListBox = uilistbox(app.SequencesPanel);
            app.SequencesListBox.ValueChangedFcn = createCallbackFcn(app, @SequencesListBoxValueChanged, true);
            app.SequencesListBox.Position = [82 34 1238 298];

            % Create LoadtoCalculatorButton
            app.LoadtoCalculatorButton = uibutton(app.SequencesPanel, 'push');
            app.LoadtoCalculatorButton.ButtonPushedFcn = createCallbackFcn(app, @LoadtoCalculatorButtonPushed, true);
            app.LoadtoCalculatorButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.LoadtoCalculatorButton.Position = [1345 25 138 26];
            app.LoadtoCalculatorButton.Text = 'Load to Calculator';

            % Create LoadtoViewerButton
            app.LoadtoViewerButton = uibutton(app.SequencesPanel, 'push');
            app.LoadtoViewerButton.ButtonPushedFcn = createCallbackFcn(app, @LoadtoViewerButtonPushed, true);
            app.LoadtoViewerButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.LoadtoViewerButton.Position = [1357 111 120 26];
            app.LoadtoViewerButton.Text = 'Load to Viewer';

            % Create FramesEditFieldLabel
            app.FramesEditFieldLabel = uilabel(app.SequencesPanel);
            app.FramesEditFieldLabel.HorizontalAlignment = 'right';
            app.FramesEditFieldLabel.Position = [1340 282 46 22];
            app.FramesEditFieldLabel.Text = 'Frames';

            % Create FramesEditField
            app.FramesEditField = uieditfield(app.SequencesPanel, 'numeric');
            app.FramesEditField.Position = [1401 277 93 32];

            % Create LoadToSegmenterButton
            app.LoadToSegmenterButton = uibutton(app.SequencesPanel, 'push');
            app.LoadToSegmenterButton.ButtonPushedFcn = createCallbackFcn(app, @LoadToSegmenterButtonPushed, true);
            app.LoadToSegmenterButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.LoadToSegmenterButton.Position = [1345 66 138 26];
            app.LoadToSegmenterButton.Text = 'Load To Segmenter';

            % Create LoadtoBOLDPcsButton
            app.LoadtoBOLDPcsButton = uibutton(app.SequencesPanel, 'push');
            app.LoadtoBOLDPcsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadtoBOLDPcsButtonPushed, true);
            app.LoadtoBOLDPcsButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.LoadtoBOLDPcsButton.Position = [1344 144 149 24];
            app.LoadtoBOLDPcsButton.Text = 'Load to BOLD Pcs.';

            % Create PatientsPanel
            app.PatientsPanel = uipanel(app.DicomTab);
            app.PatientsPanel.Title = 'Patients';
            app.PatientsPanel.Position = [1 466 1512 360];

            % Create PatientsListBoxLabel
            app.PatientsListBoxLabel = uilabel(app.PatientsPanel);
            app.PatientsListBoxLabel.HorizontalAlignment = 'right';
            app.PatientsListBoxLabel.Position = [13 307 48 22];
            app.PatientsListBoxLabel.Text = 'Patients';

            % Create PatientsListBox
            app.PatientsListBox = uilistbox(app.PatientsPanel);
            app.PatientsListBox.ValueChangedFcn = createCallbackFcn(app, @PatientsListBoxValueChanged, true);
            app.PatientsListBox.Position = [76 12 1418 319];

            % Create DicomState
            app.DicomState = uieditfield(app.DicomTab, 'text');
            app.DicomState.Position = [15 21 491 21];

            % Create ViewerTab
            app.ViewerTab = uitab(app.TabGroup);
            app.ViewerTab.Title = 'Viewer';

            % Create ImageViewPanel_2
            app.ImageViewPanel_2 = uipanel(app.ViewerTab);
            app.ImageViewPanel_2.Title = 'Image View';
            app.ImageViewPanel_2.Position = [1 35 1525 789];

            % Create UIAxes8
            app.UIAxes8 = uiaxes(app.ImageViewPanel_2);
            title(app.UIAxes8, 'IMAGE')
            xlabel(app.UIAxes8, 'X')
            ylabel(app.UIAxes8, 'Y')
            zlabel(app.UIAxes8, 'Z')
            app.UIAxes8.Position = [82 134 1353 611];

            % Create ViewSliderLabel
            app.ViewSliderLabel = uilabel(app.ImageViewPanel_2);
            app.ViewSliderLabel.HorizontalAlignment = 'right';
            app.ViewSliderLabel.Position = [421 76 62 22];
            app.ViewSliderLabel.Text = 'ViewSlider';

            % Create ViewSlider
            app.ViewSlider = uislider(app.ImageViewPanel_2);
            app.ViewSlider.Limits = [1 100];
            app.ViewSlider.ValueChangingFcn = createCallbackFcn(app, @ViewSliderValueChanging, true);
            app.ViewSlider.Position = [504 85 672 3];
            app.ViewSlider.Value = 1;

            % Create InstanceNumberEditFieldLabel
            app.InstanceNumberEditFieldLabel = uilabel(app.ImageViewPanel_2);
            app.InstanceNumberEditFieldLabel.HorizontalAlignment = 'right';
            app.InstanceNumberEditFieldLabel.Position = [469 10 96 22];
            app.InstanceNumberEditFieldLabel.Text = 'Instance Number';

            % Create InstanceNumberEditField
            app.InstanceNumberEditField = uieditfield(app.ImageViewPanel_2, 'numeric');
            app.InstanceNumberEditField.Position = [580 9 41 23];

            % Create ColormapDropDownLabel
            app.ColormapDropDownLabel = uilabel(app.ImageViewPanel_2);
            app.ColormapDropDownLabel.HorizontalAlignment = 'right';
            app.ColormapDropDownLabel.Position = [1340 62 57 22];
            app.ColormapDropDownLabel.Text = 'Colormap';

            % Create ColormapDropDown
            app.ColormapDropDown = uidropdown(app.ImageViewPanel_2);
            app.ColormapDropDown.ValueChangedFcn = createCallbackFcn(app, @ColormapDropDownValueChanged, true);
            app.ColormapDropDown.Position = [1412 62 100 22];

            % Create ContrastSliderLabel
            app.ContrastSliderLabel = uilabel(app.ImageViewPanel_2);
            app.ContrastSliderLabel.HorizontalAlignment = 'right';
            app.ContrastSliderLabel.Position = [1455 109 50 22];
            app.ContrastSliderLabel.Text = 'Contrast';

            % Create ContrastSlider
            app.ContrastSlider = uislider(app.ImageViewPanel_2, 'range');
            app.ContrastSlider.Limits = [0 5];
            app.ContrastSlider.Orientation = 'vertical';
            app.ContrastSlider.ValueChangingFcn = createCallbackFcn(app, @ContrastSliderValueChanging, true);
            app.ContrastSlider.Position = [1480 152 3 569];
            app.ContrastSlider.Value = [0 2.5];

            % Create BritghtnessSliderLabel
            app.BritghtnessSliderLabel = uilabel(app.ImageViewPanel_2);
            app.BritghtnessSliderLabel.HorizontalAlignment = 'right';
            app.BritghtnessSliderLabel.Position = [5 110 65 22];
            app.BritghtnessSliderLabel.Text = 'Britghtness';

            % Create BritghtnessSlider
            app.BritghtnessSlider = uislider(app.ImageViewPanel_2);
            app.BritghtnessSlider.Limits = [-1 1];
            app.BritghtnessSlider.Orientation = 'vertical';
            app.BritghtnessSlider.ValueChangingFcn = createCallbackFcn(app, @BritghtnessSliderValueChanging, true);
            app.BritghtnessSlider.Position = [33 151 3 580];

            % Create MapsTab
            app.MapsTab = uitab(app.TabGroup);
            app.MapsTab.Title = 'Maps';

            % Create DWIMAPSPanel
            app.DWIMAPSPanel = uipanel(app.MapsTab);
            app.DWIMAPSPanel.Title = 'DWI MAPS';
            app.DWIMAPSPanel.Position = [1 35 1525 790];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3, 'DP Free')
            xlabel(app.UIAxes3, 'X')
            ylabel(app.UIAxes3, 'Y')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Position = [22 532 409 228];

            % Create UIAxes3_2
            app.UIAxes3_2 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_2, 'DT Free')
            xlabel(app.UIAxes3_2, 'X')
            ylabel(app.UIAxes3_2, 'Y')
            zlabel(app.UIAxes3_2, 'Z')
            app.UIAxes3_2.Position = [535 532 409 228];

            % Create UIAxes3_3
            app.UIAxes3_3 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_3, 'pF Free')
            xlabel(app.UIAxes3_3, 'X')
            ylabel(app.UIAxes3_3, 'Y')
            zlabel(app.UIAxes3_3, 'Z')
            app.UIAxes3_3.Position = [1055 533 409 228];

            % Create UIAxes3_4
            app.UIAxes3_4 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_4, 'DP Seg-Fit')
            xlabel(app.UIAxes3_4, 'X')
            ylabel(app.UIAxes3_4, 'Y')
            zlabel(app.UIAxes3_4, 'Z')
            app.UIAxes3_4.Position = [22 271 409 228];

            % Create UIAxes3_5
            app.UIAxes3_5 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_5, 'DT Seg-Fit')
            xlabel(app.UIAxes3_5, 'X')
            ylabel(app.UIAxes3_5, 'Y')
            zlabel(app.UIAxes3_5, 'Z')
            app.UIAxes3_5.Position = [535 270 409 228];

            % Create UIAxes3_6
            app.UIAxes3_6 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_6, 'pF Seg-Fit')
            xlabel(app.UIAxes3_6, 'X')
            ylabel(app.UIAxes3_6, 'Y')
            zlabel(app.UIAxes3_6, 'Z')
            app.UIAxes3_6.Position = [1055 261 409 228];

            % Create UIAxes3_7
            app.UIAxes3_7 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_7, 'ADC')
            xlabel(app.UIAxes3_7, 'X')
            ylabel(app.UIAxes3_7, 'Y')
            zlabel(app.UIAxes3_7, 'Z')
            app.UIAxes3_7.Position = [535 8 409 228];

            % Create ColorMapsDropDownLabel
            app.ColorMapsDropDownLabel = uilabel(app.DWIMAPSPanel);
            app.ColorMapsDropDownLabel.HorizontalAlignment = 'right';
            app.ColorMapsDropDownLabel.Position = [1328 22 66 22];
            app.ColorMapsDropDownLabel.Text = 'Color Maps';

            % Create ColorMapsDropDown
            app.ColorMapsDropDown = uidropdown(app.DWIMAPSPanel);
            app.ColorMapsDropDown.ValueChangedFcn = createCallbackFcn(app, @ColorMapsDropDownValueChanged, true);
            app.ColorMapsDropDown.Position = [1409 22 100 22];

            % Create SegmenterBetaTab
            app.SegmenterBetaTab = uitab(app.TabGroup);
            app.SegmenterBetaTab.Title = 'Segmenter(Beta)';

            % Create UIAxes9
            app.UIAxes9 = uiaxes(app.SegmenterBetaTab);
            title(app.UIAxes9, 'Title')
            xlabel(app.UIAxes9, 'X')
            ylabel(app.UIAxes9, 'Y')
            zlabel(app.UIAxes9, 'Z')
            app.UIAxes9.Position = [36 165 689 633];

            % Create UIAxes9_2
            app.UIAxes9_2 = uiaxes(app.SegmenterBetaTab);
            title(app.UIAxes9_2, 'Title')
            xlabel(app.UIAxes9_2, 'X')
            ylabel(app.UIAxes9_2, 'Y')
            zlabel(app.UIAxes9_2, 'Z')
            app.UIAxes9_2.Position = [818 165 689 633];

            % Create SegmenterSliderLabel
            app.SegmenterSliderLabel = uilabel(app.SegmenterBetaTab);
            app.SegmenterSliderLabel.HorizontalAlignment = 'right';
            app.SegmenterSliderLabel.Position = [37 88 98 22];
            app.SegmenterSliderLabel.Text = 'Segmenter Slider';

            % Create SegmenterSlider
            app.SegmenterSlider = uislider(app.SegmenterBetaTab);
            app.SegmenterSlider.ValueChangingFcn = createCallbackFcn(app, @SegmenterSliderValueChanging, true);
            app.SegmenterSlider.Position = [156 97 557 3];

            % Create SegmentButton
            app.SegmentButton = uibutton(app.SegmenterBetaTab, 'push');
            app.SegmentButton.ButtonPushedFcn = createCallbackFcn(app, @SegmentButtonPushed, true);
            app.SegmentButton.Icon = fullfile(pathToMLAPP, 'icons', 'segment.png');
            app.SegmentButton.Position = [734 73 124 35];
            app.SegmentButton.Text = 'Segment';

            % Create SegmentROIButton
            app.SegmentROIButton = uibutton(app.SegmenterBetaTab, 'push');
            app.SegmentROIButton.ButtonPushedFcn = createCallbackFcn(app, @SegmentROIButtonPushed, true);
            app.SegmentROIButton.Icon = fullfile(pathToMLAPP, 'icons', 'draw.png');
            app.SegmentROIButton.Position = [737 22 121 36];
            app.SegmentROIButton.Text = 'Segment ROI';

            % Create SegmenterStateField
            app.SegmenterStateField = uieditfield(app.SegmenterBetaTab, 'text');
            app.SegmenterStateField.Position = [15 9 469 27];

            % Create AddPointButton
            app.AddPointButton = uibutton(app.SegmenterBetaTab, 'push');
            app.AddPointButton.ButtonPushedFcn = createCallbackFcn(app, @AddPointButtonPushed, true);
            app.AddPointButton.BackgroundColor = [0.3922 0.8314 0.0745];
            app.AddPointButton.Position = [734 773 82 23];
            app.AddPointButton.Text = 'Add Point';

            % Create DropPointButton
            app.DropPointButton = uibutton(app.SegmenterBetaTab, 'push');
            app.DropPointButton.ButtonPushedFcn = createCallbackFcn(app, @DropPointButtonPushed, true);
            app.DropPointButton.BackgroundColor = [1 0 0];
            app.DropPointButton.Position = [734 743 82 23];
            app.DropPointButton.Text = 'Drop Point';

            % Create RegionStatisticsButton
            app.RegionStatisticsButton = uibutton(app.SegmenterBetaTab, 'push');
            app.RegionStatisticsButton.ButtonPushedFcn = createCallbackFcn(app, @RegionStatisticsButtonPushed, true);
            app.RegionStatisticsButton.Icon = fullfile(pathToMLAPP, 'icons', 'stats.png');
            app.RegionStatisticsButton.Position = [912 77 154 31];
            app.RegionStatisticsButton.Text = 'Region Statistics';

            % Create AreaEditFieldLabel
            app.AreaEditFieldLabel = uilabel(app.SegmenterBetaTab);
            app.AreaEditFieldLabel.HorizontalAlignment = 'right';
            app.AreaEditFieldLabel.Position = [1129 86 30 22];
            app.AreaEditFieldLabel.Text = 'Area';

            % Create AreaEditField
            app.AreaEditField = uieditfield(app.SegmenterBetaTab, 'numeric');
            app.AreaEditField.Position = [1174 86 100 22];

            % Create CircularityEditFieldLabel
            app.CircularityEditFieldLabel = uilabel(app.SegmenterBetaTab);
            app.CircularityEditFieldLabel.HorizontalAlignment = 'right';
            app.CircularityEditFieldLabel.Position = [1102 52 58 22];
            app.CircularityEditFieldLabel.Text = 'Circularity';

            % Create CircularityEditField
            app.CircularityEditField = uieditfield(app.SegmenterBetaTab, 'numeric');
            app.CircularityEditField.Position = [1175 52 100 22];

            % Create BOLD_Proceesor
            app.BOLD_Proceesor = uitab(app.TabGroup);
            app.BOLD_Proceesor.Title = 'BOLD Processor';

            % Create R2ImagesPanel
            app.R2ImagesPanel = uipanel(app.BOLD_Proceesor);
            app.R2ImagesPanel.Title = 'R2*Images';
            app.R2ImagesPanel.Position = [36 308 626 493];

            % Create UIAxes10
            app.UIAxes10 = uiaxes(app.R2ImagesPanel);
            title(app.UIAxes10, 'Title')
            xlabel(app.UIAxes10, 'X')
            ylabel(app.UIAxes10, 'Y')
            zlabel(app.UIAxes10, 'Z')
            app.UIAxes10.Position = [6 18 607 453];

            % Create R2MapandFittingPanel
            app.R2MapandFittingPanel = uipanel(app.BOLD_Proceesor);
            app.R2MapandFittingPanel.Title = 'R2* Map and Fitting';
            app.R2MapandFittingPanel.Position = [703 23 779 778];

            % Create UIAxes11
            app.UIAxes11 = uiaxes(app.R2MapandFittingPanel);
            title(app.UIAxes11, 'Title')
            xlabel(app.UIAxes11, 'X')
            ylabel(app.UIAxes11, 'Y')
            zlabel(app.UIAxes11, 'Z')
            app.UIAxes11.Position = [20 229 734 512];

            % Create UIAxes12
            app.UIAxes12 = uiaxes(app.R2MapandFittingPanel);
            title(app.UIAxes12, 'Title')
            xlabel(app.UIAxes12, 'X')
            ylabel(app.UIAxes12, 'Y')
            zlabel(app.UIAxes12, 'Z')
            app.UIAxes12.Position = [268 23 280 201];

            % Create MeanR2S1EditFieldLabel
            app.MeanR2S1EditFieldLabel = uilabel(app.R2MapandFittingPanel);
            app.MeanR2S1EditFieldLabel.HorizontalAlignment = 'right';
            app.MeanR2S1EditFieldLabel.Position = [554 13 105 22];
            app.MeanR2S1EditFieldLabel.Text = 'Mean R2*(S^-1)';

            % Create MeanR2S1EditField
            app.MeanR2S1EditField = uieditfield(app.R2MapandFittingPanel, 'numeric');
            app.MeanR2S1EditField.Position = [668 13 99 22];

            % Create TEKnobLabel
            app.TEKnobLabel = uilabel(app.BOLD_Proceesor);
            app.TEKnobLabel.HorizontalAlignment = 'center';
            app.TEKnobLabel.Position = [180 88 25 22];
            app.TEKnobLabel.Text = 'TE';

            % Create TEKnob
            app.TEKnob = uiknob(app.BOLD_Proceesor, 'discrete');
            app.TEKnob.ValueChangedFcn = createCallbackFcn(app, @TEKnobValueChanged, true);
            app.TEKnob.Position = [141 124 104 104];

            % Create BoldSliceLabel
            app.BoldSliceLabel = uilabel(app.BOLD_Proceesor);
            app.BoldSliceLabel.HorizontalAlignment = 'center';
            app.BoldSliceLabel.Position = [424 93 58 22];
            app.BoldSliceLabel.Text = 'Bold Slice';

            % Create BoldSliceKnob
            app.BoldSliceKnob = uiknob(app.BOLD_Proceesor, 'discrete');
            app.BoldSliceKnob.ValueChangedFcn = createCallbackFcn(app, @BoldSliceKnobValueChanged, true);
            app.BoldSliceKnob.Position = [401 130 101 101];

            % Create GenerateR2MapButton
            app.GenerateR2MapButton = uibutton(app.BOLD_Proceesor, 'push');
            app.GenerateR2MapButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateR2MapButtonPushed, true);
            app.GenerateR2MapButton.Position = [574 170 115 23];
            app.GenerateR2MapButton.Text = 'Generate R2* Map';

            % Create BoldState
            app.BoldState = uieditfield(app.BOLD_Proceesor, 'text');
            app.BoldState.Position = [20 15 355 26];

            % Create DrawROIforFitButton
            app.DrawROIforFitButton = uibutton(app.BOLD_Proceesor, 'state');
            app.DrawROIforFitButton.ValueChangedFcn = createCallbackFcn(app, @DrawROIforFitButtonValueChanged, true);
            app.DrawROIforFitButton.Text = 'Draw ROI for Fit';
            app.DrawROIforFitButton.Position = [581 127 102 23];

            % Show the figure after all components are created
            app.DiffusonCalculatorforSIEMENSUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DiffAppV02_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.DiffusonCalculatorforSIEMENSUIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.DiffusonCalculatorforSIEMENSUIFigure)
        end
    end
end