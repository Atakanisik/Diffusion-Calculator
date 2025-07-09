classdef DiffAppV02reduced_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        DiffusonCalculatorforSIEMENSUIFigure  matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        CalculatorTab                  matlab.ui.container.Tab
        CopyDataButton                 matlab.ui.control.Button
        UseMedullaROIFromSegmenterButton  matlab.ui.control.Button
        GenerateMapsWithSegmenterROIButton  matlab.ui.control.Button
        UseCortexROIFromSegmenterButton  matlab.ui.control.Button
        LocationEditField              matlab.ui.control.NumericEditField
        LocationEditFieldLabel         matlab.ui.control.Label
        SmoothingDegreeEditField       matlab.ui.control.NumericEditField
        SmoothingDegreeEditFieldLabel  matlab.ui.control.Label
        mapRoiSwitch                   matlab.ui.control.Switch
        FilterDWIBetaButton            matlab.ui.control.Button
        StateField                     matlab.ui.control.EditField
        KernelSizeEditField            matlab.ui.control.NumericEditField
        KernelSizeEditFieldLabel       matlab.ui.control.Label
        SliceKnob                      matlab.ui.control.DiscreteKnob
        SliceKnobLabel                 matlab.ui.control.Label
        BValueKnob                     matlab.ui.control.DiscreteKnob
        BValueKnobLabel                matlab.ui.control.Label
        GenerateMapsButton             matlab.ui.control.Button
        EditField_2                    matlab.ui.control.NumericEditField
        EditField                      matlab.ui.control.NumericEditField
        AlgorithmButtonGroup           matlab.ui.container.ButtonGroup
        TriExpButton                   matlab.ui.control.RadioButton
        BayesianButton                 matlab.ui.control.RadioButton
        BiExpSegButton                 matlab.ui.control.RadioButton
        BiexpFreeButton                matlab.ui.control.RadioButton
        MonoExpButton                  matlab.ui.control.RadioButton
        LoadDICOMButton                matlab.ui.control.Button
        OptionsPanel                   matlab.ui.container.Panel
        pFslowEditField                matlab.ui.control.NumericEditField
        pFslowEditFieldLabel           matlab.ui.control.Label
        pFinterEditField               matlab.ui.control.NumericEditField
        pFinterEditFieldLabel          matlab.ui.control.Label
        pFfastEditField                matlab.ui.control.NumericEditField
        pFfastEditFieldLabel           matlab.ui.control.Label
        ImportFromLastDatabaseButton   matlab.ui.control.Button
        DslowEditField                 matlab.ui.control.NumericEditField
        DslowEditField_5Label          matlab.ui.control.Label
        DinterEditField                matlab.ui.control.NumericEditField
        DinterEditFieldLabel           matlab.ui.control.Label
        DfastEditField                 matlab.ui.control.NumericEditField
        DfastEditFieldLabel            matlab.ui.control.Label
        MultiRoiSwitch                 matlab.ui.control.Switch
        DrawMultipleROIButton          matlab.ui.control.Button
        DrawROIButton                  matlab.ui.control.StateButton
        ROIStyleDropDown               matlab.ui.control.DropDown
        ROIStyleDropDownLabel          matlab.ui.control.Label
        CalculateButton                matlab.ui.control.Button
        ADCEditField                   matlab.ui.control.NumericEditField
        ADCEditFieldLabel              matlab.ui.control.Label
        pFEditField                    matlab.ui.control.NumericEditField
        pFEditFieldLabel               matlab.ui.control.Label
        DpEditField                    matlab.ui.control.NumericEditField
        DpEditFieldLabel               matlab.ui.control.Label
        DtEditField                    matlab.ui.control.NumericEditField
        DtEditFieldLabel               matlab.ui.control.Label
        ImageViewPanel                 matlab.ui.container.Panel
        UIAxes                         matlab.ui.control.UIAxes
        FittingGraphicsPanel           matlab.ui.container.Panel
        RSquareEditField               matlab.ui.control.NumericEditField
        RSquareEditFieldLabel          matlab.ui.control.Label
        UIAxes2                        matlab.ui.control.UIAxes
        DicomTab                       matlab.ui.container.Tab
        ImageNumberSpinner             matlab.ui.control.Spinner
        ImageNumberSpinnerLabel        matlab.ui.control.Label
        DicomState                     matlab.ui.control.EditField
        PatientsPanel                  matlab.ui.container.Panel
        PatientsListBox                matlab.ui.control.ListBox
        PatientsListBoxLabel           matlab.ui.control.Label
        SequencesPanel                 matlab.ui.container.Panel
        LoadDICOMDatabaseButton        matlab.ui.control.Button
        SaveDICOMDatabaseButton        matlab.ui.control.Button
        TimeEditField                  matlab.ui.control.EditField
        TimeEditFieldLabel             matlab.ui.control.Label
        LoadtoBOLDPcsButton            matlab.ui.control.Button
        FramesEditField                matlab.ui.control.NumericEditField
        FramesEditFieldLabel           matlab.ui.control.Label
        LoadtoViewerButton             matlab.ui.control.Button
        LoadtoCalculatorButton         matlab.ui.control.Button
        SequencesListBox               matlab.ui.control.ListBox
        SequencesListBoxLabel          matlab.ui.control.Label
        PreviewAxes                    matlab.ui.control.UIAxes
        ViewerTab                      matlab.ui.container.Tab
        ImageViewPanel_2               matlab.ui.container.Panel
        BritghtnessSlider              matlab.ui.control.Slider
        BritghtnessSliderLabel         matlab.ui.control.Label
        ContrastSlider                 matlab.ui.control.RangeSlider
        ContrastSliderLabel            matlab.ui.control.Label
        ColormapDropDown               matlab.ui.control.DropDown
        ColormapDropDownLabel          matlab.ui.control.Label
        InstanceNumberEditField        matlab.ui.control.NumericEditField
        InstanceNumberEditFieldLabel   matlab.ui.control.Label
        ViewSlider                     matlab.ui.control.Slider
        ViewSliderLabel                matlab.ui.control.Label
        UIAxes8                        matlab.ui.control.UIAxes
        MapsTab                        matlab.ui.container.Tab
        DWIMAPSPanel                   matlab.ui.container.Panel
        SaveMapsButton                 matlab.ui.control.Button
        ColorMapsDropDown              matlab.ui.control.DropDown
        ColorMapsDropDownLabel         matlab.ui.control.Label
        UIAxes3_10                     matlab.ui.control.UIAxes
        UIAxes3_9                      matlab.ui.control.UIAxes
        UIAxes3_8                      matlab.ui.control.UIAxes
        UIAxes3_7                      matlab.ui.control.UIAxes
        UIAxes3_6                      matlab.ui.control.UIAxes
        UIAxes3_5                      matlab.ui.control.UIAxes
        UIAxes3_4                      matlab.ui.control.UIAxes
        UIAxes3_3                      matlab.ui.control.UIAxes
        UIAxes3_2                      matlab.ui.control.UIAxes
        UIAxes3                        matlab.ui.control.UIAxes
        BOLD_Proceesor                 matlab.ui.container.Tab
        DrawROIforFitButton            matlab.ui.control.StateButton
        BoldState                      matlab.ui.control.EditField
        GenerateR2MapButton            matlab.ui.control.Button
        BoldSliceKnob                  matlab.ui.control.DiscreteKnob
        BoldSliceLabel                 matlab.ui.control.Label
        TEKnob                         matlab.ui.control.DiscreteKnob
        TEKnobLabel                    matlab.ui.control.Label
        R2MapandFittingPanel           matlab.ui.container.Panel
        MeanR2msec1EditField           matlab.ui.control.NumericEditField
        MeanR2msec1EditFieldLabel      matlab.ui.control.Label
        UIAxes12                       matlab.ui.control.UIAxes
        UIAxes11                       matlab.ui.control.UIAxes
        R2ImagesPanel                  matlab.ui.container.Panel
        UIAxes10                       matlab.ui.control.UIAxes
        RegisterTab                    matlab.ui.container.Tab
        Panel                          matlab.ui.container.Panel
        SegmentGuidedRegisterButton    matlab.ui.control.Button
        MovingSegmentButton            matlab.ui.control.Button
        FixedSegmentButton             matlab.ui.control.Button
        DrawROIButton_3                matlab.ui.control.Button
        DrawROIButton_2                matlab.ui.control.Button
        RegisterStatusEditField        matlab.ui.control.EditField
        StatusEditFieldLabel           matlab.ui.control.Label
        FilterButton                   matlab.ui.control.Button
        CoveredDepthEditField_2        matlab.ui.control.NumericEditField
        CoveredDepthEditField_2Label   matlab.ui.control.Label
        CoveredDepthEditField          matlab.ui.control.NumericEditField
        CoveredDepthEditFieldLabel     matlab.ui.control.Label
        LocationEditField_3            matlab.ui.control.NumericEditField
        LocationEditField_3Label       matlab.ui.control.Label
        LocationEditField_2            matlab.ui.control.NumericEditField
        LocationEditField_2Label       matlab.ui.control.Label
        ShowPairButton                 matlab.ui.control.Button
        LoadButton_2                   matlab.ui.control.Button
        LoadButton                     matlab.ui.control.Button
        SendRegisteredImagetoSegmenterButton  matlab.ui.control.Button
        RegisterButton                 matlab.ui.control.Button
        MovingSliceSlider              matlab.ui.control.Slider
        MovingSliceSliderLabel         matlab.ui.control.Label
        FixedSliceSlider               matlab.ui.control.Slider
        FixedSliceSliderLabel          matlab.ui.control.Label
        MovingListBox                  matlab.ui.control.ListBox
        MovingListBoxLabel             matlab.ui.control.Label
        MovingImageButton              matlab.ui.control.Button
        FixedImageButton               matlab.ui.control.Button
        FixedListBox                   matlab.ui.control.ListBox
        FixedListBoxLabel              matlab.ui.control.Label
        UIAxes13_3                     matlab.ui.control.UIAxes
        UIAxes13_2                     matlab.ui.control.UIAxes
        UIAxes13                       matlab.ui.control.UIAxes
        SegmenterTab                   matlab.ui.container.Tab
        Segmenter                      matlab.ui.container.Panel
        CalculateSimilarityCoefficientsButton  matlab.ui.control.Button
        DICEMedullaEditField           matlab.ui.control.NumericEditField
        DICEMedullaEditFieldLabel      matlab.ui.control.Label
        JaccardMedullaEditField        matlab.ui.control.NumericEditField
        JaccardMedullaEditFieldLabel   matlab.ui.control.Label
        DICECortexEditField            matlab.ui.control.NumericEditField
        DICECortexEditFieldLabel       matlab.ui.control.Label
        JaccardCortexEditField         matlab.ui.control.NumericEditField
        JaccardCortexEditFieldLabel    matlab.ui.control.Label
        DrawCortexGroundTruthButton    matlab.ui.control.Button
        SegmentButton                  matlab.ui.control.Button
        UIAxes14_2                     matlab.ui.control.UIAxes
        UIAxes14                       matlab.ui.control.UIAxes
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
        fo_triexp;
        triexptype;
        mymodel;
        DpSegmap;
        DtSegmap;
        pFSegmap;
        kernel;
        DpUpperLimitfree = 0.05;
        DtUpperLimitfree = 0.005;
        DpUpperLimitseg = 0.05;
        DtUpperLimitseg = 0.005;
        DpUpperLimitBayes = 0.05;
        DtUpperLimitBayes = 0.005;
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
        seriestime;
        smoothDegree = 0.5;
        triexp_result;
        triexp_f1;
        triexp_f2;
        triexp_f3;
        triexp_d1;
        triexp_d2;
        triexp_d3;
        
        
        
        
        contrastLowLimit=0;
        contrastUpLimit=2.5;
        contrastValue;
        brightnessValue=0;
        viewImage;
        viewAdjustedImage;
        
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
        OEFvol;
        dPbayesMap;
        dTbayesMap;
        pFbayesMap;
        DpMCMCmap;
        DtMCMCmap;
        PfMCMCmap;
        movingValue;
        fixedValue;
        fixedVol;
        moveVol;
        fixedImage;
        moveImage;
        fixedM;
        fixedN;
        fixedZ;
        moveM;
        moveN;
        moveZ;
        fixedSliderValue;
        movingSliderValue;
        regFix;
        regMove;
        registerStruct;
        moveLoc;
        fixLoc;
        calcLoc;

      
        fixThick;
        moveThick;

        SegmenterImage;
        segmentRegion;
        segmentForegroundPoints;
        segmentBackgroundPoints;
        segmentMasks;
        medsam;
        sam;
        samembedds;
        medembedds;
        imageSegmented; % Description
        segmentMask;

        Prevol;
        preM;
        preN;
        preZ;
        segmentFixedRegion;
        segmentMovingRegion;
        regSegFix;
        regSegMove;
        SegFixEmbd;
        SegMoveEmbd;
        segmentFixedMask;
        segmentMoveMask;
        segmentedFixed;
        segmentedMove;

        Cortex;
        Medulla;
        CortexMask;
        MedullaMask;

        CortexGtRoi;
        CortexGt;
        MedullaGt;
        drawSegmented;



        JaccardCortex;
        DiceCortex;
        JaccardMedulla;
        DiceMedulla;
        
    end
    
    methods (Access = public)
        
        function [posterior_mean, posterior_std, chain,Rsquare,fitSignal] = bayesian(app,b_values, signal_data)

    % Define the IVIM model

    model_fun = @(params, b) params(1) * exp(-b * params(2)) + (1-params(1)) * exp(-b * params(3));


    % Define the prior distributions for the parameters

    prior_f = @(f) unifpdf(f, 0, 1); % uniform prior for perfusion fraction (f)

    prior_Dstar = @(Dstar) unifpdf(Dstar, 0.01, 0.5); % lognormal prior for pseudo-diffusion coefficient (D*)

    prior_D = @(D) unifpdf(D, 0.001, 0.5); % lognormal prior for true diffusion coefficient (D)


    % Define the likelihood function

    likelihood_fun = @(params, b, signal) prod(normpdf(signal, model_fun(params, b), 0.01));


    % Define the Bayesian model

    bayes_model = @(params) prior_f(params(1)) * prior_Dstar(params(2)) * prior_D(params(3)) * likelihood_fun(params, b_values, signal_data);


    % Initialize the MCMC chain

    n_samples = 1000;

    params_init = [0.5 0.01 0.001]; % initial values for f, D*, and D

    chain = zeros(n_samples, 3);


    % Run the MCMC chain

    for i = 1:n_samples

        params_prop = params_init + randn(1, 3) .* [0.1 0.01 0.001]; % propose new parameters

        params_prop(1) = max(0, min(1, params_prop(1))); % ensure f is between 0 and 1

        params_prop(2) = max(0, params_prop(2)); % ensure D* is positive

        params_prop(3) = max(0, params_prop(3)); % ensure D is positive

        posterior_prop = bayes_model(params_prop);

        posterior_init = bayes_model(params_init);

        alpha = min(1, posterior_prop / posterior_init);

        if rand < alpha

            params_init = params_prop;

        end

        chain(i, :) = params_init;

    end


    % Compute the posterior means and standard deviations

    posterior_mean = mean(chain, 1);

    posterior_std = std(chain, 0, 1);


    % Display the results

    fprintf('Posterior mean (f, D*, D): [%f, %f, %f]\n', posterior_mean);

    fprintf('Posterior standard deviation (f, D*, D): [%f, %f, %f]\n', posterior_std);


    % Plot the fitted IVIM model and the data

    fitSignal = model_fun(posterior_mean, b_values);
    residual = signal_data - fitSignal;
    S_res_ss = sum(residual.^2);

    S_tot_ss = sum((signal_data - mean(signal_data)).^2);

    Rsquare = 1 - (S_res_ss / S_tot_ss);

    % plot(b_values, signal_data, 'o', b_values, fitSignal, '-','Parent',axes);
    % figure;
    % 
    % subplot(1, 3, 1);
    % 
    % histogram(chain(:, 1), 20);
    % 
    % xlabel('Perfusion Fraction (f)');
    % 
    % ylabel('Frequency');
    % 
    % title('Posterior Distribution of f');
    % 
    % 
    % subplot(1, 3, 2);
    % 
    % histogram(chain(:, 2), 20);
    % 
    % xlabel('Pseudo-Diffusion Coefficient (D*)');
    % 
    % ylabel('Frequency');
    % 
    % title('Posterior Distribution of D*');
    % 
    % 
    % subplot(1, 3, 3);
    % 
    % histogram(chain(:, 3), 20);
    % 
    % xlabel('True Diffusion Coefficient (D)');
    % 
    % ylabel('Frequency');
    % 
    % title('Posterior Distribution of D');
    % 
        end
        
        
        
        function [MOVINGREG] = registerImages(app,MOVING,FIXED)
            fixedRefObj = imref2d(size(FIXED));
movingRefObj = imref2d(size(MOVING));

% Intensity-based registration
[optimizer, metric] = imregconfig('multimodal');
metric.NumberOfSpatialSamples = 500;
metric.NumberOfHistogramBins = 50;
metric.UseAllPixels = true;
optimizer.GrowthFactor = 1.050000;
optimizer.Epsilon = 1.50000e-06;
optimizer.InitialRadius = 5.25000e-03;
optimizer.MaximumIterations = 100;

% Align centers
fixedCenterXWorld = mean(fixedRefObj.XWorldLimits);
fixedCenterYWorld = mean(fixedRefObj.YWorldLimits);
movingCenterXWorld = mean(movingRefObj.XWorldLimits);
movingCenterYWorld = mean(movingRefObj.YWorldLimits);
translationX = fixedCenterXWorld - movingCenterXWorld;
translationY = fixedCenterYWorld - movingCenterYWorld;

% Coarse alignment
initTform = affinetform2d();
initTform.A(1:2,3) = [translationX ; translationY];

% Apply Gaussian blur
fixedInit = imgaussfilt(FIXED,1.000000);
movingInit = imgaussfilt(MOVING,1.000000);

% Normalize images
movingInit = mat2gray(movingInit);
fixedInit = mat2gray(fixedInit);

% Apply transformation
tform = imregtform(movingInit,movingRefObj,fixedInit,fixedRefObj,'similarity',optimizer,metric,'PyramidLevels',3,'InitialTransformation',initTform);
MOVINGREG.Transformation = tform;
MOVINGREG.RegisteredImage = imwarp(MOVING, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);

% Nonrigid registration
[MOVINGREG.DisplacementField,MOVINGREG.RegisteredImage] = imregdemons(MOVINGREG.RegisteredImage,FIXED,100,'AccumulatedFieldSmoothing',1.0,'PyramidLevels',3);

% Store spatial referencing object
MOVINGREG.SpatialRefObj = fixedRefObj;
            
        end
        function [cortex,cortexBW,medulla,medullaBW] = segmentCortexMedulla(app,RegisteredImage)
            RegisteredImage = imadjust(RegisteredImage);

            % Threshold image with adaptive threshold
            BW = imbinarize(im2gray(RegisteredImage), 'adaptive', 'Sensitivity', 0.500000, 'ForegroundPolarity', 'bright');
            cortexBW = BW;
            cortex = RegisteredImage;
            cortex(~cortexBW) = 0;
            medulla = imsubtract(RegisteredImage,cortex);
            medullaBW = imbinarize(medulla);
            medullaBW = bwareaopen(medullaBW,20);
            medullaBW = imfill(medullaBW,"holes");
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: LoadDICOMButton
        function LoadDICOMButtonPushed(app, event)
            
            folder = uigetdir;
            try
            
            app.collection = dicomCollection(folder,IncludeSubfolders=true);
            path = pwd;
            foldername = 'database';
            if ~exist(fullfile(path,foldername),'dir')
                mkdir(path,foldername)
            else
                savepath = [path '\database'];
                savefile = app.collection
                save(fullfile(savepath,'Data.mat'),'savefile');
            end
            app.PatientsListBox.Items = unique(app.collection.PatientName);
            % app.PatientsListBox.ItemsData = 1:length(app.PatientsListBox.Items);
            app.patientvalue = app.PatientsListBox.Value;
            app.SequencesListBox.Items = app.collection.SeriesDescription(find(app.collection.PatientName == app.patientvalue));
            app.SequencesListBox.ItemsData = 1:length(app.SequencesListBox.Items);
            app.FixedListBox.Items = app.SequencesListBox.Items;
            app.FixedListBox.ItemsData = app.SequencesListBox.ItemsData;
            app.MovingListBox.Items = app.SequencesListBox.Items;
            app.MovingListBox.ItemsData = app.SequencesListBox.ItemsData;
            app.sequencevalue = app.SequencesListBox.Items(app.SequencesListBox.Value);
            app.ColormapDropDown.Items = app.cmaps;
            app.ColorMapsDropDown.Items = app.cmaps;
            app.SmoothingDegreeEditField.Value = app.smoothDegree;
            app.StateField.Value = "DICOM Collection Installed. See Dicom Tab.";
            app.StateField.BackgroundColor = [1.0,0.48,0.51];
            savepath = [path '\database'];
            savefile = app.collection
            save(fullfile(savepath,'Data.mat'),'savefile');
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
             app.triexptype = fittype( ...
    'S0 * (f1 * exp(-b * D1) + f2 * exp(-b * D2) + (1 - f1 - f2) * exp(-b * D3))', ...
    'independent', 'b', ...
    'coefficients', {'f1', 'f2', 'D1', 'D2', 'D3', 'S0'});
            app.fo_triexp = fitoptions('Method','NonLinearLeastSquares', ...
    'Algorithm','Levenberg-Marquardt', ...
    'StartPoint',[0.1, 0.1, 0.01, 0.003, 0.0005, 1], ...
    'Lower',[0, 0, 0.005, 0.0005, 0.0001, 0.5], ...
    'Upper',[0.8, 0.8, 0.1, 0.02, 0.005, 2]);
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
             [app.triexp_result,gtriexp] = fit(app.bvals_ivim.',app.s_s0_ivim.',app.triexptype,app.fo_triexp);
            app.triexp_f1 = app.triexp_result.f1;
            app.triexp_f2 = app.triexp_result.f2;
            app.triexp_f3 = 1 - app.triexp_f1 - app.triexp_f2;
            app.triexp_d1 = app.triexp_result.D1;
            app.triexp_d2 = app.triexp_result.D2;
            app.triexp_d3 = app.triexp_result.D3;
            if (app.MonoExpButton.Value == 1)
                app.rsqr = gadc.rsquare;
                app.RSquareEditField.Value = app.rsqr;
                app.ADCEditField.Value = app.ADC;
                app.DpEditField.Value = 0;
                app.DtEditField.Value = 0;
                app.pFEditField.Value = 0;
                app.DfastEditField.Value = 0;
                app.DinterEditField.Value = 0;
                app.DslowEditField.Value = 0;
                app.pFfastEditField.Value = 0;
                app.pFinterEditField.Value = 0;
                app.pFslowEditField.Value = 0;
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
                app.DfastEditField.Value = 0;
                app.DinterEditField.Value = 0;
                app.DslowEditField.Value = 0;
                app.pFfastEditField.Value = 0;
                app.pFinterEditField.Value = 0;
                app.pFslowEditField.Value = 0;
                app.DtUpperLimitfree = app.DtEditField.Value*2;
                app.DpUpperLimitfree = app.DpEditField.Value*2;
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
                app.DfastEditField.Value = 0;
                app.DinterEditField.Value = 0;
                app.DslowEditField.Value = 0;
                app.pFfastEditField.Value = 0;
                app.pFinterEditField.Value = 0;
                app.pFslowEditField.Value = 0;
                app.DtUpperLimitseg = app.DtEditField.Value*2;
                app.DpUpperLimitseg = app.DpEditField.Value*2;
                plot(app.UIAxes2,app.segFitAll,app.bvals_ivim,app.s_s0_ivim);
                ylim(app.UIAxes2,[0 1]);
                title(app.UIAxes2,'Fitting Graphic');
                xlabel(app.UIAxes2,'B-Value');
                ylabel(app.UIAxes2,'S / S0');
            end
            if (app.BayesianButton.Value == 1)
                [p,l,c,rsqrb,fits] = bayesian(app,app.bvals_ivim,app.s_s0_ivim);
                app.rsqr = rsqrb;
                app.RSquareEditField.Value = app.rsqr;
                app.ADCEditField.Value = 0;
                app.DpEditField.Value = p(2);
                app.DtEditField.Value = p(3);
                app.pFEditField.Value = p(1);
                app.DfastEditField.Value = 0;
                app.DinterEditField.Value = 0;
                app.DslowEditField.Value = 0;
                app.pFfastEditField.Value = 0;
                app.pFinterEditField.Value = 0;
                app.pFslowEditField.Value = 0;
                app.DtUpperLimitBayes = app.DtEditField.Value*2;
                app.DpUpperLimitBayes = app.DpEditField.Value*2;
                plot(app.UIAxes2,app.bvals_ivim,app.s_s0_ivim,'o',app.bvals_ivim,fits,'-');
                ylim(app.UIAxes2,[0 1]);
                title(app.UIAxes2,'Fitting Graphic');
                xlabel(app.UIAxes2,'B-Value');
                ylabel(app.UIAxes2,'S / S0');
            end
            if (app.TriExpButton.Value == 1)
                app.rsqr = gtriexp.rsquare;
                app.RSquareEditField.Value = app.rsqr;
                app.ADCEditField.Value = 0;
                app.DpEditField.Value = 0;
                app.DtEditField.Value = 0;
                app.pFEditField.Value = 0;
                app.DfastEditField.Value = app.triexp_d1;
                app.DinterEditField.Value = app.triexp_d2;
                app.DslowEditField.Value = app.triexp_d3;
                app.pFfastEditField.Value = app.triexp_f1;
                app.pFinterEditField.Value = app.triexp_f2;
                app.pFslowEditField.Value = app.triexp_f3;
                plot(app.UIAxes2,app.triexp_result,app.bvals_ivim,app.s_s0_ivim);
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
            app.LocationEditField.Value = app.calcLoc(app.sliceNum);
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
                        
                        pixelADCfit = fit(pixelBvals.',pixelSigMatrix.',app.adctype,app.foadc);
                        app.ADCmap(xc,yc) = pixelADCfit.b;
                        catch
                            app.ADCmap(xc,yc) = 0;
                        end
                        try
                        pixelFreefit = fit(pixelBvals.',pixelSigMatrix',app.fotype,app.fo);
                        pixelParamFree = [pixelFreefit.D pixelFreefit.D_star pixelFreefit.f];

                        app.Dpmap(xc,yc) = pixelParamFree(2);
                        app.Dtmap(xc,yc) = pixelParamFree(1);
                        app.pFmap(xc,yc) = pixelParamFree(3);
                        catch
                        app.Dpmap(xc,yc) = 0;
                        app.Dtmap(xc,yc) = 0;
                        app.pFmap(xc,yc) = 0;
                        end
                        try
                        pixelSigMatrixHigh = pixelSigMatrix(find(app.b_map_ivim > 200));
                        pixelsegFitHigh = fit(pixelBvalHigh.',pixelSigMatrixHigh.','exp1',app.segfo);
                        DT = pixelsegFitHigh.b*-1;
                        SI = pixelsegFitHigh.a * pixelSigMatrix(1);
                        PF = (pixelSigMatrix(1)-SI)/pixelSigMatrix(1);
                        pixelmodel = @(DP,x) (PF*(exp(-x*DP)) + (1-PF)*exp(-x*DT));
                        pixelSegAllFit = fit(pixelBvals.',pixelSigMatrix.',pixelmodel,app.fo2);
                        DP = pixelSegAllFit.DP;
                        app.DpSegmap(xc,yc) =DP;
                        app.DtSegmap(xc,yc) = DT;
                        app.pFSegmap(xc,yc) = PF;
                        catch
                        
                        app.DpSegmap(xc,yc) =0;
                        app.DtSegmap(xc,yc) = 0;
                        app.pFSegmap(xc,yc) = 0;
                        end
                        try
                        [pmap,~,~,~,~] = bayesian(app,pixelBvals,pixelSigMatrix);
                        
    
                       
                        
                       
                        app.dPbayesMap(xc,yc) = pmap(2);
                        app.dTbayesMap(xc,yc) = pmap(3);
                        app.pFbayesMap(xc,yc) = pmap(1);
                       
                        catch
                        
                        
                        
                        app.dPbayesMap(xc,yc) = 0;
                        app.dTbayesMap(xc,yc) = 0;
                        app.pFbayesMap(xc,yc) = 0;
                      
                        end
                    else
                        app.ADCmap(xc,yc) = 0;
                        app.Dpmap(xc,yc) = 0;
                        app.Dtmap(xc,yc) = 0;
                        app.pFmap(xc,yc) = 0;
                        app.DpSegmap(xc,yc) =0;
                        app.DtSegmap(xc,yc) = 0;
                        app.pFSegmap(xc,yc) = 0;
                        app.dPbayesMap(xc,yc) = 0;
                        app.dTbayesMap(xc,yc) = 0;
                        app.pFbayesMap(xc,yc) = 0;
                       
                    end
                    
                end
            end
                    
                    
                    

                     
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

                     imagesc(app.UIAxes3_8,app.dPbayesMap,[0,app.DpUpperLimitBayes]);
                     colormap(app.UIAxes3_8,"jet");
                     colorbar(app.UIAxes3_8);
                     xlim(app.UIAxes3_8,[0 ,app.N]);
                     ylim(app.UIAxes3_8,[0, app.M]);

                     imagesc(app.UIAxes3_9,app.dTbayesMap,[0,app.DtUpperLimitBayes]);
                     colormap(app.UIAxes3_9,"jet");
                     colorbar(app.UIAxes3_9);
                     xlim(app.UIAxes3_9,[0 ,app.N]);
                     ylim(app.UIAxes3_9,[0, app.M]);

                     imagesc(app.UIAxes3_10,app.pFbayesMap,[0,1]);
                     colormap(app.UIAxes3_10,"jet");
                     colorbar(app.UIAxes3_10);
                     xlim(app.UIAxes3_10,[0 ,app.N]);
                     ylim(app.UIAxes3_10,[0, app.M]);

                    


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
            app.calcLoc=[];
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                info = dicominfo(f);
                seq = info.SequenceName;
                app.calcLoc=horzcat(app.calcLoc,info.SliceLocation);
                b_val = horzcat(b_val,str2double(extract(seq, digitsPattern)));
                svol = dicomread(f);
                app.vol = cat(3,app.vol,svol);

            
            end
            app.calcLoc = unique(app.calcLoc);
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
            app.seriestime = app.patienttable.SeriesDateTime;
            app.TimeEditField.Value = string(app.seriestime(I));
            app.FramesEditField.Value = app.frames(I);
            fileNames = app.patienttable.Filenames(I);
            app.ImageNumberSpinner.Value = 1;
            rowPre = app.patienttable.Rows(I);
            colPre = app.patienttable.Columns(I);
            app.Prevol = zeros(rowPre , colPre);
            fileNames = fileNames{1,1};
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                
                
                
                svol = dicomread(f);
                app.Prevol = cat(3,app.Prevol,svol);

            
            end
            [app.preM,app.preN,app.preZ] = size(app.Prevol);
            app.Prevol = app.Prevol(:,:,2:end);
            imagesc(app.PreviewAxes,app.Prevol(:,:,1));
            colormap(app.PreviewAxes,'gray');
            xlim(app.PreviewAxes,[0,app.preN]);
            ylim(app.PreviewAxes,[0,app.preM]);
            app.ImageNumberSpinner.Limits=[1,app.preZ];

        end

        % Value changed function: PatientsListBox
        function PatientsListBoxValueChanged(app, event)
            app.patientvalue = app.PatientsListBox.Value;
            app.SequencesListBox.Items = app.collection.SeriesDescription(find(app.collection.PatientName == app.patientvalue));
            app.FixedListBox.Items = app.SequencesListBox.Items;
            app.MovingListBox.Items = app.SequencesListBox.Items;
            
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
            colormap(app.UIAxes3_8,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_9,app.ColorMapsDropDown.Value);
            colormap(app.UIAxes3_10,app.ColorMapsDropDown.Value);
            
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
            app.MeanR2msec1EditField.Value = roifit.b*1000;
           
            plot(app.UIAxes12,roifit,app.TEvals,boldFitsig);
        end

        % Button pushed function: SaveMapsButton
        function SaveMapsButtonPushed(app, event)
               adc = app.ADCmap;
               dTfree = app.Dtmap;
               dPfree = app.Dpmap;
               pFfree = app.pFmap;
               dtSegmented = app.DtSegmap;
               dpSegmented = app.DpSegmap;
               pfSegmented = app.pFSegmap;
               dtBayesian = app.dTbayesMap;
               dpBayesian = app.dPbayesMap;
               pfBayesian = app.pFbayesMap;
               uisave({'adc','dTfree','dPfree','pFfree','dtSegmented','dpSegmented','pfSegmented','dtBayesian','dpBayesian','pfBayesian'},'Maps');
        end

        % Button pushed function: ImportFromLastDatabaseButton
        function ImportFromLastDatabaseButtonPushed(app, event)
            loadpath = pwd;
            loadpath = [loadpath,'\database','\Data.mat'];
            app.collection = load(loadpath);
            app.collection = app.collection.savefile;
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

        end

        % Button pushed function: SaveDICOMDatabaseButton
        function SaveDICOMDatabaseButtonPushed(app, event)
            data = app.collection;
            uisave({'data'},'DicomDataBase');
        end

        % Button pushed function: LoadDICOMDatabaseButton
        function LoadDICOMDatabaseButtonPushed(app, event)
            [fileSelected ,pathSelected] = uigetfile;
            app.collection = load(fullfile(pathSelected,fileSelected));
            app.collection = app.collection.data;
            app.PatientsListBox.Items = unique(app.collection.PatientName);
            % app.PatientsListBox.ItemsData = 1:length(app.PatientsListBox.Items);
            app.patientvalue = app.PatientsListBox.Value;
            app.SequencesListBox.Items = app.collection.SeriesDescription(find(app.collection.PatientName == app.patientvalue));
            app.SequencesListBox.ItemsData = 1:length(app.SequencesListBox.Items);
            app.FixedListBox.Items = app.SequencesListBox.Items;
            app.FixedListBox.ItemsData = app.SequencesListBox.ItemsData;
            app.MovingListBox.Items = app.SequencesListBox.Items;
            app.MovingListBox.ItemsData = app.SequencesListBox.ItemsData;
            app.sequencevalue = app.SequencesListBox.Items(app.SequencesListBox.Value);
            app.ColormapDropDown.Items = app.cmaps;
            app.ColorMapsDropDown.Items = app.cmaps;
            app.SmoothingDegreeEditField.Value = app.smoothDegree;
            app.StateField.Value = "DICOM Collection Installed. See Dicom Tab.";
            app.StateField.BackgroundColor = [1.0,0.48,0.51];
        end

        % Value changed function: FixedListBox
        function FixedListBoxValueChanged(app, event)
            app.fixedValue = app.FixedListBox.Items(app.FixedListBox.Value);
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            Ifixed = find(app.patienttable.SeriesDescription == app.fixedValue);
            Ifixed = Ifixed(find(Ifixed == app.FixedListBox.Value));
            
        end

        % Value changed function: MovingListBox
        function MovingListBoxValueChanged(app, event)
            app.movingValue = app.MovingListBox.Items(app.MovingListBox.Value);
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            Imoving = find(app.patienttable.SeriesDescription == app.movingValue);
            Imoving = Imoving(find(Imoving == app.MovingListBox.Value));
            
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            app.RegisterStatusEditField.Value = "Fixed Image Lodaing";
            pause(0.001);
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            Ifixed = find(app.patienttable.SeriesDescription == app.fixedValue);
            Ifixed = Ifixed(find(Ifixed == app.FixedListBox.Value));
            
            fileNames = app.patienttable.Filenames(Ifixed);
            fileNames = fileNames{1,1};
             % app.rows = app.patienttable.Rows(I);
            % app.cols = app.patienttable.Columns(I);
            app.fixedVol = [];
            app.fixLoc = [];
            
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                svol = dicomread(f);
                info = dicominfo(f);
                app.fixThick = info.SliceThickness;
                app.fixLoc = horzcat(app.fixLoc,info.SliceLocation);
                app.fixedVol = cat(3,app.fixedVol,svol);
            end
            [app.fixedM, app.fixedN, app.fixedZ] = size(app.fixedVol);
            
            app.fixedVol = app.fixedVol(:,:,1:end);
            
            imagesc(app.UIAxes13,imadjust(app.fixedVol(:,:,1)));
            colormap(app.UIAxes13,'gray');
            
            xlim(app.UIAxes13,[0 ,app.fixedN]);
            ylim(app.UIAxes13,[0, app.fixedM]);
            if app.fixedZ > 1
            app.FixedSliceSlider.Limits = [1,app.fixedZ];
            app.FixedSliceSlider.MajorTicks = 1:round(app.fixedZ/10):app.fixedZ;
            app.FixedSliceSlider.MajorTickLabels = string(1:round(app.fixedZ/10):app.fixedZ);
            app.FixedSliceSlider.MinorTicks = 1:1:app.fixedZ;
            else
                app.FixedSliceSlider.Limits = [0,1];
                app.FixedSliceSlider.MajorTicks = 0:1;
                app.FixedSliceSlider.MajorTickLabels = string(0:1);
                app.FixedSliceSlider.MinorTicks = [];
            end
            app.FixedSliceSlider.Value = 1;
            app.RegisterStatusEditField.Value = "Fixed Image Loaded";
        end

        % Value changing function: FixedSliceSlider
        function FixedSliceSliderValueChanging(app, event)
            app.fixedSliderValue = round(event.Value);
            app.LocationEditField_2.Value = app.fixLoc(app.fixedSliderValue);
            app.CoveredDepthEditField.Value = app.fixLoc(app.fixedSliderValue) + app.fixThick;
            app.fixedImage = imadjust(app.fixedVol(:,:,app.fixedSliderValue));
            imagesc(app.UIAxes13,app.fixedImage);
            colormap(app.UIAxes13,'gray');
            xlim(app.UIAxes13,[0 ,app.fixedN]);
            ylim(app.UIAxes13,[0, app.fixedM]);
            
        end

        % Button pushed function: LoadButton_2
        function LoadButton_2Pushed(app, event)
           app.RegisterStatusEditField.Value = "Moving Image Loading";
            pause(0.001);
            app.patienttable = app.collection(find(app.collection.PatientName == app.patientvalue),1:end);
            Imoving = find(app.patienttable.SeriesDescription == app.movingValue);
            Imoving = Imoving(find(Imoving == app.MovingListBox.Value));
            
            fileNames = app.patienttable.Filenames(Imoving);
            fileNames = fileNames{1,1};
             % app.rows = app.patienttable.Rows(I);
            % app.cols = app.patienttable.Columns(I);
            app.moveVol = [];
            app.moveLoc = [];
            for i = 1:1:length(fileNames)
                f = fileNames(i);
                svol = dicomread(f);
                info = dicominfo(f);
                app.moveThick=info.SliceThickness;
                app.moveLoc = horzcat(app.moveLoc,info.SliceLocation);
                app.moveVol = cat(3,app.moveVol,svol);
            end
            [app.moveM, app.moveN, app.moveZ] = size(app.moveVol);
            
            app.moveVol = app.moveVol(:,:,1:end);
            
            imagesc(app.UIAxes13_2,imadjust(app.moveVol(:,:,1)));
            colormap(app.UIAxes13_2,'gray');
            
            xlim(app.UIAxes13_2,[0 ,app.moveN]);
            ylim(app.UIAxes13_2,[0, app.moveM]);
            if app.moveZ > 1
            app.MovingSliceSlider.Limits = [1,app.moveZ];
            app.MovingSliceSlider.MajorTicks = 1:round(app.moveZ/10):app.moveZ;
            app.MovingSliceSlider.MajorTickLabels = string(1:round(app.moveZ/10):app.moveZ);
            app.MovingSliceSlider.MinorTicks = 1:1:app.moveZ;
            else
                app.MovingSliceSlider.Limits = [0,1];
                app.MovingSliceSlider.MajorTicks = 0:1;
                app.MovingSliceSlider.MajorTickLabels = string(0:1);
                app.MovingSliceSlider.MinorTicks = [];
            end
            app.MovingSliceSlider.Value = 1;
            app.RegisterStatusEditField.Value = "Moving Image Loaded";
        end

        % Value changing function: MovingSliceSlider
        function MovingSliceSliderValueChanging(app, event)
            app.movingSliderValue = round(event.Value);
            app.LocationEditField_3.Value = app.moveLoc(app.movingSliderValue);
            app.CoveredDepthEditField_2.Value = app.moveLoc(app.movingSliderValue) + app.moveThick;
            app.moveImage = imadjust(app.moveVol(:,:,app.movingSliderValue));
            imagesc(app.UIAxes13_2,app.moveImage);
            colormap(app.UIAxes13_2,'gray');
            xlim(app.UIAxes13_2,[0 ,app.moveN]);
            ylim(app.UIAxes13_2,[0, app.moveM]);
            
        end

        % Button pushed function: FixedImageButton
        function FixedImageButtonPushed(app, event)
            app.regFix = imadjust(app.fixedVol(:,:,app.fixedSliderValue));
            app.FixedImageButton.BackgroundColor = [0 1 0];
        end

        % Button pushed function: MovingImageButton
        function MovingImageButtonPushed(app, event)
            app.regMove  = imadjust(app.moveVol(:,:,app.movingSliderValue));
            app.regMove = imresize(app.regMove,[app.fixedM,app.fixedN]);
            app.MovingImageButton.BackgroundColor = [0 1 0];
        end

        % Button pushed function: ShowPairButton
        function ShowPairButtonPushed(app, event)
            imshowpair(app.regSegFix,app.regSegMove,'falsecolor','Scaling','independent','Parent',app.UIAxes13_3);
            app.FixedImageButton.BackgroundColor = [0 1 1];
            app.MovingImageButton.BackgroundColor = [0 1 1];
        end

        % Button pushed function: RegisterButton
        function RegisterButtonPushed(app, event)
            app.RegisterStatusEditField.Value = "Registering Images";
            pause(0.001);
            app.registerStruct = registerImages(app,app.regMove,app.regFix);
            imshowpair(app.regFix,app.registerStruct.RegisteredImage,'falsecolor','ColorChannel',[1 2 0],'Parent',app.UIAxes13_3);
            app.RegisterStatusEditField.Value = "Images Registered";
        end

        % Button pushed function: FilterButton
        function FilterButtonPushed(app, event)
            for idx = 1:1:app.fixedZ

                patch = std2(app.fixedVol(:,:,idx));
                app.fixedVol(:,:,idx) = imnlmfilt(app.fixedVol(:,:,idx),DegreeOfSmoothing=0.5*patch);
            end
        end

        % Button pushed function: SendRegisteredImagetoSegmenterButton
        function SendRegisteredImagetoSegmenterButtonPushed(app, event)
            app.RegisterStatusEditField.Value = "Image is Sending Segmenter";
            pause(0.001);
            app.SegmenterImage = app.registerStruct.RegisteredImage;
            imagesc(app.UIAxes14,app.SegmenterImage);
            xlim(app.UIAxes14,[0 ,app.fixedN]);
            ylim(app.UIAxes14,[0, app.fixedM]);
            colormap(app.UIAxes14,'gray');
            app.RegisterStatusEditField.Value = "Creating Model and Embeddings";
            pause(0.001);

            app.RegisterStatusEditField.Value = "Image Sent";
        end

        % Button pushed function: SegmentButton
        function SegmentButtonPushed(app, event)
          
            [app.Cortex,app.CortexMask,app.Medulla,app.MedullaMask] = segmentCortexMedulla(app,app.SegmenterImage);
            app.imageSegmented = insertObjectMask(app.SegmenterImage,app.CortexMask,"MaskColor","yellow","Opacity",0.3);
            app.imageSegmented=insertObjectMask(app.imageSegmented,app.MedullaMask,"MaskColor","red","Opacity",0.3);
            
            imagesc(app.UIAxes14_2,app.imageSegmented);
            app.UseCortexROIFromSegmenterButton.Enable = "on";
            app.UseMedullaROIFromSegmenterButton.Enable = "on";
            app.DrawCortexGroundTruthButton.Enable = "on";


        end

        % Button pushed function: UseCortexROIFromSegmenterButton
        function UseCortexROIFromSegmenterButtonPushed(app, event)
            app.GenerateMapsWithSegmenterROIButton.Enable = "on";
            mask = app.CortexMask;
            app.segmentMask = app.CortexMask;
            app.sigData = [0];

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
            overlayimage = app.vol(:,:,(app.showBval)+((app.sliceNum-1)*length(app.bvals)));
            overlayimage = imadjust(overlayimage);
            overlayimage = insertObjectMask(overlayimage,mask,"MaskColor","magenta","Opacity",0.3);
            imagesc(app.UIAxes,overlayimage);
            xlim(app.UIAxes,[0 ,app.N]);
            ylim(app.UIAxes,[0, app.M]);
            title(app.UIAxes,'IMAGE');
            
            
       

        end

        % Button pushed function: GenerateMapsWithSegmenterROIButton
        function GenerateMapsWithSegmenterROIButtonPushed(app, event)
            
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
                app.mapMask = app.segmentMask;
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
                        
                        pixelADCfit = fit(pixelBvals.',pixelSigMatrix.',app.adctype,app.foadc);
                        app.ADCmap(xc,yc) = pixelADCfit.b;
                        catch
                            app.ADCmap(xc,yc) = 0;
                        end
                        try
                        pixelFreefit = fit(pixelBvals.',pixelSigMatrix',app.fotype,app.fo);
                        pixelParamFree = [pixelFreefit.D pixelFreefit.D_star pixelFreefit.f];

                        app.Dpmap(xc,yc) = pixelParamFree(2);
                        app.Dtmap(xc,yc) = pixelParamFree(1);
                        app.pFmap(xc,yc) = pixelParamFree(3);
                        catch
                        app.Dpmap(xc,yc) = 0;
                        app.Dtmap(xc,yc) = 0;
                        app.pFmap(xc,yc) = 0;
                        end
                        try
                        pixelSigMatrixHigh = pixelSigMatrix(find(app.b_map_ivim > 200));
                        pixelsegFitHigh = fit(pixelBvalHigh.',pixelSigMatrixHigh.','exp1',app.segfo);
                        DT = pixelsegFitHigh.b*-1;
                        SI = pixelsegFitHigh.a * pixelSigMatrix(1);
                        PF = (pixelSigMatrix(1)-SI)/pixelSigMatrix(1);
                        pixelmodel = @(DP,x) (PF*(exp(-x*DP)) + (1-PF)*exp(-x*DT));
                        pixelSegAllFit = fit(pixelBvals.',pixelSigMatrix.',pixelmodel,app.fo2);
                        DP = pixelSegAllFit.DP;
                        app.DpSegmap(xc,yc) =DP;
                        app.DtSegmap(xc,yc) = DT;
                        app.pFSegmap(xc,yc) = PF;
                        catch
                        
                        app.DpSegmap(xc,yc) =0;
                        app.DtSegmap(xc,yc) = 0;
                        app.pFSegmap(xc,yc) = 0;
                        end
                        try
                        [pmap,~,~,~,~] = bayesian(app,pixelBvals,pixelSigMatrix);
                        
    
                       
                        
                       
                        app.dPbayesMap(xc,yc) = pmap(2);
                        app.dTbayesMap(xc,yc) = pmap(3);
                        app.pFbayesMap(xc,yc) = pmap(1);
                       
                        catch
                        
                        
                        
                        app.dPbayesMap(xc,yc) = 0;
                        app.dTbayesMap(xc,yc) = 0;
                        app.pFbayesMap(xc,yc) = 0;
                      
                        end
                    else
                        app.ADCmap(xc,yc) = 0;
                        app.Dpmap(xc,yc) = 0;
                        app.Dtmap(xc,yc) = 0;
                        app.pFmap(xc,yc) = 0;
                        app.DpSegmap(xc,yc) =0;
                        app.DtSegmap(xc,yc) = 0;
                        app.pFSegmap(xc,yc) = 0;
                        app.dPbayesMap(xc,yc) = 0;
                        app.dTbayesMap(xc,yc) = 0;
                        app.pFbayesMap(xc,yc) = 0;
                       
                    end
                    
                end
            end
                    
                    
                    

                     
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

                     imagesc(app.UIAxes3_8,app.dPbayesMap,[0,app.DpUpperLimitBayes]);
                     colormap(app.UIAxes3_8,"jet");
                     colorbar(app.UIAxes3_8);
                     xlim(app.UIAxes3_8,[0 ,app.N]);
                     ylim(app.UIAxes3_8,[0, app.M]);

                     imagesc(app.UIAxes3_9,app.dTbayesMap,[0,app.DtUpperLimitBayes]);
                     colormap(app.UIAxes3_9,"jet");
                     colorbar(app.UIAxes3_9);
                     xlim(app.UIAxes3_9,[0 ,app.N]);
                     ylim(app.UIAxes3_9,[0, app.M]);

                     imagesc(app.UIAxes3_10,app.pFbayesMap,[0,1]);
                     colormap(app.UIAxes3_10,"jet");
                     colorbar(app.UIAxes3_10);
                     xlim(app.UIAxes3_10,[0 ,app.N]);
                     ylim(app.UIAxes3_10,[0, app.M]);

                    


                    app.StateField.Value = "Map Calculation is Done";
                    app.StateField.BackgroundColor = [0.29,0.84,0.16];
                    
                 catch
                    app.StateField.Value = "Maps Could Not be Generated, Try Changing Kernel Size";
                 end



        end

        % Value changing function: ImageNumberSpinner
        function ImageNumberSpinnerValueChanging(app, event)
            Value = event.Value;
            imagesc(app.PreviewAxes,app.Prevol(:,:,Value));
            colormap(app.PreviewAxes,'gray');
            xlim(app.PreviewAxes,[0,app.preN]);
            ylim(app.PreviewAxes,[0,app.preM]);
            
        end

        % Button pushed function: DrawROIButton_2
        function DrawROIButton_2Pushed(app, event)
             delete(findall(app.UIAxes13,'Type','images.roi'));
         
            app.segmentFixedRegion = drawrectangle(app.UIAxes13);
            isempty(app.medsam)
            
            if isempty(app.medsam) 
                app.medsam = medicalSegmentAnythingModel;
            end
            app.regSegFix = imadjust(app.fixedVol(:,:,app.fixedSliderValue));
            app.SegFixEmbd = extractEmbeddings(app.medsam,app.regSegFix);
            app.segmentFixedMask = segmentObjectsFromEmbeddings(app.medsam,app.SegFixEmbd,size(app.regSegFix),BoundingBox = app.segmentFixedRegion.Position);
            app.segmentedFixed = insertObjectMask(app.regSegFix,app.segmentFixedMask);
            imagesc(app.UIAxes13,app.segmentedFixed);
            colormap(app.UIAxes13,"gray");

        end

        % Button pushed function: DrawROIButton_3
        function DrawROIButton_3Pushed(app, event)
             delete(findall(app.UIAxes13_2,'Type','images.roi'));
            app.regSegMove  = imadjust(app.moveVol(:,:,app.movingSliderValue));
            app.regSegMove = imresize(app.regSegMove,[app.fixedM,app.fixedN]);
            imagesc(app.UIAxes13_2,app.regSegMove);
            xlim(app.UIAxes13_2,[0,app.fixedN]);
            ylim(app.UIAxes13_2,[0,app.fixedM]);
            app.segmentMovingRegion = drawrectangle(app.UIAxes13_2);
            if isempty(app.medsam)
                app.medsam = medicalSegmentAnythingModel;
            end
            
            app.SegMoveEmbd = extractEmbeddings(app.medsam,app.regSegMove);
            app.segmentMoveMask = segmentObjectsFromEmbeddings(app.medsam,app.SegMoveEmbd,size(app.regSegMove),BoundingBox = app.segmentMovingRegion.Position);
            app.segmentedMove = insertObjectMask(app.regSegMove,app.segmentMoveMask);
            imagesc(app.UIAxes13_2,app.segmentedMove);
            colormap(app.UIAxes13_2,"gray");


        end

        % Button pushed function: FixedSegmentButton
        function FixedSegmentButtonPushed(app, event)
            app.regSegFix(~app.segmentFixedMask) = 0;
            imagesc(app.UIAxes13_3,app.regSegFix);
        end

        % Button pushed function: MovingSegmentButton
        function MovingSegmentButtonPushed(app, event)
            app.regSegMove(~app.segmentMoveMask) = 0;
            imagesc(app.UIAxes13_3,app.regSegMove);
        end

        % Button pushed function: SegmentGuidedRegisterButton
        function SegmentGuidedRegisterButtonPushed(app, event)
            app.RegisterStatusEditField.Value = "Registering Images";
            pause(0.001);
            app.registerStruct = registerImages(app,app.regSegMove,app.regSegFix);
            imshowpair(app.regSegFix,app.registerStruct.RegisteredImage,'falsecolor','ColorChannel',[1 2 0],'Parent',app.UIAxes13_3);
            app.RegisterStatusEditField.Value = "Images Registered";
        end

        % Button pushed function: UseMedullaROIFromSegmenterButton
        function UseMedullaROIFromSegmenterButtonPushed(app, event)
            app.GenerateMapsWithSegmenterROIButton.Enable = "on";
            mask = app.MedullaMask;
            app.segmentMask = app.MedullaMask;
            app.sigData = [0];

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
            overlayimage = app.vol(:,:,(app.showBval)+((app.sliceNum-1)*length(app.bvals)));
            overlayimage = imadjust(overlayimage);
            overlayimage = insertObjectMask(overlayimage,mask,"MaskColor","magenta","Opacity",0.3);
            imagesc(app.UIAxes,overlayimage);
            xlim(app.UIAxes,[0 ,app.N]);
            ylim(app.UIAxes,[0, app.M]);
            title(app.UIAxes,'IMAGE');
        end

        % Button pushed function: DrawCortexGroundTruthButton
        function DrawCortexGroundTruthButtonPushed(app, event)
            app.CortexGtRoi = drawassisted(findobj(app.UIAxes14,'Type','Image'));
            app.CortexGt = createMask(app.CortexGtRoi);
            temp = app.SegmenterImage;
            tempDraw = app.SegmenterImage;
            temp(~app.CortexGt) = 0;
            app.MedullaGt = imsubtract(app.SegmenterImage,temp);
            app.MedullaGt = imbinarize(app.MedullaGt);
            app.MedullaGt = bwareaopen(app.MedullaGt,20);
            app.MedullaGt = imfill(app.MedullaGt,"holes");
            app.drawSegmented = insertObjectMask(tempDraw,app.CortexGt,"MaskColor","magenta","Opacity",0.3);
            app.drawSegmented=insertObjectMask(app.drawSegmented,app.MedullaMask,"MaskColor","green","Opacity",0.3);
            imagesc(app.UIAxes14,app.drawSegmented);
        end

        % Button pushed function: CalculateSimilarityCoefficientsButton
        function CalculateSimilarityCoefficientsButtonPushed(app, event)
            app.JaccardCortex = jaccard(app.CortexMask,app.CortexGt);
            app.DiceCortex = dice(app.CortexMask,app.CortexGt);
            app.JaccardMedulla = jaccard(app.MedullaMask,app.MedullaGt);
            app.DiceMedulla = dice(app.MedullaMask,app.MedullaGt);
            app.JaccardCortexEditField.Value = app.JaccardCortex;
            app.DICECortexEditField.Value = app.DiceCortex;
            app.JaccardMedullaEditField.Value = app.JaccardMedulla;
            app.DICEMedullaEditField.Value = app.DiceMedulla;
            figure(1);
            imshowpair(app.MedullaMask,app.MedullaGt,'falsecolor','Scaling','independent');
            figure(2);
            imshowpair(app.CortexMask,app.CortexGt,'falsecolor','Scaling','independent');
        end

        % Button pushed function: CopyDataButton
        function CopyDataButtonPushed(app, event)
            fields = {
                app.ADCEditField;
                app.DpEditField;
                app.DtEditField;
                app.pFEditField;
                app.DfastEditField;
                app.DinterEditField;
                app.DslowEditField;
                app.pFfastEditField;
                app.pFinterEditField;
                app.pFslowEditField;
                }
    values = [];

    for i = 1:numel(fields)
        fieldVal = fields{i}.Value;

        % Değer string ise sayıya çevir
        if ischar(fieldVal) || isstring(fieldVal)
            val = str2double(fieldVal);
        else
            val = fieldVal;
        end

        % Geçerli ve sıfır değilse kaydet
        if ~isnan(val) && val ~= 0
            values(end+1) = val; %#ok<AGROW>
        end
    end

    % 3. Panoya tek satır, tab ayraclı metin olarak kopyala
    if ~isempty(values)
        clipboardStr = strjoin(arrayfun(@(x) sprintf('%.5f', x), values, 'UniformOutput', false), '\t');
        clipboard('copy', clipboardStr);
    end
       

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
            app.OptionsPanel.Position = [1015 1 518 307];

            % Create DtEditFieldLabel
            app.DtEditFieldLabel = uilabel(app.OptionsPanel);
            app.DtEditFieldLabel.HorizontalAlignment = 'right';
            app.DtEditFieldLabel.Position = [149 260 25 22];
            app.DtEditFieldLabel.Text = 'Dt';

            % Create DtEditField
            app.DtEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DtEditField.ValueDisplayFormat = '%.8f';
            app.DtEditField.Position = [185 262 59 17];

            % Create DpEditFieldLabel
            app.DpEditFieldLabel = uilabel(app.OptionsPanel);
            app.DpEditFieldLabel.HorizontalAlignment = 'right';
            app.DpEditFieldLabel.Position = [149 231 25 22];
            app.DpEditFieldLabel.Text = 'Dp';

            % Create DpEditField
            app.DpEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DpEditField.ValueDisplayFormat = '%.8f';
            app.DpEditField.Position = [185 233 59 18];

            % Create pFEditFieldLabel
            app.pFEditFieldLabel = uilabel(app.OptionsPanel);
            app.pFEditFieldLabel.HorizontalAlignment = 'right';
            app.pFEditFieldLabel.Position = [149 203 25 22];
            app.pFEditFieldLabel.Text = 'pF';

            % Create pFEditField
            app.pFEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.pFEditField.ValueDisplayFormat = '%.8f';
            app.pFEditField.Position = [185 204 59 19];

            % Create ADCEditFieldLabel
            app.ADCEditFieldLabel = uilabel(app.OptionsPanel);
            app.ADCEditFieldLabel.HorizontalAlignment = 'right';
            app.ADCEditFieldLabel.Position = [144 177 30 22];
            app.ADCEditFieldLabel.Text = 'ADC';

            % Create ADCEditField
            app.ADCEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.ADCEditField.ValueDisplayFormat = '%.8f';
            app.ADCEditField.Position = [185 179 59 18];

            % Create CalculateButton
            app.CalculateButton = uibutton(app.OptionsPanel, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.Icon = fullfile(pathToMLAPP, 'icons', 'calculate.png');
            app.CalculateButton.Position = [414 165 98 35];
            app.CalculateButton.Text = 'Calculate';

            % Create ROIStyleDropDownLabel
            app.ROIStyleDropDownLabel = uilabel(app.OptionsPanel);
            app.ROIStyleDropDownLabel.HorizontalAlignment = 'right';
            app.ROIStyleDropDownLabel.Position = [354 252 56 22];
            app.ROIStyleDropDownLabel.Text = 'ROI Style';

            % Create ROIStyleDropDown
            app.ROIStyleDropDown = uidropdown(app.OptionsPanel);
            app.ROIStyleDropDown.Items = {'Rectangle', 'Polygon', 'Freehand', 'Assisted'};
            app.ROIStyleDropDown.Position = [425 252 84 20];
            app.ROIStyleDropDown.Value = 'Rectangle';

            % Create DrawROIButton
            app.DrawROIButton = uibutton(app.OptionsPanel, 'state');
            app.DrawROIButton.ValueChangedFcn = createCallbackFcn(app, @DrawROIButtonValueChanged, true);
            app.DrawROIButton.Icon = fullfile(pathToMLAPP, 'icons', 'draw.png');
            app.DrawROIButton.Text = 'DrawROI';
            app.DrawROIButton.Position = [420 212 94 35];

            % Create DrawMultipleROIButton
            app.DrawMultipleROIButton = uibutton(app.OptionsPanel, 'push');
            app.DrawMultipleROIButton.ButtonPushedFcn = createCallbackFcn(app, @DrawMultipleROIButtonPushed, true);
            app.DrawMultipleROIButton.Icon = fullfile(pathToMLAPP, 'icons', 'multidraw.png');
            app.DrawMultipleROIButton.Position = [256 211 158 35];
            app.DrawMultipleROIButton.Text = 'Draw Multiple ROI';

            % Create MultiRoiSwitch
            app.MultiRoiSwitch = uiswitch(app.OptionsPanel, 'slider');
            app.MultiRoiSwitch.Position = [325 183 45 20];

            % Create DfastEditFieldLabel
            app.DfastEditFieldLabel = uilabel(app.OptionsPanel);
            app.DfastEditFieldLabel.HorizontalAlignment = 'right';
            app.DfastEditFieldLabel.Position = [146 150 33 22];
            app.DfastEditFieldLabel.Text = 'Dfast';

            % Create DfastEditField
            app.DfastEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DfastEditField.ValueDisplayFormat = '%.8f';
            app.DfastEditField.Position = [190 152 59 18];

            % Create DinterEditFieldLabel
            app.DinterEditFieldLabel = uilabel(app.OptionsPanel);
            app.DinterEditFieldLabel.HorizontalAlignment = 'right';
            app.DinterEditFieldLabel.Position = [142 122 37 22];
            app.DinterEditFieldLabel.Text = 'Dinter';

            % Create DinterEditField
            app.DinterEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DinterEditField.ValueDisplayFormat = '%.8f';
            app.DinterEditField.Position = [190 124 59 18];

            % Create DslowEditField_5Label
            app.DslowEditField_5Label = uilabel(app.OptionsPanel);
            app.DslowEditField_5Label.HorizontalAlignment = 'right';
            app.DslowEditField_5Label.Position = [143 98 38 22];
            app.DslowEditField_5Label.Text = 'Dslow';

            % Create DslowEditField
            app.DslowEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.DslowEditField.ValueDisplayFormat = '%.8f';
            app.DslowEditField.Position = [192 100 59 18];

            % Create ImportFromLastDatabaseButton
            app.ImportFromLastDatabaseButton = uibutton(app.OptionsPanel, 'push');
            app.ImportFromLastDatabaseButton.ButtonPushedFcn = createCallbackFcn(app, @ImportFromLastDatabaseButtonPushed, true);
            app.ImportFromLastDatabaseButton.Position = [187 1 161 26];
            app.ImportFromLastDatabaseButton.Text = 'Import From Last Database';

            % Create pFfastEditFieldLabel
            app.pFfastEditFieldLabel = uilabel(app.OptionsPanel);
            app.pFfastEditFieldLabel.HorizontalAlignment = 'right';
            app.pFfastEditFieldLabel.Position = [143 73 38 22];
            app.pFfastEditFieldLabel.Text = 'pFfast';

            % Create pFfastEditField
            app.pFfastEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.pFfastEditField.ValueDisplayFormat = '%.8f';
            app.pFfastEditField.Position = [192 75 59 18];

            % Create pFinterEditFieldLabel
            app.pFinterEditFieldLabel = uilabel(app.OptionsPanel);
            app.pFinterEditFieldLabel.HorizontalAlignment = 'right';
            app.pFinterEditFieldLabel.Position = [140 52 42 22];
            app.pFinterEditFieldLabel.Text = 'pFinter';

            % Create pFinterEditField
            app.pFinterEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.pFinterEditField.ValueDisplayFormat = '%.8f';
            app.pFinterEditField.Position = [193 54 59 18];

            % Create pFslowEditFieldLabel
            app.pFslowEditFieldLabel = uilabel(app.OptionsPanel);
            app.pFslowEditFieldLabel.HorizontalAlignment = 'right';
            app.pFslowEditFieldLabel.Position = [139 32 43 22];
            app.pFslowEditFieldLabel.Text = 'pFslow';

            % Create pFslowEditField
            app.pFslowEditField = uieditfield(app.OptionsPanel, 'numeric');
            app.pFslowEditField.ValueDisplayFormat = '%.8f';
            app.pFslowEditField.Position = [193 34 59 18];

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
            app.AlgorithmButtonGroup.Position = [1015 149 100 136];

            % Create MonoExpButton
            app.MonoExpButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.MonoExpButton.Text = 'MonoExp';
            app.MonoExpButton.Position = [8 92 73 22];
            app.MonoExpButton.Value = true;

            % Create BiexpFreeButton
            app.BiexpFreeButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.BiexpFreeButton.Text = 'BiexpFree';
            app.BiexpFreeButton.Position = [8 71 77 22];

            % Create BiExpSegButton
            app.BiExpSegButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.BiExpSegButton.Text = 'BiExpSeg';
            app.BiExpSegButton.Position = [8 48 75 22];

            % Create BayesianButton
            app.BayesianButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.BayesianButton.Text = 'Bayesian';
            app.BayesianButton.Position = [7 26 71 22];

            % Create TriExpButton
            app.TriExpButton = uiradiobutton(app.AlgorithmButtonGroup);
            app.TriExpButton.Text = 'TriExp';
            app.TriExpButton.Position = [7 5 56 22];

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
            app.GenerateMapsButton.Position = [986 57 153 51];
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

            % Create SmoothingDegreeEditFieldLabel
            app.SmoothingDegreeEditFieldLabel = uilabel(app.CalculatorTab);
            app.SmoothingDegreeEditFieldLabel.HorizontalAlignment = 'right';
            app.SmoothingDegreeEditFieldLabel.Position = [1346 131 105 22];
            app.SmoothingDegreeEditFieldLabel.Text = 'Smoothing Degree';

            % Create SmoothingDegreeEditField
            app.SmoothingDegreeEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.SmoothingDegreeEditField.Position = [1475 131 51 22];

            % Create LocationEditFieldLabel
            app.LocationEditFieldLabel = uilabel(app.CalculatorTab);
            app.LocationEditFieldLabel.HorizontalAlignment = 'right';
            app.LocationEditFieldLabel.Position = [745 265 50 22];
            app.LocationEditFieldLabel.Text = 'Location';

            % Create LocationEditField
            app.LocationEditField = uieditfield(app.CalculatorTab, 'numeric');
            app.LocationEditField.Position = [810 262 62 27];

            % Create UseCortexROIFromSegmenterButton
            app.UseCortexROIFromSegmenterButton = uibutton(app.CalculatorTab, 'push');
            app.UseCortexROIFromSegmenterButton.ButtonPushedFcn = createCallbackFcn(app, @UseCortexROIFromSegmenterButtonPushed, true);
            app.UseCortexROIFromSegmenterButton.Enable = 'off';
            app.UseCortexROIFromSegmenterButton.Position = [18 49 194 32];
            app.UseCortexROIFromSegmenterButton.Text = 'Use Cortex ROI From Segmenter';

            % Create GenerateMapsWithSegmenterROIButton
            app.GenerateMapsWithSegmenterROIButton = uibutton(app.CalculatorTab, 'push');
            app.GenerateMapsWithSegmenterROIButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateMapsWithSegmenterROIButtonPushed, true);
            app.GenerateMapsWithSegmenterROIButton.Enable = 'off';
            app.GenerateMapsWithSegmenterROIButton.Position = [6 144 219 35];
            app.GenerateMapsWithSegmenterROIButton.Text = 'Generate Maps With Segmenter ROI';

            % Create UseMedullaROIFromSegmenterButton
            app.UseMedullaROIFromSegmenterButton = uibutton(app.CalculatorTab, 'push');
            app.UseMedullaROIFromSegmenterButton.ButtonPushedFcn = createCallbackFcn(app, @UseMedullaROIFromSegmenterButtonPushed, true);
            app.UseMedullaROIFromSegmenterButton.Enable = 'off';
            app.UseMedullaROIFromSegmenterButton.Position = [15 89 201 32];
            app.UseMedullaROIFromSegmenterButton.Text = 'Use Medulla ROI From Segmenter';

            % Create CopyDataButton
            app.CopyDataButton = uibutton(app.CalculatorTab, 'push');
            app.CopyDataButton.ButtonPushedFcn = createCallbackFcn(app, @CopyDataButtonPushed, true);
            app.CopyDataButton.Position = [1041 313 140 31];
            app.CopyDataButton.Text = 'Copy Data';

            % Create DicomTab
            app.DicomTab = uitab(app.TabGroup);
            app.DicomTab.Title = 'Dicom';

            % Create PreviewAxes
            app.PreviewAxes = uiaxes(app.DicomTab);
            title(app.PreviewAxes, 'Title')
            xlabel(app.PreviewAxes, 'X')
            ylabel(app.PreviewAxes, 'Y')
            zlabel(app.PreviewAxes, 'Z')
            app.PreviewAxes.Position = [1066 494 470 329];

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
            app.LoadtoCalculatorButton.Position = [1346 65 138 26];
            app.LoadtoCalculatorButton.Text = 'Load to Calculator';

            % Create LoadtoViewerButton
            app.LoadtoViewerButton = uibutton(app.SequencesPanel, 'push');
            app.LoadtoViewerButton.ButtonPushedFcn = createCallbackFcn(app, @LoadtoViewerButtonPushed, true);
            app.LoadtoViewerButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.LoadtoViewerButton.Position = [1338 97 156 26];
            app.LoadtoViewerButton.Text = 'Load to Viewer';

            % Create FramesEditFieldLabel
            app.FramesEditFieldLabel = uilabel(app.SequencesPanel);
            app.FramesEditFieldLabel.HorizontalAlignment = 'right';
            app.FramesEditFieldLabel.Position = [1340 282 46 22];
            app.FramesEditFieldLabel.Text = 'Frames';

            % Create FramesEditField
            app.FramesEditField = uieditfield(app.SequencesPanel, 'numeric');
            app.FramesEditField.Position = [1401 277 93 32];

            % Create LoadtoBOLDPcsButton
            app.LoadtoBOLDPcsButton = uibutton(app.SequencesPanel, 'push');
            app.LoadtoBOLDPcsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadtoBOLDPcsButtonPushed, true);
            app.LoadtoBOLDPcsButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.LoadtoBOLDPcsButton.Position = [1328 35 177 24];
            app.LoadtoBOLDPcsButton.Text = 'Load to BOLD Pcs.';

            % Create TimeEditFieldLabel
            app.TimeEditFieldLabel = uilabel(app.SequencesPanel);
            app.TimeEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeEditFieldLabel.Position = [1322 231 32 22];
            app.TimeEditFieldLabel.Text = 'Time';

            % Create TimeEditField
            app.TimeEditField = uieditfield(app.SequencesPanel, 'text');
            app.TimeEditField.Position = [1368 231 138 22];

            % Create SaveDICOMDatabaseButton
            app.SaveDICOMDatabaseButton = uibutton(app.SequencesPanel, 'push');
            app.SaveDICOMDatabaseButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDICOMDatabaseButtonPushed, true);
            app.SaveDICOMDatabaseButton.Icon = fullfile(pathToMLAPP, 'icons', 'import.png');
            app.SaveDICOMDatabaseButton.Position = [1328 179 166 28];
            app.SaveDICOMDatabaseButton.Text = 'Save DICOM Database';

            % Create LoadDICOMDatabaseButton
            app.LoadDICOMDatabaseButton = uibutton(app.SequencesPanel, 'push');
            app.LoadDICOMDatabaseButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDICOMDatabaseButtonPushed, true);
            app.LoadDICOMDatabaseButton.Icon = fullfile(pathToMLAPP, 'icons', 'import.png');
            app.LoadDICOMDatabaseButton.Position = [1328 152 166 23];
            app.LoadDICOMDatabaseButton.Text = 'Load DICOM Database';

            % Create PatientsPanel
            app.PatientsPanel = uipanel(app.DicomTab);
            app.PatientsPanel.Title = 'Patients';
            app.PatientsPanel.Position = [1 466 1056 360];

            % Create PatientsListBoxLabel
            app.PatientsListBoxLabel = uilabel(app.PatientsPanel);
            app.PatientsListBoxLabel.HorizontalAlignment = 'right';
            app.PatientsListBoxLabel.Position = [13 307 48 22];
            app.PatientsListBoxLabel.Text = 'Patients';

            % Create PatientsListBox
            app.PatientsListBox = uilistbox(app.PatientsPanel);
            app.PatientsListBox.ValueChangedFcn = createCallbackFcn(app, @PatientsListBoxValueChanged, true);
            app.PatientsListBox.Position = [76 12 965 319];

            % Create DicomState
            app.DicomState = uieditfield(app.DicomTab, 'text');
            app.DicomState.Position = [15 21 491 21];

            % Create ImageNumberSpinnerLabel
            app.ImageNumberSpinnerLabel = uilabel(app.DicomTab);
            app.ImageNumberSpinnerLabel.HorizontalAlignment = 'right';
            app.ImageNumberSpinnerLabel.Position = [1199 460 84 22];
            app.ImageNumberSpinnerLabel.Text = 'Image Number';

            % Create ImageNumberSpinner
            app.ImageNumberSpinner = uispinner(app.DicomTab);
            app.ImageNumberSpinner.ValueChangingFcn = createCallbackFcn(app, @ImageNumberSpinnerValueChanging, true);
            app.ImageNumberSpinner.Position = [1298 460 100 22];

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
            app.UIAxes3.Position = [22 532 371 234];

            % Create UIAxes3_2
            app.UIAxes3_2 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_2, 'DT Free')
            xlabel(app.UIAxes3_2, 'X')
            ylabel(app.UIAxes3_2, 'Y')
            zlabel(app.UIAxes3_2, 'Z')
            app.UIAxes3_2.Position = [22 283 371 234];

            % Create UIAxes3_3
            app.UIAxes3_3 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_3, 'pF Free')
            xlabel(app.UIAxes3_3, 'X')
            ylabel(app.UIAxes3_3, 'Y')
            zlabel(app.UIAxes3_3, 'Z')
            app.UIAxes3_3.Position = [22 29 371 234];

            % Create UIAxes3_4
            app.UIAxes3_4 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_4, 'DP Seg-Fit')
            xlabel(app.UIAxes3_4, 'X')
            ylabel(app.UIAxes3_4, 'Y')
            zlabel(app.UIAxes3_4, 'Z')
            app.UIAxes3_4.Position = [431 532 371 234];

            % Create UIAxes3_5
            app.UIAxes3_5 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_5, 'DT Seg-Fit')
            xlabel(app.UIAxes3_5, 'X')
            ylabel(app.UIAxes3_5, 'Y')
            zlabel(app.UIAxes3_5, 'Z')
            app.UIAxes3_5.Position = [431 283 371 234];

            % Create UIAxes3_6
            app.UIAxes3_6 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_6, 'pF Seg-Fit')
            xlabel(app.UIAxes3_6, 'X')
            ylabel(app.UIAxes3_6, 'Y')
            zlabel(app.UIAxes3_6, 'Z')
            app.UIAxes3_6.Position = [431 29 371 234];

            % Create UIAxes3_7
            app.UIAxes3_7 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_7, 'ADC')
            xlabel(app.UIAxes3_7, 'X')
            ylabel(app.UIAxes3_7, 'Y')
            zlabel(app.UIAxes3_7, 'Z')
            app.UIAxes3_7.Position = [1220 283 293 250];

            % Create UIAxes3_8
            app.UIAxes3_8 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_8, 'DP Bayes Fit')
            xlabel(app.UIAxes3_8, 'X')
            ylabel(app.UIAxes3_8, 'Y')
            zlabel(app.UIAxes3_8, 'Z')
            app.UIAxes3_8.Position = [828 532 371 234];

            % Create UIAxes3_9
            app.UIAxes3_9 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_9, 'DT Bayes Fit')
            xlabel(app.UIAxes3_9, 'X')
            ylabel(app.UIAxes3_9, 'Y')
            zlabel(app.UIAxes3_9, 'Z')
            app.UIAxes3_9.Position = [828 283 371 234];

            % Create UIAxes3_10
            app.UIAxes3_10 = uiaxes(app.DWIMAPSPanel);
            title(app.UIAxes3_10, 'pF Bayes Fit')
            xlabel(app.UIAxes3_10, 'X')
            ylabel(app.UIAxes3_10, 'Y')
            zlabel(app.UIAxes3_10, 'Z')
            app.UIAxes3_10.Position = [828 29 371 234];

            % Create ColorMapsDropDownLabel
            app.ColorMapsDropDownLabel = uilabel(app.DWIMAPSPanel);
            app.ColorMapsDropDownLabel.HorizontalAlignment = 'right';
            app.ColorMapsDropDownLabel.Position = [1328 22 66 22];
            app.ColorMapsDropDownLabel.Text = 'Color Maps';

            % Create ColorMapsDropDown
            app.ColorMapsDropDown = uidropdown(app.DWIMAPSPanel);
            app.ColorMapsDropDown.ValueChangedFcn = createCallbackFcn(app, @ColorMapsDropDownValueChanged, true);
            app.ColorMapsDropDown.Position = [1409 22 100 22];

            % Create SaveMapsButton
            app.SaveMapsButton = uibutton(app.DWIMAPSPanel, 'push');
            app.SaveMapsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveMapsButtonPushed, true);
            app.SaveMapsButton.Position = [1319 85 190 36];
            app.SaveMapsButton.Text = 'Save Maps';

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
            app.UIAxes12.Position = [275 19 280 201];

            % Create MeanR2msec1EditFieldLabel
            app.MeanR2msec1EditFieldLabel = uilabel(app.R2MapandFittingPanel);
            app.MeanR2msec1EditFieldLabel.HorizontalAlignment = 'right';
            app.MeanR2msec1EditFieldLabel.Position = [548 13 111 22];
            app.MeanR2msec1EditFieldLabel.Text = 'Mean R2*(msec^-1)';

            % Create MeanR2msec1EditField
            app.MeanR2msec1EditField = uieditfield(app.R2MapandFittingPanel, 'numeric');
            app.MeanR2msec1EditField.Position = [668 13 99 22];

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

            % Create RegisterTab
            app.RegisterTab = uitab(app.TabGroup);
            app.RegisterTab.Title = 'Register';

            % Create Panel
            app.Panel = uipanel(app.RegisterTab);
            app.Panel.Title = 'Panel';
            app.Panel.Position = [1 1 1537 823];

            % Create UIAxes13
            app.UIAxes13 = uiaxes(app.Panel);
            title(app.UIAxes13, 'Title')
            xlabel(app.UIAxes13, 'X')
            ylabel(app.UIAxes13, 'Y')
            zlabel(app.UIAxes13, 'Z')
            app.UIAxes13.Position = [42 357 399 432];

            % Create UIAxes13_2
            app.UIAxes13_2 = uiaxes(app.Panel);
            title(app.UIAxes13_2, 'Title')
            xlabel(app.UIAxes13_2, 'X')
            ylabel(app.UIAxes13_2, 'Y')
            zlabel(app.UIAxes13_2, 'Z')
            app.UIAxes13_2.Position = [566 357 399 432];

            % Create UIAxes13_3
            app.UIAxes13_3 = uiaxes(app.Panel);
            title(app.UIAxes13_3, 'Title')
            xlabel(app.UIAxes13_3, 'X')
            ylabel(app.UIAxes13_3, 'Y')
            zlabel(app.UIAxes13_3, 'Z')
            app.UIAxes13_3.Position = [1085 355 399 432];

            % Create FixedListBoxLabel
            app.FixedListBoxLabel = uilabel(app.Panel);
            app.FixedListBoxLabel.HorizontalAlignment = 'right';
            app.FixedListBoxLabel.Position = [47 170 34 22];
            app.FixedListBoxLabel.Text = 'Fixed';

            % Create FixedListBox
            app.FixedListBox = uilistbox(app.Panel);
            app.FixedListBox.ValueChangedFcn = createCallbackFcn(app, @FixedListBoxValueChanged, true);
            app.FixedListBox.Position = [96 7 339 187];

            % Create FixedImageButton
            app.FixedImageButton = uibutton(app.Panel, 'push');
            app.FixedImageButton.ButtonPushedFcn = createCallbackFcn(app, @FixedImageButtonPushed, true);
            app.FixedImageButton.Position = [460 756 88 31];
            app.FixedImageButton.Text = ' Fixed Image';

            % Create MovingImageButton
            app.MovingImageButton = uibutton(app.Panel, 'push');
            app.MovingImageButton.ButtonPushedFcn = createCallbackFcn(app, @MovingImageButtonPushed, true);
            app.MovingImageButton.Position = [978 755 95 31];
            app.MovingImageButton.Text = 'Moving Image';

            % Create MovingListBoxLabel
            app.MovingListBoxLabel = uilabel(app.Panel);
            app.MovingListBoxLabel.HorizontalAlignment = 'right';
            app.MovingListBoxLabel.Position = [568 172 47 22];
            app.MovingListBoxLabel.Text = 'Moving ';

            % Create MovingListBox
            app.MovingListBox = uilistbox(app.Panel);
            app.MovingListBox.ValueChangedFcn = createCallbackFcn(app, @MovingListBoxValueChanged, true);
            app.MovingListBox.Position = [630 9 339 187];

            % Create FixedSliceSliderLabel
            app.FixedSliceSliderLabel = uilabel(app.Panel);
            app.FixedSliceSliderLabel.HorizontalAlignment = 'right';
            app.FixedSliceSliderLabel.Position = [7 325 64 22];
            app.FixedSliceSliderLabel.Text = 'Fixed Slice';

            % Create FixedSliceSlider
            app.FixedSliceSlider = uislider(app.Panel);
            app.FixedSliceSlider.Limits = [1 100];
            app.FixedSliceSlider.ValueChangingFcn = createCallbackFcn(app, @FixedSliceSliderValueChanging, true);
            app.FixedSliceSlider.Position = [92 334 298 3];
            app.FixedSliceSlider.Value = 1;

            % Create MovingSliceSliderLabel
            app.MovingSliceSliderLabel = uilabel(app.Panel);
            app.MovingSliceSliderLabel.HorizontalAlignment = 'right';
            app.MovingSliceSliderLabel.Position = [559 316 73 22];
            app.MovingSliceSliderLabel.Text = 'Moving Slice';

            % Create MovingSliceSlider
            app.MovingSliceSlider = uislider(app.Panel);
            app.MovingSliceSlider.Limits = [1 100];
            app.MovingSliceSlider.ValueChangingFcn = createCallbackFcn(app, @MovingSliceSliderValueChanging, true);
            app.MovingSliceSlider.Position = [653 325 298 3];
            app.MovingSliceSlider.Value = 1;

            % Create RegisterButton
            app.RegisterButton = uibutton(app.Panel, 'push');
            app.RegisterButton.ButtonPushedFcn = createCallbackFcn(app, @RegisterButtonPushed, true);
            app.RegisterButton.Icon = fullfile(pathToMLAPP, 'icons', 'register.png');
            app.RegisterButton.Position = [1296 58 144 30];
            app.RegisterButton.Text = 'Register';

            % Create SendRegisteredImagetoSegmenterButton
            app.SendRegisteredImagetoSegmenterButton = uibutton(app.Panel, 'push');
            app.SendRegisteredImagetoSegmenterButton.ButtonPushedFcn = createCallbackFcn(app, @SendRegisteredImagetoSegmenterButtonPushed, true);
            app.SendRegisteredImagetoSegmenterButton.Icon = fullfile(pathToMLAPP, 'icons', 'send.png');
            app.SendRegisteredImagetoSegmenterButton.Position = [1269 158 217 30];
            app.SendRegisteredImagetoSegmenterButton.Text = 'Send Registered Image to Segmenter';

            % Create LoadButton
            app.LoadButton = uibutton(app.Panel, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Icon = fullfile(pathToMLAPP, 'icons', 'import.png');
            app.LoadButton.Position = [441 156 78 29];
            app.LoadButton.Text = 'Load';

            % Create LoadButton_2
            app.LoadButton_2 = uibutton(app.Panel, 'push');
            app.LoadButton_2.ButtonPushedFcn = createCallbackFcn(app, @LoadButton_2Pushed, true);
            app.LoadButton_2.Icon = fullfile(pathToMLAPP, 'icons', 'import.png');
            app.LoadButton_2.Position = [985 157 78 29];
            app.LoadButton_2.Text = 'Load';

            % Create ShowPairButton
            app.ShowPairButton = uibutton(app.Panel, 'push');
            app.ShowPairButton.ButtonPushedFcn = createCallbackFcn(app, @ShowPairButtonPushed, true);
            app.ShowPairButton.Icon = fullfile(pathToMLAPP, 'icons', 'pair.png');
            app.ShowPairButton.Position = [1393 312 116 26];
            app.ShowPairButton.Text = 'Show Pair';

            % Create LocationEditField_2Label
            app.LocationEditField_2Label = uilabel(app.Panel);
            app.LocationEditField_2Label.HorizontalAlignment = 'right';
            app.LocationEditField_2Label.Position = [332 258 50 22];
            app.LocationEditField_2Label.Text = 'Location';

            % Create LocationEditField_2
            app.LocationEditField_2 = uieditfield(app.Panel, 'numeric');
            app.LocationEditField_2.Position = [397 251 64 36];

            % Create LocationEditField_3Label
            app.LocationEditField_3Label = uilabel(app.Panel);
            app.LocationEditField_3Label.HorizontalAlignment = 'right';
            app.LocationEditField_3Label.Position = [828 257 50 22];
            app.LocationEditField_3Label.Text = 'Location';

            % Create LocationEditField_3
            app.LocationEditField_3 = uieditfield(app.Panel, 'numeric');
            app.LocationEditField_3.Position = [893 250 64 36];

            % Create CoveredDepthEditFieldLabel
            app.CoveredDepthEditFieldLabel = uilabel(app.Panel);
            app.CoveredDepthEditFieldLabel.HorizontalAlignment = 'right';
            app.CoveredDepthEditFieldLabel.Position = [47 258 86 22];
            app.CoveredDepthEditFieldLabel.Text = 'Covered Depth';

            % Create CoveredDepthEditField
            app.CoveredDepthEditField = uieditfield(app.Panel, 'numeric');
            app.CoveredDepthEditField.Position = [148 257 47 23];

            % Create CoveredDepthEditField_2Label
            app.CoveredDepthEditField_2Label = uilabel(app.Panel);
            app.CoveredDepthEditField_2Label.HorizontalAlignment = 'right';
            app.CoveredDepthEditField_2Label.Position = [580 258 86 22];
            app.CoveredDepthEditField_2Label.Text = 'Covered Depth';

            % Create CoveredDepthEditField_2
            app.CoveredDepthEditField_2 = uieditfield(app.Panel, 'numeric');
            app.CoveredDepthEditField_2.Position = [681 257 47 23];

            % Create FilterButton
            app.FilterButton = uibutton(app.Panel, 'push');
            app.FilterButton.ButtonPushedFcn = createCallbackFcn(app, @FilterButtonPushed, true);
            app.FilterButton.Position = [222 215 105 24];
            app.FilterButton.Text = 'Filter ';

            % Create StatusEditFieldLabel
            app.StatusEditFieldLabel = uilabel(app.Panel);
            app.StatusEditFieldLabel.HorizontalAlignment = 'right';
            app.StatusEditFieldLabel.Position = [1188 10 39 22];
            app.StatusEditFieldLabel.Text = 'Status';

            % Create RegisterStatusEditField
            app.RegisterStatusEditField = uieditfield(app.Panel, 'text');
            app.RegisterStatusEditField.Position = [1242 9 285 24];

            % Create DrawROIButton_2
            app.DrawROIButton_2 = uibutton(app.Panel, 'push');
            app.DrawROIButton_2.ButtonPushedFcn = createCallbackFcn(app, @DrawROIButton_2Pushed, true);
            app.DrawROIButton_2.Position = [455 664 97 25];
            app.DrawROIButton_2.Text = 'Draw ROI';

            % Create DrawROIButton_3
            app.DrawROIButton_3 = uibutton(app.Panel, 'push');
            app.DrawROIButton_3.ButtonPushedFcn = createCallbackFcn(app, @DrawROIButton_3Pushed, true);
            app.DrawROIButton_3.Position = [974 664 97 25];
            app.DrawROIButton_3.Text = 'Draw ROI';

            % Create FixedSegmentButton
            app.FixedSegmentButton = uibutton(app.Panel, 'push');
            app.FixedSegmentButton.ButtonPushedFcn = createCallbackFcn(app, @FixedSegmentButtonPushed, true);
            app.FixedSegmentButton.Position = [456 704 98 24];
            app.FixedSegmentButton.Text = 'Fixed Segment';

            % Create MovingSegmentButton
            app.MovingSegmentButton = uibutton(app.Panel, 'push');
            app.MovingSegmentButton.ButtonPushedFcn = createCallbackFcn(app, @MovingSegmentButtonPushed, true);
            app.MovingSegmentButton.Position = [973 704 105 24];
            app.MovingSegmentButton.Text = 'Moving Segment';

            % Create SegmentGuidedRegisterButton
            app.SegmentGuidedRegisterButton = uibutton(app.Panel, 'push');
            app.SegmentGuidedRegisterButton.ButtonPushedFcn = createCallbackFcn(app, @SegmentGuidedRegisterButtonPushed, true);
            app.SegmentGuidedRegisterButton.Icon = fullfile(pathToMLAPP, 'icons', 'register.png');
            app.SegmentGuidedRegisterButton.Position = [1256 93 253 30];
            app.SegmentGuidedRegisterButton.Text = 'Segment Guided Register';

            % Create SegmenterTab
            app.SegmenterTab = uitab(app.TabGroup);
            app.SegmenterTab.Title = 'Segmenter';

            % Create Segmenter
            app.Segmenter = uipanel(app.SegmenterTab);
            app.Segmenter.Title = 'Segmenter';
            app.Segmenter.Position = [1 1 1537 823];

            % Create UIAxes14
            app.UIAxes14 = uiaxes(app.Segmenter);
            title(app.UIAxes14, 'Title')
            xlabel(app.UIAxes14, 'X')
            ylabel(app.UIAxes14, 'Y')
            zlabel(app.UIAxes14, 'Z')
            app.UIAxes14.Position = [32 221 634 570];

            % Create UIAxes14_2
            app.UIAxes14_2 = uiaxes(app.Segmenter);
            title(app.UIAxes14_2, 'Title')
            xlabel(app.UIAxes14_2, 'X')
            ylabel(app.UIAxes14_2, 'Y')
            zlabel(app.UIAxes14_2, 'Z')
            app.UIAxes14_2.Position = [891 221 634 570];

            % Create SegmentButton
            app.SegmentButton = uibutton(app.Segmenter, 'push');
            app.SegmentButton.ButtonPushedFcn = createCallbackFcn(app, @SegmentButtonPushed, true);
            app.SegmentButton.Icon = fullfile(pathToMLAPP, 'icons', 'segment.png');
            app.SegmentButton.Position = [682 167 190 40];
            app.SegmentButton.Text = 'Segment';

            % Create DrawCortexGroundTruthButton
            app.DrawCortexGroundTruthButton = uibutton(app.Segmenter, 'push');
            app.DrawCortexGroundTruthButton.ButtonPushedFcn = createCallbackFcn(app, @DrawCortexGroundTruthButtonPushed, true);
            app.DrawCortexGroundTruthButton.Enable = 'off';
            app.DrawCortexGroundTruthButton.Position = [233 147 223 41];
            app.DrawCortexGroundTruthButton.Text = 'Draw Cortex Ground Truth';

            % Create JaccardCortexEditFieldLabel
            app.JaccardCortexEditFieldLabel = uilabel(app.Segmenter);
            app.JaccardCortexEditFieldLabel.HorizontalAlignment = 'right';
            app.JaccardCortexEditFieldLabel.Position = [923 156 86 22];
            app.JaccardCortexEditFieldLabel.Text = 'Jaccard Cortex';

            % Create JaccardCortexEditField
            app.JaccardCortexEditField = uieditfield(app.Segmenter, 'numeric');
            app.JaccardCortexEditField.Position = [1024 156 61 22];

            % Create DICECortexEditFieldLabel
            app.DICECortexEditFieldLabel = uilabel(app.Segmenter);
            app.DICECortexEditFieldLabel.HorizontalAlignment = 'right';
            app.DICECortexEditFieldLabel.Position = [1161 156 72 22];
            app.DICECortexEditFieldLabel.Text = 'DICE Cortex';

            % Create DICECortexEditField
            app.DICECortexEditField = uieditfield(app.Segmenter, 'numeric');
            app.DICECortexEditField.Position = [1248 156 61 22];

            % Create JaccardMedullaEditFieldLabel
            app.JaccardMedullaEditFieldLabel = uilabel(app.Segmenter);
            app.JaccardMedullaEditFieldLabel.HorizontalAlignment = 'right';
            app.JaccardMedullaEditFieldLabel.Position = [917 108 92 22];
            app.JaccardMedullaEditFieldLabel.Text = 'Jaccard Medulla';

            % Create JaccardMedullaEditField
            app.JaccardMedullaEditField = uieditfield(app.Segmenter, 'numeric');
            app.JaccardMedullaEditField.Position = [1024 108 61 22];

            % Create DICEMedullaEditFieldLabel
            app.DICEMedullaEditFieldLabel = uilabel(app.Segmenter);
            app.DICEMedullaEditFieldLabel.HorizontalAlignment = 'right';
            app.DICEMedullaEditFieldLabel.Position = [1154 106 79 22];
            app.DICEMedullaEditFieldLabel.Text = 'DICE Medulla';

            % Create DICEMedullaEditField
            app.DICEMedullaEditField = uieditfield(app.Segmenter, 'numeric');
            app.DICEMedullaEditField.Position = [1248 106 61 22];

            % Create CalculateSimilarityCoefficientsButton
            app.CalculateSimilarityCoefficientsButton = uibutton(app.Segmenter, 'push');
            app.CalculateSimilarityCoefficientsButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateSimilarityCoefficientsButtonPushed, true);
            app.CalculateSimilarityCoefficientsButton.Position = [249 79 183 23];
            app.CalculateSimilarityCoefficientsButton.Text = 'Calculate Similarity Coefficients';

            % Show the figure after all components are created
            app.DiffusonCalculatorforSIEMENSUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DiffAppV02reduced_exported

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