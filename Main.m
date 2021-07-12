%% Locomotor Resilience Codes %%
% A method to calculate recovery time of COM movement after a large perturbation
% Public Version: 1, Date: 12.07.2021
% Contributed authors: Ravi Deepak, Marc Bartholet, Caroline Heimhofer

% Ravi DK et al. 2021, "Rhythmic auditory stimuli modulate movement recovery
% in response to perturbation during locomotion" DOI: 10.1242/jeb.237073
% or
% Ravi DK et al. 2021, "Adapting footfall rhythmicity to auditory perturbations
% affects resilience of locomotor behavior: a proof-of-concept study" DOI: 10.3389/fnins.2021.678965

% Please cite one of the above work if you are using this code!

% For clarification, suggestion and/or collaboration contact:
% Deepak Ravi (depakroshanblu@gmail.com) or (deepak.ravi@hest.ethz.ch)

% Deepak Ravi's current affiliation:
% PhD student, Laboratory of Movement Biomechanics, ETH Zurich, Switzerland


%%
clear all
clc
close all

% Add path to functions
addpath(genpath('...\functions'));
% Add path to data
addpath(genpath('...\example data'));
% Add path to other info
addpath(genpath('...\other info'));

% Path to the data folder
MainDir = '...\example data';
SubjectName = dir(MainDir);
SubjectName = SubjectName(3:end);

% Load data from other info
Tau = xlsread('General Data_V01.xlsx', 'B2:B21'); % Pre calculated time delay for state space reconstruction
Dim = xlsread('General Data_V01.xlsx', 'F2:F21'); % Pre calculated dimension for state space reconstruction

%%

for pp = 2 %Participant AMC002
    %% --------------------------------------first part - stability boundaries construction ------------------------------ %%
    %%
    pp_dir = strcat(MainDir,'\', SubjectName(pp).name, '\S01\mat');
    cd(fullfile(pp_dir));
    
    % Determination of steady state COM using baseline walking 1 data
    [XoUncued, YoUncued, ZoUncued, NewxpUncued,NewypUncued,NewzpUncued,...
        ellipse3d_01, ellipse3d_02, ellipse3d_03, normalVector] = ssr_uncued(Tau(pp), Dim(pp));
    
    % The vertical displacement time series of the sacrum marker from cued walking data
    
    load('cue_pert.mat')
    VD.COM_scalar = VD.SACR(:,3);
    remove = 5 * VD.SF; % remove start-up effects
    
    Frame_Start_Perturb = 3*60*VD.SF; % perturbation at 3min into trial
    VD.COM_scalar = VD.COM_scalar - mean(VD.COM_scalar);% demeaning
    
    %% Low pass filtering of data using a 4th order Butterworth filter
    
    forder = 4;
    cutfreq = 5;
    srate  = VD.SF;
    
    beforefilter = VD.COM_scalar;
    [b,a] = butter(forder, cutfreq / srate);
    afterfilter = filtfilt(b, a, beforefilter);
    VD.COM_scalar = afterfilter;
    
    %% Splitting data into pre- and post-perturbation
    
    data_cued = VD.COM_scalar(remove+1:end,:);
    Frame_Start_Perturb = Frame_Start_Perturb - remove;
    
    FirstNormalWalking = data_cued(1:Frame_Start_Perturb); % data_cued from before perturbation
    PerturbedWalking = data_cued(Frame_Start_Perturb:end); % all data_cued from start of perturbation to the end
    
    % filling nans with empty spaces
    FirstNormalWalking(isnan(FirstNormalWalking)) = [];
    PerturbedWalking(isnan(PerturbedWalking)) = [];
    data_cued(isnan(data_cued)) = [];
    
    %% State Space Reconstruction - MATLAB
    
    [XR,tau,dim] = phaseSpaceReconstruction(data_cued, Tau(pp), Dim(pp));
    
    if tau ~= Tau(pp) || dim ~= Dim(pp)
        fprintf('something went wrong with Tau or Dim!')
    end
    
    % data_cued sorted corresponding to phaseSpaceReconstruction divided into "x,y,z-axis"
    Original = XR(:,1); % x-axis
    Delayed1 = XR(:,2); % y-axis
    Delayed2 = XR(:,3); % z-axis
    
    % Data of Pre-Perturbation Walking
    PreOriginal = Original(1: length(FirstNormalWalking) - tau);
    PreDelayed1 = Delayed1(1: length(FirstNormalWalking) - tau);
    PreDelayed2 = Delayed2(1: length(FirstNormalWalking) - tau);
    
    
    %% Reference Trajectory
    
    xAxis = [1,0,0];
    yAxis = [0,1,0];
    zAxis = [0,0,1];
    
    % Fit Plane in Pre-Perturbation Attractor
    [fitresult4, ~] = fct_createFit1(PreOriginal, PreDelayed1, PreDelayed2);
    normalFitresult4 = [-fitresult4.p10,-fitresult4.p01,1]; % direction of plane normal vector
    
    % Rotation of Pre-Perturbation Attractor to be parallel to xy-Plane
    rotationVector5 = vrrotvec(zAxis,normalFitresult4);
    rotationMatrix5 = vrrotvec2mat(rotationVector5);
    
    ProjData1 = zeros(length(PreOriginal),3);
    ProjData2 = zeros(length(Original),3);
    
    for i = 1:length(PreOriginal)
        [ProjData1(i,:)] = [PreOriginal(i)-XoUncued,PreDelayed1(i)-YoUncued,PreDelayed2(i)-ZoUncued] * rotationMatrix5;
    end
    
    Projxp1 = ProjData1(1:end,1); % Pre-Perturbation data projected to be parallel to xy-Plane and centered around Centroid
    Projyp1 = ProjData1(1:end,2);
    Projzp1 = ProjData1(1:end,3);
    
    for i = 1:length(Original)
        [ProjData2(i,:)] = [Original(i)-XoUncued,Delayed1(i)-YoUncued,Delayed2(i)-ZoUncued] * rotationMatrix5;
    end
    
    Projxp2 = ProjData2(1:end,1); % Full data projected to be parallel to xy-Plane and centered around Centroid
    Projyp2 = ProjData2(1:end,2);
    Projzp2 = ProjData2(1:end,3);
    
    % Calculation of Phaseangles of Pre-Perturbation Attractor
    phaseAnglePre = zeros(1,length(PreOriginal));
    for i = 1: length(PreOriginal)
        phaseAnglePre(i) = atan2d(Projyp1(i),Projxp1(i));
    end
    
    phaseAnglePre = phaseAnglePre + 180;
    phaseAnglePre = phaseAnglePre.';
    
    % Calculation of Phaseangles of Attractor
    phaseAngle = zeros(1,length(Original));
    for i = 1: length(Original)
        phaseAngle(i) = atan2d(Projyp2(i),Projxp2(i));
    end
    
    phaseAngle = phaseAngle + 180;
    phaseAngle = phaseAngle.';
    
    resolution1 = 0; % or 1 or 2, resolution for rounding
    resolution2 = 360; % or 3600 or 36000; resolution for splitting the phaseangles
    resolution3 = 50; % any integer greater 3; resolution for the ellipses
    
    % Sort Pre-Perturbation Data according to Phaseangles
    % Round Phaseangles according to 'resolution1'
    
    [phaseAnglePreRS, Idx1] = sort(phaseAnglePre); % sorts phase angles (for every data point) and gives index
    phaseAnglePreRS = round(phaseAnglePreRS,resolution1);
    phaseAnglePreR = round(phaseAnglePre);
    
    PreOriginalS = PreOriginal(Idx1);
    PreDelayed1S = PreDelayed1(Idx1);
    PreDelayed2S = PreDelayed2(Idx1);
    
    [phaseAngleRS, Idx4] = sort(phaseAngle);
    phaseAngleRS = round(phaseAngleRS,resolution1);
    phaseAngleR = round(phaseAngle,resolution1);
    
    OriginalS = Original(Idx4);
    Delayed1S = Delayed1(Idx4);
    Delayed2S = Delayed2(Idx4);
    
    phaseAngleRef = linspace(0,360-360/resolution2,resolution2); % creates a vector with resolution2 points.
    phaseAngleRef = round(phaseAngleRef, resolution1);
    
    PreDataPhaseAngleRS = [PreOriginalS PreDelayed1S PreDelayed2S phaseAnglePreRS];
    PreDataPhaseAngleRSstd = [PreOriginalS PreDelayed1S PreDelayed2S phaseAnglePreRS];
    for i = 1:length(PreOriginalS)
        if  PreDataPhaseAngleRSstd(i,4) == resolution2
            PreDataPhaseAngleRSstd(i,4) = 0;
        end
    end
    [PreDataPhaseAngleRSstd(:,4), Idx] = sort(PreDataPhaseAngleRSstd(:,4));
    PreDataPhaseAngleRSstd(:,3) = PreDataPhaseAngleRSstd(Idx,3);
    PreDataPhaseAngleRSstd(:,2) = PreDataPhaseAngleRSstd(Idx,2);
    PreDataPhaseAngleRSstd(:,1) = PreDataPhaseAngleRSstd(Idx,1);
    
    [stdx,stdy,stdz] = fct_std_RefTraj(PreOriginalS,phaseAngleRef,PreDataPhaseAngleRSstd,NewxpUncued,NewypUncued,NewzpUncued,resolution2);
    
    
    %% Ellipse Torus around Pre-Perturbation Attractor
    
    % Find maximal and second-maximal Deviation at each Phaseangle from Reference Trajectory and Pre-Perturbation Attractor
    
    [maxDeviation] = fct_findEllipseDimension(phaseAngleRef,PreDataPhaseAngleRS,PreOriginalS,PreDelayed1S,PreDelayed2S,...
        NewxpUncued,NewypUncued,NewzpUncued,resolution2);
    
    for i = 1:length(phaseAngleRef) % creating the ellipses with max std as one and secondmax std as second axis
        maxDeviation(i,4) = max([stdx(i),stdy(i),stdz(i)]);
    end
    
    % Rotation of Ellipse around Rference Trajectory to align the maximal Deviations
    [~,maxDeviationellipse02,maxDeviationProj02,~] = fct_ellipseOrientation(ellipse3d_02,NewxpUncued,NewypUncued,NewzpUncued,...
        resolution2,resolution3,normalVector,maxDeviation);
    [~,maxDeviationellipse01,maxDeviationProj01,~] = fct_ellipseOrientation(ellipse3d_01,NewxpUncued,NewypUncued,NewzpUncued,...
        resolution2,resolution3,normalVector,maxDeviation);
    [~,maxDeviationellipse03,maxDeviationProj03,~] = fct_ellipseOrientation(ellipse3d_03,NewxpUncued,NewypUncued,NewzpUncued,...
        resolution2,resolution3,normalVector,maxDeviation);
    
    %Plot Ellipse Torus around Pre-Perturbation Attractor
    % 2nd std Ellipses
    j = 1;
    f101 = figure(101);
    hold on;
    scatter3(XoUncued,YoUncued,ZoUncued,100,'.','r');
    for i =1:resolution2
        if j > (resolution2*3)
            break;
        else
            %plot3(ellipse3d_01(:,j), ellipse3d_01(:,j+1), ellipse3d_01(:,j+2),'Color',[0, 0, 1],'LineWidth',0.2);
            plot3(ellipse3d_02(:,j), ellipse3d_02(:,j+1), ellipse3d_02(:,j+2),'Color',[237/255, 179/255, 20/255],'LineWidth',0.8);
            %plot3(ellipse3d_03(:,j), ellipse3d_03(:,j+1), ellipse3d_03(:,j+2),'Color',[1, 0.7, 0],'LineWidth',0.2);
            j = j + 3;
        end
    end
    %plot3(PreOriginal,PreDelayed1,PreDelayed2,'Color',[0, 1, 0]);
    %plot3(Original(Frame_Start_Perturb:end), Delayed1(Frame_Start_Perturb:end), Delayed2(Frame_Start_Perturb:end), 'Color',[0, 0.9, 0]) ;
    plot3(Original(17500:17700), Delayed1(17500:17700), Delayed2(17500:17700), 'Color',[0, 0, 0]) ;
    plot3(NewxpUncued,NewypUncued,NewzpUncued,'Color', [164/255, 91/255, 170/255], 'LineWidth', 2);
    
    %     title('Stability Tube','Fontsize', 12);
    zlabel({'COM Position delayed by 2*Tau'});
    ylabel({'COM Position delayed by Tau'});
    xlabel({'COM Position'});
    xlim([-25 25])
    ylim([-25 25])
    zlim([-15 20])
    %     axis equal
    %     legend({'Reference Trajectory and Centroid','Pre-Perturbation Attractor', 'Post-Perturbation Attractor'},'Location','northeast','Fontsize', 8)
    
    view(125,-10);
    hold off;
    
    %% --------------------------------------second part - recovery time analysis ------------------------------
    %% Create Stability Margin
    
    % Find closest Point on Ellipse to Attractor Point
    fitresult7 = fct_createFit1(NewxpUncued,NewypUncued,NewzpUncued);
    normalFitresult7 = [-fitresult7.p10,-fitresult7.p01,1];
    
    gamma = zeros(1,360);
    for k = 1:360
        maxDevPoint = [maxDeviationellipse02(k,1)-NewxpUncued(k),maxDeviationellipse02(k,2)-NewypUncued(k),maxDeviationellipse02(k,3)-NewzpUncued(k)];
        normmaxDevPoint = norm(maxDevPoint);
        if acosd(dot(maxDevPoint,normalFitresult7)/(normmaxDevPoint*norm(normalFitresult7))) < 80 % why choose 80?
            gamma(k) = 1;
        else
            gamma(k) = 0;
        end
    end
    
    % Find closest Point on Ellipse to Attractor Point
    DataPhaseAngleRS = [OriginalS Delayed1S Delayed2S phaseAngleRS zeros(length(OriginalS),7)];
    [minNormVector01,DataPhaseAngleRS01] = fct_ellipse2dataNorm(phaseAngleRef,DataPhaseAngleRS,Original,OriginalS,Delayed1S,Delayed2S,...
        ellipse3d_01,NewxpUncued,NewypUncued,NewzpUncued,resolution2,resolution3,gamma,maxDeviationellipse01,normalVector);
    [minNormVector02,DataPhaseAngleRS02] = fct_ellipse2dataNorm(phaseAngleRef,DataPhaseAngleRS,Original,OriginalS,Delayed1S,Delayed2S,...
        ellipse3d_02,NewxpUncued,NewypUncued,NewzpUncued,resolution2,resolution3,gamma,maxDeviationellipse02,normalVector);
    [minNormVector03,~] = fct_ellipse2dataNorm(phaseAngleRef,DataPhaseAngleRS,Original,OriginalS,Delayed1S,Delayed2S,...
        ellipse3d_03,NewxpUncued,NewypUncued,NewzpUncued,resolution2,resolution3,gamma,maxDeviationellipse03,normalVector);
    
    %% Stability Analysis
    
    % Calulation correction factor so that Pre-Perturbation Attractor never leaves Ellipse Torus (i.e. instable)
    resolution4 = 1.0; % constant
    resolution5 = 0.1; % or 0.1 or 0.001 or 0.0001 and so on
    
    % Calculation of
    % 1. Vector which determines if Point is outside or inside of Ellipse Torus
    % 2. Vector containing the Differences squared between the Reference
    % Trajectory and the Attractor Point
    % 3. Vector containing the Differences squared between the Ellipse Torus
    % and the Attractor Point
    % 4. Vector containing the perctentual Deviation between the Ellipse Torus
    % and the Attractor Point
    
    for i = 1:1/resolution5
        [stabilityVector_01,stabilityVector_01_02,instabilityVector_02_03,instabilityVector_03,stabilityVector,instabilityVectorTube,...
            instabilityVector,instabilityVectorPct,stabilityVector1,instabilityVector1,stabilityVector1Tube,instabilityVector1Tube,instabilityVectorTubePct,stabilityVectorSum] =...
            fct_stability02(Idx4,phaseAngleRef,DataPhaseAngleRS,PreOriginal,Original,OriginalS,Delayed1S,Delayed2S,NewxpUncued,NewypUncued,NewzpUncued,...
            minNormVector02,minNormVector01,minNormVector03,resolution2,resolution4);
        % with fct_stability02 -> stable if inside 2nd std ellipse
        % with fct_stability03 -> stable if inside 3rd std ellipse
        
        if stabilityVectorSum == length(PreOriginal) % this goes to max
            break;
        else
            resolution4 = resolution4 + resolution5;
        end
    end
    
    instabilitySum(pp,1) = length(stabilityVector)-stabilityVectorSum;
    
    
    %% Instabilities in Gait Cycle
    
    k = 1;
    n = 1;
    gaitcycleExit = zeros(length(Original),1); % Vector marking points (indices of original) where gaitcycle becomes stable again
    gaitcycleEntry = zeros(length(Original),1); % Vector marking points (indices of original) where gaitcycle becomes instable
    for i = 1:length(Original)-1
        if isnan(instabilityVector_03(i)) == 0 && isnan(instabilityVector_03(i+1)) == 1 % end instability
            gaitcycleExit(k) = i+1;
            k = k + 1;
        elseif isnan(instabilityVector_03(i)) == 1 && isnan(instabilityVector_03(i+1)) == 0 % start instability
            gaitcycleEntry(n) = i;
            n = n + 1;
        end
    end
    
    gaitcycleExit = gaitcycleExit(1:find(gaitcycleExit, 1, 'last'));
    gaitcycleEntry = gaitcycleEntry(1:find(gaitcycleEntry, 1, 'last'));
    
    %% find sequence of stability
    
    seqstability{pp,1} = findseq(stabilityVector); %Vector with 1 if "stable"(i.e. inside 2nd std) and 0 if "instable"
    % Char_resilience{pp,1} = char_resilience(Frame_Start_Perturb, seqstability{pp,1}, VD.SF, instabilityVector_03); % characterises perturbation
    
    
    %% import cue here and find first beat after perturbation to start recovery time analysis
    pp_dir = strcat(MainDir,'\', SubjectName(pp).name, '\S01\cue');
    cd(fullfile(pp_dir));
    
    [cue, Fs_cue]= audioread('cue_pert.wav');
    
    % pre-processing of cue
    cue100 = downsample(cue, 441); % downsampling to 100 Hz
    for i = 1:length(cue100) % rectify and mark points of cues
        if cue100(i) ~=0
            cue100(i) = 10;
        end
    end
    for i=2:length(cue100)-1 % if cue crossed 0 line at the exact 100Hz sampling point, it resulted in a 0 here -> i.e. makes those points also = 10
        if cue100(i-1)~=0 && cue100(i+1)~=0
            cue100(i) = 10;
        end
    end
    cue100_pks=cue100;
    for i=2:length(cue100)
        if cue100(i-1) == 0 && cue100(i)~= 0
            cue100_pks(i) = cue100(i)+1;
        end
    end
    
    cue100 = cue100(remove:end); % remove first part of cue
    cue100_pks = cue100_pks(remove:end);
    
    [~, pklocs] = findpeaks(cue100_pks); % find peaks of the cue
    pkpert = find(pklocs >=Frame_Start_Perturb); % find first peak during perturbation
    Frame_Stop_Perturb = pklocs(pkpert(5)); % 4 peaks later is the first peak after the perturbation
    
    
    %% Error Between Reference and Attractor Trajectory
    
    % figure;
    f130= figure(130);
    plot(instabilityVector,'Color',[0.403921568627451 0.537254901960784 0.631372549019608]);
    hold on;
    plot(stabilityVector_01,'.','Color',[229/255 78/255 209/255]); % makes points
    plot(stabilityVector_01_02,'.','Color',[237/255 179/255 20/255]);
    plot(instabilityVector_02_03,'.','Color',[191/255 216/255 52/255]);
    plot(instabilityVector_03,'.','Color',[0/255 7/255 97/255]);
    
    plot(stabilityVector_01,'LineWidth',1.5,'Color',[229/255 78/255 209/255]); %'Color',[0 0 0.6] % connects points
    plot(stabilityVector_01_02,'LineWidth',1.5,'Color',[237/255 179/255 20/255]); %'Color',[0 1 0]
    plot(instabilityVector_02_03,'LineWidth',1.5,'Color',[191/255 216/255 52/255]); %'Color',[0 0 1]
    plot(instabilityVector_03,'LineWidth',1.5,'Color',[0/255 7/255 97/255]); %'Color',[1 0 0]
    
    
    ylabel({'Distance between Reference Trajectory and Attractor Trajectory'},...
        'FontSize',11);
    xlabel('Frames [100*sec]','FontSize',11);
    axis([-100 length(stabilityVector_01)+100 -0.1 max(instabilityVector_03)+0.3])
    hold on;
    
    %% Start Recovery from max peak after perturbation
    
    for i = 1:length(instabilityVectorTubePct) % remove the 0s in the beginning
        if instabilityVectorTubePct(i) == 0
            instabilityVectorTubePct(i) = NaN;
        end
    end
    
    %% changed from instabilityVectorTubePct to instabilityVector_03 for paper 2 rhythm perturbations
    [~,Frame_startRecovery] = max(instabilityVector_03(Frame_Start_Perturb:(Frame_Start_Perturb+1000),1)); % max deviation after perturbation to start recovery
    Frame_startRecovery = Frame_Start_Perturb + Frame_startRecovery;
    
    %% Calculation of Gait Speed
    
    leftheel = VD.LHEE (:,3);
    rightheel = VD.RHEE (:,3);
    Fs = VD.SF;
    leftheelShort = leftheel(30*Fs:90*Fs);
    %     findpeaks(leftheelShort,'MinPeakHeight',0.2);
    [pks,locs] = findpeaks(leftheelShort,'MinPeakHeight',150);
    
    %     findpeaks(leftheelShort);
    %     [pks,locs] = findpeaks(leftheelShort);
    
    k = 1;
    gaitspeedCalc = zeros(1,length(locs));
    for i = 1:length(locs)-1
        gaitspeedCalc(k) = locs(i+1)-locs(i);
        k = k + 1;
    end
    
    gaitspeed = mean(gaitspeedCalc(1:length(gaitspeedCalc)-1)); % number of frames
    %take peaks of leftheel trajectory, difference between peaks gives you
    %speed for one gaitcycle (divided by frames)
    
    %% Calculation of Recovery Time, based on 5 gait cycles
    
    clear recoveryTime;
    gaitspeed5 = round(gaitspeed*5); %5 includes five gait cycles (in frames)
    
    stabilityVectorPert = stabilityVector(Frame_startRecovery:end); % recovery Time starts withing 10s after perturbation
    for i = 1:length(stabilityVectorPert)-gaitspeed5
        if  stabilityVectorPert(i) == 1    % Vector with 1 if "stable"(i.e. inside 2nd std) and 0 if "instable"
            k=0;
            for j= 1: gaitspeed5
                if stabilityVectorPert(i+j) ~= 1
                    
                    k=k+1;
                end
            end
            if k<4
                recoveryTime(pp) = i/VD.SF;%(i+gaitspeed5)/VD.SF;
                break;
            end
        end
    end
    
    
    X = ['The recovery time is ',num2str(recoveryTime(pp)), ' sec'];
    disp(X)
    
    
end





