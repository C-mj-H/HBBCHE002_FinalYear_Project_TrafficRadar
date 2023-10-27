%%Ground Truths Dataset 4

%%Clear workspace and all figures

clear all;clc;close all

%% Load Data
[FileName,Path2RadarData,filter_index]=uigetfile('','Select Radar Dataset');
RadarData = load([Path2RadarData filesep FileName]);

%% Extract Range Profiles before, after Equalisation and after Notch filtering

RangeProfiles_BeforeEq = RadarData.RangeLines_BeforeEq;
RangeProfiles_AfterEq = RadarData.RangeLines_AfterEq;
RangeProfiles_AfterEqNotch = RadarData.RangeLines_AfterEQ_Notch;

%% Extract other radar parameters

 PRF_Hz = RadarData.Info.PRF_Hz;
 Bandwidth_Hz = RadarData.Info.Bandwidth_Hz;
 RangeStart_m = RadarData.Info.RangeStart_m;
 BlindRange_m = RadarData.Info.BlindRange_m;
 
[NumOfPulses,NumOfRangeBins]=size(RangeProfiles_AfterEqNotch);

%% Plot Range Profiles

fontsize1 = 12;
clims = [-40 0];

% Normalise data to have a peak of 0dB or 1 in linear scale
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));
MaxRangeLine = (1/6)*MaxRangeLine
% Plot range lines
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Range lines: after Eq Notch','fontsize',fontsize1);

%% Detection Point Calculations for Dataset 5
%final variable of detection points
Target_locations_DS5 = zeros(NumOfPulses, NumOfRangeBins);

% Parameters for sections of data to screenshot
stepSize = 50;
numIterations = ceil(NumOfPulses/stepSize)
segment1end = 1035;
Iterations_segment1 = ceil(segment1end/stepSize);

%loops for first 2 target 

for z = 1:numIterations
    if z <numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):stepSize*z,418:480);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);

    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS5(1+stepSize*(z-1):stepSize*z,418:480) = finalTargetPnts;

elseif z == numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):end,418:480);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);
    
    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS5(1+stepSize*(z-1):end,418:480) = finalTargetPnts; 
    end
end


%loop for end last target

for z = 1:numIterations
    if z <numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):stepSize*z,792:798);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);

    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS5(1+stepSize*(z-1):stepSize*z,792:798) = finalTargetPnts;

elseif z == numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):end,792:798);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);
    
    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS5(1+stepSize*(z-1):end,792:798) = finalTargetPnts; 
    end
end


% Manual removal of noise data points registered as targets 
Target_locations_DS5(823, 420) = 0;
Target_locations_DS5(1066:1452, 464:480) = 0;
Target_locations_DS5(2472:2950, 460:480) = 0;
Target_locations_DS5(1:129, 419:443) = 0;
Target_locations_DS5(2001, 419) = 0;
Target_locations_DS5(3063:end, 445) = 0;
Target_locations_DS5(3110:end, 444) = 0;
Target_locations_DS5(2883:end, 448) = 0;
Target_locations_DS5(3080:end, 449) = 0;
Target_locations_DS5(2120:end, 454) = 0;
Target_locations_DS5(2410:end, 452) = 0;
Target_locations_DS5(2180:end, 459) = 0;
Target_locations_DS5(2050:end, 457) = 0;
Target_locations_DS5(2030:end, 455) = 0;
Target_locations_DS5(2241, 451) = 0;
Target_locations_DS5(2876, 448) = 0;
Target_locations_DS5(2787, 441) = 0;
Target_locations_DS5(2754:end, 442) = 0;

% % Save the dataset to a .mat file
save('Target_locations_DS5.mat', 'Target_locations_DS5');

%plotting
fontsize1 = 12;
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));
clims = [-40 0];
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Detections','fontsize',fontsize1);
hold on

[pulseCoord, rangeCoord] = find(Target_locations_DS5 == 1);

% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');


function finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine) % function to determine target locations
fontsize1 = 12;
    clims = [-40 0];
    % Normalise data to have a peak of 0dB or 1 in linear scale
    %[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles)));
    
    % Plot range lines
    figure; axes('fontsize',fontsize1);
    imagesc(20*log10(abs(RangeProfiles)./MaxRangeLine),clims);
    colormap('hot');
    colorbar;
    xlabel('Range (bins)','fontsize',fontsize1);
    ylabel('Number of pulses','fontsize',fontsize1);
    title('Range lines: after Eq Notch','fontsize',fontsize1);
    
    % Get the current figure handl
    fig = gcf;
    frame = getframe(fig);
    % Convert the captured frame to an RGB image:
    rgb_image = frame.cdata;
    % Convert RGB image to HSV color space:
    hsv_image = rgb2hsv(rgb_image);
    hsv_image = hsv_image(41:466, 83:562, :);
    
    
    
    hue_threshold = [0, 0.9]; % Adjust as needed
    saturation_threshold = [0, 0.8]; % Adjust as needed

    % Create binary mask based on hue and saturation thresholds:
    target_pnts = (hsv_image(:,:,1) >= hue_threshold(1) & hsv_image(:,:,1) <= hue_threshold(2)) & ...
              (hsv_image(:,:,2) >= saturation_threshold(1) & hsv_image(:,:,2) <= saturation_threshold(2));
    
    %imshow(target_pnts)

    target_pnts = double(target_pnts);
    [numtargetRows,numtargetCols] = size(target_pnts);
    
    target_pnts_resized = zeros(numtargetRows, NumOfRangeBins);
    for i = 1:numtargetRows
        targetsOriginalRowsize = 1:numel(target_pnts(i,:));
        targetsNewRowSize = linspace(1, numtargetCols, NumOfRangeBins);
        target_pnts_resized(i, :) = interp1(targetsOriginalRowsize, target_pnts(i,:), targetsNewRowSize, 'linear');
    end
    
    target_pnts_resized = abs(target_pnts_resized) ==1;
    target_pnts_resized = double(target_pnts_resized);
    finalTargetPnts = zeros(NumOfPulses, NumOfRangeBins); 
    
    for i = 1:NumOfRangeBins
        targetsOriginalColsize = 1:numel(target_pnts_resized(:,i));
        targetsNewRowSize = linspace(1, numtargetRows, NumOfPulses);
        finalTargetPnts(:, i) = interp1(targetsOriginalColsize, target_pnts_resized(:,i), targetsNewRowSize, 'linear');
    end
    
    finalTargetPnts = abs(finalTargetPnts) ==1;
    finalTargetPnts = double(finalTargetPnts);
end