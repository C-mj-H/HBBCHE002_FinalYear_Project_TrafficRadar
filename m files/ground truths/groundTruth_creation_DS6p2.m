%% Ground Truths Dataset 6p2
% ground truths for dataset 6 is split into three parts because there are
% more targets and running the entire thing in one script caused MATLAB to
% crash

%% Clear workspace and all figures

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
MaxRangeLine = (1/3.5)*MaxRangeLine
% Plot range lines
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Range lines: after Eq Notch','fontsize',fontsize1);

%% Detection Point Calculations for Dataset 6
%final variable of detection points
Target_locations_DS6_p2 = zeros(NumOfPulses, NumOfRangeBins);

% Parameters for sections of data to screenshot
stepSize = 100;
numIterations = ceil(NumOfPulses/stepSize)
segment1end = 1035;
Iterations_segment1 = ceil(segment1end/stepSize);

%only finds ground truths in first third of the dataset

for z = 1:numIterations
    if z <numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):stepSize*z,101:200);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);

    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS6_p2(1+stepSize*(z-1):stepSize*z,101:200) = finalTargetPnts;

elseif z == numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):end,101:200);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);
    
    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS6_p2(1+stepSize*(z-1):end,101:200) = finalTargetPnts; 
    end
end

% Manual removal of noise data points registered as targets 
Target_locations_DS6_p2(4213:end, 118:155) = 0;
Target_locations_DS6_p2(1:3275, 158:177) = 0;
Target_locations_DS6_p2(3496:end, 172:end) = 0;
Target_locations_DS6_p2(1:1042, 121:136) = 0;
Target_locations_DS6_p2(1925:3275, 138:183) = 0;
Target_locations_DS6_p2(3173:4168, 131:152) = 0;
Target_locations_DS6_p2(1525:2416, 110:124) = 0;
Target_locations_DS6_p2(3408:4001, 101:109) = 0;
Target_locations_DS6_p2(2712:3197, 107:115) = 0;
Target_locations_DS6_p2(1137:1645, 122:130) = 0;
Target_locations_DS6_p2(32:339, 137:140) = 0;
Target_locations_DS6_p2(1:651, 184:end) = 0;
Target_locations_DS6_p2(1:2592, 198:end) = 0;
Target_locations_DS6_p2(2713:3379, 184:194) = 0;
Target_locations_DS6_p2(1940:2468, 184:189) = 0;
Target_locations_DS6_p2(5232:end, 109:116) = 0;
Target_locations_DS6_p2(3563:4123, 124:130) = 0;
Target_locations_DS6_p2(2811:3069, 133:137) = 0;
Target_locations_DS6_p2(2294:2571, 136) = 0;
Target_locations_DS6_p2(1578:1919, 142:153) = 0;
Target_locations_DS6_p2(1090:1470, 149:153) = 0;
Target_locations_DS6_p2(626:801, 154:157) = 0;
Target_locations_DS6_p2(3279:4209, 162:174) = 0;

% % Save the dataset to a .mat file
save('Target_locations_DS6_p2.mat', 'Target_locations_DS6_p2');

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

[pulseCoord, rangeCoord] = find(Target_locations_DS6_p2 == 1);

% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');


function finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine)
fontsize1 = 12;
    clims = [-40 0];
    
    
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
    hsv_image = hsv_image(41:466, 82:556, :);
    
    
   
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