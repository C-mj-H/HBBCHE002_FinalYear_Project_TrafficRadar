%% Ground Truths Dataset 6p1
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
MaxRangeLine = (1/3.6)*MaxRangeLine
%MaxRangeLine = (1/4)*MaxRangeLine
% Plot range lines
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Range lines: after Eq Notch','fontsize',fontsize1);

%% Detection Point Calculations for Dataset 6
%final variable of detection points
Target_locations_DS6_p1 = zeros(NumOfPulses, NumOfRangeBins);

% Parameters for sections of data to screenshot
stepSize = 100;
numIterations = ceil(NumOfPulses/stepSize)
segment1end = 1035;
Iterations_segment1 = ceil(segment1end/stepSize);


for z = 1:numIterations
    if z <numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):stepSize*z,1:100);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);

    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS6_p1(1+stepSize*(z-1):stepSize*z,1:100) = finalTargetPnts;

elseif z == numIterations
    RangeProfiles = RangeProfiles_AfterEqNotch(1+stepSize*(z-1):end,1:100);
    [NumOfPulses,NumOfRangeBins]=size(RangeProfiles);
    
    finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine);
    Target_locations_DS6_p1(1+stepSize*(z-1):end,1:100) = finalTargetPnts; 
    end
end


% Manual removal of noise data points registered as targets 
Target_locations_DS6_p1(5000:end, 1:35) = 0;
Target_locations_DS6_p1(:, 275:end) = 0;
Target_locations_DS6_p1(:, 49:71) = 0;
Target_locations_DS6_p1(1:2250, 29:98) = 0;
Target_locations_DS6_p1(3500:end , 126:152) = 0;
Target_locations_DS6_p1(5400:5650, 89:93) = 0;
Target_locations_DS6_p1(4886:5110, 92:97) = 0;
Target_locations_DS6_p1(4205:4600, 90:100) = 0;
Target_locations_DS6_p1(1:3368, 35:90) = 0;
Target_locations_DS6_p1(1:3197, 34:91) = 0;
Target_locations_DS6_p1(1:4500, 40:85) = 0;
Target_locations_DS6_p1(5400:5700, 88:93) = 0;
Target_locations_DS6_p1(4887:5109, 92:94) = 0;
Target_locations_DS6_p1(4288:4599, 96:98) = 0;
Target_locations_DS6_p1(3947:4444, 98:100) = 0;
Target_locations_DS6_p1(4595:5109, 92:98) = 0;
Target_locations_DS6_p1(5245:5648, 89:90) = 0;
Target_locations_DS6_p1(1:4597, 42:84) = 0;
Target_locations_DS6_p1(1:3987, 37:47) = 0;
Target_locations_DS6_p1(4798, 47) = 0;
Target_locations_DS6_p1(2379, 95) = 0;
Target_locations_DS6_p1(2620, 94) = 0;
Target_locations_DS6_p1(2418, 94) = 0;
Target_locations_DS6_p1(5197, 72) = 0;
Target_locations_DS6_p1(4927, 79) = 0;
Target_locations_DS6_p1(1447, 99) = 0;
Target_locations_DS6_p1(1:1870, 98:101) = 0;
Target_locations_DS6_p1(2960:end, 1:26) = 0;


% % Save the dataset to a .mat file
save('Target_locations_DS6_p1.mat', 'Target_locations_DS6_p1');

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

[pulseCoord, rangeCoord] = find(Target_locations_DS6_p1 == 1);

% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');


function finalTargetPnts = detectionFinder(RangeProfiles,NumOfPulses,NumOfRangeBins, MaxRangeLine) % function to determine target locations
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

%     hue_threshold = [0.1, 0.8]; % Adjust as needed
% saturation_threshold = [0.5, 1]; % Adjust as needed

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