%% CA-CFAR Porcessor
clear all;
close all;

%% Parameter Selection

PFA = 10^-3 ;
RefWindow = 34; %reference window length
GuardCells = 2; % Number of guard cells on each 

%% CA-CFAR Detection Verification in Time Domain

%Generation of simulated dataset
Length = 100000;  %sample length
y_complex = randn(1,Length) + 1i*randn(1,Length); %dataset

%CA CFAR function

%number of cells needed to wrap-around the dataset
wrapFactor = RefWindow/2 +GuardCells;

RangeProfiles=  horzcat(y_complex(:, end-wrapFactor+1:end), y_complex(:, 1:end), y_complex(:,1:wrapFactor)); %extension of dataset by wrapping around
[detections, false_alarm_count] = CA_CFAR_2D(RangeProfiles, PFA, RefWindow, GuardCells); %CA-CFAR algorithm run on dataset
detections = detections(:, wrapFactor+1:end-wrapFactor); %the results of the CA-CFAR algorthm converted back to the original dataset size

%PFA error verification
PFA_expected = PFA*Length;
PFA_simulation = false_alarm_count;
PFA_errorTime = (PFA_expected - PFA_simulation)/PFA_expected*100;
disp("Question 1")
disp("PFA expected is " + PFA_expected)
disp("PFA simulation is " + PFA_simulation)
disp("PFA error is " + PFA_errorTime + "%")


%% OS-CFAR Detections in the Time Domain


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

%%CA CFAR of each range bin
Length = NumOfRangeBins;
PFA = 10^-3;
%RefWindow = 32; %reference window length
%GuardCells = 2; % Number of guard cells on each
Detection_matrix1 = zeros(NumOfPulses, NumOfRangeBins);
detection_count_matrix1 = zeros(NumOfPulses, 1);

%number of cells needed to wrap-around the dataset
wrapFactor = RefWindow/2 +GuardCells;

RangeProfiles_AfterEq=  horzcat(RangeProfiles_AfterEqNotch(:, end-wrapFactor+1:end), RangeProfiles_AfterEqNotch(:, 1:end), RangeProfiles_AfterEqNotch(:,1:wrapFactor)); %extension of dataset by wrapping around
[Detection_matrix1, detection_count_matrix1] = CA_CFAR_2D(RangeProfiles_AfterEq, PFA, RefWindow, GuardCells);
Detection_matrix1 = Detection_matrix1(:, wrapFactor+1:end-wrapFactor);

%Plotting of results
fontsize1 = 12;
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));
clims = [-40 0];
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Time Domain Detections of CA-CFAR','fontsize',fontsize1);
hold on
[pulseCoord, rangeCoord] = find(Detection_matrix1 == 1);
% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');

hold off


%% CA-CFAR Detection Verification in Frequency Domain


%creation of simulated data
NumberOfPulses = 2000; 
NumberOfRangeBins = 10; 
RangeProfilesSimulated = randn(NumberOfPulses, NumberOfRangeBins) + 1i*randn(NumberOfPulses, NumberOfRangeBins);

%CA-CFAR algorthm along slow time 
Length = NumberOfPulses;
PFA = 10^-3;
% RefWindow = 32 ; %reference window length
% GuardCells = 2; % Number of guard cells on each
detection_matrix = zeros(NumberOfPulses, NumberOfRangeBins);
false_alarm_count_matrix = zeros(1, NumberOfRangeBins);

%Length of CPI
frameLength = 150;
frameNum = round(NumberOfPulses/frameLength); %Number of CPIs

%loop that runs through each CPI and implements the CA-CFAR algorithm along
%the slow-time dimesion
for n = 1:frameNum
    %wrap around required for first CPI
    if n == 1
        StartIdx = 1 + (n-1)*frameLength;
        StopIdx = frameLength + (n-1)*frameLength+ wrapFactor; %wraps around into the next CPI
        frame = RangeProfilesSimulated(StartIdx:StopIdx, :);
        frame = vertcat(RangeProfilesSimulated( end-wrapFactor+1:end, :), frame(1:end, :));%Wraps around into the last CPI
    %wrap-around for standard CPIs     
    elseif n == frameNum
        StartIdx = 1 + (n-1)*frameLength - wrapFactor +1; %wraps around into the previous CPI
        StopIdx = frameLength + (n-1)*frameLength;
        frame = RangeProfilesSimulated(StartIdx:StopIdx, :);
        frame = vertcat(frame(1:end, :), RangeProfilesSimulated(1:wrapFactor,:) ); %wraps around back into the first CPI
   %standard wrap-around for CPIs 
    else
        StartIdx = 1 + (n-1)*frameLength - wrapFactor +1;%wraps around into the previous CPI
        StopIdx = frameLength + (n-1)*frameLength+ wrapFactor;%wraps around into the next CPI
        frame = RangeProfilesSimulated(StartIdx:StopIdx, :);
    end 

    RangeDopplerSimulated = (fft(frame, [ ], 1)); %coversion of CPI to Doppler-range
    temp_detect_matrix_transposed = zeros(NumOfRangeBins, frameLength);
    temp_detect_count_transposed = zeros(NumOfRangeBins, 1);
    [temp_detect_matrix_transposed, temp_detect_count_transposed] = CA_CFAR_2D(transpose(RangeDopplerSimulated), PFA, RefWindow, GuardCells); %CA-CFAR along the slow dimension
    temp_detect_matrix = transpose(temp_detect_matrix_transposed); %coversion back to correct dimension. Initial transpose was needed so algorithm would run along slow dimension
    temp_detect_count = transpose(temp_detect_count_transposed);

    temp_detect_matrix = temp_detect_matrix(wrapFactor+1:end-wrapFactor, :);
    
    %Single detection in a range bin means entire range bin is designated
    %as target locations:
    columns_with_ones = any(temp_detect_matrix, 1); 
    detection_matrix(StartIdx:StopIdx, columns_with_ones) = 1;

    false_alarm_count_matrix= false_alarm_count_matrix + temp_detect_count;
    
end

%PFA error calculation
PFA_expected = PFA*Length;
PFA_simulation_avg = (sum(false_alarm_count_matrix))/length(false_alarm_count_matrix);
PFA_errorTime = (PFA_expected - PFA_simulation_avg)/PFA_expected*100;
disp("Question 3")
disp("PFA expected is " + PFA_expected)
disp("Average PFA simulation is " + PFA_simulation_avg)
disp("PFA error is " + PFA_errorTime + "%")



%% OS-CFAR Detections in the Frquency Domain




%CA-CFAR algorthm along slow time 
Length = NumOfPulses;
PFA = 10^-3;
%RefWindow = 16; %reference window length
%GuardCells = 2; % Number of guard cells on each
detection_matrix2 = zeros(NumOfPulses, NumOfRangeBins);
false_alarm_count_matrix2 = zeros(1, NumOfRangeBins);

%Length of CPI
frameLength = 128;
frameNum = floor(NumOfPulses/frameLength);

%loop that runs through each CPI and implements the CA-CFAR algorithm along
%the slow-time dimesion - SAME LOOP AS SIMULATED DATASET
for n = 1:frameNum
    if n == 1
        StartIdx = 1 + (n-1)*frameLength;
        StopIdx = frameLength + (n-1)*frameLength+ wrapFactor;
        frame = RangeProfiles_AfterEqNotch(StartIdx:StopIdx, :);
        frame = vertcat(RangeProfiles_AfterEqNotch( end-wrapFactor+1:end, :), frame(1:end, :));
    elseif n == frameNum
        StartIdx = 1 + (n-1)*frameLength - wrapFactor +1;
        StopIdx = frameLength + (n-1)*frameLength;
        frame = RangeProfiles_AfterEqNotch(StartIdx:StopIdx, :);
        frame = vertcat(frame(1:end, :), RangeProfiles_AfterEqNotch(1:wrapFactor,:) );
    else
        StartIdx = 1 + (n-1)*frameLength - wrapFactor +1;
        StopIdx = frameLength + (n-1)*frameLength+ wrapFactor;
        frame = RangeProfiles_AfterEqNotch(StartIdx:StopIdx, :);
    end 

    rangeDoppler = (fft(frame, [ ], 1));
    temp_detect_matrix_transposed = zeros(NumOfRangeBins, frameLength);
    temp_detect_count_transposed = zeros(NumOfRangeBins, 1);
    [temp_detect_matrix_transposed, temp_detect_count_transposed] = CA_CFAR_2D(transpose(rangeDoppler), PFA, RefWindow, GuardCells);
    temp_detect_matrix = transpose(temp_detect_matrix_transposed);
    temp_detect_count = transpose(temp_detect_count_transposed);

    temp_detect_matrix = temp_detect_matrix(wrapFactor+1:end-wrapFactor, :);

    columns_with_ones = any(temp_detect_matrix, 1);
    detection_matrix2(StartIdx:StopIdx, columns_with_ones) = 1;

    false_alarm_count_matrix2 = false_alarm_count_matrix2 + temp_detect_count;
    
end

%Plotting
fontsize1 = 12;
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));
clims = [-40 0];
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Frequency Detections of CA-CFAR','fontsize',fontsize1);
hold on

[pulseCoord, rangeCoord] = find(detection_matrix2 == 1);

% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');
hold off

%PFA error calculation
PFA_expected = PFA*Length;
PFA_simulation_avg = (sum(false_alarm_count_matrix2))/length(false_alarm_count_matrix2);
PFA_errorTime = (PFA_expected - PFA_simulation_avg)/PFA_expected*100;
disp("Question 4")
disp("PFA expected is " + PFA_expected)
disp("Average PFA simulation is " + PFA_simulation_avg)
disp("PFA error is " + PFA_errorTime + "%")

%% AND gate combination

%The combination of the detection algotthms to reduce false alarmsa:
% AND gate attempt
Detections = zeros(NumOfPulses, NumOfRangeBins);
for j = 1:NumOfRangeBins
    for i = 1: NumOfPulses
        if (Detection_matrix1(i,j) == detection_matrix2(i,j)) && (Detection_matrix1(i,j) == 1) 
            Detections(i,j) = 1;
        else
            Detections(i,j) = 0;
        end
    end
end

%Plotting
fontsize1 = 12;
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));
clims = [-40 0];
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('AND Gate Detections of CA-CFAR','fontsize',fontsize1);
hold on

[pulseCoord, rangeCoord] = find(Detections == 1);

% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');
hold off
%% Detection Performance Test
dataset = FileName(1:2)
dataset = str2num(dataset)

[numRows,numCols] = size(Detections);
Target_locations = zeros(numRows, numCols);

%loading of required dataset
if dataset == 1
    load("Target_locations_DS1.mat");
    Target_locations = Target_locations_DS1;
elseif dataset == 2
    load("Target_locations_DS2.mat");
    Target_locations = Target_locations_DS2;
elseif dataset == 3
    load("Target_locations_DS3.mat");
    Target_locations = Target_locations_DS3;
elseif dataset == 4
    load("Target_locations_DS4.mat");
    Target_locations = Target_locations_DS4;
elseif dataset == 5
    load("Target_locations_DS5.mat");
    Target_locations = Target_locations_DS5;
else 
    load("Target_locations_DS6.mat");
    Target_locations = Target_locations_DS6;
end


%detection performance test run on AND gate processing
[DetectionAccPerCellAND, DetectionAccPerTargetAND, false_alarms_perCellAND, missed_detections_perCellAND, detections_perTargetAND, missed_detections_perTargetAND, false_alarms_perTargetAND] = DetectionAccTest(Target_locations, Detections );
display("Detection accuracy per cell (AND): "+DetectionAccPerCellAND )
display("Detection accuracy per target (AND): "+DetectionAccPerTargetAND)

%detection performance test run on time domain processing
[DetectionAccPerCellTime, DetectionAccPerTargetTime, false_alarms_perCellTime, missed_detections_perCellTime, detections_perTargetTime, missed_detections_perTargetTime, false_alarms_perTargetTime] = DetectionAccTest(Target_locations, Detection_matrix1 );
display("Detection accuracy per cell (time domain): "+DetectionAccPerCellTime )
display("Detection accuracy per target (time domain): "+DetectionAccPerTargetTime )

%detection performance test run on freq domain processing
[DetectionAccPerCellFreq, DetectionAccPerTargetFreq, false_alarms_perCellFreq, missed_detections_perCellFreq, detections_perTargetFreq, missed_detections_perTargetFreq, false_alarms_perTargetFreq] = DetectionAccTest(Target_locations, detection_matrix2);
display("Detection accuracy per cell(freq domain): "+DetectionAccPerCellFreq)
display("Detection accuracy per target (freq domain): "+DetectionAccPerTargetFreq)



%trimming datasets to focus on target locations
if dataset == 1
    RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch(:,10:25);
    detections_perTargetAND = detections_perTargetAND(:,10:25);
    missed_detections_perTargetAND = missed_detections_perTargetAND(:,10:25);
end
if dataset == 2
    RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch(:,11:22);
    detections_perTargetAND = detections_perTargetAND(:,11:22);
    missed_detections_perTargetAND = missed_detections_perTargetAND(:,11:22);
end
if dataset == 3
    RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch(:,680:730);
    detections_perTargetAND = detections_perTargetAND(:,680:730);
    missed_detections_perTargetAND = missed_detections_perTargetAND(:,680:730);
end
if dataset == 4
    RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch(:,131:226);
    detections_perTargetAND = detections_perTargetAND(:,131:226);
    missed_detections_perTargetAND = missed_detections_perTargetAND(:,131:226);
end
if dataset == 5
    RangeProfiles_AfterEqNotch = horzcat(RangeProfiles_AfterEqNotch(:,415:476), RangeProfiles_AfterEqNotch(:,790:810));
    detections_perTargetAND = horzcat(detections_perTargetAND(:,415:476), detections_perTargetAND(:,790:810));
    missed_detections_perTargetAND = horzcat(missed_detections_perTargetAND(:,415:476), missed_detections_perTargetAND(:,790:810));
end
if dataset == 6
    RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch(:,1:280);
    detections_perTargetAND = detections_perTargetAND(:,1:280);
    missed_detections_perTargetAND = missed_detections_perTargetAND(:,1:280);
end

% plotting results
fontsize1 = 12;
clims = [-40 0];
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Detections and Missed Detections (AND Processing)','fontsize',fontsize1);
hold on

%detected locations are red
[pulseCoord, rangeCoord] = find(detections_perTargetAND == 1);
% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'rx');

%missed locations are white
[pulseCoord, rangeCoord] = find(missed_detections_perTargetAND == 1);
% Plot 'x' at positions where matrix value is 1
scatter(rangeCoord, pulseCoord, 'wx');


%% False Alarm Check
PFA_expected = PFA*NumOfPulses*NumOfRangeBins;
PFA_datasetTime = sum(false_alarms_perTargetTime(:));
PFA_errorTime = (PFA_expected - PFA_datasetTime)/PFA_expected*100;
PFA_datasetFreq = sum(false_alarms_perTargetFreq(:));
PFA_errorFreq = (PFA_expected - PFA_datasetFreq)/PFA_expected*100;
PFA_datasetAND = sum(false_alarms_perTargetAND(:));
PFA_errorAND = (PFA_expected - PFA_datasetAND)/PFA_expected*100;

disp("PFA Estimates")
disp("PFA expected is " + PFA_expected)
disp("*****************************************************")
disp("PFA of Dataset from time domain processing is " + PFA_datasetTime)
disp("PFA error from time domain processing is is " + PFA_errorTime + "%")
disp("*****************************************************")
disp("PFA of Dataset from frequency domain processing is " + PFA_datasetFreq)
disp("PFA error from frequency domain processing is is " + PFA_errorFreq + "%")
disp("*****************************************************")
disp("PFA of Dataset from AND gate processing is " + PFA_datasetAND)
disp("PFA error from AND processing is is " + PFA_errorAND + "%")


%% functions

function [detections, false_alarm_count] = CA_CFAR(Length, y_complex, PFA, RefWindow, GuardCells)%un-optimised CA-CFAR algorithm - one range line at a time
TrainCells = (RefWindow)/2; % Number of training cells in reference window
threshold = zeros(1, Length);
detections = zeros(1, Length);
false_alarm_count = 0;

%the CA-CFAR constant (α) 
alpha = RefWindow*(PFA^(-1/RefWindow) - 1);

%square law detector.
DataAfterPowerLawDetector = abs(y_complex).^2 ;

%the CA-CFAR constant (α) 
alpha = RefWindow*(PFA^(-1/RefWindow) - 1);
startIndx = RefWindow/2 + GuardCells+ 1;
stopIndx = Length - RefWindow/2 - GuardCells;

for i = (startIndx):(stopIndx)
    % Calculate average of reference cells in the leading and lagging windows
    reference_avg_leading = sum(DataAfterPowerLawDetector(i + GuardCells + 1 : i + GuardCells + RefWindow/2)) / TrainCells;
    reference_avg_lagging = sum(DataAfterPowerLawDetector(i - RefWindow/2 - GuardCells: i - GuardCells - 1)) / TrainCells;

    % Calculate threshold using scaling factor
    threshold(i) = alpha*((reference_avg_leading + reference_avg_lagging) / 2);
    %false alarms calculated
    if DataAfterPowerLawDetector(i) > threshold(i)
        false_alarm_count = false_alarm_count +1;
        detections(i) = 1;
    end
end

end

%optimised CA-CFAR algorithm - all range lines at a time
function [Detections, detection_count] = CA_CFAR_2D(input, PFA, RefWindowSize, GuardCells)
TrainCells = (RefWindowSize)/2; % Number of training cells in reference window
[numRows,numCols] = size(input);
threshold = zeros(numRows, numCols);
Detections = zeros(numRows, numCols);
detection_count = zeros(numRows, 1);

%the CA-CFAR constant (α) 
alpha = RefWindowSize*(PFA^(-1/RefWindowSize) - 1);

%square law detector.
DataAfterPowerLawDetector = abs(input).^2 ;

%the CA-CFAR constant (α) 
startIndx = RefWindowSize/2 + GuardCells+ 1;
stopIndx = numCols - RefWindowSize/2 - GuardCells;

for i = (startIndx):(stopIndx)
    % Calculate average of reference cells in the leading and lagging windows
    reference_avg_leading = sum(DataAfterPowerLawDetector(:, i + GuardCells + 1 : i + GuardCells + RefWindowSize/2),2) ./ TrainCells;
    reference_avg_lagging = sum(DataAfterPowerLawDetector(:, i - RefWindowSize/2 - GuardCells: i - GuardCells - 1),2) ./ TrainCells;

    % Calculate threshold using scaling factor
    threshold(:, i) = alpha.*((reference_avg_leading + reference_avg_lagging) ./ 2);

    %Detection calculated
    Detections(:,i) = DataAfterPowerLawDetector(:, i) > threshold(:, i);
    detection_count = detection_count + Detections(:,i);

 
end

end


%detection performance test
function [DetectionAccPerCell, DetectionAccPerTarget, false_alarms_perCell, missed_detections_perCell, detections_perTarget, missed_detections_perTarget, false_alarms_perTarget] = DetectionAccTest(Target_locations, Detections )

%variable creation
[numRows,numCols] = size(Target_locations);
num_detections_perCell = 0;
detections_perCell = Detections;
Target_locations_perCell = Target_locations;
missed_detections_perCell = zeros(numRows,numCols);
false_alarms_perCell = Detections;


num_detections_perTarget = 0;
detections_perTarget = Detections;
target_locations_perTarget = Target_locations;
missed_detections_perTarget = zeros(numRows,numCols);
false_alarms_perTarget = Detections;

%variable used to acount for target sidelobes so they are not assumed as
%missed detections
targetSize = 10;
% DS1 = 3
% DS2 = 3
% DS3 =10
% DS4 = 10
% DS5 = 10
% DS6 = 10

%this method does not account for target sidelobes
for i = 1:numRows
    for j = 1:numCols
        if (detections_perCell(i,j) == Target_locations_perCell(i,j))  && (detections_perCell(i,j) == 1)
            num_detections_perCell = num_detections_perCell +1;
            false_alarms_perCell(i,j) = 0;
            
        elseif Target_locations_perCell(i,j) == 1 && (detections_perCell(i,j) == 0)
            missed_detections_perCell(i,j) = 1;
        end
        
    end
   
end

%this method accounts for target sidelobes 
for i = 1:numRows
    for j = 1:numCols
        if (detections_perTarget(i,j) == target_locations_perTarget(i,j))  && (detections_perTarget(i,j) == 1)
            num_detections_perTarget = num_detections_perTarget +1;
            detections_perTarget(i, j+1:j+3) = 0; %removes target points in vicinity if the target as already been detected so it is not measured more than once
            target_locations_perTarget(i, j+1:j+targetSize) = 0; %removes ground truths in vicinity as this target has already been detected
            target_locations_perTarget(i, j-targetSize:j-1) = 0;
            false_alarms_perTarget(i, j-targetSize:j+targetSize) = 0; %removes false alarms as this target has been detected
        end
        
    end
end
for i = 1:numRows
 for j = 1:numCols
        if target_locations_perTarget(i,j) ==1 && (detections_perTarget(i,j) == 0)
            missed_detections_perTarget(i,j) = 1;
            target_locations_perTarget(i, j+1:j+targetSize) = 0;
            target_locations_perTarget(i, j-targetSize:j-1) = 0;
        end
 end
end



total_detections_perCell = sum(Target_locations(:) == 1);
total_detections_perTarget = sum(target_locations_perTarget(:) == 1);


DetectionAccPerCell = num_detections_perCell/total_detections_perCell*100;
DetectionAccPerTarget = num_detections_perTarget/total_detections_perTarget*100;

end
