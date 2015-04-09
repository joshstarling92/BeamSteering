% Initialize Workspace
clear all; close all; clc

%++++++++++++++++++++++++++  Algorithm Summary  +++++++++++++++++++++++++++
% Inputs required:
%   - Vehicle 1 GPS position (Lat, Lon, Height(m))
%   - Vehicle 2 GPS position (Lat, Lon, Height(m))
%   - Radar produced vector of relative position between vehicles(m)
%
%
%   *** ***  Detection Method  *** ***
%   The beauty of this detection method is its simplicity.  Given the GPS
%   postions of both vehicles and the relative position vector (RPV) 
%   between the two, a simple cost function comparison of the GPS RPV to 
%   the known Radar RPV will provide a spoofing metric.  This
%   metric can be thresholded based on future testing to provide certainty.
%   Case 1: Both vehicles spoofed either with jammer position or with false
%           generated position (Possibly the position of one node)
%   Case 2:One vehicle is spoofed, the other is unaffected.
%          This would result in a variance in the Radar RPV compared to the
%          GPS RPV but no jump in position or convergence to a single
%          position.
%   
%   *** ***  Localization Method  *** ***
%   The localization method

%% ######################### Run Profile Paramters ##################### %%
spoofing.init = 1;                   % Initiate spoofing?  1 = yes, 0 = no
spoofing.traj = 'linear';        % Select spoofing trajectory type
spoofing.mag = -5;                    % Set spoofing magnitude
spoofing.startTime = 175;                 % Set start of spoofing (s)
spoofing.endTime = 200;
scale = 10;
%R_Error_Thresh = 2;             % Threshold for Spoofing Detection (m)


% Data information -> Select Data to Use.
FilePath = 'C:/Users/Nath/Nathaniel/Auburn/Research/Jamming and Spoofing Detection/Code/Data/';
MatfileNameV1 = 'MagnoliaLoopV1.mat';
MatfileNameV2 = 'MagnoliaLoopV2.mat';
Frequency = 1;                    % Hz (Frequency of GPS Data Collection)

%% Load Data from Mat files
FullPathNameV1 = strcat(FilePath,MatfileNameV1);
FullPathNameV2 = strcat(FilePath,MatfileNameV2);
if exist(FullPathNameV1,'file') == 2 
    load(FullPathNameV1); 
    
    % Retreive Position
    POS_1_LLA = [out.Bestpos.Lat_pos', out.Bestpos.Long_pos',...
                                                  out.Bestpos.Height_pos'];
    % Retreive time vector
    GpsTimeV1 = out.Bestpos.gpstime; 
    
    % Retreive Standard Deviations
    StdDevLat1 = out.Bestpos.standardevLat';
    StdDevLong1 = out.Bestpos.standardevLong';
    StdDevHeight1 = out.Bestpos.standardevHeight';
   
else
    display('File does not exist or is not on path.'); return
end

if exist(FullPathNameV2,'file') == 2 
    load(FullPathNameV2); 
    POS_2_LLA = [out.Bestpos.Lat_pos', out.Bestpos.Long_pos',...
                                                  out.Bestpos.Height_pos'];
    % Retreive time vector
    GpsTimeV2 = out.Bestpos.gpstime;
    
    % Retreive Standard Deviations
    StdDevLat2 = out.Bestpos.standardevLat';
    StdDevLong2 = out.Bestpos.standardevLong';
    StdDevHeight2 = out.Bestpos.standardevHeight';
    
else
    display('File does not exist or is not on path.'); return
end

%% Align Data in Time

if GpsTimeV1(1,1) > GpsTimeV2(1,1)
    % Use Vehicle one as start of vectors
    StartTime = GpsTimeV1(1,1);
    Comp = GpsTimeV2 == StartTime;
    StartIndex = find(Comp);
    
    POS_2_LLA = POS_2_LLA(StartIndex:end,:);
    
    GpsTimeV2 = GpsTimeV2(StartIndex:end);
    
    StdDevLat2 = StdDevLat2(StartIndex:end,:);
    StdDevLong2 = StdDevLong2(StartIndex:end,:);
    StdDevHeight2 = StdDevHeight2(StartIndex:end,:);
    StdDevPos2 = sqrt(StdDevLat2.^2 + StdDevLong2.^2 + StdDevHeight2.^2);
    StdDevHorizontal2 = sqrt(StdDevLat2.^2 + StdDevLong2.^2);
else 
    % Use Vehicle two as start of vectors
    StartTime = GpsTimeV2(1,1);
    Comp = GpsTimeV1 == StartTime;
    StartIndex = find(Comp);
    
    POS_1_LLA = POS_1_LLA(StartIndex:end,:);
    
    GpsTimeV1 = GpsTimeV1(StartIndex:end);
    
    StdDevLat1 = StdDevLat1(StartIndex:end,:);
    StdDevLong1 = StdDevLong1(StartIndex:end,:);
    StdDevHeight1 = StdDevHeight1(StartIndex:end,:);
    StdDevPos1 = sqrt(StdDevLat1.^2 + StdDevLong1.^2 + StdDevHeight1.^2);
    StdDevHorizontal1 = sqrt(StdDevLat1.^2 + StdDevLong1.^2);
    
end

if length(GpsTimeV1) > length(GpsTimeV2)
    
    GpsTimeV1 = GpsTimeV1(1:length(GpsTimeV2));
    
    POS_1_LLA = POS_1_LLA(1:length(GpsTimeV2),:);
    
    StdDevLat1 = StdDevLat1(1:length(GpsTimeV2),:);
    StdDevLong1 = StdDevLong1(1:length(GpsTimeV2),:);
    StdDevHeight1 = StdDevHeight1(1:length(GpsTimeV2),:);
    StdDevPos1 = sqrt(StdDevLat1.^2 + StdDevLong1.^2 + StdDevHeight1.^2);
    StdDevHorizontal1 = sqrt(StdDevLat1.^2 + StdDevLong1.^2);
else
    GpsTimeV2 = GpsTimeV2(1:length(GpsTimeV1));
    
    POS_2_LLA = POS_2_LLA(1:length(GpsTimeV1),:);
    
    StdDevLat2 = StdDevLat2(1:length(GpsTimeV1),:);
    StdDevLong2 = StdDevLong2(1:length(GpsTimeV1),:);
    StdDevHeight2 = StdDevHeight2(1:length(GpsTimeV1),:);
    StdDevPos2 = sqrt(StdDevLat2.^2 + StdDevLong2.^2 + StdDevHeight2.^2);
    StdDevHorizontal2 = sqrt(StdDevLat2.^2 + StdDevLong2.^2);
    
end


%% Perform Conversions and Evaluate Positions

% Set ORIGIN at the initial Vehicle 1 position
ORIGIN = POS_1_LLA(1,:);

%  POS_1_LLA should be matrices where each row contains [lat, long, alt]
%  for a particular time step. Time progresses top to bottom of the matrix.

POS_1_ENU = lla2enu(POS_1_LLA, ORIGIN);
POS_2_ENU = lla2enu(POS_2_LLA, ORIGIN);

% Calculated GPS RPV and Range (Magnitude of RPV)
True_RPV_GPS = POS_1_ENU - POS_2_ENU;    % Matrix of differences
True_Range_GPS = sqrt(sum(True_RPV_GPS.^2,2));    % Vertical column of ranges at each step

%% Generate Simulated Radar Ranges
RadarNoise = .25;
Noise = RadarNoise * randn(length(True_Range_GPS),1);
for ii = 1:length(POS_2_ENU)
    rP1(ii,1:2) = POS_1_ENU(ii,1:2) + [StdDevLat1(ii,1) * randn(1) / scale, StdDevLong1(ii,1) * randn(1) / scale];
    rP2(ii,1:2) = POS_2_ENU(ii,1:2) + [StdDevLat2(ii,1) * randn(1) / scale, StdDevLong2(ii,1) * randn(1) / scale];
    Simulated_RPV_GPS(ii,1:2) = rP1(ii,1:2) - rP2(ii,1:2);
    Simulated_Range_GPS(ii,1) = sqrt(sum(Simulated_RPV_GPS(ii,1:2).^2));
end
Range_radar = Simulated_Range_GPS + Noise;

%% Introduce Spoofing
if spoofing.init ==1 
    POS_2_Spoofed_ENU = spoofer(spoofing.mag,spoofing.traj,...
                                        spoofing.startTime,spoofing.endTime, Frequency,POS_2_ENU);
    Spoofed_RPV_GPS = POS_1_ENU - POS_2_Spoofed_ENU;
    Spoofed_Range_GPS = sqrt(sum(Spoofed_RPV_GPS.^2,2));
else 
    Spoofed_Range_GPS = True_Range_GPS;
end



% Check difference vector for Spoofing indication
True_Range_diff = abs(True_Range_GPS - Range_radar);
Measured_Range_diff = abs(Spoofed_Range_GPS - Range_radar);
%R_threshold = R_Error_Thresh * ones(length(Measured_Range_diff),1);

R_threshold = 3*sqrt((StdDevHorizontal2 / scale).^2 + (StdDevHorizontal1 / scale) .^ 2 +...
    (RadarNoise * ones(length(True_Range_GPS),1)) .^ 2); 


Comparison = bsxfun(@gt,Measured_Range_diff,R_threshold);
Index = find(Comparison);      % Extracts index of spoofed iterations
for jj = 1:(length(Comparison)-3);
    comp(jj) = Comparison(jj) + Comparison(jj+1) + Comparison(jj + 2) + Comparison(jj + 3);
end
    Comparison2 = bsxfun(@gt,comp,3.5*ones(1,length(comp)));
    Index2 = find(Comparison2)
if abs(Index) > 0
    spoof_RPVdiff = 1;
    display('## Your GPS Position has been Spoofed ##')
    spoofedMeas = length(Index);  % Number of Suspicious Measurements
    display(['## Number of Suspicious Measurements: '  num2str(spoofedMeas);])
else 
    spoof_RPVdiff = 0;
end 

%% Plotting

% Actual Positions of Both Vehicles
unspoofedPos = figure; hold on
plot((POS_1_ENU(:,1)),POS_1_ENU(:,2),'x')
plot((POS_2_ENU(:,1)),POS_2_ENU(:,2),'ro')
legend('Vehicle 1', 'Vehicle 2','FontSize', 12)
xlabel('East (m)','FontSize', 12), ylabel('North (m)','FontSize', 12) 
title('Unspoofed Positions','FontSize', 12)

% Spoofed Positions of Both Vehicles
if spoofing.init ==1
    spoofedPos = figure; hold on
    plot((POS_1_ENU(:,1)),POS_1_ENU(:,2),'x','linewidth',1.5)
    plot(POS_2_Spoofed_ENU(:,1),POS_2_Spoofed_ENU(:,2),'ro','linewidth',1.5)
    legend('Vehicle 1', 'Vehicle 2','FontSize', 12)
    xlabel('East (m)','FontSize', 12), ylabel('North (m)','FontSize', 12) 
    title('Spoofed Positions','FontSize', 12)
    
end

% Relative Position Vectors
RPVs = figure;
TimeScale = 1/Frequency:1/Frequency:(length(Range_radar))/Frequency;
subplot(2,1,1); hold on
bar(TimeScale,Range_radar,1,'b')
bar(TimeScale,Spoofed_Range_GPS,1,'m')
axis([-inf length(TimeScale) -inf inf])
alpha(.5)
legend('Radar Reported Range','GPS Reported Range','FontSize', 12)
xlabel('Time (s)','FontSize', 12), ylabel('Range (m)','FontSize', 12)
title('Reported Ranges between Vehicles','FontSize', 12)
hold off

% RPV Difference (Indication of Spoofing)
subplot(2,1,2), hold on
bar(TimeScale,Measured_Range_diff)
plot(R_threshold,'g','linewidth', 1.5)
xlabel('Time (s)','FontSize', 12), ylabel('Range Difference (m)','FontSize', 12)
title('Difference Between GPS and Radar Measurements','FontSize', 12)
legend('Difference Between RPVs','Spoofing Threshold')
axis([-inf length(TimeScale) -inf inf])
hold off

% Standard Deviations of GPS Positions
stdDev = figure; 
axisMax = max([ceil(max(StdDevPos1)), ceil(max(StdDevPos2))])/scale;
%subplot(2,1,1), 
hold on
RPVStdDev = sqrt((StdDevHorizontal1/scale).^2 + (StdDevHorizontal2/scale) .^ 2);
plot(TimeScale,StdDevHorizontal1/scale,'linewidth',2)
plot(TimeScale,StdDevHorizontal2/scale,'r','linewidth',2)
plot(TimeScale,RPVStdDev,'g','linewidth',2)
ylabel('1 Sigma Deviation (m)','FontSize', 12), xlabel('Time (s)','FontSize', 12)
title('Estimated Horizontal Standard Deviation','FontSize', 12)
legend('Vehicle 1','Vehicle 2','Combined Deviation','FontSize', 12)
axis([-inf inf -inf axisMax])
hold off
figure 
%subplot(2,1,2), 
hold on
plot(TimeScale,StdDevPos1,'linewidth',2)
plot(TimeScale,StdDevPos2,'r','linewidth',2)
ylabel('1 Sigma Deviation','FontSize', 12), xlabel('Time (s)','FontSize', 12)
title('Estimated Standard Deviation of Position','FontSize', 12)
legend('Vehicle 1','Vehicle 2','FontSize', 12)
axis([-inf inf -inf axisMax])
















