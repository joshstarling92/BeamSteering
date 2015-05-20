clear all; close all; clc;

%% Filtered and upconverted data
 % filename(1,:) = 'April_13_Sim_sat125_int99_Ant1.dat';
 % filename(2,:) = 'April_13_Sim_sat125_int99_Ant2.dat';
 % filename(3,:) = 'April_13_Sim_sat125_int99_Ant3.dat';
 filename(1,:) = '/Volumes/StolenDrive/NateCode/April_13_Sim_Ant1_Interference_S125_I99_GPSDO_Upconverted.dat';
 filename(2,:) = '/Volumes/StolenDrive/NateCode/April_13_Sim_Ant2_Interference_S125_I99_GPSDO_Upconverted.dat';
 filename(3,:) = '/Volumes/StolenDrive/NateCode/April_13_Sim_Ant3_Interference_S125_I99_GPSDO_Upconverted.dat';
 %filename(1,:) = '/Volumes/StolenDrive/NateCode/April_9_Sim_No_Int_GPSDO_Upconverted.dat';


%filename = 'April_9_Sim_Test.dat';

%% Adjustable parameters
init.sattId = [1 9 17 14 27 11 22 32]; %1 9 17 14 27 11 22 32 3 10 25
init.dopp_bin_size = 100;       % Hertz
init.run = 'nordnav';               % nordnav or wis
init.spoof = 0;                 % Add spoofed Peak?         1 = yes, 0 = no
init.plotAcqPlanes = 1;         % Plot correlation planes?  1 = yes, 0 = no
init.plotAcqResults = 1;        % Plot Acquisition Results? 1 = yes, 0 = no
init.trackAlgorithm = 0;         % Enter tracking loops?     1 = yes, 0 = no 
init.addNoise = 0;              % Add Noise to signal?      1 = yes, 0 = no 
init.cleanSig = 0;
init.BeamForming = 1;
init.run_type = 8; 				%0: W=1, 1: User defiend weights, 2: PM, 3: LMS, 4: Applebaum 5: NS, 6: BS

init.spoofCodePhaseChip = 22000;      % C/A Chips
init.spoofDoppShift = 1200;           % Doppler shift in Hz
init.spoofMagnitude = .17;


% Initialize the Settings
Settings = set_settings(init);
Settings.number_of_antennas = 3;
%% Perform Acquisition 
%for loop for skipping through files
for ii = 0
    Settings.pos = floor(ii*Settings.Sf);
    [Corr_results,Incoming,acqChannel,W] = Acquisition(Settings, init, filename); 
end


%% ENTER TRACKING LOOPS
if init.trackAlgorithm == 0; 
	return 
end 
[trackingResults, channel] = MultiFileTracking(filename, acqChannel, Settings, init,W); 
%% Reacquisition

% clearvars -except trackingResults init filename ii Corr_results
% % Data file


% % Set Settings
% filename = 'April_9_Sim_Test.dat';
% Settings = set_settings(init);
% Settings.pos = (354*Settings.Sf);

% ReacqChannel = Reacquisition(Settings, init, filename, trackingResults);
