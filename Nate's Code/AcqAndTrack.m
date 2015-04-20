clear all; close all; clc;

%% Filtered and upconverted data
 % filename(1,:) = 'April_13_Sim_sat125_int99_Ant1.dat';
 % filename(2,:) = 'April_13_Sim_sat125_int99_Ant2.dat';
 % filename(3,:) = 'April_13_Sim_sat125_int99_Ant3.dat';
 filename(1,:) = 'April_13_Sim_sat125_int99_Skip354_Ant1.dat';
 filename(2,:) = 'April_13_Sim_sat125_int99_Skip354_Ant2.dat';
 filename(3,:) = 'April_13_Sim_sat125_int99_Skip354_Ant3.dat';
%filename = 'April_9_Sim_Test.dat';

%% Adjustable parameters
init.sattId = [1 9 17 14 27 11 22 32 3 10 25]; %1 9 17 14 27 11 22 32 3 10 25
init.dopp_bin_size = 100;       % Hertz
init.run = 'nordnav';               % nordnav or wis
init.spoof = 0;                 % Add spoofed Peak?         1 = yes, 0 = no
init.plotAcqPlanes = 0;         % Plot correlation planes?  1 = yes, 0 = no
init.plotAcqResults = 1;        % Plot Acquisition Results? 1 = yes, 0 = no
init.trackAlgorithm = 1;         % Enter tracking loops?     1 = yes, 0 = no 
init.addNoise = 0;              % Add Noise to signal?      1 = yes, 0 = no 
init.cleanSig = 0;
init.BeamForming = 1;
init.run_type = 3;

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
	return; 
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
