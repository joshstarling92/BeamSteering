clear all; close all; clc;

% Data file and path to file
%filename = 'C:\Users\Nath\Nathaniel\Auburn\Research\Jamming and Spoofing Detection\Code\Data\NathanielNordnav.sim';
%filename = 'C:\Users\Nath\Nathaniel\Auburn\Class Work\GPS Class\GPS_Data_NordNav1e.sim';
%filename = 'C:\Users\Nath\Nathaniel\Auburn\Research\Jamming and Spoofing Detection\Code\Sarah\4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_2_et1.dat';
%filename = 'C:\Users\Nath\Nathaniel\Auburn\Research\Jamming and Spoofing Detection\Code\Sarah\4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_2_et2.dat';
%filename = 'C:\Users\Nath\Nathaniel\Auburn\Research\Jamming and Spoofing Detection\Code\Data\USRP Convertor Batch File\Converted Data\Feb_6_5M_int16_CSAC_USRP11_Powered.dat';

%filename = 'C:\Users\Nath\Nathaniel\Auburn\Research\Jamming and Spoofing Detection\Code\Data\WZMR Data\Dec_12_2013_05_00_00_et1.sim';
%filename = 'C:\Users\Nath\Nathaniel\Auburn\Research\Jamming and Spoofing Detection\Code\Data\WZMR Data\Dec_12_2013_05_00_00_et3.sim';

%% Fil     tered and upconverted data
filename(1,:) = 'NotGoingToWork.dat';
%filename(2,:) = 'Upconvert_Feb_23_LS_Ant_2_R.dat';
% filename(3,:) = 'Upconvert_April_Ant3.dat';
% Strong: 3 16 23 27
% Medium: 8 32
% Weak  : 31


%% Adjustable parameters
init.sattId = [1:32];
init.dopp_bin_size = 100;       % Hertz
init.run = 'nordnav';               % nordnav or wis
init.spoof = 0;                 % Add spoofed Peak?         1 = yes, 0 = no
init.plotAcqPlanes = 0;         % Plot correlation planes?  1 = yes, 0 = no
init.plotAcqResults = 1;        % Plot Acquisition Results? 1 = yes, 0 = no
init.trackAlgorithm = 1;         % Enter tracking loops?     1 = yes, 0 = no 
init.addNoise = 0;              % Add Noise to signal?      1 = yes, 0 = no 
init.cleanSig = 0;

init.spoofCodePhaseChip = 22000;      % C/A Chips
init.spoofDoppShift = 1200;           % Doppler shift in Hz
init.spoofMagnitude = .17;

% Initialize the Settings
Settings = set_settings(init);
number_of_antennas = 1;
%% Perform Acquisition 
for ii = 60*1
    Settings.pos = floor(ii*Settings.Sf);
    [Corr_results,Incoming,acqChannel] = Acquisition(Settings, init, filename,number_of_antennas); 
end


%% ENTER TRACKING LOOPS

if init.trackAlgorithm == 0; return; end 

[trackingResults, channel] = tracking(filename, acqChannel, Settings, init); 

%% Reacquisition

clearvars -except trackingResults init filename ii Corr_results
% Data file


% Set Settings
Settings = set_settings(init);
Settings.pos = (0*Settings.Sf);

ReacqChannel = Reacquisition(Settings, init, filename(i,:), trackingResults);
