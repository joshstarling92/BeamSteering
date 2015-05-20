%Writen by: Joshua Starling
%Email: jds0033@tigermail.auburn.edu
%The purpose of this program is to read data taken from synchronized USRP
%using acquisition.m and apply a weighting to each singal. The weighting
%can be determined from beam steering, null steering, an adaptive process
%that throws out signals of higher power, or predefined weights. The
%signals are then summed and mod_acquisition.m shows what satellites were
%acquired. 
clc
clear
clf
close all
tic
%% Desired Run Variables
Settings.NUM_MS=2;       % Number of MS used in the ingration period
Settings.skip_seconds=354;%4*60;    % skip initial # seconds (initial data sometimes not good).
Settings.sample_frequency = 5e6;
Settings.SV_array = [1 9 17 14 27 11 22 32 3];   %list of satellites that will have a correlation graph generated

%Run variables
Settings.NumberOfAntennas = 3;
Settings.run_type = 8; %0: W=1, 1: User defiend weights, 2: PM, 3: LMS, 4: Applebaum 5: NS, 6: BS
Settings.Add_noise = 0; % 1 adds noise to signal 0 does not

%Plots
Settings.unmod_acq_graph = 0; %0 to not plot acquisition graph of SVs, 1 to plot acquisition graph of SVs
Settings.mod_acq_graph = 1; %0 to not plot mod_acquisition graph of SVs, 1 to plot mod_acquisition graph of SVs
Settings.gain_graph = 0; %0 to not plot gain graph, 1 to plot graph

%% Settings/contants
Settings.NUM_DATA_SETS=Settings.NUM_MS+1;
Settings.freqL2 = 1227.6e6;
Settings.freqL1 = 154*10.23e6;
Settings.intermediate_frequency = 0;
Settings.integration_period = Settings.NUM_MS*(1e-3);
Settings.gpsPi = 3.14159265359; 
Settings.bytes_per_sample = 2;
Settings.c = 299792458; %m/s
Settings.lamda = Settings.c/Settings.freqL1; %wavelength (m)
Settings.unmod_done = 0;

%% Initialize Signals Structure
Signals.CleanSignal = zeros(Settings.NumberOfAntennas,20000);
Signals.CleanSignalComp = zeros(2,20000,Settings.NumberOfAntennas);
Signals.I = zeros(82,10000,Settings.NumberOfAntennas);
Signals.Q = zeros(82,10000,Settings.NumberOfAntennas);
Signals.CW_signal = zeros(82,10000,Settings.NumberOfAntennas);
Signals.modified_CW_signal = zeros(82,10000);

% Open figures to reduce processing time
for figCnt = 1:(length(Settings.SV_array)*Settings.unmod_acq_graph+length(Settings.SV_array)*Settings.mod_acq_graph+Settings.gain_graph)
    figure(figCnt)
end
  
%% Files
file1 ='/Volumes/StolenDrive/May_9_5M_int16_Sim03042015_Ant1_sat125_int99_Multi_Interference_GPSDOLinked.dat';
file2 ='/Volumes/StolenDrive/May_9_5M_int16_Sim03042015_Ant2_sat125_int99_Multi_Interference_GPSDOLinked.dat';
file3 ='/Volumes/StolenDrive/May_9_5M_int16_Sim03042015_Ant3_sat125_int99_Multi_Interference_GPSDOLinked.dat';
file4 ='/Volumes/StolenDrive/May_9_5M_int16_Sim03042015_Ant4_sat125_int99_Multi_Interference_GPSDOLinked.dat';
file5 ='/Volumes/StolenDrive/April_13_5M_Sim03042015_int16_Simulator_Ant1_Interference_sat125_int99_GSPDOLink.dat';
file6 ='/Volumes/StolenDrive/April_13_5M_Sim03042015_int16_Simulator_Ant2_Interference_sat125_int99_GSPDOLink.dat';
file7 ='/Volumes/StolenDrive/April_13_5M_Sim03042015_int16_Simulator_Ant3_Interference_sat125_int99_GSPDOLink.dat';

%% Reading Data file
[Signals.CleanSignal(1,:), Signals.CleanSignalComp(:,:,1)] = read_USRP_data(file5,Settings);
[Signals.CleanSignal(2,:), Signals.CleanSignalComp(:,:,2)] = read_USRP_data(file6,Settings);
[Signals.CleanSignal(3,:), Signals.CleanSignalComp(:,:,3)] = read_USRP_data(file7,Settings);
%[Signals.CleanSignal(4,:), Signals.CleanSignalComp(:,:,4)] = read_USRP_data(file4,Settings);

%% Adding noise
if Settings.Add_noise == 1
    for counter = 1:Settings.NumberOfAntennas
        sigma = 10000;
        noise = sigma * randn(1,length(Signals.CleanSignal(1,:)));
        Signals.CleanSignal(counter,:) = Signals.CleanSignal(counter,:) + noise;
    end
end

%% Acquisition of Unmodified signal to show satellites in view
for CurrentAntenna = 1:Settings.NumberOfAntennas
    [Signals] = Acquisition(Settings,Signals,CurrentAntenna);
    Signals.CW_signal(:,:,CurrentAntenna) = Signals.I(:,:,CurrentAntenna) + 1i*Signals.Q(:,:,CurrentAntenna);
end

for mu = 5*10^-12;
    Settings.mu = mu;
    %% Beam Forming
    [Signals,W] = WeightCalculation(Settings,Signals);
    % if max(real(Signals.ErrorResults))<2000
    %     break
    % end
end

fprintf('Weights Calculated\n') 
t = 1:10000;
figure()
%plot(t,real(W(1,:)),t,real(W(2,:)),t,real(W(3,:)))
%fprintf('%6.2f\n',W(end,:))
%% Acquisition of modified signal
Mod_Acquisition(Settings,Signals);
toc

%% Graphing of gain pattern
if Settings.gain_graph == 1
    d = 0.09514; %(m)
    %d = 0.19029; %(m)
    %Plot of The Output Radiation Pattern
    t = 0:0.05:2*pi;
    tplot = 2*pi:-0.05:0;
    M = 0;
    for I = 1:Settings.NumberOfAntennas
       H = exp (i*(2*pi*d*(I-1)*sin(t))/lamda)  ;
       M = M + (H*W(I)) ;
    end
    M = abs(M) ;
    figure()
    polar(tplot,M,'-r') , title ('The Generalized Null Steering Beam Former Output Radiation Pattern') , grid on ;
    view([90 -90])
end