%Joshua Starling
%The purpose of this program is to read data taken from the USRP
%then use the CRPA determanistic weighting method to  "block out" sections 
%of the sky and not be able to acquisition sattelites from those sections
%Functions used: read_USRP_Data,cacode_bevly,Acquisition,Mod_Acquisition
clc
clear
clf
close all
%% Set Up
tic
NUM_MS=1;       % Number of MS used in the ingration period
NUM_DATA_SETS=NUM_MS+1;
skip_seconds=6*60;    % skip initial # seconds (initial data sometimes not good).
freqL2 = 1227.6e6;
freqL1 = 154*10.23e6;
sample_frequency = 5e6;
intermediate_frequency = 0;
integration_period = NUM_MS*(1e-3);
gpsPi = 3.14159265359; 
bytes_per_sample = 2;
c = 299792458; %m/s
lamda = c/freqL1; %wavelength

%% User Defined Set UP
K = 3; %Number of Antennas
run_type = 3; %0 for null steering, 1 for beam steering, 2 for user chosen weights, 3 for adaptive PI
file1 ='Feb_6_5M_int16_CSAC_USRP10_GPS_TestTwo.dat';
file2 = 'Feb_6_5M_int16_CSAC_USRP11_Powered_Test2.dat';
file3 ='Feb_10_5M_int16_CSAC_USRP10_Powered_Long_LowGain.dat';
file4 ='Feb_10_5M_int16_CSAC_USRP11_Powered_Short.dat';
file5 ='Feb_23_5M_int16_CSAC_USRP10_Bias_Tee_CRPA_Apex.dat';
file6 ='Feb_23_5M_int16_CSAC_USRP11_Bias_Tee_Right_CRPA.dat';
file7 ='Feb_23_5M_int16_CSAC_USRP13_GPS_CRPA_Left.dat';
%% Reading Data file
%Functions reads the USRP data and converts it into a 2Xn array, the first
%being the real component and the second row being the imaginary
[signal(1,:), signal1_comp] = read_USRP_data(file5,sample_frequency,integration_period,skip_seconds,bytes_per_sample);
[signal(2,:), signal2_comp] = read_USRP_data(file6,sample_frequency,integration_period,skip_seconds,bytes_per_sample);
[signal(3,:), signal3_comp] = read_USRP_data(file7,sample_frequency,integration_period,skip_seconds,bytes_per_sample);

%% Acquisition of Unmodified signal to show sattelites in view
SV_array = [3 16 23];   %list of satellites that will have a correlation graph generated
[I1,Q1] = Acquisition(signal(1,:),SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
[I2,Q2] = Acquisition(signal(2,:),SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
[I3,Q3] = Acquisition(signal(3,:),SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
CW_signal1 = I1 + j*Q1;
CW_signal2 = I2 + j*Q2;
CW_signal2 = I3 + j*Q3;

%% Beam Forming
switch run_type
    case  0
        fprintf('Preforming Null Steering\n')
        %%Null Steering - can block K-1 number of direction
        elevation = [54*pi/180 54*pi/180]; %Elevation of signal(s) that are to be blocked (rads)
        azimuth = [316*pi/180 316*pi/180]; %Azimuth of signal(s) that are to be blocked (rads)
        p(1,:) = [0 0 0];
        p(2,:)= [lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from second antenna to reference antenna (m)
        p(3,:)= [-lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from third antenna to reference antenna (m)
        r(1,:) = [sin(azimuth(1))*cos(elevation(1)) sin(azimuth(1))*sin(elevation(1)) cos(azimuth(1))]; %bore sight vector to interference
        r(2,:) = [sin(azimuth(2))*cos(elevation(2)) sin(azimuth(2))*sin(elevation(2)) cos(azimuth(2))]; %bore sight vector to interference
        b = [1;-1/(K-1);-1/(K-1)];
        a = [exp(j*(2*pi*dot(p(1,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(2,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(3,:),r(2,:))/lamda))];
        W = a.*b; %weight to be applied to signal
        modified_CW_signal = W(1)*CW_signal1 + W(2)*CW_signal2;
    case 1
        %%Beam Steering
        fprintf('Preforming Beam Steering\n')
        elevation = [54*pi/180 43*pi/180]; %Elevation of signal(s) that are to be strengthened (rads)
        azimuth = [316*pi/180 236*pi/180]; %Azimuth of signal(s) that are to be strengthened (rads)
        p(1,:) = [0 0 0];
        p(2,:)= [lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from second antenna to reference antenna (m)
        p(3,:)= [-lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from third antenna to reference antenna (m)
        r(1,:) = [cos(azimuth(1))*sin(elevation(1)) sin(azimuth(1))*sin(elevation(1)) cos(elevation(1))]; %bore sight vector to satellite
        r(2,:) = [cos(azimuth(2))*sin(elevation(2)) sin(azimuth(2))*sin(elevation(2)) cos(elevation(2))]; %bore sight vector to satellite
        W = [exp(j*(2*pi*dot(p(1,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(2,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(3,:),r(2,:))/lamda))];
        modified_CW_signal = W(1)*CW_signal1 + W(2)*CW_signal2;

    case 2
        %%User Chosen Weights
        W = [1 0.5+j*1.7648 -1.1710+j*1.4118];
        modified_CW_signal = W(1)*CW_signal1 + W(2)*CW_signal2;
    case 3
        %Adaptive PI Beam Steering
        delta = [1;0;0];
        Rxx = (signal*ctranspose(signal));
        W1 = 1/(delta'*inv(Rxx)*delta)*inv(Rxx)*delta;
        
        elevation = [54*pi/180 86*pi/180 43*pi/180]; %Elevation of satellite(s) (rads)
        azimuth = [316*pi/180 357*pi/180 236*pi/180]; %Azimuth of satellite(s) (rads)
        r(1,:) = [cos(azimuth(1))*sin(elevation(1)) sin(azimuth(1))*sin(elevation(1)) cos(elevation(1))]; %bore sight vector to satellite
        r(2,:) = [cos(azimuth(2))*sin(elevation(2)) sin(azimuth(2))*sin(elevation(2)) cos(elevation(2))]; %bore sight vector to satellite
        r(3,:) = [cos(azimuth(3))*sin(elevation(3)) sin(azimuth(3))*sin(elevation(3)) cos(elevation(3))]; %bore sight vector to satellite

        p(1,:) = [0 0 0];
        p(2,:)= [lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from second antenna to reference antenna (m)
        p(3,:)= [-lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from third antenna to reference antenna (m)
        A(1,:) = [1 exp(-j*2*pi*lamda/2*cos(azimuth(1))*sin(elevation(1))) exp(-j*2*pi*lamda/2*cos(azimuth(1))*sin(elevation(1)))];
        A(2,:) = [1 exp(-j*2*pi*lamda/2*cos(azimuth(2))*sin(elevation(2))) exp(-j*2*pi*lamda/2*cos(azimuth(2))*sin(elevation(2)))];
        A(3,:) = [1 exp(-j*2*pi*lamda/2*cos(azimuth(3))*sin(elevation(3))) exp(-j*2*pi*lamda/2*cos(azimuth(3))*sin(elevation(3)))];
        OneJ = ones(3,1);
        W = inv(Rxx)*A*inv(A*inv(Rxx)*A)*OneJ;
        modified_CW_signal = W(1)*CW_signal1 + W(2)*CW_signal2;
        
end
%% Acquisition of modified signal
peak = Mod_Acquisition(modified_CW_signal,SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
toc

%% Graphing of gain pattern
d = 0.09514; %(m)
%d = 0.19029; %(m)
%Plot of The Output Radiation Pattern
t = 0:0.05:2*pi;
M = 0;
for I = 1:K
    H = exp (i*(2*pi*d*(I-1)*cos(t))/lamda)  ;
    M = M + (H*W(I)) ;
end
M = abs(M) ;
figure()
polar(t,M,'-r') , title ('The Generalized Null Steering Beam Former Output Radiation Pattern') , grid on ;