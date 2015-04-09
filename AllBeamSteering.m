%Joshua Starling
%The purpose of this program is to read data taken from the USRP
%then use Frost's algorithm to  "block out" sections of the sky
%and not be able to acquisition sattelites from those sections

clc
clear
clf
close all
%% Initial Set Up
tic
NUM_MS=1;       % Number of MS used in the ingration period
NUM_DATA_SETS=NUM_MS+1;
skip_seconds=1*60;    % skip initial # seconds (initial data sometimes not good).
freqL2 = 1227.6e6;
freqL1 = 154*10.23e6;
sample_frequency = 5e6;
intermediate_frequency = 0;
integration_period = NUM_MS*(1e-3);
gpsPi = 3.14159265359; 
bytes_per_sample = 2;



file1 ='Feb_6_5M_int16_CSAC_USRP10_GPS_TestTwo.dat';
file2 = 'Feb_6_5M_int16_CSAC_USRP11_Powered_Test2.dat';
file3 ='Feb_6_5M_int16_CSAC_USRP10_GPD_Attn.dat';
file4 ='Feb_6_5M_int16_CSAC_USRP11_Powered.dat';

%% Reading Data file
%Functions reads the USRP data and converts it into a 2Xn array, the first
%being the real component and the second row being the imaginary
[signal(1,:), signal1_comp] = read_USRP_data(file3,sample_frequency,integration_period,skip_seconds,bytes_per_sample);
[signal(2,:), signal2_comp] = read_USRP_data(file4,sample_frequency,integration_period,skip_seconds,bytes_per_sample);

%% Acquisition of Unmodified signal to show sattelites in view
SV_array = [22];%list of satellites that will have a correlation graph generated
%SV_array = 1:32;
[I1,Q1] = Acquisition(signal(1,:),SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
[I2,Q2] = Acquisition(signal(2,:),SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
% output = Mod_Acquisition(I1,Q1,SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
% output = Mod_Acquisition(I2,Q2,SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
CW_signal1 = I1 + j*Q1;
CW_signal2 = I2 + j*Q2;
%% Combining Two singals into one
%USRP's must be correctly synchronized for this to work

%Code from STAP_for_sponsor. Not sure if this is right.
%initialize files

%Done in read_USRP_Data
% fin1 = fopen(fileid(1,:));
% fin2 = fopen(fileid(2,:));

% fseek(fin1, 300*fs,'bof');
% fseek(fin2, (300-2.25)*fs,'bof');
%fseek(fin3, (300-4.5)*fs, 'bof'); %why are all of these numbers different
%fseek(fin4, (300-6.75)*fs,'bof'); %might deal with synchronizing data?
% 
% delta = [1 0]'; %weighting vector??
% sig1 = signal1_comp(1,:)+ 1j*signal1_comp(2,:);
% sig2 = signal2_comp(1,:)+ 1j*signal2_comp(2,:);
% x = [signal(1,:); signal(2,:)];
% Rxx = (x*ctranspose(x));
% weight_opt = (1/(delta'*inv(Rxx)*delta))*inv(Rxx)*delta; 
% modified_signal = ctranspose(weight_opt)*x;
% [I2,Q2] = Acquisition(modified_signal,SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);

%% Beam Steering using AppleBaum
%NOTE: Matlab's beam forming/steering using metic system

% c = 299792458; %m/s
% antenna_distance = .18; %has to be meters
% theta =45*pi/180; %Desired angle of arival
% number_of_antennas = 2;
% lamda = c/freqL1; %wavelength
% s = zeros(1,number_of_antennas);
% weight = zeros(1,number_of_antennas);
% for k = 1:number_of_antennas
%     s(k) = exp(-1i*2*pi*k*antenna_distance/lamda*sin(theta));
%     weight(k) = s(k);
%     mod_signal(k,:) = weight(k)*signal(k,:);
% end
% modified_signal = mod_signal(1,:)+mod_signal(2,:);

%% Beam steering Take Two
% syms b
% 
% d = .18; %antenna distance (meters)
% AOA = 0*pi/180;
% c = 299792458; %m/s
% lamda = c/freqL1; %wavelength
% theta1 = (2*pi*d*cos(AOA))/lamda;
% a_s = [1 exp(-j*theta1)]';
% w = [1; b];
% k = 2; %number of antennas
% b = solve(w'*a_s==k,b);
% b = eval(b);
% weight = [1; b];
% modified_CW = weight(1)*CW_signal(1,:)+weight(2)*CW_signal(2,:);

%% Beam Steering Take Three
d = .18; %antenna distance (meters)
theta = 33*pi/180;
elevation = 45*pi/180;
azimuth = 2*pi/180;
c = 299792458; %m/s
lamda = c/freqL1; %wavelength
% phi = cov(CW_signal(1,:),CW_signal(2,:));
%     y = [CW_signal(1,:); CW_signal(2,:)];
%     Ryy = (y*ctranspose(y));
%     %Sarah's Method, similar values different magnitude
%     x = [signal(1,:); signal(2,:)];
%     Rxx = (x*ctranspose(x));
% mu = 1;
p = [d*cos(theta) d*sin(theta) 0];
r = [cos(azimuth)*cos(elevation) cos(azimuth)*sin(elevation) sin(azimuth)];
T = [1; exp(-j*(2*pi*dot(p,r)/lamda))];
modified_CW_signal = T(1)*CW_signal1 + T(2)*CW_signal2;

%% Acquisition of modified signal
output = Mod_Acquisition(modified_CW_signal,SV_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS);
toc
%% Power of signal
% for i = 1:length(signal(1,:))
%     p_signal(1,i) = abs(signal(1,i));
%     p_mod_signal(i) = abs(modified_signal(i));
% end
% figure
% plot(p_signal,'k')
% hold on
% plot(p_mod_signal,'r')