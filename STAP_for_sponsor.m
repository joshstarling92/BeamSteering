%   Sarah Preston
%   STAP coding

clc
clear 
clf
close all

%%  Part 1: Synchronize the times between antennas by finding the Doppler...
...shift and code shift of the PRN from one SV for all four antennas. Then...
...perform a PLL and DLL for each antenna to get the data bits. The timing data bits for each antenna should be the same, so doing this should...
...allow for synchronization of time across all 4 antennas. Timing is critical for STAP operation. First few minutes of data files should be jam-free
%Initialize intermediate freq and files to be used
tic
f_if = 0;
samp_freq = 5e6;
integrate_period = 1e-3;    %was having issues because forgot to change the Doppler step size
f_L1 = 1575.42e6;
init_delay = 40;                %initial delay of about 5 min because testing didn't always start with "clear" sky         
chip_freq = 1.023e6;

% addpath('
% fileid(1,:) = '/Users/sarahpreston/Documents/Work for Sponsor/Data Files/untitled folder/4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_2_0115_et1.dat';
% fileid(2,:) = '/Users/sarahpreston/Documents/Work for Sponsor/Data Files/untitled folder/4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_3_0115_et1.dat';
% fileid(3,:) = '/Users/sarahpreston/Documents/Work for Sponsor/Data Files/untitled folder/4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_4_0115_et1.dat';
% fileid(4,:) = '/Users/sarahpreston/Documents/Work for Sponsor/Data Files/untitled folder/4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_5_0115_et1.dat';
% 
% file = '/Users/sarahpreston/Desktop/Sept_23_45db_5M_int16_Octoclock_plus_CSAC_DSBRXboard.dat' ; 
%file = '/Users/sarahpreston/Desktop/USRP_Roof_Data/Sept_30_5M_int16_CSAC.bin';
file ='Jan_7_5M_int16_no_CSAC.dat';
% file = '/Users/sarahpreston/Desktop/May_5_2014_22_decim_rate_int16_wireformat_int8_PPS_in_Ext_ref_clock';
% file = '../MATLAB/GPS_Data_NordNav1e.sim';

 delay_sec = init_delay;  %- 1*(n-1);
%  samp_freq = 16.3676e6;
signal1(1,:) = read_USRP_data(file,samp_freq,integrate_period,delay_sec);

%intermed_freq = 4.1304e6;
intermed_freq = 0;

[SV_inview,freq_shift(1,:),code_phase(1,:)] = parallel_search(file,samp_freq,integrate_period,delay_sec,signal1(1,:),intermed_freq);

% parfor n = 1:4

% parfor n = 1:4
%     signal1(n,:) = read_USRP_data(fileid(n,:),samp_freq,integrate_period,delay_sec);
% end
%     [SV_inview,freq_shift(1,:),code_phase(1,:)] = parallel_search(fileid(1,:),samp_freq,integrate_period,delay_sec,signal1(1,:), intermed_freq); 
% % 
%     for n=2:4
%        [freq_shift(end+1,:), code_phase(end+1,:)] = parallel_search_second(fileid(n,:),samp_freq, integrate_period, signal1(n,:),SV_inview);
%     end
%     
% %     %     fclose(fileid(n,:))
% % % end
% fprintf(['Acquisition completed \n'])
% fclose('all');
% fid = fopen('USRP_data.txt','w');
% dlmwrite('USRP_data.txt', [SV_inview; freq_shift; code_phase],'delimiter','\t');
% clear freq_shift code_phase SV_inview signal1

% % read the file that was generated during acquisition. The file will contain
% % the top six SVs with their corresponding Doppler frequencies and code
% % shifts
% fid = fopen('USRP_data.txt','r');
% 
% 
% x = dlmread('USRP_data.txt');
% SV_inview = x(1,:);
% freq = x(2:5,:);
% code = x(6:end,:);
% fclose(fid);
% init_delay=300;         %initial delay to skip the first 30-40 sec of USRP data which is generally not very good
fclose('all');
% %get first file and first SV data
% init_tracking_vars
% initialize_vars_STAP

% fprintf(['Starting tracking \n'])

% fs=round(100e6/22);
fs = 5e6;
i = 2;
% while ~feof(fid)

fclose('all');
%  for prn = 1:6
for prn = 1:6;
bw_num = 15;    % bw=8, k= 0.25 gets lock for PRN2      
% init_tracking_vars
% initialize_vars_STAP  

for i=2:3e5
    
%     if i<5000      % i<10000 works pretty well
%         
%         fll_dll_tracking_loop
%     else
        
        %pll_dll_tracking_loop
%         if parity ~=1 
%          
%         else
%             frame_decode
%         end
%         if data_tracking~=1
%             fprintf([num2str(i),'\n'])
%             clc
%         
%         end

%     end
%         get first file and first SV data

        
    
%  progressbar(i/(3*3e5));
  
 end  
 

% i = i+1;


% save(bw_work)
fclose(fid);
toc
figure; subplot(221)
plot(Ip,'b'); hold on; plot(Qp,'r')
title(sprintf('IP for PRN %i',prn))
subplot(222)
plot(f)
title(sprintf('Freq for PRN %i',prn))

subplot(223)
plot(th*180/pi)
title(sprintf(' \theta for PRN %i in degrees',prn)) 

subplot(224)
plot(cf-1.023e6)
title(sprintf('Normalized Chipping Freq'))

figure; plot([sum(Ie) sum(Ip) sum(Il)])
title(sprintf('EPL for In-Phase arm for PRN %i',prn))




end

% prn_5_data.f = f;
% prn_5_data.pll_params.bw = bw_carr;
% prn_5_data.pll_params.k = k_carr;
% prn_5_data.code_params.bw = bw_code;
% prn_5_data.code_params.k = k_code;
% prn_5_data.subframe_data = subframe_data;
% prn_5_data.data_bits = db_from_locked_tracking;
% prn_5_data.almanac = almanac;
% prn_5_data.clock = clock;
% prn_5_data.eph = ephemeris;
% prn_5_data.data = data;
% prn_5_data.I = Ip;
% prn_5_data.Q = Qp;
% 
% save('Decoded_PRN_data_5','prn_5_data');


% f_if = 0;
% f_dopp = freq(4,1);
% f(1)=f_if+f_dopp;       %initial f_dopp
% SV = SV_inview(4,1);           %SV 19, which appears for all four USRPs
% carrFreq(1) = f_if + f_dopp;t
% code_shift(1) = code(4,1);
% 
%  
% fs=samp_freq;           %sampling frequency
% samp_per_chip = fs/1.023e6;
% cf(1)=1.023e6;          %chipping frequency
% % phaseRem(1) = 0;
% % 
% fid = fopen(fileid(4,:));
% skip = fseek(fid,init_delay*samp_freq,'bof');
% theta(1) = 0; theta_hat(1) = 0;
% w(1) = 0; w_hat(1) = 0;
% alpha(1) = 0; alpha_hat(1) = 0;
% 
% lambda = 299792458/1575.42e6;
% % d = 
% % 
% % cf = 1.023e6;
% 
%% Beamforming Algorithm
% 
% %initialize files
% fin1 = fopen(fileid(1,:));
% fin2 = fopen(fileid(2,:));
% fin3 = fopen(fileid(3,:));
% fin4 = fopen(fileid(4,:));
% 
% fseek(fin1, 300*fs,'bof');
% fseek(fin2, (300-2.25)*fs,'bof');
% fseek(fin3, (300-4.5)*fs, 'bof');
% fseek(fin4, (300-6.75)*fs,'bof'); 
% % 
% % a = [1; exp(1j*pi/lambda*sin(theta(k)));...
% %     exp(1j*pi/lambda*sin(theta(k))*d);
% %     exp(1j*sqrt(2)*pi/lambda*sin(theta(k)))]; 
% 
% delta = [1 0 0 0]';
% t1 = fread(fin1, [2 fs*1e-3],'short');
% t2 = fread(fin2, [2 fs*1e-3],'short');
% t3 = fread(fin3, [2 fs*1e-3],'short');
% t4 = fread(fin4, [2 fs*1e-3],'short');
% 
% sig1 = t1(1,:)+ 1j*t1(2,:);
% sig2 = t2(1,:)+ 1j*t2(2,:);
% sig3 = t3(1,:)+ 1j*t3(2,:);
% sig4 = t4(1,:)+ 1j*t4(2,:);
% x = [sig1; sig2; sig3; sig4];
% 
% Rxx = (x*ctranspose(x));
% 
% weight_opt = (1/(delta'*inv(Rxx)*delta))*inv(Rxx)*delta; 
% 
% y_sig = ctranspose(weight_opt)*x;
% 
% [SV_inview, freqshift, codephase] = parallel_search(fin2, fs, integrate_period,0,y_sig);

% fseek(fid,300*fs,'bof');

% theta_hat = squeeze(x_hat_carr(1,1,:));                       %squeeze(squeeze(z_hat_v(1,:,k)));
% f_dopp_hat = squeeze(x_hat_carr(2,1,:));
% theta = atan(Qp./Ip); 
% figure; subplot(321); plot(Ip); title('IP')
% subplot(322); plot(Qp); title('QP')
% subplot(323); plot(f_dopp_hat); title('Doppler Freq from KF')
% subplot(324); plot(theta_hat); title('\theta from KF')
% subplot(325); plot(atan(Qp./Ip)); title('Residual \theta from QP/IP')
% % subplot(326); plot(Ahat);   title('Ahat') 
% subplot(326); plot(squeeze(y(1,:,:))); title('\theta')


% save('Tracking output',Ip);

% toc
