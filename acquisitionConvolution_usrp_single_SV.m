% This code performs aquisition.  Can specify an integration period.  
% The code will search a window of the integration period in 1 ms
% increments to find the maximum value (i.e. looking for a window where no
% data bit transitions occur)  Recomended maximum window would be 20 ms.

% ***** initial conditions *****
clc 
close all
clear

% SV to look fr
for SV_NUM = [1:32]
NUM_MS=1;       % Number of MS used in the ingration period
SKIP_SECONDS= 2.25*60;    % skip initial # seconds (initial data sometimes not good).
% open data file
%filename = 'USRP_Data/4_5_M_complex_16_8int_cable_with_clock_unknownsync_PPS_2_cleartest.dat';
%filename = 'USRP_Data/6_27_13_test.dat';
filename ='April_13_5M_Sim03042015_int16_Simulator_Ant1_Interference_sat128_Jam109_Test.dat';
% data files are complex int16
% 2 byte I followed by 2 byte Q for each sample (4 bytes per sample)

freqL2 = 1227.6e6;
freqL1 = 154*10.23e6;
sample_frequency = 5e6;
intermediate_frequency = 0;
integration_period = NUM_MS*(1e-3);


fid = fopen(sprintf('%s',filename),'rb');

gpsPi = 3.14159265359; 

bytes_per_sample = 2;

% Checks over double the integration period (in 1 ms shifts) in case of
% Data Bit Transition
samples_to_read = 2*round(sample_frequency*integration_period);  
NUM_DATA_SETS=NUM_MS+1;


samples_to_skip=SKIP_SECONDS*sample_frequency;    % how many seconds to skip (Initial USRP data is not good)
fseek(fid, samples_to_skip*bytes_per_sample, 'bof');

t = fread (fid, [2,samples_to_read], 'int16');
x_all = t(1,:) + t(2,:)*1i;
%[r, c] = size (v);
%complex_vector_1 = reshape (v, c, r);
fclose(fid);

%disp(x1(1:20))
    
ts=1/sample_frequency;	% sampling time
n=sample_frequency*integration_period;	% data pt in integration period (rounded off to ms to match code length).
nn=[0:n-1];	% total no. of pts % changed from n-1 to n
samples_per_chip=sample_frequency*(1e-3/1023);

fc0=intermediate_frequency;	% center freq without Doppler

code_1ms=cacode_bevly(SV_NUM,0,samples_per_chip);
code=code_1ms;
for k=1:(integration_period/1e-3)-1
    code=[code  code_1ms];
end

% ***** DFT of C/A code
%codefreq = conj(fft(code));
best_result=0;

for dataset=1:NUM_DATA_SETS;
    START_Dataset=(n/NUM_MS)*(dataset-1)+1;
    END_Dataset=START_Dataset+n-1;
    x=x_all(START_Dataset:END_Dataset);
    %x=x_all(1:5000);
    for i=1:41
      fc(i) = fc0 + 0.0005e6*(i-21);
      expfreq=exp(j*2*pi*fc(i)*ts*nn); %reference freq
      signal=x.*expfreq;
      I = real(signal);
      Q = imag(signal);
      IQfreq = fft(I+j*Q);
      codefreq = conj(fft(code));  
      convcodeIQ = IQfreq .* codefreq;
      %result(i,:,dataset) = abs(ifft(convcodeIQ)).^2;
      result(i,:) = abs(ifft(convcodeIQ)).^2;   % modified due to memory limits
    end
    %datasetresult=squeeze(result(:,:,dataset));
    datasetresult=result;
    [peak(dataset) codephase(dataset)]=max(max(datasetresult));
    meanresult=mean(mean(datasetresult));
    [peak(dataset) frequency(dataset)]=max(max(datasetresult'));
    if (peak(dataset)/meanresult>best_result)
        saveresult=result;
        saveddataset=dataset;
        best_result=peak(dataset)/meanresult;
    end;
end

% figure(1)
% plot(peak,'*')
% grid
% xlabel('DataSet #');
% ylabel('Correlation');
frequency = fc(frequency);

codephaseChips = round(1023 - (codephase/11999)*1023);
gold_rate = 1.023e6;			% Gold code clock rate in Hz
ts=1/sample_frequency;
tc=1/gold_rate;
b=[1:n/NUM_MS];      % BEVLY HAD TO ADD ONE HERE - NOT SURE WHY
c=ceil((ts*b)/tc);
x_axis=c; %code phase axis
y_axis=fc/1e6; %dopler shift axis

%[SortedPeaks,SortedDataSets]=sort(peak,2,'descend');
    figure
    %     dataset=SortedDataSets(1);  % modified due to memory
    %     dataset=3;
    %     datasetresult=squeeze(result(:,:,dataset));
    datasetresult=saveresult;
    s=surf(x_axis,y_axis,datasetresult(:,1:n/NUM_MS));
    set(s,'EdgeColor','none','Facecolor','interp');
    axis([min(x_axis) max(x_axis) min(y_axis) max(y_axis) min(min(datasetresult)) max(max(datasetresult))]);
    caxis([0 max(max(datasetresult))]);
    xlabel('Code Phase [chips]');
    ylabel('Frequency [MHz]');
    zlabel('Magnitude');
    text=sprintf('SV Number: %i',SV_NUM);
    title(text);


% figure(3)   
%     plot(x_axis,max(datasetresult(:,1:n/NUM_MS)))

end