%%  Function to read binary file
%Inputs are the fileid, which is where the file is located,
%samp_freq--> sampling frequency used at time of test(should be in the
%filename)
%ingegrate_period is the time over which to integrate the parallel search

%delay_sec is the number of seconds to not sample from. USRPs do not have
%good data in the first 30ish seconds
%this reads data that has been recorded in complex int16 format


function [signal,signal_comp]  = read_USRP_data(file,sample_frequency,integration_period,skip_seconds,bytes_per_sample)
fid = fopen(sprintf('%s',file),'rb');
% Checks over double the integration period (in 1 ms shifts) in case of
% Data Bit Transition
samples_to_read = 2*round(sample_frequency*integration_period);  

samples_to_skip=skip_seconds*sample_frequency;   
% how many seconds to skip (Initial USRP data is not good)

fseek(fid, samples_to_skip*bytes_per_sample, 'bof');

t = fread (fid, [2, samples_to_read], 'short');
signal = t(1,:) + t(2,:)*1i;
signal_comp(1,:) = t(1,:);
signal_comp(2,:) = t(2,:)*1i;
fclose(fid);


