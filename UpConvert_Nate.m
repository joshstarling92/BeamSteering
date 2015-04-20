% int8modulatebyfs4.m
% Upmodulate by fs/4
% Create a real signal by upconverting by a carrier of freq fs/4 using the %formula:
% mixer output = Real(data)*( cos2*pi*n*fm/fs)-Imag(data)*(sin2*pi*n*fm/fs)
tic
clc; clear all; close all;
%disp('Trims off first block of 2000 samples')
% Number of complex samples in a block
blksize = 5000000*.1;
processLength = 100/.1;

% samplingFreq = 5.0e6 [Hz]
% carrierOffset = 4e3  [Hz]
% phase = 0:(blksize-1);
% angleVal = 2*pi.*phase/blksize;
% modulation = exp(j*angleVal); %
% Where the cos values are [ 1 0 -1 0 ] repeating pattern
% Where the sin values are [ 0 1 0 -1 ] repeating pattern
cosv4 =  [ 1; 0; -1; 0]; % [ 1; cos(pi/4); 0; cos(3*pi/4); -1; cos(5*pi/4); 0; cos(7*pi/4)];%
sinv4 = [ 0; 1; 0; -1];  % [ 0; sin(pi/4); 1; sin(3*pi/4); 0; sin(5*pi/4); -1; sin(7*pi/4)];% 
maximumV = 0;
minimumV = 0;

% Replicate sine & cosine cycles
repmat    = ones(1,blksize/length(cosv4));
modsine   = sinv4 * repmat;
modcosine = cosv4 * repmat;
modsine   = modsine(:)';
modcosine = modcosine(:)';

fileInput = ...
 'April_13_5M_Sim03042015_int16_Simulator_Ant3_Interference_sat125_int99_GSPDOLink.dat';

fileOutput = ...
 'April_13_Sim_sat125_int99_Skip354_Ant3.dat';
% if exist(fileOutput,'file') ~= 0
%     fprintf('\n##### WARNING: This file exists. Would you like to overwrite it? ##### \n')
%     disp(fileOutput);
%     cont = input('... Enter 1 to continue. Enter 0 to exit: ');
%     if cont == 0; return; end
% end
disp( [ 'Processing input file ' fileInput ] )
% Open file handle
fin = fopen(fileInput,'rb');
% Open file handle
fout = fopen(fileOutput,'wb');

Sf =  5e6;
SkipPos = Sf*354;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in first block, drop samples, and gets rid of transients
fseek(fin, SkipPos, 'bof');
[rawSignalx2, samplesReadx2] = fread(fin, 2*blksize, 'int16');
if (samplesReadx2 ~= 2*blksize)
 disp('Cannot read first block of data - exiting program!')
 fclose(fin);
 fclose(fout);
 return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads input file and finds the maximum and minimum values
a = waitbar(0,'First While Loop');
count = 0;
while count < processLength % ( (samplesReadx2 == 2*blksize) )
     count = count + 1;
     waitbar(count/processLength)
     [rawSignalx2, samplesReadx2] = fread(fin, 2*blksize, 'int16');
     if (samplesReadx2 ~= 2*blksize)
         disp('No more full blocks of data')
         disp('Breaking out of first "while" loop!')
         disp('Minimum')
         disp(minimumV)
         disp('Maximum')
         disp(maximumV)
         break
     end
     samplesRead = samplesReadx2/2;
     % Extract even samples as I
     rawSignalI = rawSignalx2(1:2:(2*blksize));
     % Extract odd samples as Q
     rawSignalQ = rawSignalx2(2:2:(2*blksize));
     % Transpose vector ( Note do not take complex conjugate )
     rawSignalI = rawSignalI.';
     % Transpose vector ( Note do not take complex conjugate )
     rawSignalQ = rawSignalQ.';
     % Complex mix up by fs/4 - Create a real waveform of blksize samples
     outputWaveform = rawSignalI .* modcosine - rawSignalQ .* modsine;
     localmax = max (outputWaveform);
     maximumV = max([localmax maximumV]);
     localmin = min (outputWaveform);
     minimumV = min([localmin minimumV]);
end
close(a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reopens and reloads file from the beginning
% Read in first block, drop samples, and gets rid of transients again
disp('Loading File')
disp('Triming first block of 2000 samples')
disp( [ 'Creating output file ' fileOutput ] )
fseek(fin,5e6*60,-1);
[rawSignalx2, samplesReadx2] = fread(fin, 4000, 'int16');
if (samplesReadx2 ~= 4000)
    disp('Cannot read first block of data - exiting program!')
    fclose(fin);  fclose(fout);
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads input file and writes the modulated I and Q values into output file
cnt = 0;
b = waitbar(0,'Processing Data...');
while cnt < processLength %( (samplesReadx2 == 2*blksize) )
     cnt = cnt + 1;
     waitbar(cnt/processLength)
     [rawSignalx2, samplesReadx2] = fread(fin, 2*blksize, 'int16');
     if (samplesReadx2 ~= 2*blksize)
         disp('No more full blocks of data - exiting program!')
         fclose(fin);
         fclose(fout);
         disp('Minimum')
         disp(minimumV)
         disp('Maximum')
         disp(maximumV)
         return
     end
     samplesRead = samplesReadx2/2;
     % Extract even samples as I
     rawSignalI = rawSignalx2(1:2:(2*blksize));
     % Extract odd samples as Q
     rawSignalQ = rawSignalx2(2:2:(2*blksize));
     % Transpose vector ( Note do not take complex conjugate )
     rawSignalI = rawSignalI.';
     % Transpose vector ( Note do not take complex conjugate )
     rawSignalQ = rawSignalQ.';
     % Complex mix up by fs/4 - Create a real waveform of blksize samples
     outputWaveform = rawSignalI .* modcosine - rawSignalQ .* modsine;
     % Scales the output waveform in order to fit 8 bit data
     [count ] = fwrite(fout, outputWaveform/((max([maximumV abs(minimumV)]))/127) ,'int8');
     % Check if 2*blksize samples written
     if (count ~= blksize )
         disp('Not able to write the specified number of samples for tracking - exiting!')
         fclose(fin);
         fclose(fout);
         disp('Minimum')
         disp(minimumV)
         disp('Maximum')
         disp(maximumV)
         return
     end
end
close(b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Closes files and displays the maximum and minimum
fclose(fin);
fclose(fout);
disp('Minimum')
disp(minimumV)
disp('Maximum')
disp(maximumV)
toc