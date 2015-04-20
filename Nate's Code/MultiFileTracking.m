function [trackResults, channel]= tracking(filename, channel, settings, init,W)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute starting 
%                       positions of spreading codes, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.32 2007/01/30 09:45:12 dpl Exp $

%% Initialize result structure ============================================

for count = 1:settings.number_of_antennas
    fid(count) = fopen(sprintf('%s',filename(count,:)));
    fid
end
tic
% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock
% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

trackResults.remCarrPhase   = zeros(1,settings.msToProcess);
trackResults.samplesRead    = zeros(1, settings.msToProcess);

if init.cleanSig == 1
    cleanedSig = zeros(settings.msToProcess, ceil(settings.Sf * 1e-3));
end
    %trackResults.cleanedSig = zeros(settings.msToProcess * ceil(settings.Sf * 1e-3));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels
    qBasebandSignalTemp = zeros(1500,5001,settings.number_of_antennas);
    iBasebandSignalTemp = zeros(1500,5001,settings.number_of_antennas);
    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;
        
        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample) 
        for count = 1:settings.number_of_antennas
            fseek(fid(count), ...
              settings.pos + channel(channelNr).codePhase-1, ...
              'bof');
        end
        channel(channelNr).PRN
        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1023) caCode caCode(1)];

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;
        
        hwb = waitbar(0,['Tracking... PRN ', num2str(trackResults(channelNr).PRN)]);
        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods


%% Read next block of data ------------------------------------------------            
            waitbar(loopCnt/codePeriods)
            % Find the size of a "block" or code period in whole samples
            
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.Sf;
            
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
            % Read in the appropriate number of samples to process this
            % interation
            for count = 1:settings.number_of_antennas
                switch init.run
                    case 'nordnav'
                        [rawSignal(:,count), samplesRead(count)] = fread(fid(count), blksize, 'int8');
                          %transpose vector
                    case 'wis'
                        [IQmixed, samplesRead] = fread(fid,[2, blksize] ,'int16');
                        samplesRead = samplesRead / 2;
                        rawSignal = IQmixed(1,:) + IQmixed(2,:) * 1i;
                end
            end
            rawSignal = rawSignal';
            % If did not read in enough samples, then could be out of 
            % data - better exit 
            if (samplesRead ~= blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);
                return
            end
            
            % Store samples read to determine where to remove spoofing
            %trackResults(channelNr).samplesRead(loopCnt) = samplesRead;
            

%% Set up all the code phase tracking information -------------------------
            % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);
            
            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);
            
            remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;

%% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.Sf;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Store the remainder carrier phase 
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;
            
            % Finally compute the signal to mix the collected data to bandband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));

%% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            for count = 1:settings.number_of_antennas
                switch init.run
                    case 'nordnav'
                        qBasebandSignaltemp(count,:) = carrCos .* rawSignal(count,:);
                        iBasebandSignaltemp(count,:) = carrSin .* rawSignal(count,:);
                        % qBasebandSignalTemp(loopCnt,1:length(qBasebandSignal),count) = qBasebandSignal(count,:);
                        % iBasebandSignalTemp(loopCnt,1:length(iBasebandSignal),count) = iBasebandSignal(count,:);
                        % qBasebandSignalOut(:,:,loopCnt,channelNr) = qBasebandSignal(:,1:4999);
                        % iBasebandSignalOut(:,:,loopCnt,channelNr) = iBasebandSignal(:,1:4999);
                    case 'wis'
                        qBasebandSignal = carrCos .* IQmixed(1,:);
                        iBasebandSignal = carrSin .* IQmixed(1,:);
                end
            end
            IQ = iBasebandSignaltemp(:,:) + 1i*qBasebandSignaltemp(:,:);
            modIQ = IQ(1,:)*W(1)+IQ(2,:)*W(2)+IQ(3,:)*W(3);
            qBasebandSignal = imag(modIQ);
            iBasebandSignal = real(modIQ);
            clear qBasebandSignaltemp iBasebandSignaltemp
            % Now get early, late, and prompt values for each
            %I_p = abs(sum((promptCode .* iBasebandSignal).^2));
            %Q_e = abs(sum((earlyCode  .* qBasebandSignal).^2));
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
%% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);
            
            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

%% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            
            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;
            
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

%% Record various measures to show in postprocessing ----------------------
            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid(1));

            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
            
%% Remove Spoofing from IF signal =========================================
            
            if init.cleanSig == 1
                A = sqrt(2*sqrt(I_P .* I_P - Q_P .* Q_P)/length(iBasebandSignal));
                powerI = (1/length(iBasebandSignal)) * I_p;
                powerQ = (1/length(iBasebandSignal)) * Q_e;
                A = sqrt(2*powerI)-sqrt(2*powerQ);%/(2*powerQ);
%                 disp(I_P)
%                 disp('I^')
%                 disp(Q_L)
%                 disp('Q^')
                replicaSpoofSig = 50*sign(I_P)* A * promptCode .* carrSin;
               
                time = toc;
                if time < 50
                   cleanedSig(loopCnt,1:blksize) = rawSignal - replicaSpoofSig;  
                end
            else 
                %cleanedSig(loopCnt,1:blksize) = rawSignal;
            end 
            clear rawSignaltemp rawSignal IQ

        end % for loopCnt
        close(hwb)
        % Convert signal to row vector and eliminate 0 values
        
        
        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        % trackResults(channelNr).status  = channel(channelNr).status;        
        
    end % if a PRN is assigned
    
    if init.cleanSig == 1
        cleanedSig = cleanedSig';
        cleanedSig = cleanedSig(:)';
        cleanedSig(cleanedSig == 0) = [];
        trackResults(channelNr).cleanedSig = cleanedSig;
    end
    clear qBasebandSignaltemp iBasebandSignaltemp
end % for channelNr 
  % Close the waitbar
  toc
return

