function ReacqChannel = Reacquisition(Settings, init, filename, trackingResults)
%% Initialize Acquisition Structure =======================================
ReacqChannel.PRN = 0;
ReacqChannel.acquiredFreq = 0;
ReacqChannel.codePhase = 0;
Incoming.a = 0;
ReacqChannel = repmat(ReacqChannel, 1, Settings.numberOfChannels);

A = 800;


%% Read in Data and create time vector
switch init.run
    case 'nordnav'
        Incoming = readIF(Incoming,Settings,filename,1);
    case 'wis'
        Incoming = readIF_Wis(Settings, filename);
end

Incoming.signal = trackingResults(1).cleanedSig(1,:);
% length(Incoming.signal)
% Incoming.SigSize
Incoming.correlationSignal = Incoming.signal(A*16367:2*Incoming.SigSize(2)+A*16367);
SigSize = Incoming.SigSize;
timeVec = 0:1/Settings.Sf:(SigSize(2) - 1)/Settings.Sf;

% Add random noise to the signal if desired
if init.addNoise == 1
    sigma = 10;
    noise = sigma * randn(1,length(fullsignal));
    Incoming.signal = fullsignal + noise;
end

%% Create and resample CA code

% Initialize matrix to hold ca codes
caMatrix = zeros(33,SigSize(2));
if length(init.sattId) > 5; h = waitbar(0,'Generating C/A Code'); end

for sattelite = 1:length(init.sattId)
    if exist('h','var'); waitbar(sattelite/length(init.sattId)); end
    ca = cacodeA(Settings, init.sattId(1,sattelite),...
                 Settings.Sf/(1.023*10^6),Settings.acq_int_MS);
             
    ca(ca == 0) = -1;   ca(ca == 1) = 1;   % Reset CA values
    caMatrix(init.sattId(1,sattelite),:) = ca;  % Store CA in row of matrix
end
if exist('h','var'); close(h); end
clear sattelite h

% Set doppler search window
dopp_range = -5000:init.dopp_bin_size:5000;



%% Introduce Spoofed C/A ==================================================
% if init.spoof == 1;
%     
%     % Determine code phase in terms of samples
%     spoofCodePhase = floor(init.spoofCodePhaseChip * Settings.Sf/(1e3*1023));
%     
%     % Create spoofed CA
%     spoofed_ca = spoofMagnitude*circshift(ca, [0  spoofCodePhase])...
%           .*cos(2*pi*(Settings.IF+init.spoofDoppShift).*timeVec-3.5/(pi));
%       
%     % Introduce spoofed peak to correlation signal
%     Incoming.correlationSignal = Incoming.correlationSignal...
%                                    + [spoofed_ca, spoofed_ca];
% end

%% FFT Search =============================================================
% Initialize matrices, vectors, and counters
channelCount = 0;
[Y, X] = meshgrid(dopp_range,1:floor(Settings.Sf*10^-3));
S1 = zeros(floor(Settings.Sf * 1e-3),length(dopp_range),33);
PeaksVec = zeros(3,33);   % Row 1 = peak val, Row 2 = dopp shift, row 3 = code shift
AcqrdSatts = zeros(33,1);

% for figCnt = 34:length(init.sattId) + 34
%     figure(figCnt)
% end

% Initialize waitbar if acquisition will take long
if Settings.acq_int_MS * length(dopp_range) < 500; 
   h = waitbar(0,'Correlating') ; end

for sattelite = 1:length(init.sattId) 
    
    if exist('h','var')
        waitbar(sattelite/length(init.sattId),...
        h,['Correlating Sattelite PRN: ' num2str(init.sattId(sattelite))]); 
    end
 
    Corr_results = fftAcqCorrelation(Settings,... 
                   Incoming, init, dopp_range, timeVec,...
                   caMatrix(init.sattId(1,sattelite),:));
    size(S1)
    size(Corr_results.corr')
    S1(:,:,init.sattId(1,sattelite)) = Corr_results.corr';

    if init.plotAcqPlanes == 1
        figure(sattelite + 33+A)
        surf(Y,X,S1(1:floor(Settings.Sf*1e-3),:,...
             init.sattId(1,sattelite)),'EdgeColor','none');
          title(['Correlation Plane: SV ' num2str(init.sattId(1,sattelite))])
          axis([-Inf Inf -Inf Inf -Inf Inf])
          ylabel('CA code shift')
          xlabel('Doppler shift (Hz)**')
          zlabel('Correlation Value')
    end
 
    CorrPlaneSv = S1(:,:,init.sattId(sattelite));
    [peakVal, idx] = max(CorrPlaneSv(:));
    [codePhase, dopp_idx] = ind2sub(size(CorrPlaneSv),idx);
    doppShift = dopp_range(dopp_idx);
     % Store doppler shifts % Store code shift % Store peak values
    PeaksVec(1:3,init.sattId(sattelite)) = [peakVal; doppShift; codePhase];
    
    
    thresh = 5 * mean(CorrPlaneSv(:));
    if max(CorrPlaneSv(:)) > thresh
        channelCount = channelCount  + 1;
        ReacqChannel(channelCount).PRN = init.sattId(sattelite);
        ReacqChannel(channelCount).acquiredFreq = Settings.IF + doppShift;
        ReacqChannel(channelCount).codePhase = codePhase;
    end  
   
end

if exist('h','var'); close(h); clear h; end
AcqrdSatts = AcqrdSatts(1:channelCount);

if init.plotAcqResults == 1
    b = figure;
    bar(1:33,PeaksVec(1,:));
    AcqBarX = 1:33;
    AcqBarX(1,AcqrdSatts) = AcqrdSatts;
    AcqBarY = zeros(1,33);
    AcqBarY(1, AcqrdSatts) = PeaksVec(1,AcqrdSatts);
    hold on; bar(AcqBarX,AcqBarY,'g')
    xlabel('PRN Number (no bar - Sattelite PRN not Searched)')
    ylabel('Acquisition Magnitude')
    title('Acquisition Results')
end