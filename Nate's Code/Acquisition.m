function [Corr_results,Incoming, acqChannel] = Acquisition(Settings, init, filename,number_of_antennas)


%% Initialize Acquisition Structure =======================================
acqChannel.PRN = 0;
acqChannel.acquiredFreq = 0;
acqChannel.codePhase = 0;
Incoming.a = 0;
Corr_results.a = 0;
Mod_corr_results.a = 0;
acqChannel = repmat(acqChannel, 1, Settings.numberOfChannels);

%% Read in Data and create time vector

switch init.run
    case 'nordnav'
        for i = 1:number_of_antennas;
            Incoming = readIF(Incoming,Settings, filename(i,:),i);
        end
    case 'wis'
        Incoming = readIF_Wis(Settings, filename);
end
SigSize = Incoming.SigSize(i,:);
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
if init.spoof == 1;
    
    % Determine code phase in terms of samples
    spoofCodePhase = floor(init.spoofCodePhaseChip * Settings.Sf/(1e3*1023));
    
    % Create spoofed CA
    spoofed_ca = init.spoofMagnitude*circshift(ca, [0  spoofCodePhase])...
          .*cos(2*pi*(Settings.IF+init.spoofDoppShift).*timeVec-3.5/(pi));
      
    % Introduce spoofed peak to correlation signal
    Incoming.correlationSignal = Incoming.correlationSignal...
                                   + [spoofed_ca, spoofed_ca];
end

 
%% FFT Search =============================================================
% Initialize matrices, vectors, figures, and counters
channelCount = 0;
[Y, X] = meshgrid(dopp_range,1:floor(Settings.Sf*10^-3));
S1 = zeros(Settings.code_block,length(dopp_range),33);
PeaksVec = zeros(3,33);   % Row 1 = peak val, Row 2 = dopp shift, row 3 = code shift
AcqrdSatts = zeros(33,1);

% Open figures to reduce processing time
for figCnt = 1:length(init.sattId)
    %figure(figCnt)
end


% Initialize waitbar if acquisition will take long
if Settings.acq_int_MS * length(dopp_range) < 500; 
   h = waitbar(0,'Correlating') ; end

for sattelite = 1:length(init.sattId) 
    if exist('h','var')
        waitbar(sattelite/length(init.sattId),...
        h,['Correlating Sattelite PRN: ' num2str(init.sattId(sattelite))]); 
    end

    for i = 1:number_of_antennas
        Corr_results = fftAcqCorrelation(Corr_results,Settings,... 
                        Incoming, init, dopp_range, timeVec,...
                       caMatrix(init.sattId(1,sattelite),:),i,sattelite);  
        unmod_signal(:,:,i,sattelite) = Corr_results.I(:,:,i,sattelite)+1i*Corr_results.Q(:,:,i,sattelite);
    end
    %Beam Stearing using power minimization
    delta = [1;0;0];
    %Rxx = (Incoming.signal*ctranspose(Incoming.signal));
    %W = 1/(delta'*inv(Rxx)*delta)*inv(Rxx)*delta;
    W = [1; -0.4360 + 0.4268i];
    modified_CW_signal = W(1)*unmod_signal(:,:,1);%+ W(2)*unmod_signal(:,:,2) ;%+ W(3)*unmod_signal(:,: ,3);
    Mod_corr_results = modfftAcqCorrelation(Mod_corr_results,Settings,... 
                        Incoming, init, dopp_range, timeVec,...
                       caMatrix(init.sattId(1,sattelite),:),modified_CW_signal,sattelite);
    S1(:,:,init.sattId(1,sattelite)) = Mod_corr_results.corr(1:101,1:5000)';
    %S1(:,:,init.sattId(1,sattelite)) = Corr_results.corr(1:101,1:5000)';

   
    if init.plotAcqPlanes == 1
          figure%(sattelite)
          surf(Y,X,S1(1:floor(Settings.Sf*1e-3),:,...
                            init.sattId(1,sattelite)),'EdgeColor','none')
          title(['Correlation Plane: SV ' num2str(init.sattId(1,sattelite))])
          axis([-Inf Inf -Inf Inf -Inf inf])
          ylabel('CA code shift')
          xlabel('Doppler shift (Hz)')
          zlabel('Correlation Value')
          
          %{
          figure(sattelite+1)
          dopp_range = -2000:init.dopp_bin_size:5000;
          [Y, X] = meshgrid(dopp_range,11000:14000);
          surf(Y,X,S1(11000:14000,51:121,...
             init.sattId(1,sattelite)),'EdgeColor','none');
          title(['Correlation Plane: SV ' num2str(init.sattId(1,sattelite))])
          axis([-Inf Inf -Inf Inf -Inf 1])
          ylabel('CA code shift')
          xlabel('Doppler shift (Hz)')
          zlabel('Correlation Value')
          %}
    end
 
    CorrPlaneSv = S1(:,:,init.sattId(sattelite));
    [peakVal, idx] = max(CorrPlaneSv(:));
    [codePhase, dopp_idx] = ind2sub(size(CorrPlaneSv),idx);
    doppShift = dopp_range(dopp_idx);
     % Store doppler shifts % Store code shift % Store peak values
    PeaksVec(1:3,init.sattId(sattelite)) = [peakVal; doppShift; codePhase];
    
    
    thresh = 25*mean(CorrPlaneSv(:));
    if max(CorrPlaneSv(:)) > thresh
        channelCount = channelCount  + 1;
        AcqrdSatts(channelCount) = init.sattId(sattelite);
        acqChannel(channelCount).PRN = init.sattId(sattelite);
        acqChannel(channelCount).acquiredFreq = Settings.IF + doppShift;
        acqChannel(channelCount).codePhase = codePhase;
    end  
   
end

% Update number of channels to the number of acquired sattelites
Settings.numberOfChannels = channelCount;

% Remove trailing zeros from acquired sattelite vector
AcqrdSatts(AcqrdSatts == 0) = [];

if exist('h','var'); close(h); clear h; end

if init.plotAcqResults == 1

    bar(1:33,PeaksVec(1,:));
    if isempty(AcqrdSatts) == 0;
        AcqBarX = 1:33;
        AcqBarY = zeros(1,33);
        AcqBarY(1, AcqrdSatts) = PeaksVec(1,AcqrdSatts);
        hold on; bar(AcqBarX,AcqBarY,'g')
    end
    xlabel('PRN Number (no bar - Sattelite PRN not Searched)')
    ylabel('Acquisition Magnitude')
    title('Acquisition Results')
end
return