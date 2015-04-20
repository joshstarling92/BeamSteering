function set = set_settings(init)

switch init.run
    case 'nordnav'
        set.IF = 1.25e6;    %4.1304e6; %             % Hz for Nordnav
        set.Sf =  5e6;  %16.3676e6; %             % Hz for Nordnav
    case 'wis'
        set.IF = 0;                        % For Wismer
        set.Sf = 5e6;            % For Wismer 
end

set.L1 = 154*10.23e6;                      % Hz
set.GPSpi = 3.14159265359;
set.c = 299792458; %m/s
set.lamda = 3/4*set.L1/set.c;
% Integration Periods
posS = 354;                                 % Position in data file in Sec
set.pos = floor(posS * set.Sf);           % Position in data file in bytes
set.acq_int_MS = 3;
set.acq_int_Sec = set.acq_int_MS*1e-3;
set.track_int_MS = 1;
set.track_int_Sec = set.track_int_MS*1e-3;

% Processing Settings
set.msToProcess = 1500;                         % ms
set.acq_block = floor(set.Sf*set.acq_int_Sec);  % Acq block in samples
set.code_block = floor(set.Sf*1e-3);            % Length of C/A in samples
set.codeFreqBasis = 1.023e6;                    %[Hz]            
set.codeLength          = 1023;                 % Length of C/A in chips

%% Tracking loops settings ================================================
% Code tracking loop parameters
set.dllDampingRatio         = 0.7;
set.dllNoiseBandwidth       = 2;       %[Hz]
set.dllCorrelatorSpacing    = .3;     %[chips]
set.numberOfChannels = 2;

% Carrier tracking loop parameters
set.pllDampingRatio         = 0.7;
set.pllNoiseBandwidth       = 10;      %[Hz]

return