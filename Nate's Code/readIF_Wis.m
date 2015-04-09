function Incoming = readIF_Wis(Settings, filename)

%settings.posS  = 700;

% open file
fid = fopen(sprintf('%s',filename));

% Position in the file
fseek(fid, Settings.pos, 'bof');

% number of bytes
bytes_to_read = Settings.acq_block;
Incoming.signalmatrix = zeros(Settings.msToProcess,Settings.acq_block);

for count = 1:Settings.msToProcess
 
    % read  chunk of data
    Incoming.IQmixed = fread(fid,[2, Settings.acq_block] ,'int16');
    signal = Incoming.IQmixed(1,:)+ Incoming.IQmixed(2,:) * 1i;
    
    Incoming.signalmatrix(count,:) = signal;
end

Incoming.SigSize = size(Incoming.signalmatrix);
intermediatestep = Incoming.signalmatrix';
Incoming.signal = intermediatestep(:)';
Incoming.correlationSignal = Incoming.signal(1:2*Incoming.SigSize(2));

% close file when done
fclose(fid);

