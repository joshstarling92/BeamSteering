function Incoming = readIF(Incoming,Settings, filename,antenna)


% open file
fid = fopen(sprintf('%s',filename))
filename
% Position in the file
fseek(fid, Settings.pos, 'bof');

% number of bytes
bytes_to_read = Settings.acq_block;

% read  chunk of data
Incoming.signalmatrix(:,:,antenna) = fread(fid,[bytes_to_read, Settings.msToProcess/Settings.acq_int_MS] ,'int8')';
Incoming.SigSize(antenna,:) = size(Incoming.signalmatrix(:,:,antenna));
intermediatestep = Incoming.signalmatrix(:,:,antenna)';
Incoming.signal(antenna,:) = intermediatestep(:)';
Incoming.correlationSignal(antenna,:) = Incoming.signal(antenna,1:2*Incoming.SigSize(antenna,2));
% close file when done
fclose(fid);

