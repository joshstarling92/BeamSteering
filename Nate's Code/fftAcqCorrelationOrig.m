function corr_results = fftAcqCorrelationOrig( settings, incoming, init, dopp_range, t, ca )
%This function takes in IF data in 'signal' a range of doppler frequencies
%to search, the IF values, a time vector, and ca code for a specific
%sattelite.  It generates as its main output, a correlation matrix
%containing correlation peaks.

IF = settings.IF;

SigSize = incoming.SigSize;
correlationSignal1 = incoming.correlationSignal(1:SigSize(2));
corr_results.corr = zeros(length(dopp_range),settings.code_block);
corrMat = zeros(length(dopp_range),settings.code_block);
t = t(1:settings.acq_block); %Cut timevec down to 1 ms length

if settings.acq_int_MS * length(dopp_range) > 500; h = waitbar(0,'Checking for Correlation') ; end
for tCount = 1:1%settings.acq_int_MS;
    doppCount = 0;
    oneMsSig = correlationSignal1((tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block);
    oneMsCa = ca((tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block);
    for doppler = dopp_range
        doppCount = doppCount + 1;
        if exist('h','var'); waitbar(tCount/settings.acq_int_MS,h); end

        switch init.run
            case 'nordnav'
                corr_results.I = cos(2*pi*(IF + doppler)*t).*(oneMsSig);
                corr_results.Q = sin(2*pi*(IF + doppler)*t).*(oneMsSig);

                X1 = fft(corr_results.I+1i*corr_results.Q);
            case 'wis'
                mixedExp = oneMsSig .* exp(1i*2*pi*doppler.*t);
                I = real(mixedExp);
                Q = imag(mixedExp);
                X1 = fft(I + 1i * Q);
        end
        F_CA = conj(fft(oneMsCa));%/length(X1);

        corr = (abs(ifft(X1.*F_CA))).^2;
        corrMat(doppCount,:) = squeeze(corr(1,1:settings.code_block)); % trim matrix
        
    end
    corr_results.corr(:,:) = corr_results.corr(:,:) + corrMat(:,:);
    %corr_results.corr(:,:) = corrMat(:,:);
end
corr_results.corr(:,:) = abs(corr_results.corr(:,:));
    
if exist('h','var') == 1; close(h); end


if settings.acq_int_MS <= 10
   
    if max(corr_results.corr(:)) <= 5 * mean(corr_results.corr(:))
        correlationSignal2 = incoming.correlationSignal(SigSize(2)+1:2*SigSize(2));
        clear corr_results.corr X1 F_CA corr corrMat
        corr_results.corr = zeros(length(dopp_range),settings.acq_block);
        corrMat = zeros(length(dopp_range),settings.acq_block);
        for tCount = 1:settings.acq_int_MS;
            doppCount = 0;
            oneMsSig = correlationSignal2((tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block);
            
            for doppler = dopp_range
                doppCount = doppCount + 1;
                if exist('h','var'); waitbar(tCount/settings.acq_int_MS,h); end

                corr_results.I = cos(2*pi*(IF + doppler)*t).*(oneMsSig);
                corr_results.Q = sin(2*pi*(IF + doppler)*t).*(oneMsSig);

                X1 = fft(corr_results.I+1i*corr_results.Q);
                F_CA = conj(fft(oneMsCa));%/length(corr_results.I);

                corr = (abs(ifft(X1.*F_CA))).^2;
                corrMat(doppCount,:) = corr(1,1:settings.acq_block)'; % trim matrix
        
            end
            corr_results.corr(:,:) = corr_results.corr(:,:) + corrMat(:,:);
        end 
    end
end
corr_results.corr(:,:) = abs(corr_results.corr(:,:));
return    
if exist('h','var') == 1; close(h); end

return

