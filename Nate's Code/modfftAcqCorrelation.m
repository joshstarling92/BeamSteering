function Mod_corr_results = fftAcqCorrelation(Mod_corr_results,settings, incoming, init, dopp_range, t, ca,mod_signal,sattelite)
%This function takes in IF data in 'signal' a range of doppler frequencies
%to search, the IF values, a time vector, and ca code for a specific
%sattelite.  It generates as its main output, a correlation matrix
%containing correlation peaks.
% I = real(signal);
% Q = imag(signal);
Mod_corr_results.I = real(mod_signal);
Mod_corr_results.Q = imag(mod_signal);

IF = settings.IF;
SigSize = incoming.SigSize(1,:);
Mod_corr_results.corrshort(:,:) = zeros(length(dopp_range),settings.code_block);
corrMat = zeros(length(dopp_range),settings.code_block);
t = t(1:settings.acq_block); % Cut timevec down to 1 ms length

if settings.acq_int_MS * length(dopp_range) > 500; 
    h = waitbar(0,'Checking for Correlation'); 
end

for tCount = 1:1%settings.acq_int_MS;
    state = 0;
    doppCount = 0;

    oneMsCa = ca((tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block);
    for doppler = dopp_range
        doppCount = doppCount + 1;
        if exist('h','var'); waitbar(tCount/settings.acq_int_MS,h); end

        switch init.run
            case 'nordnav'
                X1 = fft(Mod_corr_results.I(doppCount,:)+1i*Mod_corr_results.Q(doppCount,:));
            case 'wis'
                X1 = fft(I + 1i * Q);
                
        end
        F_CA = conj(fft(oneMsCa));%/length(X1);
        % Mod_corr_results.value(doppCount,:,antenna) = Mod_corr_results.I(doppCount,:,antenna);
        corrshort = (abs(ifft(X1.*F_CA))).^2;
        Mod_corr_results.corrMat(doppCount,:) = squeeze(corrshort(1,1:5000)); % trim matrix
    end
    Mod_corr_results.corrshort(:,:) = Mod_corr_results.corrshort(:,:) + Mod_corr_results.corrMat(:,:);
    %Mod_corr_results.corr(:,:) = corrMat(:,:);
end
Mod_corr_results.corrshort(:,:) = abs(Mod_corr_results.corrshort(:,:));
    
if exist('h','var') == 1;
    close(h); 
end


if settings.acq_int_MS <= 10
    if max(Mod_corr_results.corrshort(:,:)) <= 5 * mean(Mod_corr_results.corrshort(:,:))
    	fprintf('Mod using long\n')
        clear Mod_corr_results.corr X1 F_CA corr corrMat
        Mod_corr_results.corrlong(:,:) = zeros(length(dopp_range),settings.acq_block);
        corrMatlong = zeros(length(dopp_range),settings.acq_block);
        for tCount = 1:settings.acq_int_MS;
            doppCount = 0;
                                          
            for doppler = dopp_range
                doppCount = doppCount + 1;
                if exist('h','var'); waitbar(tCount/settings.acq_int_MS,h); end

                X1 = fft(Mod_corr_results.I(doppCount,:)+1i*Mod_corr_results.Q(doppCount,:));
                F_CA = conj(fft(oneMsCa));%/length(Mod_corr_results.I);

                corr = (abs(ifft(X1.*F_CA))).^2;
                corrMatlong(doppCount,:) = corr(1,1:settings.acq_block)'; % trim matrix
        
            end
            Mod_corr_results.corrlong(:,:) = Mod_corr_results.corrlong(:,:) + corrMatlong(:,:);
        end 
        state = 1;
    end
end
if state == 1
    Mod_corr_results.corr(:,:) = abs(Mod_corr_results.corrlong(:,1:settings.code_block));
else 
    Mod_corr_results.corr(:,:) = abs(Mod_corr_results.corrshort(:,:));
end
return    
if exist('h','var') == 1; close(h); end

return