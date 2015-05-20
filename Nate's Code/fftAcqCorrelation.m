function corr_results = fftAcqCorrelation( corr_results,settings, incoming, init, dopp_range, t, ca, antenna,sattelite)
%This function takes in IF data in 'signal' a range of doppler frequencies
%to search, the IF values, a time vector, and ca code for a specific
%sattelite.  It generates as its main output, a correlation matrix
%containing correlation peaks.

if init.run_type == 3
    %Perfrom LMS on clean signal
    load('April_9_Sim_Skip_354_2MS_DesiredCorrelationSig.mat') %Signal
    d = DesiredCorrelationSignal;

    %Optimal weight chooseing
    X = [incoming.correlationSignal(1,:); incoming.correlationSignal(2,:);incoming.correlationSignal(3,:)];%Signals.CleanSignal(4,:)];
    R = [mean(X(1,:).*X(1,:)) mean(X(1,:).*X(2,:)) mean(X(1,:).*X(3,:));% mean(X(1,:).*X(4,:));
         mean(X(2,:).*X(1,:)) mean(X(2,:).*X(2,:)) mean(X(2,:).*X(3,:));% mean(X(2,:).*X(4,:));
         mean(X(3,:).*X(1,:)) mean(X(3,:).*X(2,:)) mean(X(3,:).*X(3,:))];% mean(X(3,:).*X(4,:))];
         %mean(X(4,:).*X(1,:)) mean(X(4,:).*X(2,:)) mean(X(4,:).*X(3,:)) mean(X(4,:).*X(4,:))];
    P = [mean(X(1,:).*d); mean(X(2,:).*d); mean(X(3,:).*d)];
    Wopt = inv(R)*P;

    %Weight Searching via gradient estimation
    WLMS = Wopt;
    mu = 4*10^-8;
    count = 1;
    while count < 30000
        y(count) = incoming.correlationSignal(1,count)*WLMS(1,count)+incoming.correlationSignal(2,count)*WLMS(2,count)...
            +incoming.correlationSignal(3,count)*WLMS(3,count);%+Signals.CleanSignal(4,count)*W(4,count);

        error(count) = d(count)-y(count);

        WLMS(:,count+1) = WLMS(:,count)-2*mu*error(count)...
            *[incoming.correlationSignal(1,count);incoming.correlationSignal(2,count);...
            incoming.correlationSignal(3,count)];%incoming.correlationSignal(4,count)];

        if isfinite(WLMS(1,count)) == 0
            fprintf('Weights are not a number\n')
            WLMS(:,count+1) = WLMS(:,count-5);
            break
        end
        count = count + 1;
    end
    incoming.correlationSignal = WLMS(1,end)*incoming.correlationSignal(1,:)...
        + WLMS(2,end)*incoming.correlationSignal(2,:) + WLMS(3,end)*incoming.correlationSignal(3,:);

    if sattelite == 1
        Wopt
        WLMS(:,end)
        figure(2)
        plot(incoming.correlationSignal)
        legend('mixed')
        figure(3)
        plot(d)
        legend('Desired')

        figure(4)
        time = 1:length(WLMS);
        plot(time,WLMS(1,:),time,WLMS(2,:),time,WLMS(3,:))

        figure(5)
        plot(error)
        legend('Error')
        figure(1)
    end
end
IF = settings.IF;
SigSize = incoming.SigSize(antenna,:);
correlationSignal1 = incoming.correlationSignal(antenna,1:SigSize(2));
corr_results.corrshort(:,:) = zeros(length(dopp_range),settings.code_block);
corrMat = zeros(length(dopp_range),settings.code_block);
t = t(1:settings.acq_block); % Cut timevec down to 1 ms length  

if settings.acq_int_MS * length(dopp_range) > 500; 
    h = waitbar(0,'Checking for Correlation'); 
end

for tCount = 1:1%settings.acq_int_MS;
    state = 0;
    doppCount = 0;
    oneMsSig = correlationSignal1((tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block);

    oneMsCa = ca((tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block);
    for doppler = dopp_range
        doppCount = doppCount + 1;
        if exist('h','var'); 
            waitbar(tCount/settings.acq_int_MS,h); 
        end

        switch init.run
            case 'nordnav'
                corr_results.I(doppCount,:,antenna,sattelite) = cos(2*pi*(IF + doppler)*t).*(oneMsSig);
                corr_results.Q(doppCount,:,antenna,sattelite) = sin(2*pi*(IF + doppler)*t).*(oneMsSig);
                X1 = fft(corr_results.I(doppCount,:,antenna,sattelite)+1i*corr_results.Q(doppCount,:,antenna,sattelite));
            case 'wis'
                mixedExp = oneMsSig .* exp(1i*2*pi*doppler.*t);
                I = real(mixedExp);
                Q = imag(mixedExp);
                X1 = fft(I + 1i * Q);
    
        end
        F_CA = conj(fft(oneMsCa));%/length(X1);
        corrshort = (abs(ifft(X1.*F_CA))).^2;
        corr_results.corrMat(doppCount,:) = squeeze(corrshort(1,1:5000)); % trim matrix
    end
    corr_results.corrshort(:,:) = corr_results.corrshort(:,:) + corr_results.corrMat(:,:);
    %corr_results.corr(:,:) = corrMat(:,:);
end
corr_results.corrshort(:,:) = abs(corr_results.corrshort(:,:));
    
if exist('h','var') == 1;
    close(h); 
end


if settings.acq_int_MS <= 10
    if max(corr_results.corrshort(:,:)) <= 5 * mean(corr_results.corrshort(:,:))
        correlationSignal2 = incoming.correlationSignal(antenna,SigSize(2)+1:2*SigSize(2));
        clear corr_results.corr X1 F_CA corr corrMat
        corr_results.corrlong(:,:,antenna) = zeros(length(dopp_range),settings.acq_block);
        corrMatlong = zeros(length(dopp_range),settings.acq_block);
        for tCount = 1:settings.acq_int_MS;
            doppCount = 0;
            oneMsSig = correlationSignal2;%(1,(tCount-1) * settings.acq_block + 1 : tCount * settings.acq_block)
                                          
            for doppler = dopp_range
                doppCount = doppCount + 1;
                if exist('h','var'); waitbar(tCount/settings.acq_int_MS,h); end

                corr_results.I(doppCount,:,antenna,sattelite) = cos(2*pi*(IF + doppler)*t).*(oneMsSig);
                corr_results.Q(doppCount,:,antenna,sattelite) = sin(2*pi*(IF + doppler)*t).*(oneMsSig);
                X1 = fft(corr_results.I(doppCount,:,antenna,sattelite)+1i*corr_results.Q(doppCount,:,antenna,sattelite));
                
                F_CA = conj(fft(oneMsCa));%/length(corr_results.I);

                corr = (abs(ifft(X1.*F_CA))).^2;
                corrMatlong(doppCount,:) = corr(1,1:settings.acq_block)'; % trim matrix
        
            end
            corr_results.corrlong(:,:) = corr_results.corrlong(:,:) + corrMatlong(:,:);
        end 
        state = 1;
    end
end
if state == 1
    corr_results.corr(:,:,antenna) = abs(corr_results.corrlong(:,1:settings.code_block));
else
    corr_results.corr(:,:,antenna) = abs(corr_results.corrshort(:,:));
end
return    
if exist('h','var') == 1; close(h); end

return

