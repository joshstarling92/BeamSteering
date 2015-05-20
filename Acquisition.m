function [Signals] = Acquisition(Settings,Signals,CurrentAntenna)
%Initalize variables to improve performance 
peak = zeros(1,3);
frequency = zeros(1,3);
codephase = zeros(1,3);
result = zeros(41,10000);
fc = zeros(1,41);

if Settings.run_type == 3
	%Perfrom LMS on clean signal
	load('April_9_Sim_Skip_354_2MS_DesiredSig.mat') %Signal
	d = DesiredSignal;
	%Optimal weight chooseing
	X = [Signals.CleanSignal(1,:); Signals.CleanSignal(2,:);Signals.CleanSignal(3,:)];%Signals.CleanSignal(4,:)];
	R = [mean(X(1,:).*X(1,:)) mean(X(1,:).*X(2,:)) mean(X(1,:).*X(3,:));% mean(X(1,:).*X(4,:));
	     mean(X(2,:).*X(1,:)) mean(X(2,:).*X(2,:)) mean(X(2,:).*X(3,:));% mean(X(2,:).*X(4,:));
	     mean(X(3,:).*X(1,:)) mean(X(3,:).*X(2,:)) mean(X(3,:).*X(3,:))];% mean(X(3,:).*X(4,:))];
	     %mean(X(4,:).*X(1,:)) mean(X(4,:).*X(2,:)) mean(X(4,:).*X(3,:)) mean(X(4,:).*X(4,:))];
	P = [mean(X(1,:).*d); mean(X(2,:).*d);mean(X(3,:).*d)]*-1;
	Wopt = inv(R)*P;
    R
	%Weight Searching via gradient estimation
	WLMS = Wopt;
	mu = 1*10^-10;
	count = 1;
	while count < 10000
	    y(count) = Signals.CleanSignal(1,count)*WLMS(1,count)+Signals.CleanSignal(2,count)*WLMS(2,count)...
	        +Signals.CleanSignal(3,count)*WLMS(3,count);%+Signals.CleanSignal(4,count)*W(4,count);

	    error(count) = DesiredSignal(count)-y(count);

	    WLMS(:,count+1) = WLMS(:,count)-2*mu*error(count)...
	        *[Signals.CleanSignal(1,count);Signals.CleanSignal(2,count);...
	        Signals.CleanSignal(3,count)];%Signals.CleanSignal(4,count)];

	    if isfinite(WLMS(1,count)) == 0
	        fprintf('Weights are not a number\n')
	        WLMS(:,count+1) = WLMS(:,count-5);
	        break
	    end
	    count = count + 1;
    end
    Signals.WLMS = WLMS;
	Signals.CleanSignal = WLMS(1,end)*Signals.CleanSignal(1,:) + WLMS(2,end)*Signals.CleanSignal(2,:)...
	            + WLMS(3,end)*Signals.CleanSignal(3,:);
end

for SV_array = Settings.SV_array
	ts=1/Settings.sample_frequency;	% sampling time
	n=Settings.sample_frequency*Settings.integration_period;	
		% data pt in integration period (rounded off to ms to match code length).
	nn=(0:n-1);	% total no. of pts % changed from n-1 to n
	samples_per_chip=Settings.sample_frequency*(1e-3/1023); 

	fc0=Settings.intermediate_frequency;	% center freq without Doppler
	code_1ms=cacode_bevly(SV_array,0,samples_per_chip);
	%code_1ms=cacodeA(SV_array,samples_per_chip,NUM_MS*10^-6);
	code=code_1ms;

	for k=1:(Settings.integration_period/1e-3)-1
	    code=[code  code_1ms];
	end


	%discrete Fourier transform (DFT) of C/A code
	%codefreq = conj(fft(code));
	best_result=0; %initialize constant

	for dataset=1:Settings.NUM_DATA_SETS;
	    START_Dataset=(n/Settings.NUM_MS)*(dataset-1)+1;
	    END_Dataset=START_Dataset+n-1;
	    read_signal=Signals.CleanSignal(CurrentAntenna,START_Dataset:END_Dataset);
        if dataset == 2
            value = 41;
        else
            value = 0;
        end
	    for i=1:41
	      fc(i) = fc0 + 0.0005e6*(i-21);
	      expfreq=exp(1i*2*pi*fc(i)*ts*nn);
	      signal=read_signal.*expfreq;
          Signals.I(i+value,:,CurrentAntenna) = imag(signal);
          Signals.Q(i+value,:,CurrentAntenna) = real(signal);
	      IQfreq = fft(Signals.I(i+value,:,CurrentAntenna)+1i*Signals.Q(i+value,:,CurrentAntenna));
	      codefreq = conj(fft(code));  
	      convcodeIQ = IQfreq .* codefreq;
	      %result(i,:,dataset) = abs(ifft(convcodeIQ)).^2;
	      result(i,:) = abs(ifft(convcodeIQ)).^2;   % modified due to memory limits
%           realnum(dataset,:) = I;
%           imaginary(dataset,:) = Q;
	    end
	    datasetresult=result;
	    [peak(dataset), codephase(dataset)]=max(max(datasetresult));
	    meanresult=mean(mean(datasetresult));
	    [peak(dataset), frequency(dataset)]=max(max(datasetresult'));
	    if (peak(dataset)/meanresult>best_result)
	        saveresult=result;
	        %saveddataset=dataset;
	        best_result=peak(dataset)/meanresult;
        end
    end
	frequency = fc(frequency);
	%codephaseChips = round(1023 - (codephase/11999)*1023);
	gold_rate = 1.023e6;			% Gold code clock rate in Hz
	ts=1/Settings.sample_frequency;
	tc=1/gold_rate;
	b=1:n/Settings.NUM_MS;      % BEVLY HAD TO ADD ONE HERE - NOT SURE WHY
	c=ceil((ts*b)/tc);
	x_axis=c;%code;
	y_axis=fc/1e6;

    if Settings.unmod_acq_graph == 1 && Settings.unmod_done== 0;
        figure(find(Settings.SV_array == SV_array))
	    datasetresult=saveresult;
	    s=surf(x_axis,y_axis,datasetresult(:,1:n/Settings.NUM_MS));
	    set(s,'EdgeColor','none','Facecolor','interp');
	    axis([min(x_axis) max(x_axis) min(y_axis) max(y_axis) min(min(datasetresult)) max(max(datasetresult))]);
	    caxis([0 max(max(datasetresult))]);
	    xlabel('Code Phase [chips]');
	    ylabel('Frequency [MHz]');
	    zlabel('Magnitude');
	    text=sprintf('SV Number: %i',SV_array);
        title(text);
    end
    
end
return

