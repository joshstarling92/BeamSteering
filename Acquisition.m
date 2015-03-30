function [out_I,out_Q] = Acquisition(incoming_signal,SV_input_array,sample_frequency,integration_period,intermediate_frequency,NUM_MS,NUM_DATA_SETS)

for SV_array = SV_input_array
	ts=1/sample_frequency;	% sampling time
	n=sample_frequency*integration_period;	
		% data pt in integration period (rounded off to ms to match code length).
	nn=[0:n-1];	% total no. of pts % changed from n-1 to n
	samples_per_chip=sample_frequency*(1e-3/1023); 

	fc0=intermediate_frequency;	% center freq without Doppler

	code_1ms=cacode_bevly(SV_array,0,samples_per_chip);
	code=code_1ms;

	for k=1:(integration_period/1e-3)-1
	    code=[code  code_1ms];
	end


	%discrete Fourier transform (DFT) of C/A code
	%codefreq = conj(fft(code));
	best_result=0; %initialize constant

	for dataset=1:NUM_DATA_SETS;
	    START_Dataset=(n/NUM_MS)*(dataset-1)+1;
	    END_Dataset=START_Dataset+n-1;
	    read_signal=incoming_signal(START_Dataset:END_Dataset);
        if dataset == 2
            value = 41;
        else
            value = 0;
        end
	    for i=1:41
	      fc(i) = fc0 + 0.0005e6*(i-21);
	      expfreq=exp(j*2*pi*fc(i)*ts*nn);
	      signal=read_signal.*expfreq;
	      I = imag(signal);
	      Q = real(signal);
          temp_out_I(i+value,:) = I;
          temp_out_Q(i+value,:) = Q;
	      IQfreq = fft(I+j*Q);
	      codefreq = conj(fft(code));  
	      convcodeIQ = IQfreq .* codefreq;
	      %result(i,:,dataset) = abs(ifft(convcodeIQ)).^2;
	      result(i,:) = abs(ifft(convcodeIQ)).^2;   % modified due to memory limits
%           realnum(dataset,:) = I;
%           imaginary(dataset,:) = Q;
	    end
	    datasetresult=result;
	    [peak(dataset) codephase(dataset)]=max(max(datasetresult));
	    meanresult=mean(mean(datasetresult));
	    [peak(dataset) frequency(dataset)]=max(max(datasetresult'));
	    if (peak(dataset)/meanresult>best_result)
	        saveresult=result;
	        saveddataset=dataset;
	        best_result=peak(dataset)/meanresult;
	    end;
	end
out_I = temp_out_I;
out_Q = temp_out_Q;
	frequency = fc(frequency);
	codephaseChips = round(1023 - (codephase/11999)*1023);
	gold_rate = 1.023e6;			% Gold code clock rate in Hz
	ts=1/sample_frequency;
	tc=1/gold_rate;
	b=1:n/NUM_MS;      % BEVLY HAD TO ADD ONE HERE - NOT SURE WHY
	c=ceil((ts*b)/tc);
	x_axis=c;code;
	y_axis=fc/1e6;

% 	figure
% 	    datasetresult=saveresult;
% 	    s=surf(x_axis,y_axis,datasetresult(:,1:n/NUM_MS));
% 	    set(s,'EdgeColor','none','Facecolor','interp');
% 	    axis([min(x_axis) max(x_axis) min(y_axis) max(y_axis) min(min(datasetresult)) max(max(datasetresult))]);
% 	    caxis([0 max(max(datasetresult))]);
% 	    xlabel('Code Phase [chips]');
% 	    ylabel('Frequency [MHz]');
% 	    zlabel('Magnitude');
% 	    text=sprintf('SV Number: %i',SV_array);
% 	    title(text);
 

end

