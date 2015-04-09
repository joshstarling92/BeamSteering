clc
close all
for chnl = 1:Settings.numberOfChannels
    % Retrive PRN from structure
    PRN = trackingResults(chnl).PRN;

    % Create Figure for that PRN
    figure(chnl)
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box',...
           'off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,['\bf Satellite PRN: ',num2str(PRN)],'HorizontalAlignment' ...
                                      ,'center','VerticalAlignment', 'top')
    % In-Phase Prompt Plotting
    subplot(2,2,1)
    plot(((trackingResults(chnl).I_P)))
    title('Prompt In-Phase Arm')
    xlabel('Time (ms)')
    
    % In-Phase Power Early-Late
    subplot(2,2,2); hold on
    plot(abs(trackingResults(chnl).I_P)) 
    plot(abs(trackingResults(chnl).I_E),'r')
    plot(abs(trackingResults(chnl).I_L),'g')
    title('In-Phase Power Levels')
    xlabel('Time (ms)')
    legend('Prompt Power', 'Early Power', 'Late Power')
    
    % Quadrature Arm Plotting
    subplot(2,2,4); hold on
    plot(abs(trackingResults(chnl).Q_P)) 
    plot(abs(trackingResults(chnl).Q_E),'r')
    plot(abs(trackingResults(chnl).Q_L),'g')
    title('Quadrature Phase Power Levels')
    xlabel('Time (ms)')
    legend('Prompt Power', 'Early Power', 'Late Power')
    
    % Plot Doppler Frequency
    subplot(2,2,3)
    plot(trackingResults(chnl).carrFreq-Settings.IF,'r')
    title('Doppler Frequency')
    xlabel('Time (ms)')
    ylabel('Frequency (hz)')
%     % Title Figure
%     figure(chnl)
%     
%    
%     figure(chnl)
end
