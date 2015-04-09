function [Signals,W] = WeightCalculation(Settings,Signals)
switch Settings.run_type
    case  0
        fprintf('Preforming Null Steering\n')
        %%Null Steering - can block K-1 number of direction
        elevation = [54*pi/180 43*pi/180 43*pi/180]; %Elevation of signal(s) that are to be blocked (rads)
        azimuth = [316*pi/180 186*pi/180 196*pi/180]; %Azimuth of signal(s) that are to be blocked (rads)
        
        p(1,:) = [0 0 0];
        p(2,:)= [.5*lamda 0 0]; %vector from second antenna to reference antenna (m)
        %p(3,:)= [0 -.5*lamda 0]; %vector from third antenna to reference antenna (m)
        
        r(1,:) = [sin(azimuth(1))*cos(elevation(1)) cos(azimuth(1))*cos(elevation(1)) sin(azimuth(1))]; %bore sight vector to interference
       % r(2,:) = [sin(azimuth(2))*cos(elevation(2)) cos(azimuth(2))*cos(elevation(2)) sin(azimuth(2))]; %bore sight vector to interference
        b = [1;-1/(Settings.NumberOfAntennas-1)];
        a = [exp(j*(2*pi*dot(p(1,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(2,:),r(1,:))/lamda))];
        W = a.*b; %weight to be applied to signal
        
        Signals.modified_CW_signal = W(1)*Signals.CW_signal(:,:,1) + W(2)*Signals.CW_signal(:,:,2); %+ W(3)*Signals.CW_signal(:,:,3);
    case 1
        %%Beam Steering
        fprintf('Preforming Beam Steering\n')
        elevation = [12*pi/180 12*pi/180 ]; %Elevation of signal(s) that are to be strenghtened (rads)
        azimuth = [172*pi/180 172*pi/180]; %Azimuth of signal(s) that are to be strenghtened (rads)
        
        p(1,:) = [0 0 0];
        p(2,:)= [lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from second antenna to reference antenna (m)
        p(3,:)= [-lamda/2*cosd(60) -lamda/2*sind(60) 0]; %vector from third antenna to reference antenna (m)
        
        r(1,:) = [cos(azimuth(1))*sin(elevation(1)) sin(azimuth(1))*sin(elevation(1)) cos(elevation(1))]; %bore sight vector to satellite
        r(2,:) = [cos(azimuth(2))*sin(elevation(2)) sin(azimuth(2))*sin(elevation(2)) cos(elevation(2))]; %bore sight vector to satellite
        W = [exp(j*(2*pi*dot(p(1,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(2,:),r(1,:))/lamda)); exp(j*(2*pi*dot(p(3,:),r(2,:))/lamda))];
        
        Signals.modified_CW_signal = W(1)*Signals.CW_signal(:,:,1) + W(2)*Signals.CW_signal(:,:,2);

    case 2
        %%User Chosen Weights
        fprintf('User \n')
        W = [   0.5000 - 0.3497i; -0.4360 + 0.4268i];
        Signals.modified_CW_signal = W(1)*Signals.CW_signal(:,:,1) + W(2)*Signals.CW_signal(:,:,2);
        
    case 3
        %Adaptive PM Beam Steering
        %Power Min. without having a look direction 
        delta = [1;0;0];
        Rxx = (Signals.CleanSignal*ctranspose(Signals.CleanSignal));
        W = 1/(delta'*inv(Rxx)*delta)*inv(Rxx)*delta;
        Signals.modified_CW_signal = W(1)*Signals.CW_signal(:,:,1) + W(2)*Signals.CW_signal(:,:,2) + W(3)*Signals.CW_signal(:,:,3);
        
    case 4
        %Power Min. with a look direction 
        elevation = [46*pi/180 42*pi/180 41*pi/180 42*pi/180 32*pi/180 31*pi/180 22*pi/180 31*pi/180]; %Elevation of satellite(s) (rads)
        azimuth = [189*pi/180 351*pi/180 275*pi/180 101*pi/180 330*pi/180 168*pi/180 76*pi/180 182*pi/180]; %Azimuth of satellite(s) (rads)
        Rxx = (Signals.CleanSignal*ctranspose(Signals.CleanSignal));
        d = 3*lamda/4;
%         a(:,1) = [exp(-j*2*pi/lamda*d*cos(azimuth(1))*sin(elevation(1))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(1))*sin(elevation(1)))]'; %#ok<*IJCL>
%         a(:,2) = [exp(-j*2*pi/lamda*d*cos(azimuth(2))*sin(elevation(2))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(2))*sin(elevation(2)))]';
%         a(:,3) = [exp(-j*2*pi/lamda*d*cos(azimuth(3))*sin(elevation(3))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(3))*sin(elevation(3)))]';
%         a(:,4) = [exp(-j*2*pi/lamda*d*cos(azimuth(4))*sin(elevation(4))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(4))*sin(elevation(4)))]';
%         a(:,5) = [exp(-j*2*pi/lamda*d*cos(azimuth(5))*sin(elevation(5))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(5))*sin(elevation(5)))]';
%         a(:,6) = [exp(-j*2*pi/lamda*d*cos(azimuth(6))*sin(elevation(6))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(6))*sin(elevation(6)))]';
%         a(:,7) = [exp(-j*2*pi/lamda*d*cos(azimuth(7))*sin(elevation(7))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(7))*sin(elevation(7)))]';
%         a(:,8) = [exp(-j*2*pi/lamda*d*cos(azimuth(8))*sin(elevation(8))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(8))*sin(elevation(8)))]';
        for i = 1:8
            a(:,i) = [exp(-j*2*pi/lamda*d*cos(azimuth(i))*sin(elevation(i))) 1 exp(-j*2*pi/lamda*d*sin(azimuth(i))*sin(elevation(i)))]'; %#ok<*IJCL>
        end
        OneJ = ones(8,1);
        W = inv(Rxx)*a*inv(ctranspose(a)*inv(Rxx)*a)*OneJ;
        
        Signals.modified_CW_signal = W(1)*Signals.CW_signal(:,:,1) + W(2)*Signals.CW_signal(:,:,2) + W(3)*Signals.CW_signal(:,:,3);

    case 5
        %Real-weight control NS
        W = -cos(-(2*pi/lamda)*(3/4*lamda)*cosd(275)+(2*pi/lamda)*(3/4*lamda)*cosd(0));
        psi = (2*pi/lamda)*(3/4*lamda)*cosd(275)+(2*pi/lamda)*(3/4*lamda)*cosd(90);
        B = -2*cos(psi);
        Signals.modified_CW_signal = W(1)*Signals.CW_signal(:,:,1) + W(2)*Signals.CW_signal(:,:,2) + W(3)*Signals.CW_signal(:,:,3);   
end
end