function Trajectory = spoofer(magnitude, trajectory, startTime, endTime, Frequency, POS_2_ENU)

% Generate vector of "time" the length of the position vectors.

%t = 0:1/Frequency:(length(POS_2_ENU)- 1 - startTime * Frequency)/Frequency;
t = 0:1/Frequency:(endTime * Frequency - 1 - startTime * Frequency)/Frequency;
% Change the magnitude of the spoofing "Time Vector" to increase the
% spoofing effect.

t = (magnitude * t)';

% Initialize the spoofed position as the vehicle position and the start
% time

StartPos = POS_2_ENU(startTime * Frequency,:);
type = 'add';

switch type
    case 'append'
        switch trajectory
            case 'linear'
               x = t;
               y = t;
               z = 0*t; 

            case 'circular'
                x = t;
                y = sqrt(x.^2 + 50*ones(length(t),1));
                z = 0*t; 
            case 'sinusoidal'
                x = magnitude * sin(t);
                y = -magnitude*t;
                z = 0*t;
        end
        spoofedTrajectory = ones(length(t),length(StartPos)) * diag(StartPos) + [x y z];
        POS_2_ENU(startTime * Frequency + 1:endTime * Frequency,:) = spoofedTrajectory;
    case 'add'
        x = 25 * sin(t/(13*pi));
        y = 0.05 * t;
        z = 0 * t;
        spoofedTrajectory =  [x y z];
        POS_2_ENU(startTime * Frequency + 1:endTime * Frequency,:) = POS_2_ENU(startTime * Frequency + 1:endTime * Frequency,:) + spoofedTrajectory;
        
end

% circular






Trajectory = POS_2_ENU;
return

