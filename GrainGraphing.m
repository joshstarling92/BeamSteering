%This Program is Designed for Calculation of A Generalized Null Steering Broadside Beam former Weights Vector
%General Antenna Array Specifications
%Program assumes that all atnennas are the same legnth from reference
%antenna
a = input ('The Number of Array Elements : ') ;
d = input ('The Separation Distance Between Elements : ') ; %meters
f = input ('The Operation Frequency : ') ; % in Hz
%Input Received Signals Arrival
b = input ('The Angle of Incidence of The Desired Source Signal in Degrees: ') ; 
for I = 1 : a-1
    A(I) = input ('The Angle of Incidence of the Undesired Interference Source Signal in Degrees: ') ;
end
c = 299792458; %m/s
lamda = c/f; %meters

%Estimation of The Weight Vector of A Null Steering Beamformer
for I = 1 : a
    for J = 1 : a
        if I ==1
            S ( I , J ) = exp ( i*( 2*pi*d*( J-1 )*cosd( b ))/lamda) ;
            Z ( I ) = 1 ;
        else  
            S ( I , J ) = exp ( i*(2*pi*d*( J-1 )*cosd( A(I-1) ))/lamda)  ;
            Z ( I )= 0 ;
        end
    end
end
W = ( pinv ( S ) * Z' ) ;
%Plot of The Output Radiation Pattern
t = 0 : 0.05 : 2*pi ;
M = 0 ;
for I = 1 : a
    H = exp ( i*(2*pi*d*( I-1 )*cos( t ))/lamda)  ;
    M = M + ( H * W(I) ) ;
end
M = abs ( M ) ;
polar( t , M , '-r' ) , title ( 'The Generalized Null Steering Beam Former Output Radiation Pattern ') , grid on ;