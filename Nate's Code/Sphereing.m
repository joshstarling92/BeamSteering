clear all
close all
clc


figure
[X, Y, Z] = sphere(1000);

max = 20000;                %Kilometers
mi = 20000;
C = ones(length(X));
C1 = C * 0.2;
XA = max*X; YA = max*Y; ZA = max*Z;
XB = mi*X; YB = mi*Y; ZB = mi*Z;


surf(XA+8000,YA+2000,ZA+1000,C1,'EdgeColor','none')
alpha(.99)
hold on

surf(XA,YA,ZA,C,'EdgeColor','none')
alpha(0.19)

