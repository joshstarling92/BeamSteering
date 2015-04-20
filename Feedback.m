n=[50]; d=[1/5 1 0];
figure(1); clf; margin(n,d); grid; hold on;
ess2ramp= 1/200; Kvd=1/ess2ramp;
Kva = n(end)/d(end-1); Kzp = Kvd/Kva;
figure(2); margin(Kzp*n,d); grid; [GM,PM,wpc,wgc]=margin(Kzp*n,d);
wgcd=wgc;
% Calculation of zeta and PM difference.
Mp = 20/100; zeta =sqrt((log(Mp))^2/(pi^2+(log(Mp))^2)); PMd = zeta * 100 + 10;
phimax = (PMd-PM)*pi/180; alpha=(1-sin(phimax))/(1+sin(phimax)); z=wgcd*sqrt(alpha);
p=wgcd/sqrt(alpha);
ngc = conv(n, 1/alpha*Kzp*[1 z]); dgc = conv(d, [1 p]); figure(3); margin(tf(ngc,dgc)); grid; [ncl,dcl]=feedback(ngc,dgc,1,1);
figure(4); step(ncl,dcl); grid;
figure(5); margin(ncl*1.414,dcl); grid;