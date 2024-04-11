X=[0 0 1 1 1 0 0];
Ts=3;
[S_in, Ts_out, powerOfBlock] = TX_optical_eml(X, Ts, params);

Fs=1/Ts;
attenuation=10^(0.2/10);
Ts=3;
delta_lambda=2/Ts;
L=1e3;
D=17;%ps/nm/km
S=0.09;%ps/nm^2/km
lambda=1550;
c=3e8;
beta2=-D*lambda^2/(2*pi*c);
beta3=S*lambda^4/(4*pi^2*c^2);

S_out=opticalFiber(S_in,Fs,attenuation,beta2,beta3,L);
bins=[0:length(S_in)-1];
plot(bins,fft(S_out));
