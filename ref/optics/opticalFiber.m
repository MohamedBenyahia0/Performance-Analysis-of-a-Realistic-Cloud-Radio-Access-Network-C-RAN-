function S_out= opticalFiber(S_in,Fs,attenuation,beta2,beta3,L)
N=length(S_in);
bins=-N/2:(N-1)/2;
w=bins*Fs/N;
transferFunct=exp(-((w.^2).*(beta2/2)+(w.^3).*(beta3/6)).*L.*1i);
spec_in=fft(S_in,N);
spec_out=spec_in.*transferFunct;

S_out=ifft(spec_out,N);
S_out=S_out*exp(-attenuation/2*(L/10^3));


end