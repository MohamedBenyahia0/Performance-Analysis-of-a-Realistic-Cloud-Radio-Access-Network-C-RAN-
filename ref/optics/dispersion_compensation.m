function S_out= dispersion_compensation(S_in,Fs,beta2,beta3,L)
N=length(S_in);
bins=-N/2:(N-1)/2;
w=bins*Fs/N;
inverse_transferFunct=exp(((w.^2).*(beta2/2)+(w.^3).*(beta3/6)).*L.*1i);
spec_in=fft(S_in,N);
spec_in=fftshift(spec_in);
spec_out=spec_in.*inverse_transferFunct;
spec_out=fftshift(spec_out);
S_out=ifft(spec_out,N);
end
