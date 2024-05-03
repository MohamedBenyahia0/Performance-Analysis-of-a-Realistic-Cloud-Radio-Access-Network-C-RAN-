
T_bit=1e-9;
T_imp=1e-7;
nbits=int32(T_imp/T_bit);
t=linspace(0,T_imp,nbits);
S_in=zeros(1,nbits);
S_in(60:80)=1;

Ts=T_bit;
Fs=1/Ts;
attenuation=10^(0.2/10);
lambda=1550e-9;
delta_lambda=2*lambda^2/(Ts*c);
L=1e3;
D=17e-12;%ps/nm/km
S=0.09e-12;%ps/nm^2/km

c=3e8;
beta2=-D*lambda^2/(2*pi*c);
beta3=S*lambda^4/(4*pi^2*c^2);

S_out=opticalFiber(S_in,Fs,attenuation,beta2,beta3,L);
plot(t,S_in,'r')
hold on;
plot(t,S_out,'g');

