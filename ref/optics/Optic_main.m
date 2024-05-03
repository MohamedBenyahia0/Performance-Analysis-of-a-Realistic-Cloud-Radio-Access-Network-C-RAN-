T_imp=1e-7;
T_bit=1e-9;
nbit=int32(T_imp/T_bit);
t=linspace(0,T_imp,nbit);
S_in=zeros(1,nbit);
S_in(60:80)=1;
Ts=T_bit;
Fs=1/Ts;
attenuation=10^(0.2/10);
lambda=1550e-9;
c=3e8;
delta_lambda=(2/Ts)*lambda^2/c;
L=1e3;
D=17e-12;%ps/nm/km
S=0.09e-12;%ps/nm^2/km


beta2=-D*lambda^2/(2*pi*c);
beta3=S*lambda^4/(4*pi^2*c^2);

S_out=opticalFiber(S_in,Fs,attenuation,beta2,beta3,L);
plot(t,S_in);
hold on;
plot(t,S_out);
hold off;
