clear;
close all;
T=3e-10;
T_bit=1e-12;
nbit=int32(T/T_bit);
t=linspace(0,T,nbit);
S_in=zeros(1,nbit);
n_ones=25;
T_imp=n_ones*T_bit;
S_in(nbit/2-n_ones/2:nbit/2+n_ones/2)=1;

Ts=T_bit;
Fs=1/Ts;
attenuation=0.2*log(10)/10;
lambda=1550e-9;%m
c=3e8;%m

L=10e3;%m
D=17*10^(-6);%s/m/m
S=0.09e3;%s/m^2/m
%delta_lambda=(0.44/sigma)*lambda^2/c;
delta_lambda=(2/T_imp)*lambda^2/c;
delta_tau=L*D*delta_lambda;
beta2=-D*lambda^2/(2*pi*c);
beta3=S*lambda^4/(4*pi^2*c^2);

disp(["etalement temporel theorique en secondes",delta_tau]);
S_out=opticalFiber(S_in,Fs,attenuation,beta2,beta3,L);
plot(t,S_in);
hold on;
plot(t,abs(S_out));

xlabel('Time(s)')
ylabel('Amplitude')
title('Temporal Broadening due to Fiber Propagation'); 
legend('In','Out')
hold off;
%% BER optical system back-to-back OOK direct detection

R=1;
Fs=2.5e9;%10Gbits/s
targetBER=1e-3;
fun = @(P_out_dbm) BER_backToback_OOK(P_out_dbm,Fs,targetBER);

P_out_dbm=[-30:1:-10];
P_mean_dbm=P_out_dbm+3; % Puissance moyenne =Pmax/2
BER_array_2_5=zeros(1,length(P_out_dbm));
for i=1:length(BER_array_2_5)
    
    BER_array_2_5(i)=BER_backToback_OOK(P_out_dbm(i),R,Fs,targetBER);
end

%% BER optical system back-to-back OOK direct detection
R=1;
Fs=10e9;%10Gbits/s
targetBER=1e-3;


P_out_dbm=[-30:1:-10];
P_mean_dbm=P_out_dbm+3; % Puissance moyenne =Pmax/2
BER_array_10=zeros(1,length(P_out_dbm));
for i=1:length(BER_array_10)
    
    BER_array_10(i)=BER_backToback_OOK(P_out_dbm(i),R,Fs,targetBER);
end
figure( );
semilogy(P_mean_dbm,BER_array_2_5,"-o",'LineWidth',3.0);
hold on;
semilogy(P_mean_dbm,BER_array_10,"-o",'LineWidth',3.0);

xlabel('mean power laser')
ylabel('BER')
title('Back To back OOK BER versus mean power laser'); 
legend('2.5GHZ','10GHz')
grid on;
hold off;

%% BER optical system fiber OOK direct detection
R=1;
Fs=2.5e9;%10Gbits/s
target_BER=1e-3;
power_eml=-22;

L_array= [0:5:100].*1e3;
BER_array_2_5Gbps=zeros(1,length(L_array));
for i=1:length(BER_array_2_5Gbps)
    
    BER_array_2_5Gbps(i)=BER_Fiber_OOK(power_eml,R,Fs,attenuation,beta2,beta3,L_array(i));
end
BER_array_10Gbps=zeros(1,length(L_array));
Fs=10e9;
power_eml=-19;
for i=1:length(BER_array_10Gbps)
    
    BER_array_10Gbps(i)=BER_Fiber_OOK(power_eml,R,Fs,attenuation,beta2,beta3,L_array(i));
end
figure( );
semilogy(L_array/1e3,BER_array_2_5Gbps,"-o",'LineWidth',3.0);
hold on;
semilogy(L_array/1e3,BER_array_10Gbps,"-o",'LineWidth',3.0);

xlabel('fiber length (km)')
ylabel('BER')
title('Fiber propagation BER versus fiber length'); 
legend('2.5GHZ','10GHz')
grid on;
hold off;

%% Fiber DCF
D_DCF=-80*10^(-6);%s/m/m
lambda=1550e-9;
attenuation_DCF=0.5*log(10)/10;
beta2_DCF=-D_DCF*lambda^2/(2*pi*c);
beta3_DCF=S*lambda^4/(4*pi^2*c^2);

R=1;
target_BER=1e-3;
P_out_dbm=-19;

L_array= [0:5:100].*1e3;
L_DCF=L_array.*(-D/D_DCF);
BER_arrayDCF_10Gbps=zeros(1,length(L_array));
Fs=10e9;
for i=1:length(BER_arrayDCF_10Gbps)
    
    BER_arrayDCF_10Gbps(i)=BER_Fiber_WithCompensation_OOK(P_out_dbm,R,Fs,target_BER,attenuation,beta2,beta3,L_array(i),beta2_DCF,beta3_DCF,attenuation_DCF,L_DCF(i));
end
figure( );

semilogy(L_array/1e3,BER_arrayDCF_10Gbps,"-o",'LineWidth',3.0);
hold on;
semilogy(L_array/1e3,BER_array_10Gbps,"-o",'LineWidth',3.0);
xlabel('fiber length (km)')
ylabel('BER')
title('Fiber propagation BER at 10Gbps versus fiber length'); 
legend('with DCF','without DCF')
grid on;
hold off;

%% Zero dispersion 
D0=0;
attenuation_0=0.35*log(10)/10;
lambda0=1300e-9;
beta2_0=0;
beta3_0=S*lambda0^4/(4*pi^2*c^2);
L_array= [0:5:100].*1e3;
BER_array0_10Gbps=zeros(1,length(L_array));
Fs=10e9;
power_eml=-19;
R=1;
for i=1:length(BER_array0_10Gbps)
    
    BER_array0_10Gbps(i)=BER_Fiber_OOK(power_eml,R,Fs,attenuation_0,beta2_0,beta3_0,L_array(i));
end
figure( );
semilogy(L_array/1e3,BER_array0_10Gbps,"-o",'LineWidth',3.0);
hold on;
semilogy(L_array/1e3,BER_array_10Gbps,"-o",'LineWidth',3.0);

xlabel('fiber length (km)')
ylabel('BER')
title('Fiber propagation BER at 10Gbps versus fiber length'); 
legend('At 1300nm','At 1550nm')
grid on;
hold off;



