
T_imp=4e-11;
T_bit=1e-13;
nbit=int32(T_imp/T_bit);
t=linspace(0,T_imp,nbit);
%S_in=zeros(1,nbit);
mu=T_imp/2;
sigma=T_imp/60;
gaussian_pulse = exp(-(t - mu).^2 / (2 * sigma^2));
S_in=gaussian_pulse;
%S_in(160:240)=1;
Ts=T_bit;
Fs=1/Ts;
attenuation=0.2*log(10)/10;
lambda=1550e-9;%m
c=3e8;%m
delta_lambda=(2/T_imp)*lambda^2/c;
L=2e3;%m
D=17*10^(-6);%s/m/m
S=0.09e3;%s/m^2/m

beta2=-D*lambda^2/(2*pi*c);
beta3=S*lambda^4/(4*pi^2*c^2);

%delta_tau=(L)*(D)*(delta_lambda);

% Temporal Broadening
Delta_t0_gauss = 2 * sqrt(2 * log(2)) * sigma;
Delta_t_gauss = Delta_t0_gauss * sqrt(1 + (D * L * delta_lambda / Delta_t0_gauss)^2);

disp(["etalement temporel theorique en secondes",Delta_t_gauss]);
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

R=2;
Fs=2.5e9;%10Gbits/s
targetBER=1e-3;
fun = @(P_out_dbm) BER_backToback_OOK(P_out_dbm,Fs,targetBER);

P_out_dbm=[-30:1:-20];
P_mean_dbm=P_out_dbm-3; % Puissance moyenne =Pmax/2
BER_array=zeros(1,length(P_out_dbm));
for i=1:length(BER_array)
    
    BER_array(i)=BER_backToback_OOK(P_out_dbm(i),R,Fs,targetBER);
end
semilogy(P_mean_dbm,BER_array,"-o",'LineWidth',3.0);
xlabel('Power mean (dbm)')
ylabel('BER')
title('Back to Back OOK BER versus power emission'); 
legend('2.5GHZ')
grid on;
%% BER optical system back-to-back OOK direct detection
R=2;
Fs=10e9;%10Gbits/s
targetBER=1e-3;


P_out_dbm=[-30:1:-10];
P_mean_dbm=P_out_dbm-3; % Puissance moyenne =Pmax/2
BER_array=zeros(1,length(P_out_dbm));
for i=1:length(BER_array)
    
    BER_array(i)=BER_backToback_OOK(P_out_dbm(i),R,Fs,targetBER);
end
semilogy(P_mean_dbm,BER_array,"-o",'LineWidth',3.0);
xlabel('Power mean (dbm)')
ylabel('BER')
title('Back to Back OOK BER versus power emission'); 
legend('10GHZ')
grid on;

%% BER optical system fiber OOK direct detection
R=2;
Fs=2.5e9;%10Gbits/s
target_BER=1e-3;
P_out_dbm=-22;

L_array= [0:10:100].*1e3;
BER_array_2_5GHz=zeros(1,length(L_array));
for i=1:length(BER_array_2_5GHz)
    
    BER_array_2_5GHz(i)=BER_Fiber_OOK(P_out_dbm,R,Fs,target_BER,attenuation,beta2,beta3,L_array(i));
end
BER_array_10GHz=zeros(1,length(L_array));
Fs=10e9;
for i=1:length(BER_array_10GHz)
    
    BER_array_10GHz(i)=BER_Fiber_OOK(P_out_dbm,R,Fs,target_BER,attenuation,beta2,beta3,L_array(i));
end
figure( );
semilogy(L_array/1e3,BER_array_2_5GHz,"-o",'LineWidth',3.0);
hold on;
semilogy(L_array/1e3,BER_array_10GHz,"-o",'LineWidth',3.0);

xlabel('fiber length (km)')
ylabel('BER')
title('Fiber propagation BER versus fiber length'); 
legend('2.5GHZ','10GHz')
grid on;
hold off;

%% Fiber DCF
D_DCF=-80*10^(-6);%s/m/m
attenuation_DCF=0.5*log(10)/10;




