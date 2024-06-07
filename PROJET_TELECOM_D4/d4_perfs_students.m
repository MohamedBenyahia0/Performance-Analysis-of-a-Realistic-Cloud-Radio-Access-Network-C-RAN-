%% [rej] = d4_perfs(mode,rx)


%% Rejection rate for deliverable D4
%%
%% mode = selected data rate for all the users ('40k','400k','4M','40M')
%% rx = selected receiver for all the users ('zf','dfe')
  
%% rej = rejection rate versus K (number of users) 

%% Location : Telecom Paris
%% Author : Philippe Ciblat <ciblat@telecom-paris.fr>
%% Date   : 09/06/2023



function [K,rej] = d4_perfs_students(mode,rx)

Ts=1/(20e6);
MC=1000; %%number of Monte-Carlo simulations
FER=1e-2;
N=100;
M=2;
BER=1-(1-FER)^(1/(N*log2(M)));


if(strcmp(mode,'40k')==1) 
R=40*10^(3); 
Kmax=2000; 
end;

if(strcmp(mode,'400k')==1) 
R=400*10^(3); 
Kmax=200; 
end;

if(strcmp(mode,'4M')==1) 
R=4*10^(6); 
Kmax=20;
end;

if(strcmp(mode,'40M')==1) 
R=40*10^(6);
Kmax=2; 
end;


if(strcmp(rx,'zf')==1)
SNR_min_bpsk=[8.75,13.5,22]; %minimum Es/N0 for each channel (3-length vector);
SNR_min_8qam=[12,15.1,24]+10*log10(3); %minimum Es/N0 for each channel (3-length vector);
SNR_min_16qam=[12.145,18,24]+10*log10(4); %minimum Es/N0 for each channel (3-length vector);
end;

if(strcmp(rx,'dfe')==1)
SNR_min_bpsk=[8.5,11,14]; %minimum Es/N0 for each channel (3-length vector);
SNR_min_8qam=[11,15,17]+10*log10(3); %minimum Es/N0 for each channel (3-length vector);
SNR_min_16qam=[12,145,16,19]+10*log10(4);%minimum Es/N0 for each channel (3-length vector);
end;
%endif;
  

K=[1:(floor(Kmax/10)+1):Kmax];
rej=zeros(1,length(K));

for kk=1:length(K)
aux=0;

for mm=1:MC
xx= -1 + 2*rand(1, K(kk)); % uniformly-distributed x-axis in a square of semi-length 1 (K(kk) length vector);
yy= -1 + 2*rand(1, K(kk)); % uniformly-distributed y-axis in a square of semi-length 1 (K(kk) length vector);

d2=xx.^2+yy.^2;%% square-distance from the origin 
c=3e8;
f0=2.4e9;
a2=min(1, (c/(4*pi*f0))^2./(d2*1e6));% square magnitude attenuation in Friis equation (calculate ??) 
a2_dB=10*log10(a2);
P_max=20;%dbm
N0=-174;%dbm/Hz
B=30e6;%Hz
N_lin=(1e-3)*10^(N0/10)*B; 
N_dbm=10*log10(N_lin*1e3);
SNR_rx_max=a2_dB+ P_max-N_dbm; %calculate ??? : SNRmax at TX in dB

SNR_min=SNR_min_16qam; %%ZF or DFE/16QAM: channel 1, channel 2, channel 3
if(R*Ts*K(kk)<3)  SNR_min=SNR_min_8qam; end%%ZF or DFE/8QAM: channel 1, channel 2, channel 3
if(R*Ts*K(kk)<1) SNR_min=SNR_min_bpsk; end %%ZF or DFE/BPSK: channel 1, channel 2, channel 3
  
SNR_rx_min=SNR_min(randi(3,K(kk),1));

aux=aux+length((find(sign(SNR_rx_min-SNR_rx_max)+1)/2));

end

rej(kk)= aux/(K(kk)*MC);


end;



