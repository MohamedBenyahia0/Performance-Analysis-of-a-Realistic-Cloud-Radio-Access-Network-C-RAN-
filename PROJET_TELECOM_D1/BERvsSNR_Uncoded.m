function [SNR_bit,BER0] = BERvsSNR_Uncoded(N)
SNR_bit=[0:10]; %% Eb/No en dB

BER0 = zeros(1,length(SNR_bit));
Energy_symbol = 1;
No = Energy_symbol./10.^(SNR_bit/10);

nb_err0=zeros(1,length(SNR_bit));
for ii = 1:length(SNR_bit) 
bits0=randi([0,1],N,1); % generation des bits d'informations 
%% %%%%%%%%% BPSK modulation
symboles0= bits2symbols(bits0,'PSK',2);
symboles0=reshape(symboles0,[N,1]);
%% Application du canal AWGN
W0 =sqrt(No(ii)/2)* (randn(N,1) + i*randn(N,1));  % bruit
Z0=symboles0+W0;
Z0=reshape(Z0,[N,1]);
%% %%%%%%%%% demodulation
s_est0 = threshold_detector(Z0, 'PSK',2) ; % prise de décision
d_est0 = symbols2bits(s_est0,'PSK',2); % convertir les symboles en bits
d_est0=reshape(d_est0,[N,1]);
%%   %%%%%%%   calcul du BER
nb_err0(ii)=nb_err0(ii)+sum(abs(d_est0-bits0));% compter le nombre d'erreur
BER0(ii) = nb_err0(ii)/N;  % calcul du taux d'erreur binaire
end
end
