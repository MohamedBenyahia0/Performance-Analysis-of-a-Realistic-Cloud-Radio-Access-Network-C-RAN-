function [SNR_bit,BER0] = BERvsSNR_Uncoded(N,H,equalizer)
SNR_bit=[0:10]; %% Eb/No en dB

BER0 = zeros(1,length(SNR_bit));
Energy_symbol = 1;
No = Energy_symbol./10.^(SNR_bit/10);
nb_trames=floor(N/100);

for ii = 1:length(SNR_bit) 
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    for t=1:nb_trames
        bits0=randi([0,1],100,1); % generation des bits d'informations 
        %% %%%%%%%%% BPSK modulation
        symboles0= bits2symbols(bits0,'PSK',2);
        symboles0=reshape(symboles0,[100,1]);
        %% Application du canal AWGN
        W0 =sqrt(No(ii)/2)* (randn(100,1) + i*randn(100,1));  % bruit
        Z0=H*symboles0+W0;
        
        %% %%%%%%%%% demodulation
        if strcmp(equalizer,'ZF')
            s_est0 = ZFequalizer(Z0,H, 'PSK',2) ; % prise de décision
        
        else %DF
            s_est0 = DFequalizer(Z0,H, 'PSK',2) ; % prise de décision
        end
        
        d_est0 = symbols2bits(s_est0,'PSK',2); % convertir les symboles en bits
        d_est0=reshape(d_est0,[100,1]);
        %%   %%%%%%%   calcul du BER
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));% compter le nombre d'erreur
        BERtram(t) = nb_err_tram(t)/100;  % calcul du taux d'erreur binaire
    end
    BER0(ii)=mean(BERtram);
end
end
