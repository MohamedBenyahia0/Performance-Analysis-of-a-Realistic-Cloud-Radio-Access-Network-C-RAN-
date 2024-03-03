function [SNR_bit,BER0] = BERvsSNR_Uncoded(N,H,equalizer,mod,M,SNRbit_min,SNRbit_max,step)
SNR_bit=[SNRbit_min:step:SNRbit_max]; %% Eb/No en dB
SNR=SNR_bit+10*log10(log2(M ));
BER0 = zeros(1,length(SNR_bit));
Energy_symbol = 1*log2(M);
No = Energy_symbol./10.^(SNR/10);

nb_trames=floor(N/100);
N_bits_trame = 100*log2(M); % nombre de bits dans une trame
for ii = 1:length(SNR_bit) 
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    for t=1:nb_trames
        bits0=randi([0,1],100*log2(M),1); % generation des bits d'informations 
        %% %%%%%%%%% BPSK modulation
        symboles0= bits2symbols(bits0,mod,M);
        
        symboles0=reshape(symboles0,[100,1]);
        %% Application du canal AWGN
        W0 =sqrt(No(ii)/2)* (randn(100,1) + i*randn(100,1));  % bruit
        Z0=H*symboles0+W0;
        
        %% %%%%%%%%% demodulation
        if strcmp(equalizer,'ZF')
            s_est0 = ZFequalizer(Z0,H, mod,M) ; % prise de décision
        
        else 
            if strcmp(equalizer,'DF')
        
                s_est0 = DFequalizer(Z0,H,mod,M) ; % prise de décision
            else
                Z0=reshape(Z0,[100,1]);
                s_est0 = threshold_detector(Z0,mod, M) ; 
            end
        end
        
        d_est0 = symbols2bits(s_est0,mod,M); % convertir les symboles en bits
        d_est0=reshape(d_est0,[100*log2(M),1]);
        %%   %%%%%%%   calcul du BER
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));% compter le nombre d'erreur
        BERtram(t) = nb_err_tram(t)/N_bits_trame;  % calcul du taux d'erreur binaire
    end
    BER0(ii)=sum(BERtram)/nb_trames;
end
end
