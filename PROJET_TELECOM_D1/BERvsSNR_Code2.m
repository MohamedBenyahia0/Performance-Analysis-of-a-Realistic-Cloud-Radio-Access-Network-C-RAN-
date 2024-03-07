function [SNR_bit,BER2] = BERvsSNR_Code2(N,k2,g2,H,equalizer)

SNR_bit=[0:10]; %% Eb/No en dB
M=2;
mod='PSK';
rate=k2/31;
SNR=SNR_bit+10*log10(log2(M )*rate);
BER2= zeros(1,length(SNR_bit));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(SNR/10);
N_bits=N*log2(M);
nb_word=floor(N_bits/k2);

nb_trames=floor(nb_word/100);
nb_err2=zeros(1,length(SNR_bit));
for ii = 1:length(SNR_bit) 
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    for t=1:nb_trames
        bits2 = randi([0,1],100*k2,1); % generation des messages d'information pour  BCH2
        mess2=  reshape(bits2,[k2,100]);  
        codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:100],'UniformOutput',false); % generation des mots codes avec BCH2
        codeWords2=cell2mat(codeWords2);
        codeWords2=reshape(codeWords2,[31*100,1]);
        
        %% %%%%%%%%% BPSK modulation
        
        symboles2=bits2symbols(codeWords2,mod,M);
        symboles2=reshape(symboles2,[100,31]);
        
        %% Application du canal AWGN
        
        W2 =sqrt(No(ii)/2)* (randn(100,31) + i*randn(100,31));  % bruit
        Z2=H*symboles2+W2;
        
        %% %%%%%%%%% demodulation
        
        if strcmp(equalizer,'ZF')
            s_est2 = ZFequalizer(Z2,H, mod,M) ; % prise de décision
        else 
            Z2=reshape(Z2,[100*31,1]);
            s_est2 = threshold_detector(Z2, mod, M) ; 
            
        end
        d_est2 = symbols2bits(s_est2,mod,M); % convertir les symboles en bits
        d_est2=reshape(d_est2,[31,100]);
        m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:100],'UniformOutput',false);
        m_hat2=cell2mat(m_hat2);
        bits_hat2=reshape(m_hat2,[100*k2,1]);
        
        %%   %%%%%%%   calcul du BER
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        BERtram(t) = nb_err_tram(t)/(100*k2);  % calcul du taux d'erreur binaire
    end

BER2(ii)=sum(BERtram)/nb_trames;  % calcul du taux d'erreur binaire

end

end