function [SNR_bit,BER1] = BERvsSNR_Code1(N,k1,g1,H,equalizer)
SNR_bit=[0:10]; %% Eb/No en dB
M=2;
mod='PSK';
rate=k1/31;
SNR=SNR_bit+10*log10(log2(M )*rate);
BER1 = zeros(1,length(SNR ));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(SNR/10);
N_bits=N;
nb_word=floor(N_bits/k1);
nb_trames=floor(nb_word/100);

for ii = 1:length(SNR_bit) 
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    for t=1:nb_trames
        bits1 = randi([0,1],100*k1,1); % generation des messages d'information pour  BCH1
        mess1=reshape(bits1,[k1,100]);
        codeWords1=arrayfun(@(i) encoder(transpose(mess1(:,i)),g1),[1:100],'UniformOutput',false); % generation des mots codes avec BCH1
        codeWords1=cell2mat(codeWords1);
        
        codeWords1=reshape(codeWords1,[31*100,1]);
        symboles1=bits2symbols(codeWords1,mod,M);
        
        symboles1=reshape(symboles1,[100,31]);
        
        W1 =sqrt(No(ii)/2)* (randn(100,31) + i*randn(100,31));  % bruit
        Z1=H*symboles1+W1;
        
        if strcmp(equalizer,'ZF')
            s_est1 = ZFequalizer(Z1,H, mod,M) ; % prise de décision
        
        else 
            Z1=reshape(Z1,[100*31,1]);
            s_est1 = threshold_detector(Z1, mod, M) ; 
            
        end
        d_est1 = symbols2bits(s_est1,mod,M); % convertir les symboles en bits
        d_est1=reshape(d_est1,[31,100]);
        m_hat1=arrayfun(@(i) decoder_1error(transpose(d_est1(:,i)),g1),[1:100],'UniformOutput',false);
        m_hat1=cell2mat(m_hat1);
        bits_hat1=reshape(m_hat1,[100*k1,1]);
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat1-bits1));% compter le nombre d'erreur
        BERtram(t) = nb_err_tram(t)/(100*k1);  % calcul du taux d'erreur binaire
    end
    
    BER1(ii)=sum(BERtram)/nb_trames;


end
end
