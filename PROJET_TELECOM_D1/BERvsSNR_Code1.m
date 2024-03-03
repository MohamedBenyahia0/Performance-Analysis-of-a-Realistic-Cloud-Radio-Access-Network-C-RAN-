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

nb_err1=zeros(1,length(SNR_bit));
for ii = 1:length(SNR_bit) 
bits1 = randi([0,1],N_bits,1); % generation des messages d'information pour  BCH1
mess1=reshape(bits1,[k1,nb_word]);
codeWords1=arrayfun(@(i) encoder(transpose(mess1(:,i)),g1),[1:nb_word],'UniformOutput',false); % generation des mots codes avec BCH1
codeWords1=cell2mat(codeWords1);

codeWords1=reshape(codeWords1,[N_bits+(31-k1)*nb_word,1]);
symboles1=bits2symbols(codeWords1,mod,M);

symboles1=reshape(symboles1,[nb_word,31]);

W1 =sqrt(No(ii)/2)* (randn(nb_word,31) + i*randn(nb_word,31));  % bruit
Z1=H*symboles1+W1;

if strcmp(equalizer,'ZF')
    s_est1 = ZFequalizer(Z1,H, mod,M) ; % prise de décision

else 
    Z1=reshape(Z1,[nb_word*31,1]);
    s_est1 = threshold_detector(Z1, mod, M) ; 
    
end
d_est1 = symbols2bits(s_est1,mod,M); % convertir les symboles en bits
d_est1=reshape(d_est1,[31,nb_word]);
m_hat1=arrayfun(@(i) decoder_1error(transpose(d_est1(:,i)),g1),[1:nb_word],'UniformOutput',false);
m_hat1=cell2mat(m_hat1);
bits_hat1=reshape(m_hat1,[N_bits,1]);

nb_err1(ii)=nb_err1(ii)+sum(abs(bits_hat1-bits1),"all");
BER1(ii) = nb_err1(ii)/(N_bits);  % calcul du taux d'erreur binaire
end
end
