function [SNR_bit,BER2] = BERvsSNR_Code2(N,k2,g2,H,equalizer)

SNR_bit=[0:10]; %% Eb/No en dB
M=2;
rate=k2/31;
SNR=SNR_bit+10*log10(log2(M )*rate);
BER2= zeros(1,length(SNR_bit));
Energy_symbol = 1;
No = Energy_symbol./10.^(SNR/10);
nb_word=floor(N/k2);
nb_err2=zeros(1,length(SNR_bit));
for ii = 1:length(SNR_bit) 
bits2 = randi([0,1],N,1); % generation des messages d'information pour  BCH2
mess2=  reshape(bits2,[k2,nb_word]);  
codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:nb_word],'UniformOutput',false); % generation des mots codes avec BCH2
codeWords2=cell2mat(codeWords2);
codeWords2=reshape(codeWords2,[N+(31-k2)*nb_word,1]);

%% %%%%%%%%% BPSK modulation

symboles2=bits2symbols(codeWords2,'PSK',2);


%% Application du canal AWGN

W2 =sqrt(No(ii)/2)* (randn(N+(31-k2)*nb_word,1) + i*randn(N+(31-k2)*nb_word,1));  % bruit
Z2=H*symboles2+W2;
Z2=reshape(Z2,[N+(31-k2)*nb_word,1,1]);
%% %%%%%%%%% demodulation

if strcmp(equalizer,'ZF')
    s_est2 = ZFequalizer(Z2,H, 'PSK',2) ; % prise de décision

else %DF
    s_est2 = DFequalizer(Z2,H, 'PSK',2) ; % prise de décision
end
d_est2 = symbols2bits(s_est2,'PSK',2); % convertir les symboles en bits
d_est2=reshape(d_est2,[31,nb_word]);
m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:nb_word],'UniformOutput',false);
m_hat2=cell2mat(m_hat2);
bits_hat2=reshape(m_hat2,[N,1]);

%%   %%%%%%%   calcul du BER

nb_err2(ii)=nb_err2(ii)+sum(abs(bits_hat2-bits2),"all");

BER2(ii) = nb_err2(ii)/(N);  % calcul du taux d'erreur binaire
end

end