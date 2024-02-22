function [SNR_bit,BER1] = BERvsSNR_Code1(N)
SNR_bit=[0:10]; %% Eb/No en dB
BER1 = zeros(1,length(SNR_bit));
Energy_symbol = 1;
No = Energy_symbol./10.^(SNR_bit/10);


nb_err1=zeros(1,length(SNR_bit));
for ii = 1:length(SNR_bit) 
mess1 = randi([0,1],N,k1); % generation des messages d'information pour  BCH1
codeWords1=arrayfun(@(i) encoder(mess1(i,:),g1),[1:N],'UniformOutput',false); % generation des mots codes avec BCH1
codeWords1=cell2mat(codeWords1);
codeWords1=reshape(codeWords1,[N,31]);
symboles1=arrayfun(@(i) bits2symbols(codeWords1(i,:),'PSK',2),[1:N],'UniformOutput',false);
symboles1=cell2mat(symboles1);
symboles1=reshape(symboles1,[N,31]);
W1 =sqrt(No(ii)/2)* (randn(N,31) + i*randn(N,31));  % bruit
Z1=symboles1+W1;
Z1=reshape(Z1,[N*31,1]);
s_est1 = threshold_detector(Z1, 'PSK',2) ; % prise de décision
d_est1 = symbols2bits(s_est1,'PSK',2); % convertir les symboles en bits
d_est1=reshape(d_est1,[N,31]);
m_hat1=arrayfun(@(i) decoder_1error(d_est1(i,:),g1),[1:N],'UniformOutput',false);
m_hat1=cell2mat(m_hat1);
m_hat1=reshape(m_hat1,[N,k1]);

nb_err1(ii)=nb_err1(ii)+sum(abs(m_hat1-mess1),"all");
BER1(ii) = nb_err1(ii)/(N*k1);  % calcul du taux d'erreur binaire
end
end
