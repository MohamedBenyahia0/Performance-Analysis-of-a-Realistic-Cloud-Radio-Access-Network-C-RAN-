function [SNR_bit,BER2] = BERvsSNR_Code2(N)

SNR_bit=[0:10]; %% Eb/No en dB

BER2= zeros(1,length(SNR_bit));
Energy_symbol = 1;
No = Energy_symbol./10.^(SNR_bit/10);

nb_err2=zeros(1,length(SNR_bit));
for ii = 1:length(SNR_bit) 

mess2=  randi([0,1],N,k2);  % generation des messages d'information pour  BCH2
codeWords2=arrayfun(@(i) encoder(mess2(i,:),g2),[1:N],'UniformOutput',false); % generation des mots codes avec BCH2
codeWords2=cell2mat(codeWords2);
codeWords2=reshape(codeWords2,[N,31]);

%% %%%%%%%%% BPSK modulation

symboles2=arrayfun(@(i) bits2symbols(codeWords2(i,:),'PSK',2),[1:N],'UniformOutput',false);
symboles2=cell2mat(symboles2);
symboles2=reshape(symboles2,[N,31]);
%% Application du canal AWGN

W2 =sqrt(No(ii)/2)* (randn(N,31) + i*randn(N,31));  % bruit
Z2=symboles2+W2;
Z2=reshape(Z2,[N*31,1]);
%% %%%%%%%%% demodulation

s_est2 = threshold_detector(Z2, 'PSK',2) ; % prise de décision
d_est2 = symbols2bits(s_est2,'PSK',2); % convertir les symboles en bits
d_est2=reshape(d_est2,[N,31]);
m_hat2=arrayfun(@(i) decoder_2errors(d_est2(i,:),g2),[1:N],'UniformOutput',false);
m_hat2=cell2mat(m_hat2);
m_hat2=reshape(m_hat2,[N,k2]);

%%   %%%%%%%   calcul du BER

nb_err2(ii)=nb_err2(ii)+sum(abs(m_hat2-mess2),"all");

BER2(ii) = nb_err2(ii)/(N*k2);  % calcul du taux d'erreur binaire
end

end