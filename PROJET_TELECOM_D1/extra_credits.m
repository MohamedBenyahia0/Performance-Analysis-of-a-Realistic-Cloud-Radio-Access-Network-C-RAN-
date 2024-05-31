g1=[1,0,1,0,0,1];%g=x^5+x^2+1 de degre n-k=5
k1=26;%k=31-5
g2=[1,0,0,1,0,1,1,0,1,1,1];%g=x^10+x^9+x^8+x^6+x^5+x^3+1
k2=31-10;
Nb0=7000;
Nb1=k1*100*2;
Nb2=k2*100*3;
H0=eye(100,100);

[H1,H2,H3] =plot_channel();
%%
[SNR_bit,BER0] = BERvsSNR_Uncoded(Nb0,H1,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER1] = BERvsSNR_Code1(Nb1,k1,g1,H1,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER2] = BERvsSNR_Code2(Nb2,k2,g2,H1,"threshold detector",'PSK',2,0,10,1);

figure();


semilogy(SNR_bit,BER0,"-o",'LineWidth',3.0);
hold on
semilogy(SNR_bit,BER1,"-*",'LineWidth',3.0);
semilogy(SNR_bit,BER2,"-x",'LineWidth',3.0)

xlabel('Eb/No (in dB)')
ylabel('BER')
title('Channel 1 Empirical BER vs Eb/N0'); 
legend('Uncoded','Code BCH 1',' Code BCH 2')
grid on
hold off
%%
[SNR_bit,BER0] = BERvsSNR_Uncoded(Nb0,H2,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER1] = BERvsSNR_Code1(Nb1,k1,g1,H2,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER2] = BERvsSNR_Code2(Nb2,k2,g2,H2,"threshold detector",'PSK',2,0,10,1);

figure();


semilogy(SNR_bit,BER1,"-o",'LineWidth',3.0);
hold on

semilogy(SNR_bit,BER2,"-x",'LineWidth',3.0)

xlabel('Eb/No (in dB)')
ylabel('BER')
title('Channel 2 Empirical BER vs Eb/N0'); 
legend('Uncoded','Code BCH 1','Code BCH 2')
grid on
hold off
%%
[SNR_bit,BER0] = BERvsSNR_Uncoded(Nb0,H3,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER1] = BERvsSNR_Code1(Nb1,k1,g1,H3,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER2] = BERvsSNR_Code2(Nb2,k2,g2,H3,"threshold detector",'PSK',2,0,10,1);

figure();


semilogy(SNR_bit,BER1,"-o",'LineWidth',3.0);
hold on

semilogy(SNR_bit,BER2,"-x",'LineWidth',3.0)

xlabel('Eb/No (in dB)')
ylabel('BER')
title('Channel 3 Empirical BER vs Eb/N0'); 
legend('Uncoded','Code BCH 1','Code BCH 2')
grid on
hold off
%% ARQ retransmission BPSK BCH

EsN0_dB = 0:1:10;   
M=2;
rate=k2/31;
N=k2*100*2;
throughput = zeros(size(EsN0_dB));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);
nb_word=floor(N/k2);
nb_trames=floor(nb_word/100);
for ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
        bits2 = randi([0,1],100*k2*log2(M),1); % generation des messages d'information pour  BCH2
        mess2=  reshape(bits2,[k2,100*log2(M)]);  
        codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:100*log2(M)],'UniformOutput',false); % generation des mots codes avec BCH2
        codeWords2=cell2mat(codeWords2);
        codeWords2=reshape(codeWords2,[31*100*log2(M),1]);
        symboles2=bits2symbols(codeWords2,'PSK',M);
        symboles2=reshape(symboles2,[100,31]);
        
        W2 =sqrt(No(ii)/2)* (randn(100,31) + i*randn(100,31));  % bruit
        Z2=H0*symboles2+W2;
        
        Z2=reshape(Z2,[100*31,1]);
        s_est2 = threshold_detector(Z2, 'PSK', M) ;
        
        d_est2 = symbols2bits(s_est2,'PSK',M); % convertir les symboles en bits
        d_est2=reshape(d_est2,[31,100*log2(M)]);
        m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:100*log2(M)],'UniformOutput',false);
        m_hat2=cell2mat(m_hat2);
        bits_hat2=reshape(m_hat2,[100*k2*log2(M),1]);
        
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        if ~nb_err_tram(t)
            successful_transmissions=successful_transmissions+1;
        else
            W2 =sqrt(No(ii)/2)* (randn(100,31) + i*randn(100,31));  % bruit
            Z2=H0*symboles2+W2;
            Z2=reshape(Z2,[100*31,1]);
            s_est2 = threshold_detector(Z2, 'PSK', M) ;
            d_est2 = symbols2bits(s_est2,'PSK',M); % convertir les symboles en bits
            d_est2=reshape(d_est2,[31,100*log2(M)]);
            m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:100*log2(M)],'UniformOutput',false);
            m_hat2=cell2mat(m_hat2);
            bits_hat2=reshape(m_hat2,[100*k2*log2(M),1]);
            nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        
            if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
            end
        end

    end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput(ii) = successful_transmissions / nb_trames;
end
% Plot throughput vs Es/N0
figure;
plot(EsN0_dB, throughput, '-o');
xlabel('E_s/N_0 (dB)');
ylabel('Throughput');
title('Throughput vs E_s/N_0 for ARQ Scheme with 1 Retransmission BPSK BCH');
grid on;

%% ARQ retransmission BPSK 

EsN0_dB = 0:1:10;    
max_retransmissions = 1;  
M=2;
rate=k2/31;
N=k2*100*2;
throughput = zeros(size(EsN0_dB));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

nb_word=floor(N/k2);

nb_trames=floor(nb_word/100);
for ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
       bits0=randi([0,1],100*log2(M),1);
       symboles0= bits2symbols(bits0,'PSK',M);
       symboles0=reshape(symboles0,[100,1]);
       W0 =sqrt(No(ii)/2)* (randn(100,1) + i*randn(100,1));  % bruit
       Z0=H0*symboles0+W0;
       Z0=reshape(Z0,[100,1]);
       s_est0 = threshold_detector(Z0,'PSK', M) ;
       d_est0 = symbols2bits(s_est0,'PSK',M); % convertir les symboles en bits
       d_est0=reshape(d_est0,[100*log2(M),1]);
       nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));

       if ~nb_err_tram(t)
        successful_transmissions=successful_transmissions+1;
       else
           W0 =sqrt(No(ii)/2)* (randn(100,1) + i*randn(100,1));  % bruit
           Z0=H0*symboles0+W0;
           Z0=reshape(Z0,[100,1]);
           s_est0 = threshold_detector(Z0,'PSK', M) ;
           d_est0 = symbols2bits(s_est0,'PSK',M); % convertir les symboles en bits
           d_est0=reshape(d_est0,[100*log2(M),1]);
           nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));
            
           if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
           end
        end
           
    
     end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput(ii) = successful_transmissions / nb_trames;
end

% Plot throughput vs Es/N0
figure;
plot(EsN0_dB, throughput, '-o');
xlabel('E_s/N_0 (dB)');
ylabel('Throughput');
title('Throughput vs E_s/N_0 for ARQ Scheme with 1 Retransmission BPSK Uncoded');
grid on;
%% ARQ retransmission 8QAM BCH

EsN0_dB = 0:1:10;   
M=8;
rate=k2/31;
N=k2*100*2;
throughput = zeros(size(EsN0_dB));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);
nb_word=floor(N/k2);
nb_trames=floor(nb_word/100);
for ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
        bits2 = randi([0,1],100*k2*log2(M),1); % generation des messages d'information pour  BCH2
        mess2=  reshape(bits2,[k2,100*log2(M)]);  
        codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:100*log2(M)],'UniformOutput',false); % generation des mots codes avec BCH2
        codeWords2=cell2mat(codeWords2);
        codeWords2=reshape(codeWords2,[31*100*log2(M),1]);
        symboles2=bits2symbols(codeWords2,'QAM',M);
        symboles2=reshape(symboles2,[100,31]);
        
        W2 =sqrt(No(ii)/2)* (randn(100,31) + i*randn(100,31));  % bruit
        Z2=H0*symboles2+W2;
        
        Z2=reshape(Z2,[100*31,1]);
        s_est2 = threshold_detector(Z2, 'QAM', M) ;
        
        d_est2 = symbols2bits(s_est2,'QAM',M); % convertir les symboles en bits
        d_est2=reshape(d_est2,[31,100*log2(M)]);
        m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:100*log2(M)],'UniformOutput',false);
        m_hat2=cell2mat(m_hat2);
        bits_hat2=reshape(m_hat2,[100*k2*log2(M),1]);
        
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        if ~nb_err_tram(t)
            successful_transmissions=successful_transmissions+1;
        else
            W2 =sqrt(No(ii)/2)* (randn(100,31) + i*randn(100,31));  % bruit
            Z2=H0*symboles2+W2;
            Z2=reshape(Z2,[100*31,1]);
            s_est2 = threshold_detector(Z2, 'QAM', M) ;
            d_est2 = symbols2bits(s_est2,'QAM',M); % convertir les symboles en bits
            d_est2=reshape(d_est2,[31,100*log2(M)]);
            m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:100*log2(M)],'UniformOutput',false);
            m_hat2=cell2mat(m_hat2);
            bits_hat2=reshape(m_hat2,[100*k2*log2(M),1]);
            nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        
            if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
            end
        end

    end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput(ii) = successful_transmissions / nb_trames;
end
% Plot throughput vs Es/N0
figure;
plot(EsN0_dB, throughput, '-o');
xlabel('E_s/N_0 (dB)');
ylabel('Throughput');
title('Throughput vs E_s/N_0 for ARQ Scheme with 1 Retransmission 8QAM BCH');
grid on;



