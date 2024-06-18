g1=[1,0,1,0,0,1];%g=x^5+x^2+1 de degre n-k=5
k1=26;%k=31-5
g2=[1,0,0,1,0,1,1,0,1,1,1];%g=x^10+x^9+x^8+x^6+x^5+x^3+1
k2=31-10;
Nb0=50000;
Nb1=k1*100*20;
Nb2=k2*100*5;

%%
[H1,H2,H3] =plot_channel(62);
H0=eye(62,62);
%% ARQ retransmission BPSK BCH
EsN0_dB = 0:1:10;   
M=2;
rate=k2/31;

throughput_BCH = zeros(size(EsN0_dB));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

nb_trames=10;
parfor ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
        bits2 = randi([0,1],2*k2*log2(M),1); % generation des messages d'information pour  BCH2
        mess2=  reshape(bits2,[k2,2*log2(M)]);  
        codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false); % generation des mots codes avec BCH2
        codeWords2=cell2mat(codeWords2);
        codeWords2=reshape(codeWords2,[62*log2(M),1]);
        symboles2=bits2symbols(codeWords2,'PSK',M);
        symboles2=reshape(symboles2,[62,1]);
        
        W2 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
        Z2=H0*symboles2+W2;
        
        Z2=reshape(Z2,[62,1]);
        s_est2 = threshold_detector(Z2, 'PSK', M) ;
        
        d_est2 = symbols2bits(s_est2,'PSK',M); % convertir les symboles en bits
        d_est2=reshape(d_est2,[31,2*log2(M)]);
        m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false);
        m_hat2=cell2mat(m_hat2);
        bits_hat2=reshape(m_hat2,[2*k2*log2(M),1]);
        
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        if ~nb_err_tram(t)
            successful_transmissions=successful_transmissions+1;
        else
            W2 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
            Z2=H0*symboles2+W2;
            Z2=reshape(Z2,[62,1]);
            s_est2 = threshold_detector(Z2, 'PSK', M) ;
            d_est2 = symbols2bits(s_est2,'PSK',M); % convertir les symboles en bits
            d_est2=reshape(d_est2,[31,2*log2(M)]);
            m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false);
            m_hat2=cell2mat(m_hat2);
            bits_hat2=reshape(m_hat2,[2*k2*log2(M),1]);
            nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        
            if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
            end
        end

    end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput_BCH(ii) = rate*log2(M)*successful_transmissions / nb_trames;
end


%% ARQ retransmission Uncoded BPSK 

EsN0_dB = 0:1:10;    

M=2;


throughput = zeros(size(EsN0_dB));
Energy_symbol = 1*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

for ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
       bits0=randi([0,1],62*log2(M),1);
       symboles0= bits2symbols(bits0,'PSK',M);
       symboles0=reshape(symboles0,[62,1]);
       W0 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
       Z0=H0*symboles0+W0;
       Z0=reshape(Z0,[62,1]);
       s_est0 = threshold_detector(Z0,'PSK', M) ;
       d_est0 = symbols2bits(s_est0,'PSK',M); % convertir les symboles en bits
       d_est0=reshape(d_est0,[62*log2(M),1]);
       nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));

       if ~nb_err_tram(t)
        successful_transmissions=successful_transmissions+1;
       else
           W0 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
           Z0=H0*symboles0+W0;
           Z0=reshape(Z0,[62,1]);
           s_est0 = threshold_detector(Z0,'PSK', M) ;
           d_est0 = symbols2bits(s_est0,'PSK',M); % convertir les symboles en bits
           d_est0=reshape(d_est0,[62*log2(M),1]);
           nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));
            
           if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
           end
        end
           
    
     end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput(ii) = log2(M)*successful_transmissions / nb_trames;
end

%% Plot throughput vs Es/N0 BPSK
figure;
plot(EsN0_dB, throughput_BCH, '-o');
hold on;
plot(EsN0_dB, throughput, '-o');
xlabel('E_s/N_0 (dB)');
ylabel('Throughput');
title('Throughput vs E_s/N_0 for ARQ Scheme with 1 Retransmission BPSK ');
legend('Code BCH 2','Uncoded')
grid on;
hold off;
%% ARQ retransmission 8QAM BCH

EsN0_dB = 0:2:20;   
M=8;
rate=k2/31;

throughput_BCH = zeros(size(EsN0_dB));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

nb_trames=10;
parfor ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
        bits2 = randi([0,1],2*k2*log2(M),1); % generation des messages d'information pour  BCH2
        mess2=  reshape(bits2,[k2,2*log2(M)]);  
        codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false); % generation des mots codes avec BCH2
        codeWords2=cell2mat(codeWords2);
        codeWords2=reshape(codeWords2,[62*log2(M),1]);
        symboles2=bits2symbols(codeWords2,'QAM',M);
        symboles2=reshape(symboles2,[62,1]);
        
        W2 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
        Z2=H0*symboles2+W2;
        
        Z2=reshape(Z2,[62,1]);
        s_est2 = threshold_detector(Z2, 'QAM', M) ;
        
        d_est2 = symbols2bits(s_est2,'QAM',M); % convertir les symboles en bits
        d_est2=reshape(d_est2,[31,2*log2(M)]);
        m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false);
        m_hat2=cell2mat(m_hat2);
        bits_hat2=reshape(m_hat2,[2*k2*log2(M),1]);
        
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        if ~nb_err_tram(t)
            successful_transmissions=successful_transmissions+1;
        else
            W2 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
            Z2=H0*symboles2+W2;
            Z2=reshape(Z2,[62,1]);
            s_est2 = threshold_detector(Z2, 'QAM', M) ;
            d_est2 = symbols2bits(s_est2,'QAM',M); % convertir les symboles en bits
            d_est2=reshape(d_est2,[31,2*log2(M)]);
            m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false);
            m_hat2=cell2mat(m_hat2);
            bits_hat2=reshape(m_hat2,[2*k2*log2(M),1]);
            nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        
            if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
            end
        end

    end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput_BCH(ii) = rate*log2(M)*successful_transmissions / nb_trames;
end
%% ARQ retransmission Uncoded 8 QAM

EsN0_dB = 0:2:20;    
max_retransmissions = 1;  
M=8;


throughput = zeros(size(EsN0_dB));
Energy_symbol = 1*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

for ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
       bits0=randi([0,1],62*log2(M),1);
       symboles0= bits2symbols(bits0,'QAM',M);
       symboles0=reshape(symboles0,[62,1]);
       W0 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
       Z0=H0*symboles0+W0;
       Z0=reshape(Z0,[62,1]);
       s_est0 = threshold_detector(Z0,'QAM', M) ;
       d_est0 = symbols2bits(s_est0,'QAM',M); % convertir les symboles en bits
       d_est0=reshape(d_est0,[62*log2(M),1]);
       nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));

       if ~nb_err_tram(t)
        successful_transmissions=successful_transmissions+1;
       else
           W0 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
           Z0=H0*symboles0+W0;
           Z0=reshape(Z0,[62,1]);
           s_est0 = threshold_detector(Z0,'QAM', M) ;
           d_est0 = symbols2bits(s_est0,'QAM',M); % convertir les symboles en bits
           d_est0=reshape(d_est0,[62*log2(M),1]);
           nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));
            
           if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
           end
        end
           
    
     end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput(ii) = log2(M)*successful_transmissions / nb_trames;
end

%% Plot throughput vs Es/N0 8 QAM
figure;
plot(EsN0_dB, throughput_BCH, '-o');
hold on;
plot(EsN0_dB, throughput, '-o');
xlabel('E_s/N_0 (dB)');
ylabel('Throughput');
title('Throughput vs E_s/N_0 for ARQ Scheme with 1 Retransmission 8 QAM ');
legend('Code BCH 2','Uncoded')
grid on;
hold off;
%% ARQ retransmission BCH 16 QAM


EsN0_dB = 0:2:20;   
M=16;
rate=k2/31;

throughput_BCH = zeros(size(EsN0_dB));
Energy_symbol = 1*rate*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

nb_trames=10;
parfor ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
        bits2 = randi([0,1],2*k2*log2(M),1); % generation des messages d'information pour  BCH2
        mess2=  reshape(bits2,[k2,2*log2(M)]);  
        codeWords2=arrayfun(@(i) encoder(transpose(mess2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false); % generation des mots codes avec BCH2
        codeWords2=cell2mat(codeWords2);
        codeWords2=reshape(codeWords2,[62*log2(M),1]);
        symboles2=bits2symbols(codeWords2,'QAM',M);
        symboles2=reshape(symboles2,[62,1]);
        
        W2 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
        Z2=H0*symboles2+W2;
        
        Z2=reshape(Z2,[62,1]);
        s_est2 = threshold_detector(Z2, 'QAM', M) ;
        
        d_est2 = symbols2bits(s_est2,'QAM',M); % convertir les symboles en bits
        d_est2=reshape(d_est2,[31,2*log2(M)]);
        m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false);
        m_hat2=cell2mat(m_hat2);
        bits_hat2=reshape(m_hat2,[2*k2*log2(M),1]);
        
        nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        if ~nb_err_tram(t)
            successful_transmissions=successful_transmissions+1;
        else
            W2 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
            Z2=H0*symboles2+W2;
            Z2=reshape(Z2,[62,1]);
            s_est2 = threshold_detector(Z2, 'QAM', M) ;
            d_est2 = symbols2bits(s_est2,'QAM',M); % convertir les symboles en bits
            d_est2=reshape(d_est2,[31,2*log2(M)]);
            m_hat2=arrayfun(@(i) decoder_2errors(transpose(d_est2(:,i)),g2),[1:2*log2(M)],'UniformOutput',false);
            m_hat2=cell2mat(m_hat2);
            bits_hat2=reshape(m_hat2,[2*k2*log2(M),1]);
            nb_err_tram(t)=nb_err_tram(t)+sum(abs(bits_hat2-bits2));% compter le nombre d'erreur
        
            if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
            end
        end

    end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput_BCH(ii) = rate*log2(M)*successful_transmissions / nb_trames;
end
%% ARQ retransmission Uncoded 16 QAM

EsN0_dB = 0:2:20;    
max_retransmissions = 1;  
M=16;


throughput = zeros(size(EsN0_dB));
Energy_symbol = 1*log2(M);
No = Energy_symbol./10.^(EsN0_dB/10);

parfor ii = 1:length(EsN0_dB)
    nb_err_tram=zeros(nb_trames,1);
    BERtram=zeros(nb_trames,1);
    successful_transmissions = 0;
    for t=1:nb_trames
       bits0=randi([0,1],62*log2(M),1);
       symboles0= bits2symbols(bits0,'QAM',M);
       symboles0=reshape(symboles0,[62,1]);
       W0 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
       Z0=H0*symboles0+W0;
       Z0=reshape(Z0,[62,1]);
       s_est0 = threshold_detector(Z0,'QAM', M) ;
       d_est0 = symbols2bits(s_est0,'QAM',M); % convertir les symboles en bits
       d_est0=reshape(d_est0,[62*log2(M),1]);
       nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));

       if ~nb_err_tram(t)
        successful_transmissions=successful_transmissions+1;
       else
           W0 =sqrt(No(ii)/2)* (randn(62,1) + i*randn(62,1));  % bruit
           Z0=H0*symboles0+W0;
           Z0=reshape(Z0,[62,1]);
           s_est0 = threshold_detector(Z0,'QAM', M) ;
           d_est0 = symbols2bits(s_est0,'QAM',M); % convertir les symboles en bits
           d_est0=reshape(d_est0,[62*log2(M),1]);
           nb_err_tram(t)=nb_err_tram(t)+sum(abs(d_est0-bits0));
            
           if ~nb_err_tram(t)
                successful_transmissions = successful_transmissions + 1;
           end
        end
           
    
     end 
    
    % Calculate throughput as the ratio of successful transmissions to total transmissions
    throughput(ii) = log2(M)*successful_transmissions / nb_trames;
end
%% Plot throughput vs Es/N0 16 QAM
figure;
plot(EsN0_dB, throughput_BCH, '-o');
hold on;
plot(EsN0_dB, throughput, '-o');
xlabel('E_s/N_0 (dB)');
ylabel('Throughput');
title('Throughput vs E_s/N_0 for ARQ Scheme with 1 Retransmission 16 QAM ');
legend('Code BCH 2','Uncoded')
grid on;
hold off;

