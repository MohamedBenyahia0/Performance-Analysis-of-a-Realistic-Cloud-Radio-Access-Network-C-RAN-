g1=[1,0,1,0,0,1];%g=x^5+x^2+1 de degre n-k=5
k1=26;%k=31-5
rate1=k1/31;
m1=ones(1,k1);
c1=encoder(m1,g1);
dmin1=3;
c1_prim=c1;
c1_prim(4)=c1_prim(4)+1;
m1_hat=decoder_1error(c1_prim,g1);

%% %%%%%%%%%%%
g2=[1,0,0,1,0,1,1,0,1,1,1];%g=x^10+x^9+x^8+x^6+x^5+x^3+1
k2=31-10;
rate2=k2/31;
m2=ones(1,k2);
c2=encoder(m2,g2);
dmin2=5;
c2_prim=c2;
c2_prim(3)=c2_prim(3)+1;
c2_prim(6)=c2_prim(6)+1;
m2_hat=decoder_2errors(c2_prim,g2);
%% %%%%%
[H1,H2,H3] =plot_channel();

H0=eye(100,100);
%% %%%%%%%%%% plot BER versus Eb/No for uncoded and both BCH codes (BPSK) with AWGN channel%%%%%%%%
Nb0=3000;
Nb1=k1*100;
Nb2=k2*100;
[SNR_bit,BER0] = BERvsSNR_Uncoded(Nb0,H0,"threshold detector",'PSK',2,0,10,1);
[SNR_bit,BER1] = BERvsSNR_Code1(Nb1,k1,g1,H0,"threshold detector");
[SNR_bit,BER2] = BERvsSNR_Code2(Nb2,k2,g2,H0,"threshold detector");

EbNo=10.^(SNR_bit/10);
Pb0=0.5*erfc(sqrt(EbNo));
Pb1=0.5*erfc(sqrt(EbNo*rate1)).^dmin1;
Pb2=0.5*erfc(sqrt(EbNo*rate2)).^dmin2;
figure();
grid

semilogy(SNR_bit,BER0);
hold on
semilogy(SNR_bit,BER1)
semilogy(SNR_bit,BER2)

xlabel('Eb/No (in dB)')
ylabel('BER')
title('AWGN Channel Empirical BER vs Eb/N0'); 
legend('Uncoded','Code BCH 1','Code BCH 2')
hold off
figure();
grid
semilogy(SNR_bit,Pb0);
hold on
semilogy(SNR_bit,Pb1)
semilogy(SNR_bit,Pb2)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('AWGN Channel Theoretical BER vs Eb/N0');
legend('Uncoded','Code BCH 1','Code BCH 2')
hold off

%% %%%%%%%%%% plot BER versus Eb/No for uncoded BPSK with 4 channels%%%%%%%%
Nb0=3000;
[SNR_bit,BER_0_ZF] = BERvsSNR_Uncoded(Nb0,H0,'ZF','PSK',2,0,20,2);
[SNR_bit,BER_0_DF] = BERvsSNR_Uncoded(Nb0,H0,'DF','PSK',2,0,20,2);
[SNR_bit,BER_1_ZF] = BERvsSNR_Uncoded(Nb0,H1,'ZF','PSK',2,0,20,2);
[SNR_bit,BER_1_DF] = BERvsSNR_Uncoded(Nb0,H1,'DF','PSK',2,0,20,2);
[SNR_bit,BER_2_ZF] = BERvsSNR_Uncoded(Nb0,H2,"ZF",'PSK',2,0,20,2);
[SNR_bit,BER_2_DF] = BERvsSNR_Uncoded(Nb0,H2,"DF",'PSK',2,0,20,2);
[SNR_bit,BER_3_ZF] = BERvsSNR_Uncoded(Nb0,H3,'ZF','PSK',2,0,20,2);
[SNR_bit,BER_3_DF] = BERvsSNR_Uncoded(Nb0,H3,'DF','PSK',2,0,20,2);
figure();
grid
semilogy(SNR_bit,BER_0_ZF);
hold on
semilogy(SNR_bit,BER_0_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded BPSK AWGN');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_1_ZF);
hold on
semilogy(SNR_bit,BER_1_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded BPSK Channel 1');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_2_ZF);
hold on
semilogy(SNR_bit,BER_2_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded BPSK Channel 2');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_3_ZF);
hold on
semilogy(SNR_bit,BER_3_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded BPSK Channel 3');
legend('ZF','DF')
hold off

%% %%%%%%%%%% plot BER versus Eb/No for uncoded 8-QAM with 4 channels%%%%%%%%
Nb0=3000;
[SNR_bit,BER_0_ZF] = BERvsSNR_Uncoded(Nb0,H0,'ZF','QAM',8,0,30,3);
[SNR_bit,BER_0_DF] = BERvsSNR_Uncoded(Nb0,H0,'DF','QAM',8,0,30,3);
[SNR_bit,BER_1_ZF] = BERvsSNR_Uncoded(Nb0,H1,'ZF','QAM',8,0,30,3);
[SNR_bit,BER_1_DF] = BERvsSNR_Uncoded(Nb0,H1,'DF','QAM',8,0,30,3);
[SNR_bit,BER_2_ZF] = BERvsSNR_Uncoded(Nb0,H2,"ZF",'QAM',8,0,30,3);
[SNR_bit,BER_2_DF] = BERvsSNR_Uncoded(Nb0,H2,"DF",'QAM',8,0,30,3);
[SNR_bit,BER_3_ZF] = BERvsSNR_Uncoded(Nb0,H3,'ZF','QAM',8,0,30,3);
[SNR_bit,BER_3_DF] = BERvsSNR_Uncoded(Nb0,H3,'DF','QAM',8,0,30,3);
figure();
grid
semilogy(SNR_bit,BER_0_ZF);
hold on
semilogy(SNR_bit,BER_0_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 8-QAM AWGN');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_1_ZF);
hold on
semilogy(SNR_bit,BER_1_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 8-QAM Channel 1');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_2_ZF);
hold on
semilogy(SNR_bit,BER_2_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 8-QAM Channel 2');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_3_ZF);
hold on
semilogy(SNR_bit,BER_3_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 8-QAM Channel 3');
legend('ZF','DF')
hold off

%% %%%%%%%%%% plot BER versus Eb/No for uncoded 16-QAM with 4 channels%%%%%%%%
Nb0=3000;
[SNR_bit,BER_0_ZF] = BERvsSNR_Uncoded(Nb0,H0,'ZF','QAM',16,0,30,3);
[SNR_bit,BER_0_DF] = BERvsSNR_Uncoded(Nb0,H0,'DF','QAM',16,0,30,3);
[SNR_bit,BER_1_ZF] = BERvsSNR_Uncoded(Nb0,H1,'ZF','QAM',16,0,30,3);
[SNR_bit,BER_1_DF] = BERvsSNR_Uncoded(Nb0,H1,'DF','QAM',16,0,30,3);
[SNR_bit,BER_2_ZF] = BERvsSNR_Uncoded(Nb0,H2,"ZF",'QAM',16,0,30,3);
[SNR_bit,BER_2_DF] = BERvsSNR_Uncoded(Nb0,H2,"DF",'QAM',16,0,30,3);
[SNR_bit,BER_3_ZF] = BERvsSNR_Uncoded(Nb0,H3,'ZF','QAM',16,0,30,3);
[SNR_bit,BER_3_DF] = BERvsSNR_Uncoded(Nb0,H3,'DF','QAM',16,0,30,3);
figure();
grid
semilogy(SNR_bit,BER_0_ZF);
hold on
semilogy(SNR_bit,BER_0_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 16-QAM AWGN');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_1_ZF);
hold on
semilogy(SNR_bit,BER_1_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 16-QAM Channel 1');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_2_ZF);
hold on
semilogy(SNR_bit,BER_2_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 16-QAM Channel 2');
legend('ZF','DF')
hold off

figure();
grid
semilogy(SNR_bit,BER_3_ZF);
hold on
semilogy(SNR_bit,BER_3_DF)
xlabel('Eb/No (in dB)')
ylabel('BER')
title('BER vs Eb/N0 uncoded 16-QAM Channel 3');
legend('ZF','DF')
hold off