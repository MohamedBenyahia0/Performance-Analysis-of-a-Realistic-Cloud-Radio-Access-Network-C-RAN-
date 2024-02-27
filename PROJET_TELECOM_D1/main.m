g1=[1,0,1,0,0,1];%g=x^5+x^2+1 de degre n-k=5
k1=26;%k=31-5
rate1=k1/31;
m1=ones(1,k1);
c1=encoder(m1,g1);

c1_prim=c1;
c1_prim(4)=c1_prim(4)+1;
m1_hat=decoder_1error(c1_prim,g1);

%% %%%%%%%%%%%
g2=[1,0,0,1,0,1,1,0,1,1,1];%g=x^10+x^9+x^8+x^6+x^5+x^3+1
k2=31-10;
rate2=k2/31;
m2=ones(1,k2);
c2=encoder(m2,g2);

c2_prim=c2;
c2_prim(3)=c2_prim(3)+1;
c2_prim(6)=c2_prim(6)+1;
m2_hat=decoder_2errors(c2_prim,g2);
%% %%%%%
[H1,H2,H3] =plot_channel();

%% %%%%%%%%%% plot BER versus Eb/No %%%%%%%%
Nb0=10000;
Nb1=k1*100;
Nb2=k2*100;
[SNR_bit,BER0_DF] = BERvsSNR_Uncoded(Nb0,H1,'DF');
[SNR_bit,BER0_ZF] = BERvsSNR_Uncoded(Nb0,H1,'ZF');
%[SNR_bit,BER1] = BERvsSNR_Code1(Nb1,k1,g1,H1,'ZF');
%[SNR_bit,BER2] = BERvsSNR_Code2(Nb2,k2,g2);

EbNo=10.^(SNR_bit/10);
Pb0=0.5* erfc(sqrt(EbNo*2));
Pb1=0.5* erfc(sqrt(EbNo*2*rate1));

grid
xlabel('Eb/No (in dB)')
ylabel('BER')
semilogy(SNR_bit,Pb0);
hold on
semilogy(SNR_bit,BER0_ZF)
semilogy(SNR_bit,BER0_DF)
%semilogy(SNR_bit,BER1)
%semilogy(SNR_bit,BER2)
legend('Theoretical','Empirical ZF','DF')
% hold off