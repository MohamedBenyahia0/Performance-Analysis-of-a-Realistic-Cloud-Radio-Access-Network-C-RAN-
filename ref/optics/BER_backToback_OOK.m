function BER = BER_backToback_OOK(P_out_dbm,R,Fs,target_BER)
nb_err=0;
X=randi([0,1],1,100/target_BER); % generation des bits d'informations 
params=make_emlaser('P_opt_dBm',P_out_dbm);
[S_out, Ts_out, powerOfBlock] = TX_optical_eml(X, 1/Fs, params);
params_detector=make_photodetector('B_e',Fs,'sensitivity',R);

N_opt=0; % optical noise spectral power density
[S, Ts, powerOfBlock,SNR_elec] = RX_photodetector(S_out, 1/Fs, N_opt, params_detector);
P_out_lin=1e-3*10^(P_out_dbm/10);
S=S./(P_out_lin*R);

S_hat=zeros(1,length(S));

for i=1:length(S)
    if (S(i)>=0.5)
        S_hat(i)=1;
    end
end

nb_err=nb_err+sum(abs(S_hat-X));
BER=nb_err/length(X);
end