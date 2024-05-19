function BER = BER_Fiber_OOK(P_out_dbm,R_b,Fs,target_BER,attenuation,beta2,beta3,L)
nb_err=0;
X=randi([0,1],1,100/target_BER); % generation des bits d'informations 
params=make_emlaser('P_opt_dBm',P_out_dbm);
[S_out, Ts_out, powerOfBlock] = TX_optical_eml(X, 1/Fs, params);

S_out_fiber=opticalFiber(S_out,Fs,attenuation,beta2,beta3,L);

params_detector=make_photodetector('sensitivity',R_b);

N_opt=0; % optical noise spectral power density
[S, Ts, powerOfBlock,SNR_elec] = RX_photodetector(S_out_fiber, 1/Fs, N_opt, params_detector);
P_out_lin=10^(P_out_dbm/10);
S=S./(P_out_lin*R_b*exp(-attenuation*L/1e3));
thresh=(max(S)-min(S))/2;
S_hat=zeros(1,length(S));

for i=1:length(S)
    if (S(i)>=thresh)
        S_hat(i)=1;
    end
end
nb_err=nb_err+sum(abs(S_hat-X));
BER=nb_err/length(X);
end