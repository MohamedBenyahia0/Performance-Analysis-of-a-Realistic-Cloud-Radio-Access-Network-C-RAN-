function BER = BER_Fiber_OOK(P_out_dbm,R,Fs,target_BER,attenuation,beta2,beta3,L,beta2_DCF,beta3_DCF,attenuation_DCF,L_DCF)
nb_err=0;
X=randi([0,1],1,100/target_BER); % generation des bits d'informations 
signal = repelem(X,5);
n=length(signal);
param_eml=make_emlaser('P_opt_dBm',P_out_dbm);
[S_out, Ts_out, powerOfBlock] = TX_optical_eml(signal, 1/Fs, param_eml);

S_out_fiber=opticalFiber(S_out,1/Ts_out,attenuation,beta2,beta3,L);
S_out_fiber_DCF=opticalFiber(S_out_fiber,1/Ts_out,attenuation_DCF,beta2_DCF,beta3_DCF,L_DCF);
params_detector=make_photodetector('B_e',Fs,'sensitivity',R);

N_opt=0; % optical noise spectral power density
[S, Ts, powerOfBlock,SNR_elec] = RX_photodetector(S_out_fiber_DCF, Ts_out, N_opt, params_detector);
P_out_lin=1e-3*10^(P_out_dbm/10);
S=S./(P_out_lin*R*exp(-attenuation*L/1e3)*exp(-attenuation_DCF*L_DCF/1e3));

S_hat=zeros(1,length(S));

for i=1:length(S)
    if (S(i)>=0.5)
        S_hat(i)=1;
    end
end
nb_err=nb_err+sum(abs(S_hat-signal));
BER=nb_err/n;
end