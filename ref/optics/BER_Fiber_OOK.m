function BER = BER_Fiber_OOK(power_eml,R,Fs,attenuation,beta2,beta3,L)

att_dB=0.2*(L/1e3);
power_eml_amp = power_eml + att_dB;
nb_err=0;
X=randi([0,1],1,1e5); % generation des bits d'informations 
ovs=5;
signal = repelem(X,ovs);
Tsymb=1/Fs;
Te = Tsymb/ovs;
n=length(X);

param_eml=make_emlaser('P_opt_dBm',power_eml_amp);
[S_out, Ts_out, powerOfBlock] = TX_optical_eml(signal, Tsymb, param_eml);

S_out_fiber=opticalFiber(S_out,1/Te,attenuation,beta2,beta3,L);

params_detector=make_photodetector('B_e',Fs,'sensitivity',R);

N_opt=0; % optical noise spectral power density
[S, Ts, powerOfBlock,SNR_elec] = RX_photodetector(S_out_fiber, Ts_out, N_opt, params_detector);
P_out_lin=1e-3*10^(power_eml_amp/10);
S=S./(P_out_lin*R*exp(-attenuation*L/1e3));

S_hat=zeros(1,length(S));

for i=1:length(S)
    if (S(i)>=0.5)
        S_hat(i)=1;
    end
end
S_hat = S_hat(fix(ovs/2)+1:ovs:length(S_hat));
nb_err=biterr(X,S_hat);
BER=nb_err/n;
end