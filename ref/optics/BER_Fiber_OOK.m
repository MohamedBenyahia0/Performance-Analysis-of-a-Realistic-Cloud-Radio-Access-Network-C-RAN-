function BER = BER_Fiber_OOK(R,Fs,attenuation,beta2,beta3,L)
power_eml=-19;
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
[S_out, Ts_out, powerOfBlock] = TX_optical_eml(signal, 1/Fs, param_eml);

S_out_fiber=opticalFiber(S_out,Fs,attenuation,beta2,beta3,L);

params_detector=make_photodetector('B_e',Fs,'sensitivity',R);

N_opt=0; % optical noise spectral power density
[S, Ts, powerOfBlock,SNR_elec] = RX_photodetector(S_out_fiber, 1/Fs, N_opt, params_detector);
P_out_lin=1e-3*10^(power_eml_amp/10);
S=S./(P_out_lin*R*exp(-attenuation*L/1e3));

S_hat=zeros(1,length(S));

for i=1:length(S)
    if (S(i)>=0.5)
        S_hat(i)=1;
    end
end
S_hat = S_hat(fix(ovs/2)+1:ovs:length(S_hat));
nb_err=nb_err+sum(abs(S_hat-X));
BER=nb_err/n;
end