close all;

K            = 1.38e-23;   % Boltzmann Constant
T            = 290;        % Room temperature
R            = 50;         % Input impedance

N            = 2^13;       % Number of signal points ADC
BW           = 10e6;       % Signal bandwidth
Flo          = 0.9e9;      % Central frequency




Fs_ADC       = 15e6;       % Sampling frequency ADC
Fs_Cont      = 15e9;       % Sampling frequency Simulation
Ts_Cont      = 1/Fs_Cont;  % Sampling period Simulation
Vref         = 1;          % Reference voltage of the ADC
Nbits_ADC    = 10;         % Number of bits for the ADC



N_Sim        = floor(N*Fs_Cont/Fs_ADC);                 % Number of  Simulation points
t_Sim        = 0:Ts_Cont:(N_Sim-1)*Ts_Cont;             % Time vetor Simulation
f_Sim        = 0:Fs_Cont/N_Sim:Fs_Cont/2-Fs_Cont/N_Sim; % Frequency vector Simulation


Ts_ADC       = 1/Fs_ADC;                                % Sampling period ADC
t_ADC        = 0:Ts_ADC:(N-1)*Ts_ADC;                   % Time vetor Simulation
f_ADC        = 0:Fs_ADC/N:Fs_ADC/2-Fs_ADC/N;            % Frequency vector Simulation

fin_or       = [1.1e6,2.1e6];                           % Input sine frequency
Bin_in       = round(fin_or./Fs_ADC*N);                 % Determining the input bin
fin          = Bin_in.*Fs_ADC/N;                        % 


Pin          = -40; % Pin in dBm
Ain          = sqrt(10.^((Pin-30)/10)*2*R);

AntennaNoise = randn(1,N_Sim)*sqrt(K*T*BW*R);

G_LNA        = 30;      % Gain of the first stage in dB
NF_LNA       = 4;       % Noise figure of the first stage in dB
IIP3_LNA     = -20;     % Set very high for Simplicity

Filter_order = 7;

% Generate dualtone signal
in = 0;
for i = 1:length(fin)
  in = in+Ain*sin(2*pi*(Flo+fin(i))*t_Sim+rand()*180); % Input signal
end
in        = in+AntennaNoise;
out_LNA   = rfLNA(in,G_LNA,NF_LNA,IIP3_LNA,R,Fs_Cont/2); % output Signal of the amplifier

[outI,outQ] = downMixer(out_LNA,Flo,Fs_Cont);



[filter_num,filter_denum]   = butter(Filter_order,10*BW/Fs_Cont);
Out_Filter                  = filter(filter_num,filter_denum,outI);

ADCdelay = 0;
Out_ADC  = ADC(Out_Filter,Nbits_ADC,Vref,Fs_ADC,ADCdelay,Fs_Cont);
%SNR_out=perf_estim(Out_ADC,1,Bin_in,5,1);   


subplot(2,2,1)
plot_spectrum(out_LNA*sqrt(1000/R),1,Fs_Cont);
title('Spectrum output LNA')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')
subplot(2,2,2)

plot_spectrum(outI*sqrt(1000/R),1,Fs_Cont,1);
title('Spectrum output mixer')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')

subplot(2,2,3)
plot_spectrum(Out_Filter*sqrt(1000/R),1,Fs_Cont,1);
title('Spectrum output Filter')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')


subplot(2,2,4)
plot_spectrum(Out_ADC*sqrt(1000/R),1,Fs_ADC,1);
title('Spectrum output ADC')
xlabel('Frequency(Hz)')
ylabel('PSD(dBm/bin)')


