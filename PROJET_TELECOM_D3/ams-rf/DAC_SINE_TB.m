close all;
clear;

%%% DAC Specifications %%%
Nbits_DAC               = 18;
Vref_DAC                = 1;
Fs_DAC                  = 30e6;
basebandSamplingRate    = Fs_DAC;


%%% General %%%
continuousTimeSamplingRate    = 900e6; % Operating Frequency to emulate the behavior of a continuous time system
Ts_Cont                       = 1/continuousTimeSamplingRate; % Continuous time sampling period
Subsamp_fdac                  = round(continuousTimeSamplingRate/basebandSamplingRate);
N   = 2^13;     % Number of signal points
BW  = 10e6;     % Signal bandwidth
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
t_Cont = 0:Ts_Cont:(N*Subsamp_fdac-1)*Ts_Cont; % Time vetor
f_Cont = 0:continuousTimeSamplingRate/N:continuousTimeSamplingRate/2-continuousTimeSamplingRate/N; % Frequency vector

Ts_BB  = 1/basebandSamplingRate; % Baseband sampling period
t_BB   = 0:Ts_BB:(N-1)*Ts_BB;    % Baseband Time vetor


Rin   = 50;  % Matching impedance chosen equal to 50




%%/Input Signal Specifications %%%%
fin_or    = 1e6;    % Input sine frequency
Bin_in    = round(fin_or/basebandSamplingRate*N); % Determining the input bin
fin       = Bin_in*Fs_DAC/N;
Ain       = 1;
% Input signal
Input     = Ain*sin(2*pi*fin*t_BB+rand()*180);

mode              = 'impulse'; % "Zero Padding"
AnalogOutputZP    = DAC(Input',Nbits_DAC,Vref_DAC,mode,Fs_DAC,continuousTimeSamplingRate);

n_samples_plot=100;
figure(1)
plot(t_BB(1:n_samples_plot),Input(1:n_samples_plot));
hold all
plot(t_Cont(1:n_samples_plot*Subsamp_fdac),AnalogOutputZP(1:n_samples_plot*Subsamp_fdac),'linewidth',2);
hold off
xlabel('Time (s)')
ylabel('Amplitude (V)')
legend('Input','Output')
title('DAC signal')








