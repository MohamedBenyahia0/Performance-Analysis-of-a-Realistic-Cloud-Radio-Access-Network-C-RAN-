% This script demonstrates the use of the rfPA function
% - with a thermally noised RF sinusoidal input signal

% Rev: March 2022, Chadi+Germain
% Rev: March 2023, Germain

close all;
clear;

% Simulation parameters
continuousTimeSamplingRate = 19.98e9; % Operating Frequency to emulate the behavior of a continuous time system
Ts_Cont                    = 1/continuousTimeSamplingRate; % Continuous time sampling period
N_Cont                     = 2^18;   % Number of signal points (@continuous time rate)

% System properties
BW  = 10e6;     % Signal bandwidth
K   = 1.38e-23; % Boltzmann Constant
T   = 290;      % room temperature
Rin = 50;       % Matching impedance chosen equal to 50

% RF parameters
Flo       = 2.4e9;
BW_rf     = 2*BW;

% Simulation subproperties
t_Cont       = (0:(N_Cont-1))*Ts_Cont; % Time vector (@continuous time rate)

% Input signal
fin_or          = Flo+1.1e6; % Input sine frequency
Bin_in          = round(fin_or/continuousTimeSamplingRate*N_Cont); % Determining the input bin
fin             = Bin_in*continuousTimeSamplingRate/N_Cont;

Pin             = 0;  %Pin in dBm
Ain             = sqrt(10.^((Pin-30)/10)*2*Rin);
AntennaNoise    = randn(1,N_Cont)*sqrt(K*T*continuousTimeSamplingRate/2*Rin);
in              = Ain*sin(2*pi*fin*t_Cont+rand()*2*pi)+AntennaNoise; % Input signal

% PA parameters
PA_model    = 2; %%% Select a model between 1, 2 or 3
PA_NF       = 5;
PA_Gain     = 40;

%%%%%%%%%%%% PA %%%%%%%%%%%
switch (PA_model) 
  case 1
   PA_IIP3  = 35;
  case 2
   PA_IIP3  = 30;
  case 3
   PA_IIP3  = 20;
  otherwise
   PA_IIP3  = 100;
end

% Operate PA
PA_out = rfPA(in,PA_Gain,PA_NF, PA_IIP3,Rin,continuousTimeSamplingRate/2);

%%%%%%%%%%%% Plotting %%%%%%%%%%%%
voltsq2mwatt        = 1e3/Rin; % Conversion factor from V^2 to milliWatt
window_number       = 1;
lineSpec_index      = 1;

plot_spectrum(PA_out*sqrt(voltsq2mwatt),window_number,continuousTimeSamplingRate,lineSpec_index);
hold on
plot_spectrum(in*sqrt(voltsq2mwatt),window_number,continuousTimeSamplingRate,lineSpec_index+1);
legend('output','input')
xlabel('frequency (Hz)')
ylabel('PSD (dBm/bin)')








