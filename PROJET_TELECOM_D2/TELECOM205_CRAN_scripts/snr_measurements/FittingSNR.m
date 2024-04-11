% FittingSNR.m - a script that demonstrates how to make a fitting on SNR
% measurements
%
% This script is a simple script to teach how to perform a polynomial fitting on
% measurement data
%
% For illustration purpose, the script starts by creating some dummy noisy data
% on which the fitting will be performed. %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% YOU MUST REPLACE THE DUMMY NOISY DATA BY YOUR OWN MEASUREMENT DATA ! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Chadi Jabbour, Mar. 2022
% Last revision: Germain Pham, Jan. 2023



close all
clear
% %%
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
% if isOctave
%     c_light   = physical_constant("speed of light in vacuum");
%     k         = physical_constant("boltzmann constant");
% else
%     c_light   = physconst('LightSpeed');
%     k         = physconst('Boltzmann');
% end
% 
% 
% % RF parameters
% Fc                = 2.4e9;       % Central frequency of the band
% BW                = 20e6;        % Considered bandwidth of the system
% lambda            = c_light/Fc;  % The wavelength
% alpha             = 2;           % The path loss exponent =2 for the free space
% Pt                = 20;          % transmitted power in dBm
% T                 = 290;         % Temperature of the system
% ThermNoise_dBm    = 10*log10(k*T*BW/1e-3); % Therm Noise power in dBm at the receiver input
% 
% NFreceiver        = 8; % Receiver NF in dB
% noise_pow_dBm     = NFreceiver+ThermNoise_dBm; % Noise power at receiver output
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Generate dummy noisy data %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance          = [5:5:55];     % (m) vector the distances of measurements ; PLEASE ADAPT !
% signal_pow_dBm    = Pt - 10*alpha*log10(4*pi.*distance/lambda); % signal power at receiver
% SNR_theo_dB       = signal_pow_dBm-noise_pow_dBm;
% 
% % Adding random variations to the dummy data to make it more real
% SNR_real_dB          = SNR_theo_dB+randn(1,length(SNR_theo_dB));

%% Real Measurement data %%%
distance          = [5:5:55];     % (m) vector the distances of measurements ;
Power2_5GHz_dB=zeros(1,11);
Power600MHz_dB=zeros(1,11);
for i=1:11
    filename1='freq2.5GHZ//distance_'+string(5*i)+'m_sig_pows.csv';
    M=readmatrix(filename1);
    Power2_5GHz_dB(i)=M(8,2); % Gain 35 dB

    filename2='freq600MHZ//distance_'+string(5*i)+'m_sig_pows.csv';
    A=readmatrix(filename2);
    Power600MHz_dB(i)=A(6,2); % Gain 5 dB

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Demonstrate polynomial fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distance_dB           = 10*log10(distance); % Converting the distance to dB as the SNR_dB = 10*alpha*distance_dB + Cte
poly_ord              = 1; % Order of the polynomial to be fitted
poly_coef_fit         = polyfit(distance_dB,Power2_5GHz_dB,poly_ord);  %%% Linear fitting of the coefficient

distance_interp       = [5:5:55]; %%% distance vector of the interpolated function
distance_interp_dB    = 10*log10(distance_interp);
fittedcurve           = polyval(poly_coef_fit,distance_interp_dB);


disp(['The expression of the ' num2str(poly_ord) ' order fitted polynomial is '])
vpa(poly2sym(poly_coef_fit),3)

disp('The propagation exponent estimated from measurement is: ')
disp(num2str(-poly_coef_fit(end-1))) % Please read doc of polyfit for coeff ordering details


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Display real measured data and fitted curve %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogx(distance,Power2_5GHz_dB,'s-')   
hold all
semilogx(distance_interp,fittedcurve)
xlabel('Distance(m)')
ylabel('Power(dB)')
legend('Power real', 'fitted curve')
