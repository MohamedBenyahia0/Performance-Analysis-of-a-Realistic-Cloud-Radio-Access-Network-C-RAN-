close all

% Parameters
c = 3e8;  % Speed of light in m/s
n = 1;    % Refractive index
T_bit = 1e-9;  % Bit period in seconds
lambda=1550e-9;
v = c / lambda;  %Optical frequency
Ts = T_bit;  % Sampling time

% Initialize variables for I_DC and I_AC
I_DC_values = 1000:500:10000;  % Array of I_DC values to test (in mA)
%I_DC_values = 7; %145:1:145;
I_AC = 2;  % Fixed I_AC value (in mA)

% Initialize results storage
results = struct();

% Fiber parameters
L = 20e3; % Length in meters
D = 17e-6; % s/m^2 (chromatic dispersion)
S = 0.09e3; % s/m^2/m (dispersion slope)
attenuation=0.2; % dB/km

% Calculate fiber parameters
delta_lambda = (2 / T_bit) * lambda^2 / c;
beta2 = -D * lambda^2 / (2 * pi * c);
beta3 = S * lambda^4 / (4 * pi^2 * c^2);

% Loop over different I_DC values
for I_DC=I_DC_values
    
    % Generate OOK signal
    num_bits = 1000;
    bits = randi([0, 1], 1, num_bits);
    OOK_sig = bits * I_AC * 1e-3 + I_DC * 1e-6; % Convertir mA en A
    
    % Create laser
    laser = make_laser_simple('v', v);
    
    % Transmit signal
    [S_opt, Ts_opt, powerOfBlock_tx] = TX_optical_dml(OOK_sig, Ts, laser);
    opticalPower_dBm=10*log10(mean(abs(S_opt).^2)/1e3);

    %figure();
    %plot(OOK_sig)

     % Pass through fiber
    S_opt_fiber = opticalFiber(S_opt, 1 / Ts_opt, attenuation, beta2, beta3, L);

    % Create photodetector
    photodetector = make_photodetector('B_e',1e9);

    % Conversion optique-électrique
    [S_elec, Ts_elec, powerOfBlock_rx, SNR_elec] = RX_photodetector(S_opt_fiber, Ts_opt, photodetector.N_th, photodetector);
    figure();
    plot(S_elec)
    % Normalisation du signal électrique
    S_elec_mean = mean(S_elec);
    S_elec_std = std(S_elec);
    S_elec_norm = (S_elec - S_elec_mean) / S_elec_std; % Normalisation
%     hold on
%     plot(S_elec_norm)
%     hold off

    S_elec_norm(S_elec_norm > 0) = 1;
    S_elec_norm(S_elec_norm < 0) = 0;

    S_elec_norm_downsample = downsample(S_elec_norm, 4505,round(4505/2));


    %% Calculate BER
    BER = biterr(bits,S_elec_norm_downsample)/num_bits;

    % Calculate BER
    %BER = sum(bits ~= S_elec_norm_downsample) / num_bits;


    % Stock the results
    results.(sprintf('I_DC_%duA', I_DC)).S_opt = S_opt;
    results.(sprintf('I_DC_%duA', I_DC)).Ts_opt = Ts_opt;
    results.(sprintf('I_DC_%duA', I_DC)).powerOfBlock_tx = powerOfBlock_tx;
    results.(sprintf('I_DC_%duA', I_DC)).S_elec = S_elec;
    results.(sprintf('I_DC_%duA', I_DC)).Ts_elec = Ts_elec;
    results.(sprintf('I_DC_%duA', I_DC)).powerOfBlock_rx = powerOfBlock_rx;
    results.(sprintf('I_DC_%duA', I_DC)).SNR_elec = SNR_elec;
    results.(sprintf('I_DC_%duA', I_DC)).BER = BER;
    results.(sprintf('I_DC_%duA', I_DC)).opticalPower_dBm = opticalPower_dBm;
    results.(sprintf('I_DC_%duA', I_DC)).S_elec_norm_downsample = S_elec_norm_downsample;
end

% Plot results for different I_DC values
figure;
hold on;
for I_DC = I_DC_values
    S_opt = results.(sprintf('I_DC_%duA', I_DC)).S_opt;
    plot(abs(S_opt).^2, 'DisplayName', sprintf('I_{DC} = %d uA', I_DC));
end
xlabel('Sample Index');
ylabel('Optical Power TX(W)');
title('Optical Power TX(W) for Different I_{DC} Values');
legend;
hold off;

% Plot Power of Block for different I_DC values
figure;
% Ectract all the values of powerOfBlock_tx for each I_DC values
powerValues = arrayfun(@(x) results.(sprintf('I_DC_%duA', x)).powerOfBlock_tx, I_DC_values);
plot(I_DC_values, powerValues, '-o');
xlabel('I_{DC} (uA)');
ylabel('Power of Block TX(W)');
title('Power of Block TX vs I_{DC}');
grid on;

figure;
% Ectract all the values of powerOfBlock for each I_DC values
powerValues = arrayfun(@(x) results.(sprintf('I_DC_%duA', x)).powerOfBlock_rx, I_DC_values);
plot(I_DC_values, powerValues, '-o');
xlabel('I_{DC} (uA)');
ylabel('Power of Block RX(W)');
title('Power of Block RX vs I_{DC}');
grid on;

% Traçage du SNR électrique pour différentes valeurs de I_DC
figure;
SNRValues = arrayfun(@(x) results.(sprintf('I_DC_%duA', x)).SNR_elec, I_DC_values);
plot(I_DC_values, SNRValues, '-o');
xlabel('I_{DC} (uA)');
ylabel('SNR électrique (dB)');
title('SNR électrique en fonction de I_{DC}');
grid on;

% Traçage du BER pour différentes valeurs de I_DC
figure;
BERValues = arrayfun(@(x) results.(sprintf('I_DC_%duA', x)).BER, I_DC_values);
plot(I_DC_values, BERValues, '-o');
xlabel('I_{DC} (uA)');
ylabel('Taux d erreur binaire (BER)');
title('BER en fonction de I_{DC}');
grid on;

% 
% figure;
% opticalPower = arrayfun(@(x) results.(sprintf('I_DC_%dmA', x)).opticalPower_dBm, I_DC_values);
% semilogy(opticalPower,BERValues,'-o')
% xlabel('Optical Power')
% ylabel('BER')
% title('BER versus Optical Power')
% grid on;

% Plot BER vs Optical Power for different I_DC values
figure;
hold on;
for I_DC = I_DC_values
    opticalPower_dBm = results.(sprintf('I_DC_%duA', I_DC)).opticalPower_dBm;
    BER = results.(sprintf('I_DC_%duA', I_DC)).BER;
    semilogy(opticalPower_dBm, BER, '-o', 'DisplayName', sprintf('I_{DC} = %d uA', I_DC));
end
% Les mettre dans une liste pour le plot
xlabel('Optical Power (dBm)');
ylabel('BER');
title('BER vs Optical Power for Different I_{DC} Values');
legend;
grid on;
hold off;

% % Plot results for different I_DC values (electrical signal)
% figure;
% hold on;
% 
% plot(bits, 'o', 'Color', 'r', 'DisplayName', 'Input');
% 
% % Plot results for different I_DC values (optical signal)
% for I_DC = I_DC_values
%     S_elec_norm_downsample = results.(sprintf('I_DC_%dmA', I_DC)).S_elec_norm_downsample;
%     plot(S_elec_norm_downsample, 'o', 'Color', 'b', 'DisplayName', sprintf('I_{DC}opt = %d mA', I_DC));
% end
% 
% xlabel('Sample Index');
% ylabel('Signal Normalized');
% title('Input vs Output Normalized');
% legend;
% hold off;