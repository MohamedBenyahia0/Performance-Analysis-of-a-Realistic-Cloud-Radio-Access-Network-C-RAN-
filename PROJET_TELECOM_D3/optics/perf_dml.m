close all

% Parameters
c = 3e8;  % Speed of light in m/s
n = 1;    % Refractive index
T_bit = 1e-9;  % Bit period in seconds
lambda = 1550e-9;
v = c / lambda;  % Optical frequency
Ts = T_bit;  % Sampling time

% Initialize variables for I_DC and I_AC
I_DC_values = 3000:500:9000;  % Array of I_DC values to test (in uA)
I_AC = 2;  % Fixed I_AC value (in mA)

% Preallocate results storage
results_fiber = cell(1, length(I_DC_values));
results_BtB = cell(1, length(I_DC_values));

% Fiber parameters
L = 20e3; % Length in meters
D = 17e-6; % s/m^2 (chromatic dispersion)
S = 0.09e3; % s/m^2/m (dispersion slope)
attenuation = 0.2; % dB/km

% Calculate fiber parameters
delta_lambda = (2 / T_bit) * lambda^2 / c;
beta2 = -D * lambda^2 / (2 * pi * c);
beta3 = S * lambda^4 / (4 * pi^2 * c^2);

% Generate OOK signal
num_bits = 10000;
bits = randi([0, 1], 1, num_bits);

% Loop over different I_DC values using parfor
parfor idx = 1:length(I_DC_values)
    I_DC = I_DC_values(idx);
    
    OOK_sig = bits * I_AC * 1e-3 + I_DC * 1e-6; % Convertir µA en A
    
    % Create laser
    laser = make_laser_simple('v', v);
    
    % Transmit signal
    [S_opt, Ts_opt, powerOfBlock_tx] = TX_optical_dml(OOK_sig, Ts, laser);
    opticalPower_dBm = 10 * log10(mean(abs(S_opt).^2) / 1e-3);

    % Pass through fiber
    S_opt_fiber = opticalFiber(S_opt, 1 / Ts_opt, attenuation, beta2, beta3, L);

    % Create photodetector
    photodetector = make_photodetector('B_e', 1e9);

    % Conversion optique-électrique avec fibre
    [S_elec_fiber, Ts_elec_fiber, powerOfBlock_rx_fiber, SNR_elec_fiber] = RX_photodetector(S_opt_fiber, Ts_opt, photodetector.N_th, photodetector);

    %opticalPower_dBm = 10 * log10(mean(abs(S_elec_fiber).^2) / 1e-3);
    % Normalisation du signal électrique (fibre)
    S_elec_fiber_mean = mean(S_elec_fiber);
    S_elec_fiber_std = std(S_elec_fiber);
    S_elec_fiber_norm = (S_elec_fiber - S_elec_fiber_mean) / S_elec_fiber_std; % Normalisation

    S_elec_fiber_norm(S_elec_fiber_norm > 0) = 1;
    S_elec_fiber_norm(S_elec_fiber_norm < 0) = 0;

    S_elec_norm_fiber_downsample = downsample(S_elec_fiber_norm, 4505, round(4505 / 2));

    % Conversion optique-électrique avec BtB
    [S_elec_BtB, Ts_elec_BtB, powerOfBlock_rx_BtB, SNR_elec_BtB] = RX_photodetector(S_opt, Ts_opt, photodetector.N_th, photodetector);

    % Normalisation du signal électrique (BtB)
    S_elec_BtB_mean = mean(S_elec_BtB);
    S_elec_BtB_std = std(S_elec_BtB);
    S_elec_BtB_norm = (S_elec_BtB - S_elec_BtB_mean) / S_elec_BtB_std; % Normalisation

    S_elec_BtB_norm(S_elec_BtB_norm > 0) = 1;
    S_elec_BtB_norm(S_elec_BtB_norm < 0) = 0;

    S_elec_norm_BtB_downsample = downsample(S_elec_BtB_norm, 4505, round(4505 / 2));

    % Calculate BER
    BER_fiber = biterr(bits, S_elec_norm_fiber_downsample) / num_bits;
    if BER_fiber == 0
        BER_fiber = 1e-6;
    end

    BER_BtB = biterr(bits, S_elec_norm_BtB_downsample) / num_bits;
    if BER_BtB == 0
        BER_BtB = 1e-6;
    end

    % Stock the results for fiber
    results_fiber{idx}.S_opt = S_opt;
    results_fiber{idx}.Ts_opt = Ts_opt;
    results_fiber{idx}.powerOfBlock_tx = powerOfBlock_tx;
    results_fiber{idx}.S_elec = S_elec_fiber;
    results_fiber{idx}.Ts_elec = Ts_elec_fiber;
    results_fiber{idx}.powerOfBlock_rx = powerOfBlock_rx_fiber;
    results_fiber{idx}.SNR_elec = SNR_elec_fiber;
    results_fiber{idx}.BER = BER_fiber;
    results_fiber{idx}.opticalPower_dBm = opticalPower_dBm;
    results_fiber{idx}.S_elec_norm_downsample = S_elec_norm_fiber_downsample;

    % Stock the results for BtB
    results_BtB{idx}.S_opt = S_opt;
    results_BtB{idx}.Ts_opt = Ts_opt;
    results_BtB{idx}.powerOfBlock_tx = powerOfBlock_tx;
    results_BtB{idx}.S_elec = S_elec_BtB;
    results_BtB{idx}.Ts_elec = Ts_elec_BtB;
    results_BtB{idx}.powerOfBlock_rx = powerOfBlock_rx_BtB;
    results_BtB{idx}.SNR_elec = SNR_elec_BtB;
    results_BtB{idx}.BER = BER_BtB;
    results_BtB{idx}.opticalPower_dBm = opticalPower_dBm;
    results_BtB{idx}.S_elec_norm_downsample = S_elec_norm_BtB_downsample;
end

% Convert cell arrays to structures
results.fiber = cell2struct(results_fiber, arrayfun(@(x) sprintf('I_DC_%duA', x), I_DC_values, 'UniformOutput', false), 2);
results.BtB = cell2struct(results_BtB, arrayfun(@(x) sprintf('I_DC_%duA', x), I_DC_values, 'UniformOutput', false), 2);

% Plot results for different I_DC values (Fiber)
figure;
hold on;
for I_DC = I_DC_values
    S_opt = results.fiber.(sprintf('I_DC_%duA', I_DC)).S_opt;
    plot(abs(S_opt).^2, 'DisplayName', sprintf('Fiber I_{DC} = %d uA', I_DC));
end
xlabel('Sample Index');
ylabel('Optical Power TX(W)');
title('Optical Power TX(W) for Different I_{DC} Values (Fiber)');
legend;
hold off;

% Plot results for different I_DC values (BtB)
figure;
hold on;
for I_DC = I_DC_values
    S_opt = results.BtB.(sprintf('I_DC_%duA', I_DC)).S_opt;
    plot(abs(S_opt).^2, 'DisplayName', sprintf('BtB I_{DC} = %d uA', I_DC));
end
xlabel('Sample Index');
ylabel('Optical Power TX(W)');
title('Optical Power TX(W) for Different I_{DC} Values (BtB)');
legend;
hold off;

% Plot Power of Block for different I_DC values (Fiber)
figure;
powerValues_fiber = arrayfun(@(x) results.fiber.(sprintf('I_DC_%duA', x)).powerOfBlock_tx, I_DC_values);
plot(I_DC_values, powerValues_fiber, '-o');
xlabel('I_{DC} (uA)');
ylabel('Power of Block TX(W)');
title('Power of Block TX vs I_{DC} (Fiber)');
grid on;

% Plot Power of Block for different I_DC values (BtB)
figure;
powerValues_BtB = arrayfun(@(x) results.BtB.(sprintf('I_DC_%duA', x)).powerOfBlock_tx, I_DC_values);
plot(I_DC_values, powerValues_BtB, '-o');
xlabel('I_{DC} (uA)');
ylabel('Power of Block TX(W)');
title('Power of Block TX vs I_{DC} (BtB)');
grid on;

% Plot Power of Block RX for different I_DC values (Fiber)
figure;
powerValues_rx_fiber = arrayfun(@(x) results.fiber.(sprintf('I_DC_%duA', x)).powerOfBlock_rx, I_DC_values);
plot(I_DC_values, powerValues_rx_fiber, '-o');
xlabel('I_{DC} (uA)');
ylabel('Power of Block RX(W)');
title('Power of Block RX vs I_{DC} (Fiber)');
grid on;

% Plot Power of Block RX for different I_DC values (BtB)
figure;
powerValues_rx_BtB = arrayfun(@(x) results.BtB.(sprintf('I_DC_%duA', x)).powerOfBlock_rx, I_DC_values);
plot(I_DC_values, powerValues_rx_BtB, '-o');
xlabel('I_{DC} (uA)');
ylabel('Power of Block RX(W)');
title('Power of Block RX vs I_{DC} (BtB)');
grid on;

% Traçage du SNR électrique pour différentes valeurs de I_DC (Fiber)
figure;
SNRValues_fiber = arrayfun(@(x) results.fiber.(sprintf('I_DC_%duA', x)).SNR_elec, I_DC_values);
plot(I_DC_values, SNRValues_fiber, '-o');
xlabel('I_{DC} (uA)');
ylabel('SNR électrique (dB)');
title('SNR électrique en fonction de I_{DC} (Fiber)');
grid on;

% Traçage du SNR électrique pour différentes valeurs de I_DC (BtB)
figure;
SNRValues_BtB = arrayfun(@(x) results.BtB.(sprintf('I_DC_%duA', x)).SNR_elec, I_DC_values);
plot(I_DC_values, SNRValues_BtB, '-o');
xlabel('I_{DC} (uA)');
ylabel('SNR électrique (dB)');
title('SNR électrique en fonction de I_{DC} (BtB)');
grid on;

% Traçage du BER pour différentes valeurs de I_DC (Fiber)
figure;
BERValues_fiber = arrayfun(@(x) results.fiber.(sprintf('I_DC_%duA', x)).BER, I_DC_values);
plot(I_DC_values, BERValues_fiber, '-o');
xlabel('I_{DC} (uA)');
ylabel('Taux d erreur binaire (BER)');
title('BER en fonction de I_{DC} (Fiber)');
grid on;

% Traçage du BER pour différentes valeurs de I_DC (BtB)
figure;
BERValues_BtB = arrayfun(@(x) results.BtB.(sprintf('I_DC_%duA', x)).BER, I_DC_values);
plot(I_DC_values, BERValues_BtB, '-o');
xlabel('I_{DC} (uA)');
ylabel('Taux d erreur binaire (BER)');
title('BER en fonction de I_{DC} (BtB)');
grid on;

% Plot BER vs Optical Power for different I_DC values (Fiber)
figure;
hold on;
for I_DC = I_DC_values
    opticalPower_dBm = results.fiber.(sprintf('I_DC_%duA', I_DC)).opticalPower_dBm;
    BER = results.fiber.(sprintf('I_DC_%duA', I_DC)).BER;
    semilogy(opticalPower_dBm, BER, '-o', 'DisplayName', sprintf('Fiber I_{DC} = %d uA', I_DC));
end
xlabel('Optical Power (dBm)');
ylabel('BER');
title('BER vs Optical Power for Different I_{DC} Values (Fiber)');
legend;
grid on;
hold off;

% Plot BER vs Optical Power for different I_DC values (BtB)
figure;
hold on;
for I_DC = I_DC_values
    opticalPower_dBm = results.BtB.(sprintf('I_DC_%duA', I_DC)).opticalPower_dBm;
    BER = results.BtB.(sprintf('I_DC_%duA', I_DC)).BER;
    semilogy(opticalPower_dBm, BER, '-o', 'DisplayName', sprintf('BtB I_{DC} = %d uA', I_DC));
end
xlabel('Optical Power (dBm)');
ylabel('BER');
title('BER vs Optical Power for Different I_{DC} Values (BtB)');
legend;
grid on;
hold off;

% Prepare data for BER vs Optical Power plot
opticalPower_dBm_fiber_values = [];
BER_fiber_values = [];
opticalPower_dBm_BtB_values = [];
BER_BtB_values = [];

for I_DC = I_DC_values
    opticalPower_dBm_fiber = results.fiber.(sprintf('I_DC_%duA', I_DC)).opticalPower_dBm;
    BER_fiber = results.fiber.(sprintf('I_DC_%duA', I_DC)).BER;
    opticalPower_dBm_BtB = results.BtB.(sprintf('I_DC_%duA', I_DC)).opticalPower_dBm;
    BER_BtB = results.BtB.(sprintf('I_DC_%duA', I_DC)).BER;

    opticalPower_dBm_fiber_values = [opticalPower_dBm_fiber_values, opticalPower_dBm_fiber];
    BER_fiber_values = [BER_fiber_values, BER_fiber];
    opticalPower_dBm_BtB_values = [opticalPower_dBm_BtB_values, opticalPower_dBm_BtB];
    BER_BtB_values = [BER_BtB_values, BER_BtB];
end

% Plot BER vs I_DC and BER vs Optical Power on the same figure for Fiber and BtB

figure;

% Subplot 1: BER vs I_DC
subplot(1, 2, 1);
hold on;
BERValues_fiber = arrayfun(@(x) results.fiber.(sprintf('I_DC_%duA', x)).BER, I_DC_values);
BERValues_BtB = arrayfun(@(x) results.BtB.(sprintf('I_DC_%duA', x)).BER, I_DC_values);
semilogy(I_DC_values, BERValues_fiber, '-o', 'DisplayName', 'Fiber');
semilogy(I_DC_values, BERValues_BtB, '-o', 'DisplayName', 'BtB');
xlabel('I_{DC} (uA)');
ylabel('BER');
title('BER for Different I_{DC} Values (I_{AC} = 2 mA)');
legend;
grid on;
hold off;

% Subplot 2: BER vs Optical Power
subplot(1, 2, 2);
hold on;
semilogy(opticalPower_dBm_fiber_values, BER_fiber_values, '-o', 'DisplayName', 'Fiber');
semilogy(opticalPower_dBm_BtB_values, BER_BtB_values, '-o', 'DisplayName', 'BtB');
xlabel('Optical Power of photodetector (dBm)');
ylabel('BER');
title('BER vs Optical Power for Different I_{DC} Values (I_{AC} = 2 mA)');
legend;
grid on;
hold off;
