close all
clear

% % RF parameters
% Fc                = 2.4e9;       % Central frequency of the band
% BW                = 20e6;        % Considered bandwidth of the system
% c_light           = physconst('LightSpeed');
% lambda            = c_light/Fc;  % The wavelength
% alpha             = 2;           % The path loss exponent =2 for the free space
% Pt                = 20;          % transmitted power in dBm
% T                 = 290;         % Temperature of the system
% k                 = physconst('Boltzmann');
% ThermNoise_dBm    = 10*log10(k*T*BW/1e-3); % Therm Noise power in dBm at the receiver input
% 
% NFreceiver        = 8; % Receiver NF in dB

% Real Measurement data
distance          = [5:5:55];     % (m) vector the distances of measurements ;
PowerByGain2_5GHz = zeros(8, 11); % Store power values for 2.5GHz
PowerByGain600MHz = zeros(8, 11); % Store power values for 600MHz

% Load data for 2.5 GHz
for gain = 0:5:35
    for i = 1:11
        filename = strcat('freq2.5GHZ//distance_', string(5*i), 'm_sig_pows.csv');
        A = readmatrix(filename);
        PowerByGain2_5GHz(gain/5+1, i) = A(5+gain/5, 2) - gain; % Adjusting power based on gain
    end
end

% Load data for 600 MHz
for gain = 0:5:35
    for i = 1:11
        filename = strcat('freq600MHZ//distance_', string(5*i), 'm_sig_pows.csv');
        A = readmatrix(filename);
        PowerByGain600MHz(gain/5+1, i) = A(1+gain/5, 2) - gain; % Adjusting power based on gain
    end
end

% Demonstrate polynomial fitting for each gain
for gain = 0:5:35
    % Create a new figure for each gain
    figure;
    
    distance_dB           = 10*log10(distance); % Converting the distance to dB as the SNR_dB = 10*alpha*distance_dB + Cte
    poly_ord              = 1; % Order of the polynomial to be fitted

    % Fitting for 2.5GHz
    subplot(2, 1, 1);
    hold on;
    poly_coef_fit         = polyfit(distance_dB, PowerByGain2_5GHz(gain/5+1, :), poly_ord);  %%% Linear fitting of the coefficient
    fittedcurve_2_5GHz    = polyval(poly_coef_fit, distance_dB);
    
    % Plot for 2.5GHz
    semilogx(distance, PowerByGain2_5GHz(gain/5+1, :), 's-', 'DisplayName', 'Power real');
    semilogx(distance, fittedcurve_2_5GHz, 'DisplayName', 'Fitted curve');
    xlabel('Distance(m)');
    ylabel('Power(dB)');
    title(['Fitting for 2.5GHz, Gain ', num2str(gain), ' dB']);
    grid on;
    legend('Location', 'best');
    set(gca, 'XScale', 'log');
    
    % Display fitting info on the plot
    x_pos = min(distance) * 10^(0.08*(log10(max(distance)) - log10(min(distance))));
    y_pos = min(PowerByGain2_5GHz(gain/5+1, :)) + 0.08*(max(PowerByGain2_5GHz(gain/5+1, :)) - min(PowerByGain2_5GHz(gain/5+1, :)));
    text(x_pos, y_pos, {['Polynomial: ', char(vpa(poly2sym(poly_coef_fit), 3))], ...
                        ['Propagation Exponent: ', num2str(-poly_coef_fit(end-1))]}, 'FontSize', 8, 'HorizontalAlignment', 'left');
    
    % Fitting for 600MHz
    subplot(2, 1, 2);
    hold on;
    poly_coef_fit         = polyfit(distance_dB, PowerByGain600MHz(gain/5+1, :), poly_ord);  %%% Linear fitting of the coefficient
    fittedcurve_600MHz    = polyval(poly_coef_fit, distance_dB);
    
    % Plot for 600MHz
    semilogx(distance, PowerByGain600MHz(gain/5+1, :), 's-', 'DisplayName', 'Power real');
    semilogx(distance, fittedcurve_600MHz, 'DisplayName', 'Fitted curve');
    xlabel('Distance(m)');
    ylabel('Power(dB)');
    title(['Fitting for 600MHz, Gain ', num2str(gain), ' dB']);
    grid on;
    legend('Location', 'northeast');
    set(gca, 'XScale', 'log');
    
    % Display fitting info on the plot
    x_pos = min(distance) * 10^(0.08*(log10(max(distance)) - log10(min(distance))));
    y_pos = min(PowerByGain600MHz(gain/5+1, :)) + 0.08*(max(PowerByGain600MHz(gain/5+1, :)) - min(PowerByGain600MHz(gain/5+1, :)));
    text(x_pos, y_pos, {['Polynomial: ', char(vpa(poly2sym(poly_coef_fit), 3))], ...
                        ['Propagation Exponent: ', num2str(-poly_coef_fit(end-1))]}, 'FontSize', 8, 'HorizontalAlignment', 'left');
end
