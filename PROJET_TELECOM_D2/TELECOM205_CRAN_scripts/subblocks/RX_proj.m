function [basebandAnalog_adc_I, basebandAnalog_adc_Q] = RX_proj(rxSignal,LNA_IIP3,LNA_NF,LNA_Gain,RXBB_Filt,RXBB_Filt_NF,Flo,continuousTimeSamplingRate,adcSamplingRate,nBitADC,Vref,delay,BB_gain)
    %RX_proj - A shortcut function to emulate the RX chain - TELECOM205 version
    %   Calls the functions rfLNA, downMixer, basebandAnalogFiltFake and ADC
    %
    % Syntax:  [basebandAnalog_adc_I, basebandAnalog_adc_Q]=RX_proj(rxSignal,LNA_IIP3,LNA_NF,LNA_Gain,RXBB_Filt,RXBB_Filt_NF,Flo,continuousTimeSamplingRate,adcSamplingRate,nBitADC,Vref,delay,BB_gain)
    %
    % Inputs:
    %    rxSignal                   - RX signal at antenna port
    %    LNA_IIP3                   - LNA IIP3 (dBm)
    %    LNA_NF                     - LNA NF (dB)
    %    LNA_Gain                   - LNA Gain (dB)
    %    RXBB_Filt                  - RX baseband filter impulse response ; FIR ONLY
    %    RXBB_Filt_NF               - RX baseband filter NF (dB)
    %    Flo                        - LO frequency (Hz)
    %    continuousTimeSamplingRate - simulation step rate (Hz)
    %    adcSamplingRate            - ADC sampling rate (Hz)
    %    nBitADC                    - number of bits of the ADC
    %    Vref                       - reference voltage of the ADC
    %    delay                      - delay to compensate for the filter group delay
    %    BB_gain                    - baseband gain (dB)
    %
    % Outputs:
    %    basebandAnalog_adc_I       - baseband analog equivalent I signal after ADC
    %    basebandAnalog_adc_Q       - baseband analog equivalent Q signal after ADC
    %
    % Other m-files required: rfLNA, downMixer, basebandAnalogFiltFake, ADC
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: none
    % Author: Germain PHAM, Chadi JABBOUR
    % C2S, COMELEC, Telecom Paris, Palaiseau, France
    % email address: dpham@telecom-paris.fr
    % Website: https://c2s.telecom-paristech.fr/TODO
    % Dec. 2023
    %------------- BEGIN CODE --------------



Rin = 50;

rfLNASignal = rfLNA(rxSignal,LNA_Gain,LNA_NF,LNA_IIP3,Rin,continuousTimeSamplingRate/2);

%%% Mixing down to BB %%%
[basebandAnalog_raw_I,basebandAnalog_raw_Q] = downMixer(rfLNASignal,Flo,continuousTimeSamplingRate);

% Perform filtering
disp('Filtering the mixer baseband output - This takes a while...')
basebandAnalog_filtrx_I = basebandAnalogFiltFake(basebandAnalog_raw_I,...
                                                    RXBB_Filt,...
                                                    RXBB_Filt_NF,...
                                                    Rin,...
                                                    continuousTimeSamplingRate);
basebandAnalog_filtrx_Q = basebandAnalogFiltFake(basebandAnalog_raw_Q,...
                                                    RXBB_Filt,...
                                                    RXBB_Filt_NF,...
                                                    Rin,...
                                                    continuousTimeSamplingRate);


BB_gain_lin = 10^(BB_gain/20);

% Perform conversion
basebandAnalog_adc_I = ADC(BB_gain_lin*basebandAnalog_filtrx_I,...
                            nBitADC,...
                            Vref,...
                            adcSamplingRate,...
                            delay,...
                            continuousTimeSamplingRate);
basebandAnalog_adc_Q = ADC(BB_gain_lin*basebandAnalog_filtrx_Q,...
                            nBitADC,...
                            Vref,...
                            adcSamplingRate,...
                            delay,...
                            continuousTimeSamplingRate);


