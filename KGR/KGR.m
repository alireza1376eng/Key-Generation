% APPROVED FOR ARTICLE
clc
clear all
close all

%% Initialization
% Define total power range (logarithmically spaced) and noise variance
p_total = logspace(-1, 10, 100); % Total power in logarithmic scale
varn = 3; % Noise variance

Euler_constant = 0.5772;  %Euler constant
e = 2.71828; % Euler number

% Power allocation for Near and Far users based on NOMA principle
PowerN = 0.2 .* p_total; % 20% of total power allocated to Near user
PowerF = 0.8 .* p_total; % 80% of total power allocated to Far user

%% Calculate SNR and Capacity for Our Results
% Signal-to-Noise Ratio (SNR) for Near and Far users
snr_F = 10 .* log10(PowerF ./ varn); % SNR for Far user
snr_N = 10 .* log10(PowerN ./ varn); % SNR for Near user
snr = 10 .* log10(p_total ./ varn);  % Overall SNR

% Phase-based capacity calculations
% Capacity formula for proposed phase-based systems in paper
phase_capacity_F = (0.5) .* log2((PowerF ./ (varn + PowerN)) + 1) - ((1+Euler_constant)/2)*log2(e) + 0.5*log2(pi)+1; % Far user
phase_capacity_N = (0.5) .* log2((PowerN ./ varn) + 1) - ((1+Euler_constant)/2)*log2(e) + 0.5*log2(pi)+1;            % Near user
sum_capacity = phase_capacity_F + phase_capacity_N; % Total capacity

% Plot the results for our proposed method
plot(snr_F, phase_capacity_F, 'linewidth', 2); % Capacity for Far user
hold on;
plot(snr_N, phase_capacity_N, 'linewidth', 2); % Capacity for Near user
save("SNR_db", "snr"); % Save the calculated SNR for future use

%% Plot Results from Other Papers
% Data extracted from other papers for comparison
artBP_p = [0.6 0.79 1.03 1.25 1.5 2 2.3 2.6 2.95 3.27 3.67]; % Phase-based (Ref. [46])
artBP_ph = [0.26 0.4 0.6 0.8 1.0 1.465 1.74 2.05 2.4 2.75 3.1]; % Phase-based (Ref. [44])
artBP_g = [0.05 0.1 0.18 0.28 0.4 0.71 0.93 1.2 1.5 1.8 2.2];   % Gain-based (Ref. [45])

artBP_SNR = [0 1.8 3.6 5.4 7.2 10.8 12.6 14.4 16.2 18.1 20]; % Corresponding SNR values

% Plot the results from other papers
plot(artBP_SNR, artBP_p, 'linewidth', 2);  % Phase-based (Ref. [46])
plot(artBP_SNR, artBP_ph, 'linewidth', 2); % Phase-based from Ref. [44]
plot(artBP_SNR, artBP_g, 'linewidth', 2);  % Gain-based from Ref. [45]

%% Final Plot Customization
% Set plot title, axis labels, legend, and grid
title('Key Generation Rate', 'FontSize', 17);
xlabel('SNR_m (dB)', 'FontSize', 17); % X-axis: SNR in dB
ylabel('KGR', 'FontSize', 17);        % Y-axis: Key Generation Rate (KGR)

% Add legend
h = legend('Upper Bound-Far', 'Upper Bound-Near', 'Phased Based [46]', ...
           'Phased Based [44]', 'Gain Based [45]');
h.FontSize = 15; % Set legend font size

% Set grid and axis limits
grid on;
xlim([0, 20]); % Limit X-axis to SNR range [0, 20]