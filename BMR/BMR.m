clc
clear all

numLevels = 4; % Quantization levels
iteration = 5000; % Number of iterations for key generation

% Alice BS initial bit (Phased) Generation
bit = (pi/2)*randi([0 3],1)+ (pi/4);
teta_A = angle(exp(1j*bit));    %initial random phase

% channel initialization
Variance_h_N = 10;  %10 dB Near
Variance_h_F = 1;  %0 dB Far
P_total = 10; %10 dB Total Power Allocation
Variance_n = logspace(-3,3,500); % AWGN variance
PowerN = sqrt((0.2*10) / Variance_h_N); %power allocated Near user
PowerF = sqrt((0.8*10) / Variance_h_F); %power allocated Far user

% SNR calculation
SNRF_db = 10 .* log10( 2 ./ Variance_n);
SNRN_db = 10 .* log10( 8 ./ Variance_n);

%% Process for Key Generation and BMR Calculation
for j=1:1:length(Variance_n)
    for n=1:1:iteration

        % Generate channels and noise for Near and Far users
        H_N = sqrt(Variance_h_N/2) * (randn(1) + 1i * randn(1));   %Near User channel
        H_F = sqrt(Variance_h_F/2) * (randn(1) + 1i * randn(1));   %Far User channel
        n_bN = sqrt(Variance_n(j)/2) * (randn(1) + 1i * randn(1)); % AWGN of BOB_N
        n_bF = sqrt(Variance_n(j)/2) * (randn(1) + 1i * randn(1)); % AWGN of BOB_F

        % Alice transmits signals (phase 1)
        transmited_sig_AN =  PowerN * exp(1i*teta_A);
        transmited_sig_AF =  PowerF * exp(1i*teta_A);

        % BoB received signals (phase 1)
        Recieved_sig_BOB_N =  transmited_sig_AN * H_N + n_bN ; % Near User
        Recieved_sig_BOB_F =  transmited_sig_AF * H_F + n_bF ; % Far User

        % BOB signal processing
        teta_hat_BOB_N = angle(Recieved_sig_BOB_N);  % Estimated phase for Near User
        teta_hat_BOB_F = angle(Recieved_sig_BOB_F);  % Estimated phase for Far User

        % Random signal generation for security
        Sig_BN_rand = exp(1i*pi*randi([0 1],1)) * exp(1i * teta_hat_BOB_N);
        Sig_BF_rand = exp(1i*pi*randi([0 1],1)) * exp(1i * teta_hat_BOB_F);

        % Key generation and quantization
        teta_BN_rand  = angle(Sig_BN_rand);
        teta_BF_rand  = angle(Sig_BF_rand); 
        [KEY_BOB_N(:,n),KEY_BOB_phased_N(j,n)] = P_Quantization(teta_BN_rand);
        [KEY_BOB_F(:,n),KEY_BOB_phased_F(j,n)] = P_Quantization(teta_BF_rand);

        %% Alice and Eve Processing
        Variance_h_E = 5;  % Eve's channel variance
        H_E_N = sqrt(Variance_h_E/2) * (randn(1) + 1i * randn(1));
        H_E_F = sqrt(Variance_h_E/2) * (randn(1) + 1i * randn(1));
        n_a = sqrt(Variance_n(j)/2) * (randn(1) + 1i * randn(1));  % AWGN for Alice
        n_e = sqrt(Variance_n(j)/2) * (randn(1) + 1i * randn(1));  % AWGN for Eve

        % BoB transmits signals (phase 2)
        transmited_sig_B_N = PowerN * exp(1i*(KEY_BOB_phased_N(j,n) - teta_hat_BOB_N));
        transmited_sig_B_F = PowerF * exp(1i*(KEY_BOB_phased_F(j,n) - teta_hat_BOB_F));

        % Received signals at Alice and Eve
        Recieved_sig_A_F =  (transmited_sig_B_N * H_N) + (transmited_sig_B_F * H_F) + n_a;
        Recieved_sig_A_N =  (transmited_sig_B_N * H_N) + n_a;
        Recieved_sig_A_EVE_F =  (transmited_sig_B_N * H_E_N) + (transmited_sig_B_F * H_E_F) + n_e;
        Recieved_sig_A_EVE_N =  (transmited_sig_B_N * H_E_N) + n_e;

        % Signal processing at Alice and Eve
        teta_hat_Alice_F = angle(Recieved_sig_A_F);
        teta_hat_Alice_N = angle(Recieved_sig_A_N);
        teta_hat_EVE_F = angle(Recieved_sig_A_EVE_F);
        teta_hat_EVE_N = angle(Recieved_sig_A_EVE_N);

        % Key quantization for Alice and Eve
        [KEY_Alice_F(:,n),KEY_Alice_phased_F] = P_Quantization(angle(exp(1i*teta_A) * exp(1i * teta_hat_Alice_F)));
        [KEY_Alice_N(:,n)] = P_Quantization(angle(exp(1i*teta_A) * exp(1i * teta_hat_Alice_N)));

        [KEY_Alice_E_F(:,n),KEY_Alice_phased_E_F] = P_Quantization(teta_hat_EVE_F);
        [KEY_Alice_E_N(:,n)] = P_Quantization(teta_hat_EVE_N);

    end

    % XOR operation for BMR calculation
    XOR_F = xor(KEY_BOB_F,KEY_Alice_F);
    XOR_N = xor(KEY_BOB_N,KEY_Alice_N);
    XOR_E_F = xor(KEY_BOB_F,KEY_Alice_E_F);
    XOR_E_N = xor(KEY_BOB_N,KEY_Alice_E_N);

    % Calculate error rates
    Errore_F(j) = sum(sum(XOR_F),2)/(iteration*2);
    Errore_N(j) = sum(sum(XOR_N),2)/(iteration*2);
    Errore_E_F(j) = sum(sum(XOR_E_F),2)/(iteration*2);
    Errore_E_N(j) = sum(sum(XOR_E_N),2)/(iteration*2);
end

%% Plotting Results

% other articls results for comparison
Channel_gain = [0.34 0.28 0.25 0.2 0.17 0.135 0.105 0.090 0.07 0.06 0.05];
Channel_phased = [0.405 0.373 0.35 0.305 0.27 0.225 0.180 0.15 0.118 0.10 0.075];

Direct_Forensics = [0.29 0.262 0.23 0.19 0.169 0.13 0.11 0.085 0.063 0.056 0.047];
Relay_Forensics = [0.34 0.312 0.28 0.24 0.209 0.166 0.142 0.115 0.09 0.08 0.070];

P_db = [0 2.5 4 6.4 8 10 12.4 14 16.1 17.7 20];
P_db_Forensics = [5 7 9 11 12.5 15 17 19.2 21.5 22.6 25];

plot(P_db,Channel_gain,'linewidth',2)
hold on
plot(P_db,Channel_phased,'linewidth',2)
plot(P_db_Forensics,Direct_Forensics,'linewidth',2)
plot(P_db_Forensics,Relay_Forensics,'linewidth',2)

% Proposed methods
erroreF = movmean(Errore_F, 5);
erroreN = movmean(Errore_N, 5);
erroreEF = movmean(Errore_E_F, 3);
erroreEN = movmean(Errore_E_N, 3);
erroreEF = movmean(erroreEF, 3);
erroreEN = movmean(erroreEN, 3);


% Our results

plot(SNRF_db ,erroreF,'linewidth',2)
plot(SNRN_db ,erroreN,'linewidth',2)
plot(SNRF_db ,erroreEF,'linewidth',2)
plot(SNRN_db ,erroreEN,'linewidth',2)

xlim([0 25])
title('Bit Mismatch Rate')
xlabel('SNR_m(db)','FontSize',17)
ylabel('BMR','FontSize',17)
h=legend('phased based [44]','Relay based SKG [40]','Direct SKG [40]','gain based [45]','BMR-Far','BMR-Near','errore-EF','errore-EN');
h.FontSize=15;
grid on


%% Binary Quantization Function
function [Qbit,KEY_phased] = P_Quantization(phase)
    normalizedPhase = mod(phase, 2*pi); % Normalize phase to [0, 2*pi)
    if normalizedPhase >= 0 && normalizedPhase < pi/2 
        Qbit = [0 0];
        KEY_phased = (pi/4);
    elseif normalizedPhase >= pi/2 && normalizedPhase < pi
        Qbit = [0 1];
        KEY_phased = (3*pi/4);
    elseif normalizedPhase >= -pi && normalizedPhase < -pi/2 || normalizedPhase == pi
        Qbit = [1 1];
        KEY_phased = (-3*pi/4);
    else % normalizedPhase >= -pi/2 && normalizedPhase < 0
        Qbit = [1 0];
        KEY_phased = (-pi/4);
    end
end