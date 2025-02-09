% Clear workspace and initialize
clc
clear all

% Initialize parameters
numLevels = 4; % Number of quantization levels
iteration = 12800; % Number of iterations for key generation
total_KEY_Alice_N = zeros(1, 2 * iteration);
total_KEY_BOB_N = zeros(1, 2 * iteration);

% Alice's initial phase bit generation (phased generation)
bit = (pi / 2) * randi([0, 3], 1) + (pi / 4); % Random phase initialization
teta_A = angle(exp(1j * bit)); % Compute initial phase

% Channel parameters
Variance_h_N = 10;  % Variance for Near User (10 dB)
Variance_h_F = 1;   % Variance for Far User (0 dB)
P_total = 10;       % Total power (10 dB)
Variance_n = logspace(-3, 3, 500); % Noise variance array
PowerN = sqrt((0.2 * P_total) / Variance_h_N); % Power allocation for Near User
PowerF = sqrt((0.8 * P_total) / Variance_h_F); % Power allocation for Far User
PowerT = PowerN + PowerF; % Total power

%% Process for key generation
for j = 1:length(Variance_n)
    for n = 1:iteration

        % Generate random channels for Near and Far users
        H_N = sqrt(Variance_h_N / 2) * (randn(1) + 1i * randn(1)); % Near User channel
        H_F = sqrt(Variance_h_F / 2) * (randn(1) + 1i * randn(1)); % Far User channel

        % Additive White Gaussian Noise (AWGN) for BOB
        n_bN = sqrt(Variance_n(j) / 2) * (randn(1) + 1i * randn(1));
        n_bF = sqrt(Variance_n(j) / 2) * (randn(1) + 1i * randn(1));

        % Alice transmits signals (phase 1)
        transmited_sig_AN = PowerN * exp(1i * teta_A);
        transmited_sig_AF = PowerF * exp(1i * teta_A);

        % Signals received by BOB
        Recieved_sig_BOB_N = transmited_sig_AN * H_N + n_bN; % Near User received signal
        Recieved_sig_BOB_F = transmited_sig_AF * H_F + n_bF; % Far User received signal

        % Signal processing at BOB
        teta_hat_BOB_N = angle(Recieved_sig_BOB_N); % Estimated phase for Near User
        teta_hat_BOB_F = angle(Recieved_sig_BOB_F); % Estimated phase for Far User

        % Random signal generation for security
        Sig_BN_rand = exp(1i * pi * randi([0, 1], 1)) * exp(1i * teta_hat_BOB_N);
        Sig_BF_rand = exp(1i * pi * randi([0, 1], 1)) * exp(1i * teta_hat_BOB_F);

        % Phase extraction
        teta_BN_rand = angle(Sig_BN_rand);
        teta_BF_rand = angle(Sig_BF_rand);

        % Key generation and quantization for BOB
        [KEY_BOB_N(:, n), KEY_BOB_phased_N] = P_Quantization(teta_BN_rand);
        [KEY_BOB_F(:, n), KEY_BOB_phased_F] = P_Quantization(teta_BF_rand);


        

        %% Alice's reception and key generation (phase 2)

        n_a = sqrt(Variance_n(j) / 2) * (randn(1) + 1i * randn(1)); % AWGN for Alice

        % BOB transmits signal
        transmited_sig_B_N = PowerN * exp(1i * (KEY_BOB_phased_N - teta_hat_BOB_N));
        transmited_sig_B_F = PowerF * exp(1i * (KEY_BOB_phased_F - teta_hat_BOB_F));

        % Signals received by Alice
        Recieved_sig_A_F = (transmited_sig_B_N * H_N) + (transmited_sig_B_F * H_F) + n_a;
        Recieved_sig_A_N = (transmited_sig_B_N * H_N) + n_a;

        % Alice signal processing
        teta_hat_Alice_F = angle(Recieved_sig_A_F);
        teta_hat_Alice_N = angle(Recieved_sig_A_N);
        Sig_A_nrand_F = exp(1i * teta_A) * exp(1i * teta_hat_Alice_F);
        Sig_A_nrand_N = exp(1i * teta_A) * exp(1i * teta_hat_Alice_N);

        % Key prediction and quantization at Alice
        teta_A_nrand_F = angle(Sig_A_nrand_F);
        teta_A_nrand_N = angle(Sig_A_nrand_N);
        [KEY_Alice_F(:, n), KEY_Alice_phased_F] = P_Quantization(teta_A_nrand_F);
        [KEY_Alice_N(:, n)] = P_Quantization(teta_A_nrand_N);

        %% EVE's reception and key generation (phase 2)

        Variance_h_E = 5;  % Eve's channel variance
        H_E_N = sqrt(Variance_h_E/2) * (randn(1) + 1i * randn(1));
        H_E_F = sqrt(Variance_h_E/2) * (randn(1) + 1i * randn(1));
        n_e = sqrt(Variance_n(j)/2) * (randn(1) + 1i * randn(1));  % AWGN for Eve

        % Received signals at Eve
        Recieved_sig_A_EVE_F =  (transmited_sig_B_N * H_E_N) + (transmited_sig_B_F * H_E_F) + n_e;
        Recieved_sig_A_EVE_N =  (transmited_sig_B_N * H_E_N) + n_e;

        % Signal processing at Eve
        teta_hat_EVE_F = angle(Recieved_sig_A_EVE_F);
        teta_hat_EVE_N = angle(Recieved_sig_A_EVE_N);

        % Key quantization for Alice and Eve
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



    % Combine keys for evaluation
    preliminary_KEY_BOB_N(j, 1:2:(2 * iteration) - 1) = KEY_BOB_N(1, :);
    preliminary_KEY_BOB_N(j, 2:2:2 * iteration) = KEY_BOB_N(2, :);
    preliminary_KEY_Alice_N(j, 1:2:(2 * iteration) - 1) = KEY_Alice_N(1, :);
    preliminary_KEY_Alice_N(j, 2:2:2 * iteration) = KEY_Alice_N(2, :);
    preliminary_KEY_BOB_F(j, 1:2:(2 * iteration) - 1) = KEY_BOB_F(1, :);
    preliminary_KEY_BOB_F(j, 2:2:2 * iteration) = KEY_BOB_F(2, :);
    preliminary_KEY_Alice_F(j, 1:2:(2 * iteration) - 1) = KEY_Alice_F(1, :);
    preliminary_KEY_Alice_F(j, 2:2:2 * iteration) = KEY_Alice_F(2, :);
    preliminary_KEY_EVE_N(j, 1:2:(2 * iteration) - 1) = KEY_Alice_E_N(1, :);
    preliminary_KEY_EVE_N(j, 2:2:2 * iteration) = KEY_Alice_E_N(2, :);
    preliminary_KEY_EVE_F(j, 1:2:(2 * iteration) - 1) = KEY_Alice_E_F(1, :);
    preliminary_KEY_EVE_F(j, 2:2:2 * iteration) = KEY_Alice_E_F(2, :);
end


%% Plotting results


SNRF_db = 10 .* log10( 2 ./ Variance_n);
SNRN_db = 10 .* log10( 8 ./ Variance_n);

save("preliminary_KEY_BOB_N", "preliminary_KEY_BOB_N");
save("preliminary_KEY_Alice_N", "preliminary_KEY_Alice_N");
save("preliminary_KEY_BOB_F", "preliminary_KEY_BOB_F");
save("preliminary_KEY_Alice_F", "preliminary_KEY_Alice_F");
save("preliminary_KEY_EVE_N", "preliminary_KEY_EVE_N");
save("preliminary_KEY_EVE_F", "preliminary_KEY_EVE_F");

save("SNRN_db", "SNRN_db");
save("SNRF_db", "SNRF_db");

%% Binary Quantization Function
function [Qbit, KEY_phased] = P_Quantization(phase)
    normalizedPhase = mod(phase, 2 * pi); % Ensure phase is within [0, 2*pi)

    if normalizedPhase >= 0 && normalizedPhase < pi / 2
        Qbit = [0 0];
        KEY_phased = (pi / 4);
    elseif normalizedPhase >= pi / 2 && normalizedPhase < pi
        Qbit = [0 1];
        KEY_phased = (3 * pi / 4);
    elseif normalizedPhase >= -pi && normalizedPhase < -pi / 2 || normalizedPhase == pi
        Qbit = [1 1];
        KEY_phased = (-3 * pi / 4);
    else % normalizedPhase >= -pi / 2 && normalizedPhase < 0
        Qbit = [1 0];
        KEY_phased = (-pi / 4);
    end
end
