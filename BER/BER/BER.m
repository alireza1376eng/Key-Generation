clc
clear all

%% Loading Generated Keys
% Load pre-saved preliminary keys and SNR data for reconciliation
load("preliminary_KEY_EVE_F.mat");
load("preliminary_KEY_EVE_N.mat");
load("preliminary_KEY_BOB_N.mat");
load("preliminary_KEY_BOB_F.mat");
load("preliminary_KEY_Alice_N.mat");
load("preliminary_KEY_Alice_F.mat");
load("SNRN_db.mat");
load("SNRF_db.mat");

% save("SNRN_db", "SNRN_db");
% save("SNRF_db", "SNRF_db");


%% Initialization
iteration = 100; % Number of blocks to process in each key set
lenght = 500;    % Number of key sets to process

% Parity-check matrix for error detection
H = [1 0 1 1 1 0 0 0; 
     1 1 1 0 0 1 0 0; 
     0 1 1 1 0 0 1 0; 
     1 1 1 1 1 1 1 1];

% Predefined error patterns for single-bit error correction
erorre_tabel = [0 0 0 0 0 0 0 1; 
                0 0 0 0 0 0 1 0; 
                0 0 0 0 0 1 0 0; 
                0 0 0 0 1 0 0 0; 
                0 0 0 1 0 0 0 0; 
                0 0 1 0 0 0 0 0; 
                0 1 0 0 0 0 0 0; 
                1 0 0 0 0 0 0 0];

% Compute syndrome table for error detection and correction
syndrome_table = mod(erorre_tabel * H', 2);

%% Processing Each Key Set
for k = 1:lenght
    % Extract keys for Bob, Alice, and Eve
    K_BOB_N = preliminary_KEY_BOB_N(k, :);
    K_Alice_N = preliminary_KEY_Alice_N(k, :);
    K_BOB_F = preliminary_KEY_BOB_F(k, :);
    K_Alice_F = preliminary_KEY_Alice_F(k, :);
    K_EVE_N = preliminary_KEY_EVE_N(k, :);
    K_EVE_F = preliminary_KEY_EVE_F(k, :);

    % Process each block of 256 bits in the current key set
    for m = 1:iteration
        % Extract 256-bit blocks for each participant
        T_BOB_N = K_BOB_N(1, 256*(m-1)+1:256*m);
        T_Alice_N = K_Alice_N(1, 256*(m-1)+1:256*m);
        T_BOB_F = K_BOB_F(1, 256*(m-1)+1:256*m);
        T_Alice_F = K_Alice_F(1, 256*(m-1)+1:256*m);
        T_EVE_N = K_EVE_N(1, 256*(m-1)+1:256*m);
        T_EVE_F = K_EVE_F(1, 256*(m-1)+1:256*m);

        % Reshape blocks into 8x32 matrices for error correction
        C_BOB_N = reshape(T_BOB_N, 8, 32)';
        C_Alice_N = reshape(T_Alice_N, 8, 32)';
        C_BOB_F = reshape(T_BOB_F, 8, 32)';
        C_Alice_F = reshape(T_Alice_F, 8, 32)';
        C_EVE_N = reshape(T_EVE_N, 8, 32)';
        C_EVE_F = reshape(T_EVE_F, 8, 32)';

        % Compute syndromes for error detection
        Syndrome_BOB_N = mod(C_BOB_N * H', 2);
        Syndrome_Alice_N = mod(C_Alice_N * H', 2);
        Syndrome_BOB_F = mod(C_BOB_F * H', 2);
        Syndrome_Alice_F = mod(C_Alice_F * H', 2);
        Syndrome_EVE_N = mod(C_EVE_N * H', 2);
        Syndrome_EVE_F = mod(C_EVE_F * H', 2);

        % Compute discrepancies between syndromes
        Syndrome_discrepancy_N = xor(Syndrome_BOB_N, Syndrome_Alice_N);
        Syndrome_discrepancy_F = xor(Syndrome_BOB_F, Syndrome_Alice_F);
        Syndrome_discrepancy_EVE_N = xor(Syndrome_BOB_N, Syndrome_EVE_N);
        Syndrome_discrepancy_EVE_F = xor(Syndrome_BOB_F, Syndrome_EVE_F);

        % Initialize error patterns
        erorre_patern_N = zeros(32, 8);
        erorre_patern_F = zeros(32, 8);
        erorre_patern_EVE_N = zeros(32, 8);
        erorre_patern_EVE_F = zeros(32, 8);

        % Correct errors using syndrome table
        for t = 1:32
            for idx = 1:8
                if isequal(Syndrome_discrepancy_N(t, :), syndrome_table(idx, :))
                    erorre_patern_N(t, :) = erorre_tabel(idx, :);
                end
                if isequal(Syndrome_discrepancy_EVE_N(t, :), syndrome_table(idx, :))
                    erorre_patern_EVE_N(t, :) = erorre_tabel(idx, :);
                end
                if isequal(Syndrome_discrepancy_F(t, :), syndrome_table(idx, :))
                    erorre_patern_F(t, :) = erorre_tabel(idx, :);
                end
                if isequal(Syndrome_discrepancy_EVE_F(t, :), syndrome_table(idx, :))
                    erorre_patern_EVE_F(t, :) = erorre_tabel(idx, :);
                end
            end
        end

        % Reconcile keys using corrected patterns
        Reconciled_N = xor(C_Alice_N, erorre_patern_N);
        Reconciled_F = xor(C_Alice_F, erorre_patern_F);
        Reconciled_EVE_N = xor(C_EVE_N, erorre_patern_EVE_N);
        Reconciled_EVE_F = xor(C_EVE_F, erorre_patern_EVE_F);

        % Check for mismatches after reconciliation
        Check_N = xor(Reconciled_N, C_BOB_N);
        Check_F = xor(Reconciled_F, C_BOB_F);
        Check_EVE_N = xor(Reconciled_EVE_N, C_BOB_N);
        Check_EVE_F = xor(Reconciled_EVE_F, C_BOB_F);

        % Compute error rates
        Errore_N(m) = sum(sum(Check_N,2)) / 256;
        Errore_F(m) = sum(sum(Check_F,2)) / 256;
        Errore_EVE_N(m) = sum(Check_EVE_N(:)) / 256;
        Errore_EVE_F(m) = sum(Check_EVE_F(:)) / 256;
    end

    % Average mismatch rates over iterations
    BER_F(k) = mean(Errore_F);
    BER_N(k) = mean(Errore_N);
    BER_EVE_F(k) = mean(Errore_EVE_F);
    BER_EVE_N(k) = mean(Errore_EVE_N);
end

%% Plot Results

% our Results
plot(SNRN_db, BER_N, 'linewidth', 2);
hold on;
plot(SNRF_db, BER_F, 'linewidth', 2);
plot(SNRF_db, BER_EVE_F, 'linewidth', 2);
plot(SNRN_db, BER_EVE_N, 'linewidth', 2);


% same articles results
BER_Fornsics_Direct = [0.005 0 0 0 0 0 0 0 0 0 0];
BER_Fornsics_Relay =  [0.01 0 0 0 0 0 0 0 0 0 0];
P_db_Forensics = [5 7 9 11 12.5 15 17 19.2 21.5 22.6 25];

plot(P_db_Forensics,BER_Fornsics_Direct,'linewidth',2)
plot(P_db_Forensics,BER_Fornsics_Relay,'linewidth',2)

xlim([5 25])

% Add plot labels and legend
title('Key Mismatch Rate');
xlabel('SNR (dB)', 'FontSize', 17);
ylabel('Mismatch Rate', 'FontSize', 17);
legend('Near User', 'Far User', 'Eve (Far)', 'Eve (Near)', 'Direc_SKG[40]', 'Relay_based_SKG[40]' );
grid on;