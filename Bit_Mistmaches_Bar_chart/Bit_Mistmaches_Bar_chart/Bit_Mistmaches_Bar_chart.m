clc
clear all

%% Load the generated keys for Near User and Far User Form Preliminary Key Generation Folder
 
load("preliminary_KEY_BOB_N.mat");
load("preliminary_KEY_BOB_F.mat");
load("preliminary_KEY_Alice_N.mat");
load("preliminary_KEY_Alice_F.mat");

% Initialize matrices to store results for Near User (N) and Far User (F)
Numb_N = zeros(9,100);
Numb_F = zeros(9,100);

% Select keys corresponding to a specific Signal-to-Noise Ratio (SNR = 0 dB) point
K_BOB_N = preliminary_KEY_BOB_N(182,:);
K_Alice_N = preliminary_KEY_Alice_N(182,:);
K_BOB_F = preliminary_KEY_BOB_F(182,:);
K_Alice_F = preliminary_KEY_Alice_F(182,:);

% Iterate over 100 packets
for m = 1:1:100
    % Extract and reshape 256-bit keys into 8x32 matrices for both users
    T_BOB_N = K_BOB_N(1, 256*(m-1)+1:256*m);
    T_Alice_N = K_Alice_N(1, 256*(m-1)+1:256*m);
    C_BOB_N = reshape(T_BOB_N, 8, 32);
    C_Alice_N = reshape(T_Alice_N, 8, 32);

    T_BOB_F = K_BOB_F(1, 256*(m-1)+1:256*m);
    T_Alice_F = K_Alice_F(1, 256*(m-1)+1:256*m);
    C_BOB_F = reshape(T_BOB_F, 8, 32);
    C_Alice_F = reshape(T_Alice_F, 8, 32);

    % Compute the number of bit mismatches (syndrome discrepancy)
    Syndrome_discrepancy_N = sum(xor(C_BOB_N, C_Alice_N));
    Syndrome_discrepancy_F = sum(xor(C_BOB_F, C_Alice_F));

    % Classify mismatches for each packet into categories (0-8 mismatches)
    for t = 1:1:32
        % Near User mismatch classification
        if Syndrome_discrepancy_N(1,t) == 0
           Numb_N(1,m) =  Numb_N(1,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 1
           Numb_N(2,m) =  Numb_N(2,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 2
           Numb_N(3,m) =  Numb_N(3,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 3
           Numb_N(4,m) =  Numb_N(4,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 4
           Numb_N(5,m) =  Numb_N(5,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 5
           Numb_N(6,m) =  Numb_N(6,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 6
           Numb_N(7,m) =  Numb_N(7,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 7
           Numb_N(8,m) =  Numb_N(8,m) + 1;
        elseif Syndrome_discrepancy_N(1,t) == 8
           Numb_N(9,m) =  Numb_N(9,m) + 1;
        end

        % Far User mismatch classification
        if Syndrome_discrepancy_F(1,t) == 0
           Numb_F(1,m) =  Numb_F(1,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 1
           Numb_F(2,m) =  Numb_F(2,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 2
           Numb_F(3,m) =  Numb_F(3,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 3
           Numb_F(4,m) =  Numb_F(4,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 4
           Numb_F(5,m) =  Numb_F(5,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 5
           Numb_F(6,m) =  Numb_F(6,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 6
           Numb_F(7,m) =  Numb_F(7,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 7
           Numb_F(8,m) =  Numb_F(8,m) + 1;
        elseif Syndrome_discrepancy_F(1,t) == 8
           Numb_F(9,m) =  Numb_F(9,m) + 1;
        end
    end 
end

% Compute the mean number of mismatches for each category (0-4 mismatches)
N1 = mean(Numb_N(1,:));
N2 = mean(Numb_N(2,:));
N3 = mean(Numb_N(3,:));
N4 = mean(Numb_N(4,:));
N5 = mean(Numb_N(5,:));

F1 = mean(Numb_F(1,:));
F2 = mean(Numb_F(2,:));
F3 = mean(Numb_F(3,:));
F4 = mean(Numb_F(4,:));
F5 = mean(Numb_F(5,:));

% Data for the plot
x = categorical({'0', '1', '2', '3', '4+'}); % Symbolic x-coordinates
x = reordercats(x, {'0', '1', '2', '3', '4+'}); % Preserve order
near_user = [N1, N2, N3, N4, N5]; % Near User data
far_user = [F1, F2, F3, F4, F5];   % Far User data

% Create the bar chart
figure;
bar(x, [near_user; far_user]', 'grouped');
ylim([0 32]); % Set y-axis limits

% Customize axes labels and legend
ylabel('Number of Packets');
xlabel('Number of Bit Mismatches');
legend({'Near User', 'Far User'}, 'Location', 'northoutside', 'Orientation', 'horizontal');

% Set the figure properties
set(gca, 'FontSize', 12); % Adjust font size for better readability
title('Number of Packets with Bit Mismatches');
grid on; % Add grid for better visibility