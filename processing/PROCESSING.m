%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — Processing Script
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: PROCESSING.m
%
% Description:
% This script estimates the CSI (Channel State Information), visualises amplitude 
% variations accross subcarriers and time, and performs the following
% tasks: CSI matrix dimensionality reduction via PCA followed by multiclass
% classification using SVM; extraction of 22 CSI-based features;
% cross-correlation analysis; and relative power evaluation of received signals. 
% The goal is to evaluate the detectability of human presence based on CSI variations 
% and other statistical processess.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. OFDM Parameters and Transmitted Signal Generation
clear; close all; clc

Nfft = 208;                         % FFT Size
Ncp = 13;                            % Cyclic prefix length
ofdmSampleRate = 25e6;              % Sampling rate (Hz)
NofdmSyms = 100;                    % Number of OFDM symbols in the frame
scs = ofdmSampleRate/(1e3*Nfft);    % Subcarrier spacing (kHz)

N_zc = 199;
[txSignal, dataGrid, nullIndices] = generation_txSignalZC(Nfft,Ncp,ofdmSampleRate,N_zc,NofdmSyms);
type = 'ZC';

%N_crs = 207;
%[txSignal, dataGrid, nullIndices] = generation_txSignalCRS(Nfft,Ncp,ofdmSampleRate,N_crs,NofdmSyms);
%type = 'CRS';

%N_golay = 2^(floor(log2(Nfft)));
%[txSignal, dataGrid, nullIndices] = generation_txSignalGolay(Nfft, Ncp, ofdmSampleRate, N_golay, NofdmSyms);
%type = 'Golay';

%% 2. Loading and Processing Received Signals

folder_empty = 'empty';
folder_1subject = '1subject';
folder_2subject = '2subject';
folder_3subject = '3subject';

archivos_empty = fullfile(folder_empty, {dir(fullfile(folder_empty, '*.mat')).name});
archivos_1subject = fullfile(folder_1subject, {dir(fullfile(folder_1subject, '*.mat')).name});
archivos_2subject = fullfile(folder_2subject, {dir(fullfile(folder_2subject, '*.mat')).name});
archivos_3subject = fullfile(folder_3subject, {dir(fullfile(folder_3subject, '*.mat')).name});

CSI_est_empty = cell(1, numel(archivos_empty));
CSI_est_1subject = cell(1, numel(archivos_1subject));
CSI_est_2subject = cell(1, numel(archivos_2subject));
CSI_est_3subject = cell(1, numel(archivos_3subject));


for i = 1:numel(archivos_empty)
    rxSignal_empty = load(archivos_empty{i}).IQ;
    [appearances_empty, rxGrids_empty, ~] = demodulation(txSignal, rxSignal_empty, Nfft, Ncp, nullIndices, NofdmSyms, type);
    CSI_est_empty{i} = CSI(rxGrids_empty, dataGrid, NofdmSyms, type);
end

for i = 1:numel(archivos_1subject)
    rxSignal_1subject = load(archivos_1subject{i}).IQ;
    [appearances_1subject, rxGrids_1subject, ~] = demodulation(txSignal, rxSignal_1subject, Nfft, Ncp, nullIndices, NofdmSyms, type);
    CSI_est_1subject{i} = CSI(rxGrids_1subject, dataGrid, NofdmSyms, type);
end

for i = 1:numel(archivos_2subject)
    rxSignal_2subject = load(archivos_2subject{i}).IQ;
    [appearances_2subject, rxGrids_2subject, ~] = demodulation(txSignal, rxSignal_2subject, Nfft, Ncp, nullIndices, NofdmSyms, type);
    CSI_est_2subject{i}= CSI(rxGrids_2subject, dataGrid, NofdmSyms, type);
end

for i = 1:numel(archivos_3subject)
    rxSignal_3subject = load(archivos_3subject{i}).IQ;
    [appearances_3subject, rxGrids_3subject, ~] = demodulation(txSignal, rxSignal_3subject, Nfft, Ncp, nullIndices, NofdmSyms, type);
    CSI_est_3subject{i} = CSI(rxGrids_3subject, dataGrid, NofdmSyms, type);
end
%% 3. Construction of CSI Matrices
CSI_empty_complex = cell2mat(CSI_est_empty');  % [subcarriers x symbols]
CSI_1subject_complex = cell2mat(CSI_est_1subject');
CSI_2subject_complex = cell2mat(CSI_est_2subject');
CSI_3subject_complex = cell2mat(CSI_est_3subject');

CSI_empty = abs(cell2mat(CSI_est_empty'));  % [subcarriers x symbols]
CSI_1subject = abs(cell2mat(CSI_est_1subject'));
CSI_2subject = abs(cell2mat(CSI_est_2subject'));
CSI_3subject = abs(cell2mat(CSI_est_3subject'));

%% 4. Basic Visualisation
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize',12);

figure;
subplot(4,1,1);
plot(CSI_empty');
title('|CSI| - Empty'); xlabel('Symbol'); ylabel('|CSI|');
grid on;

subplot(4,1,2);
plot(CSI_1subject');
title('|CSI| - 1 Subject'); xlabel('Symbol'); ylabel('|CSI|');
grid on;

subplot(4,1,3);
plot(CSI_2subject');
title('|CSI| - 2 Subjects'); xlabel('Symbol'); ylabel('|CSI|');
grid on;

subplot(4,1,4);
plot(CSI_3subject');
title('|CSI| - 3 Subjects'); xlabel('Symbol'); ylabel('|CSI|');
grid on;

%% 5. Temporal Analysis (Per Subcarrier)
sc_example = round(size(CSI_empty,1)/2);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize',12);

figure;
subplot(4,1,1);
plot(CSI_empty(sc_example,:), 'b'); title(sprintf('Subcarrier %d - Empty', sc_example));
xlabel('Symbol'); ylabel('|CSI|'); grid on;

subplot(4,1,2);
plot(CSI_1subject(sc_example,:), 'r'); title(sprintf('Subcarrier %d - 1 Subject', sc_example));
xlabel('Symbol'); ylabel('|CSI|'); grid on;

subplot(4,1,3);
plot(CSI_2subject(sc_example,:), 'r'); title(sprintf('Subcarrier %d - 2 Subjects', sc_example));
xlabel('Symbol'); ylabel('|CSI|'); grid on;

subplot(4,1,4);
plot(CSI_3subject(sc_example,:), 'r'); title(sprintf('Subcarrier %d - 3 Subjects', sc_example));
xlabel('Symbol'); ylabel('|CSI|'); grid on;
%% 6. Heatmap Visualisation
figure;
imagesc([CSI_empty, CSI_1subject, CSI_2subject, CSI_3subject]); axis xy;
xlabel('Symbol'); ylabel('Subcarrier');
title('Heatmap of |CSI| (Empty | 1 Subject | 2 Subjects | 3 Subjects)');
colorbar;

%% 7. Histograms and CDFs
all_empty = CSI_empty(:);
all_1subject = CSI_1subject(:);
all_2subject = CSI_2subject(:);
all_3subject = CSI_3subject(:);

figure;
histogram(all_empty, 50, 'Normalization','pdf', 'FaceColor','b'); hold on;
histogram(all_1subject, 50, 'Normalization','pdf', 'FaceColor','r');
histogram(all_2subject, 50, 'Normalization','pdf', 'FaceColor','g');
histogram(all_3subject, 50, 'Normalization','pdf', 'FaceColor','y');
xlabel('|CSI|'); ylabel('Density'); grid on;
title('Histogram of |CSI| amplitude');
legend('Empty','1 Subject', '2 Subjects', '3 Subjects');

figure;
cdfplot(all_empty); hold on; cdfplot(all_1subject); cdfplot(all_2subject); cdfplot(all_3subject);
xlabel('|CSI|'); ylabel('Cumulative Probability');
title('CDF of |CSI| amplitude'); grid on;
legend('Empty','1 Subject', '2 Subjects', '3 Subjects');
%% 8. PCA + Clasification

response = '';

while ~strcmpi(response, 'Yes')
    response = input('Have you saved all the previous figures? Type "Yes" to continue: ', 's');
end

close all; clc;
disp('Continuing with execution...');

% ---------- Initial Setup ----------
CSI_total = [CSI_empty, CSI_1subject, CSI_2subject, CSI_3subject];
labels = [ones(1, size(CSI_empty, 2)), 2*ones(1, size(CSI_1subject, 2)), 3*ones(1, size(CSI_2subject, 2)), 4*ones(1, size(CSI_3subject, 2))];  % Labels: 1=Empty, 2=1 Subject, 3=2 Subjects, 4=3 Subjects

% ---------- PCA on Raw CSI (by subcarrier) ----------
% CSI size: [subcarriers x symbols]
CSI_norm_sub = zscore(CSI_total, 0, 2);          % Normalize along each subcarrier (rows)
[coeff1, score1, latent1] = pca(CSI_norm_sub'); % PCA on transposed matrix (samples x features)

varExplained = latent1 / sum(latent1) * 100;     % Percentage of variance explained

% Plot cumulative variance explained
figure;
plot(cumsum(varExplained), 'LineWidth', 1.5);
xlabel('Principal Components');
ylabel('Cumulative Variance Explained (%)');
title('Variance Explained by PCA (Raw CSI by Subcarrier)');
grid on;

% ---------- Multiclass SVM Classification ----------
[trainIdx, testIdx] = crossvalind('HoldOut', labels, 0.2);
% Use ECOC model for multiclass SVM
SVM_model1 = fitcecoc(score1(trainIdx, 1:2), labels(trainIdx));
predicted_labels1 = predict(SVM_model1, score1(testIdx, 1:2));
confMat1 = confusionmat(labels(testIdx), predicted_labels1);

classNames = ["Empty"; "1 Subject"; "2 Subjects"; "3 Subjects"];

% ---------- Visualization ----------
figure;
% Colors: 'rbgy' = red (Empty), blue (1 Subject), green (2 Subjects), yellow (3 Subjects)
gscatter(score1(:,1), score1(:,2), labels, 'rbgy', 'o', 8);
set(gca, 'FontName', 'Times New Roman');
xlabel(sprintf('Principal Component 1 (%.1f%% variance)', varExplained(1)), 'FontName', 'Times New Roman');
ylabel(sprintf('Principal Component 2 (%.1f%% variance)', varExplained(2)), 'FontName', 'Times New Roman');
title('\it{Class Distribution in PCA Space}', 'FontName', 'Times New Roman');
grid on;

% Plot centroids
hold on;
centroidsX = [mean(score1(labels==1,1)), mean(score1(labels==2,1)), mean(score1(labels==3,1)), mean(score1(labels==4,1))];
centroidsY = [mean(score1(labels==1,2)), mean(score1(labels==2,2)), mean(score1(labels==3,2)), mean(score1(labels==4,2))];
scatter(centroidsX, centroidsY, 100, 'kx', 'LineWidth', 2);
offset = 0.5;
for k = 1:4
    text(centroidsX(k) + offset, centroidsY(k), classNames(k), 'FontWeight', 'bold', 'FontName', 'Times New Roman');
end
hold off;

legend('Empty','1 Subject','2 Subjects','3 Subjects','Centroid','Location','best', 'FontName', 'Times New Roman');

% ---------- Numerical Summary ----------
PC1_means = [mean(score1(labels==1,1)); mean(score1(labels==2,1)); mean(score1(labels==3,1)); mean(score1(labels==4,1))];
PC2_means = [mean(score1(labels==1,2)); mean(score1(labels==2,2)); mean(score1(labels==3,2)); mean(score1(labels==4,2))];
summaryTable = table(classNames, PC1_means, PC2_means);
disp('--- Numerical Summary (Method A: Raw CSI by Subcarrier) ---');
disp(summaryTable);
disp('--- Confusion Matrix (Raw CSI by Subcarrier) ---');
disp(confMat1);

% Compute Euclidean distances between centroids
dist12 = sqrt((centroidsX(1)-centroidsX(2))^2 + (centroidsY(1)-centroidsY(2))^2);
dist13 = sqrt((centroidsX(1)-centroidsX(3))^2 + (centroidsY(1)-centroidsY(3))^2);
dist14 = sqrt((centroidsX(1)-centroidsX(4))^2 + (centroidsY(1)-centroidsY(4))^2);

dist23 = sqrt((centroidsX(2)-centroidsX(3))^2 + (centroidsY(2)-centroidsY(3))^2);
dist24 = sqrt((centroidsX(2)-centroidsX(4))^2 + (centroidsY(2)-centroidsY(4))^2);

dist34 = sqrt((centroidsX(3)-centroidsX(4))^2 + (centroidsY(3)-centroidsY(4))^2);

fprintf('\n--- Euclidean Distances Between Centroids ---\n');
fprintf('Empty - 1 Subject: %.3f\n', dist12);
fprintf('Empty - 2 Subjects: %.3f\n', dist13);
fprintf('Empty - 3 Subjects: %.3f\n', dist14);
fprintf('1 Subject - 2 Subjects: %.3f\n', dist23);
fprintf('1 Subject - 3 Subjects: %.3f\n', dist24);
fprintf('2 Subjects - 3 Subjects: %.3f\n', dist34);


%% 9. Feature Extraction via Pseudo Super-Resolution by Scenario
% From the CSI matrices (subcarriers x symbols), extract:
%   1. The top P principal eigenvalues (multipath profile and temporal variation)
%   2. Statistics of the magnitude entropy (channel energy)
%   3. Statistics of phase difference between subcarriers
%   4. Minor channel variation metrics (variance/energy ratio)
% Output: vector of 22 features per sample

response = '';

while ~strcmpi(response, 'Yes')
    response = input('Have you saved all previous figures? Type "Yes" to continue: ', 's');
end

close all; clc;
disp('Continuing execution...');

P = 5;  % Number of dominant paths to consider (expected low channel rank)

% Initialize containers for features and labels
F_all = [];      % Matrix to store all features (22 x total_samples)
labels_all = []; % Vector to store corresponding labels (1 x total_samples)

% List of datasets corresponding to 4 scenarios: empty, 1 subject, 2 subjects, 3 subjects
datasets = {CSI_empty_complex, CSI_1subject_complex, CSI_2subject_complex, CSI_3subject_complex};
labels = [1, 2, 3, 4];  % Labels for each scenario

for d = 1:numel(datasets)
    X = datasets{d};             % Complex CSI matrix [subcarriers x time_samples]
    [subc, T] = size(X);         

    % --- FEATURE EXTRACTION ---

    % 1 & 2) Top P singular values of X: capture dominant channel energy
    % Note: svds avoids explicit covariance matrix calculation
    [~, S_f, ~] = svds(X, P, 'largest');  % Top singular values in frequency domain
    eig_vals = diag(S_f).^2;             
    eig_f = eig_vals;       

    [~, S_t, ~] = svds(X.', P, 'largest');  % Top singular values in time domain
    eig_t = diag(S_t).^2;                 

    % 3) Frequency-specific features (amplitude and phase)
    amp = abs(X);                      % Magnitude of CSI
    phs = angle(X);                    % Phase of CSI

    % Normalized entropy of amplitude (uniformity of energy distribution)
    dist_amp = amp ./ sum(amp + eps,1);                   
    ent_amp = -sum(dist_amp .* log2(dist_amp + eps),1);   
    ent_stats = [mean(ent_amp), median(ent_amp), max(ent_amp), min(ent_amp), std(ent_amp, 0)]; 

    % Phase difference between consecutive subcarriers
    phs_diff = diff(phs,1,1);                              
    phs_stats = [mean(phs_diff(:)), std(phs_diff(:))];     

    % 4) Minor channel variation: stability relative to average energy
    t_vec = sum(abs(X).^2,1)/T;            
    var_X = var(X,0,1);                    
    delta = var_X ./ (t_vec + eps);        
    delta_stats = [mean(delta), median(delta), max(delta), min(delta), std(delta, 0)];

    % Concatenate all 22 features into a single column vector
    F = [eig_f; eig_t; ent_stats.'; phs_stats.'; delta_stats.'];

    %  Accumulate features and labels
    F_all = [F_all, F];  
    
end

% Define feature names (rows correspond to each feature)
FeatureNames = [ ...
    arrayfun(@(i) sprintf('Multipath_eig_%d', i), 1:P, 'UniformOutput', false), ...
    arrayfun(@(i) sprintf('Temporal_eig_%d', i), 1:P, 'UniformOutput', false), ...
    {'Ent_amp_mean','Ent_amp_median','Ent_amp_max','Ent_amp_min','Ent_amp_std'}, ...
    {'Phs_diff_mean','Phs_diff_std'}, ...
    {'Delta_mean','Delta_median','Delta_max','Delta_min','Delta_std'} ...
];

% Display feature matrix summary (22 x total samples)
FeatureTable = array2table(F_all, 'RowNames', FeatureNames, 'VariableNames', {'Empty','1 Subject','2 Subjects', '3 Subjects'});
disp('--- Matriz de características (22x4) ---');
disp(FeatureTable);

%% 10. Average Maximum Cross-Correlation Values
list_appearances = {appearances_empty, appearances_1subject, appearances_2subject, appearances_3subject};
names_scenarios = {'Empty', '1 Subject', '2 Subjects', '3 Subjects'};

max_corr_values = zeros(1, numel(list_appearances));

for esc = 1:numel(list_appearances)
    appearances = list_appearances{esc};
    N = length(appearances);
    
    max_corrs = zeros(1, N);
    
    for i = 1:N
        rx = appearances{i};
       
        
        [c, lags] = xcorr(rx, txSignal);  
        c_mag = abs(c);
        
        max_corrs(i) = max(c_mag);
    end

    max_corr_values(esc) = mean(max_corrs);

end

% ----- Numerical Summary of Cross-Correlation Values ----
T = table(names_scenarios', max_corr_values', 'VariableNames', {'Scenario', 'MaxCorrelationValue'});
disp('--- Cross-Correlation Metrics by Scenario ---');
disp(T);

%% 11. Calculation and Comparison of Relative Powers
 
power_empty = zeros(1, numel(appearances_empty));
for i = 1:numel(appearances_empty)
    data = rxGrids_empty{i}; % cada datos{i} es una matriz 199×100
    power_empty(i) = mean(abs(data(:)).^2);
end

meanPower_empty = mean(power_empty);
stdPower_empty = std(power_empty);

power_1subject = zeros(1, numel(appearances_1subject));
for i = 1:numel(appearances_1subject)
    data = rxGrids_1subject{i}; % cada datos{i} es una matriz 199×100
    power_1subject(i) = mean(abs(data(:)).^2);
end

meanPower_1subject = mean(power_1subject);
stdPower_1subject = std(power_1subject);

power_2subject = zeros(1, numel(appearances_2subject));
for i = 1:numel(appearances_2subject)
    data = rxGrids_2subject{i}; % cada datos{i} es una matriz 199×100
    power_2subject(i) = mean(abs(data(:)).^2);
end

meanPower_2subject = mean(power_2subject);
stdPower_2subject = std(power_2subject);

power_3subject = zeros(1, numel(appearances_3subject));
for i = 1:numel(appearances_3subject)
    data = rxGrids_3subject{i}; % cada datos{i} es una matriz 199×100
    power_3subject(i) = mean(abs(data(:)).^2);
end

meanPower_3subject = mean(power_3subject);
stdPower_3subject = std(power_3subject);

% ------ Nummerical Summary of Relative Powers ------
Class = ["Empty"; "1 Subject"; "2 Subjects"; "3 Subjects"];
mean_relativepowers = [meanPower_empty; meanPower_1subject; meanPower_2subject; meanPower_3subject];
T = table(Class, mean_relativepowers);
disp('--- Relative Powers ---'); disp(T);

