%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title       : Human detection by 5G radio signals — Correlation-based Peak Detection
% Author      : Sheila Moro Robledo
% Institution : Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree      : Bachelor's Degree in Telecommunications Engineering
% Date        : June, 2025
% File        : find_best_match.m
%
% Description :
% This function computes the normalized cross-correlation (NCC) between a known
% transmitted signal pattern (txSignal) and a received signal (rxSignal) to detect
% occurrences of the transmitted pattern within the received data. It returns the
% index of the best match, the array of NCC scores, and a table of top candidates
% surpassing a dynamic threshold.
%
% [bestIdx, nccScores, topCandidates] = find_best_match(txSignal, rxSignal)
%
% Input parameters:
%    txSignal : Known transmitted signal vector (complex column vector)
%    rxSignal : Received signal vector (complex column vector)
%
% Output:
%    bestIdx      : Index in rxSignal with the highest normalized cross-correlation
%    nccScores    : Vector of normalized cross-correlation scores for each valid lag
%    topCandidates: Table with indices and NCC scores above the significance threshold
%
% Notes:
%   - The NCC score is computed as the absolute value of the correlation normalized
%     by the energy of the signals.
%   - The function handles NaN values by ignoring invalid scores.
%   - Throws errors if input signals are improperly sized or if no significant
%     matches are found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestIdx, nccScores, topCandidates] = find_best_match(txSignal, rxSignal)

% Ensure column vector
txSignal = txSignal(:);  
rxSignal = rxSignal(:); 

Ltx = length(txSignal);
Lrx = length(rxSignal);

if Lrx < Ltx
    error("rxSignal must be longer than txSignal.");
end

% Normalize by energy of txSignal
txNorm = norm(txSignal);
if txNorm == 0
    error("txSignal is empty or contains only zeros.");
end

% Compute cross-correlation using convolution
corrVal = conv(rxSignal, flipud(conj(txSignal)), 'valid');

% Compute energy of each window in rxSignal
rxEnergy = conv(abs(rxSignal).^2, ones(Ltx,1), 'valid');

% Compute normalized cross-correlation (NCC) score
nccScores = abs(corrVal) ./ (sqrt(rxEnergy) * txNorm);

% Filter out NaN scores
validMask = ~isnan(nccScores);
validScores = nccScores(validMask);
validIndices = find(validMask);

if isempty(validScores)
    error("All NCC scores are NaN. Cannot determine bestIdx.");
end

% Sort valid scores descending
[sortedScores, sortOrder] = sort(validScores, 'descend');
sortedIdx = validIndices(sortOrder);

bestIdx = sortedIdx(1);  % Best valid match index

% Dynamic threshold to identify significant matches
threshold = 0.7;  
significantIdx = sortedScores > threshold;

% Extract top candidates passing threshold
topCandidates = table(sortedIdx(significantIdx), sortedScores(significantIdx), ...
    'VariableNames', {'Index', 'NCC_Score'});

if isempty(topCandidates)
    error("No significant matches detected.");
end
end
