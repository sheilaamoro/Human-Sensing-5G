%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title       : Human detection by 5G radio signals — Demodulation Function
% Author      : Sheila Moro Robledo
% Institution : Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree      : Bachelor's Degree in Telecommunications Engineering
% Date        : June, 2025
% File        : demodulation.m
%
% Description :
% This function performs synchronization and OFDM demodulation of a received signal 
% based on a known transmitted signal pattern. It uses a robust correlation-based 
% method to detect occurrences of the transmitted sequence in the received stream, 
% taking into account the specific structure of the training sequence (e.g., Golay 
% with symbol pairs). It returns the demodulated OFDM grids for each detected burst 
% and information about the number of complete appearances.
%
% [apariciones, rxGrids, topHits] = demodulation(txSignal, rxSignal,Nfft, Ncp, nullIndices, NofdmSyms, type)
%
% Input parameters:
%    txSignal: Transmitted signal (time domain) used as reference
%    rxSignal: Received signal (time domain)
%    Nfft: FFT size (number of subcarriers)
%    Ncp: Cyclic prefix length
%    nullIndices: Indices of subcarriers that are deactivated (e.g., guard bands)
%    NofdmSyms: Number of OFDM symbols expected per transmission frame
%    type: Type of training sequence ('golay' | 'zc' | 'crs')
%
% Output:
%    apariciones: Cell array containing individual bursts found in the received signal
%    rxGrids: Cell array of demodulated OFDM grids (one per burst)
%    topHits: Top correlation matches found in the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [apariciones, rxGrids, topHits] = demodulation(txSignal, rxSignal,Nfft, Ncp, nullIndices, NofdmSyms, type)

% Adjust the detection pattern depending on the type of sequence
if strcmpi(type, 'golay')
    txPattern = txSignal(1:2*(Nfft+Ncp));
    NofdmSyms = NofdmSyms/2;
else
    txPattern = txSignal(1:(Nfft+Ncp));
end


% Perform correlation-based detection
[bestIdx, scores, topHits] = find_best_match(txPattern, rxSignal);

% Extract the indices of the best matches
topIndices = topHits.Index;
patternLength = length(txPattern);

% Find the actual start index aligned with full repetitions
trueStartIdx = find_aligned_start(topIndices, patternLength,NofdmSyms);

% Locate complete appearances of the transmitted pattern
lengthTxSignal = length(txSignal);
apariciones = {}; 

i = 0;
while true
    startIdx = trueStartIdx + i * lengthTxSignal;
    endIdx = startIdx + lengthTxSignal - 1;

    if endIdx > length(rxSignal)
        break;
    end

    apariciones{end+1} = rxSignal(startIdx:endIdx);
    i = i + 1;
end

disp(['Number of complete appearances found: ', num2str(length(apariciones))]);

% OFDM demodulation for each detected burst
rxGrids = cell(1, length(apariciones));

for i = 1:length(apariciones)
    rxGrid = ofdmdemod(apariciones{i}, Nfft, Ncp,Ncp, nullIndices);
    rxGrids{i} = rxGrid;
end
end