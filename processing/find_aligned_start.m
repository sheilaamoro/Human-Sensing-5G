%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title       : Human detection by 5G radio signals — Aligned Start Index Finder
% Author      : Sheila Moro Robledo
% Institution : Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree      : Bachelor's Degree in Telecommunications Engineering
% Date        : June, 2025
% File        : find_aligned_start.m
%
% Description :
% This function searches for the smallest index within a list of candidate indices
% ('topIndices') such that a sequence of 'minHits' consecutive indices spaced by
% 'patternLength' all exist in 'topIndices'. It is used to determine the true start
% index of repeated signal patterns within a received signal.
%
% trueStartIdx = find_aligned_start(topIndices, patternLength, minHits)
%
% Input parameters:
%    topIndices    : Sorted vector of candidate indices (integer array)
%    patternLength : Length of the repeated pattern (integer)
%    minHits       : Minimum number of consecutive aligned occurrences required (integer)
%
% Output:
%    trueStartIdx  : The smallest index fulfilling the alignment condition
%
% Throws an error if no such aligned start index with sufficient repetitions is found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trueStartIdx = find_aligned_start(topIndices, patternLength, minHits)

    topIndices = sort(topIndices);  % Ensure ascending order

    for i = 1:length(topIndices)
        startIdx = topIndices(i);
        hits = 0;

        % Check if startIdx + k * patternLength are in topIndices for k = 0 to minHits-1
        for k = 0:minHits-1
            candidate = startIdx + k * patternLength;
            if ismember(candidate, topIndices)
                hits = hits + 1;
            else
                break; % Sequence broken
            end
        end

        if hits == minHits
            trueStartIdx = startIdx;
            return;
        end
    end

    error('No aligned start index found with the required number of repetitions.');
end
