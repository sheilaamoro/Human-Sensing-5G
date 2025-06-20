%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — CSI Estimation Function
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: CSI.m
%
% Description :
% This function estimates the frequency response H(f), also known as Channel State 
% Information (CSI), from received OFDM grids and known transmitted data grids. 
% The estimation method depends on the type of sensing sequence used (Golay, Zadoff-Chu,
% or CRS). For Golay sequences, the function exploits their perfect complementary 
% autocorrelation by processing symbol pairs. For ZC and CRS sequences, a direct 
% complex division is used.
%
% CSI_out = CSI(rxGrids, dataGrid, NofdmSyms, type)
%
% Input parameters:
%    rxGrids: Cell array of received grids (one per antenna), each of size [Nsc x Nsym]
%    dataGrid: Known transmitted data grid of the same size as each rxGrid
%    NofdmSyms: Number of OFDM symbols in the frame
%    type: Type of training sequence used: 'Golay' | 'ZC' | 'TSS' (CRS)
%
% Output:
%    CSI_out: Estimated CSI matrix [Nsc x (Nsym × Nantennas)], concatenated in frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CSI_out = CSI(rxGrids, dataGrid, NofdmSyms, type)
CSI_cells = cell(1,numel(rxGrids));

for i = 1:numel(rxGrids)
    Y = rxGrids{i};   % Received symbols
    X = dataGrid;     % Known transmitted symbols

    switch lower(type)
      case 'golay'
        H = zeros(size(Y));
        % Golay sequences: use pairwise processing
        for k = 1:2:NofdmSyms
          Ya = Y(:,k);
          Yb = Y(:,k+1);
          Xa = X(:,k);
          Xb = X(:,k+1);
          H_pair = (Xa.*Ya + Xb.*Yb) ./ (abs(Xa).^2 + abs(Xb).^2);
          H(:,k)   = H_pair;
          H(:,k+1) = H_pair;
        end

      otherwise
        % Zadoff-Chu or CRS: element-wise complex division
        H = Y ./ X;
    end
    CSI_cells{i} = H;
end

CSI_out = cat(2, CSI_cells{:});
end
