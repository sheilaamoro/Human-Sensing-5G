%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — Zadoff-Chu TX Sequence Generation
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: generation_txSignalGolay.m
%
% Description:
% This function generates the OFDM transmit signal when using Golay complementary 
% sequences as the sensing pattern. Specifically, odd-numbered OFDM symbols 
% carry the Golay sequence Ga, while even-numbered symbols carry Gb, in order 
% to exploit their perfect complementary autocorrelation properties. It
% maps the GCPs onto an OFDM resource
% grid, inserts cyclic prefix, defines null subcarriers (including DC), and
% returns the full transmit signal, the data grid, and the null indices. The resulting
% signal is ready to be passed to the USRP for transmission.
%
% [txSignal, dataGrid, nullIndices] = generacion_txSignalGolayAlterno(Nfft, Ncp, ofdmSampleRate, N_golay, NofdmSyms)
%
% Input parameters:
%     Nfft: Total number of subcarriers (FFT size)
%     Ncp: Lenght of the cyclic prefix
%     ofdmSampleRate: Sampling rate (defines bandwidth)
%     N_golay: Lenght of each Golay Complementary Sequence (power of 2)
%     NofdmSyms: Number of OFDM symbols to generate
%
% Output parameters:
%     txSignal: Time-domain transmit signal
%     dataGrid: CRS sequence mapped onto the OFDM symbols (N_crs x NofdmSyms)
%     nullIndices: Indices of null subcarriers (guard bands + DC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [txSignal, dataGrid, nullIndices] = generation_txSignalGolay(Nfft, Ncp, ofdmSampleRate, N_golay, NofdmSyms)

% Check the lenght
  if mod(log2(N_golay),1)~=0
    error('N_golay debe ser una potencia de 2.');
  end
  if mod(NofdmSyms,2)~=0
    error('NofdmSyms debe ser par para alternar Ga/Gb.');
  end

% ---- Generate GCPs ---- 
[Ga, Gb] = wlanGolaySequence(N_golay);

% ---- OFDM Data Grid ----
% The sequence is mapped over Nfft subcarriers, with N_zc active subcarriers. 
DC_index  = Nfft/2 + 1; % DC subcarrier index (for even Nfft)                         
leftGuard = floor((Nfft - (N_golay+1))/2);           
activeIndices = (leftGuard+1):(leftGuard+N_golay+1); 

% Index of the DC subcarrier within the active block:
dcPos = DC_index - leftGuard;

% Define null subcarriers
nullIndices = [1:leftGuard, DC_index, (leftGuard+N_golay+2):Nfft].';

% Define data grid
dataGrid = zeros(N_golay, NofdmSyms);

for symIdx = 1:NofdmSyms
    % Insert Ga if even, Gb if odd
    if mod(symIdx,2)==1
        s = Ga;
    else
        s = Gb;
    end
    % Insert the GCPs and set DC component to 0
    fullSym = [s(1:dcPos-1); 0; s(dcPos:end)];  
    % Remove DC for ofdmod
    fullSym(dcPos) = [];
    dataGrid(:,symIdx) = fullSym;
end


% ---- Resource Grid Visualisation ----
resourceGridType = zeros(Nfft, NofdmSyms);% Initialise with zeros (nulls)
resourceGridType(setdiff(activeIndices, DC_index), :) = 2; % Active data
resourceGridType(setdiff(nullIndices, DC_index), :) = 0; % Null 
resourceGridType(DC_index, :) = 1; % DC carrier

figure();
hasNulls = any(resourceGridType(:) == 0);

if hasNulls
    customColormap = [0.88 0.88 0.88;  % Gray for Guard
        0.95 0.55 0.55;  % Red for DC
        0.30 0.75 0.60]; % Green for Sequence
    imagesc(resourceGridType);
    colormap(customColormap);
    c = colorbar;
    c.Ticks = [0.33 1 1.67];
    c.TickLabels = {'Guard', 'DC', 'Sequence'};
else

    customColormap = [0.95 0.55 0.55;  % Red for DC
                      0.30 0.75 0.60]; % Green for Sequence
    imagesc(resourceGridType - 1);
    colormap(customColormap);
    c = colorbar;
    c.Ticks = [0.25 0.75];
    c.TickLabels = {'DC', 'Sequence'};
end


xlim([1 NofdmSyms]);
if NofdmSyms <= 10
    xticks(1:NofdmSyms);
    xtickangle(0);
else

    if NofdmSyms <= 50
        xticks(unique([1, 5:5:NofdmSyms]));
        xtickangle(0);
    else
        xticks(unique([1, 10:10:NofdmSyms])); 
    end
end

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize',12);
title('\it{Subcarrier Allocation in OFDM Resource Grid}');
xlabel('OFDM Symbol Index');
ylabel('Subcarrier Index');

% ---- OFDM Modulation ----
ofdmModSignal = ofdmmod(dataGrid, Nfft, Ncp, nullIndices);
txSignal = ofdmModSignal(:);
end
