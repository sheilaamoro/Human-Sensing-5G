%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — Zadoff-Chu TX Sequence Generation
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: generation_txSignalZC.m
%
% Description:
% This function generates the baseband OFDM signal to be transmitted when the 
% sensing sequence is Zadoff-Chu. It maps the ZC sequence onto an OFDM resource
% grid, inserts cyclic prefix, defines null subcarriers (including DC), and
% returns the full transmit signal, the data grid, and the null indices. The resulting
% signal is ready to be passed to the USRP for transmission.
%
% function [txSignal, dataGrid, nullIndices] = generation_txSignalZC(Nfft,Ncp, ofdmSampleRate,N_zc, NofdmSyms)
%
% Input parameters:
%     Nfft: Total number of subcarriers (FFT size)
%     Ncp: Lenght of the cyclic prefix
%     ofdmSampleRate: Sampling rate (defines bandwidth)
%     N_zc: Lenght of the Zadoff-Chu sequence
%     NofdmSyms: Number of OFDM symbols to generate
%
% Output parameters:
%     txSignal: Time-domain transmit signal
%     dataGrid: ZC sequence mapped onto the OFDM symbols (N_zc x NofdmSyms)
%     nullIndices: Indices of null subcarriers (guard bands + DC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [txSignal, dataGrid, nullIndices] = generation_txSignalZC(Nfft,Ncp, ofdmSampleRate,N_zc, NofdmSyms)

% ---- Generate Zadoff-Chu Sequence ---- 
R = 25;
zcSeq_time = zadoffChuSeq(R, N_zc);
zcSeq = fft(zcSeq_time)./norm(zcSeq_time); 

% ---- OFDM Data Grid ----
% The sequence is mapped over Nfft subcarriers, with N_zc active subcarriers. 
leftGuard = floor((Nfft - (N_zc + 1)) / 2);        
activeIndices = (leftGuard+1):(leftGuard+N_zc+1);  
DC_index = Nfft/2 + 1; % DC subcarrier index (for even Nfft)

% Index of the DC subcarrier within the active block:
dcDataIdx = DC_index - leftGuard;  

% Create one OFDM symbol vector:
% Insert the ZC sequence and set DC component to O
dataSymbol_full = [zcSeq(1:dcDataIdx-1); 0; zcSeq(dcDataIdx:end)]; 
dataSymbol = dataSymbol_full;
dataSymbol(dcDataIdx) = []; % Remove DC for ofdmod

% Repeat the same data symbol across all OFDM symbols
dataGrid = repmat(dataSymbol, 1, NofdmSyms);  % Size: [N_zc x NofdmSyms]

% Define null subcarriers
nullIndices = [1:leftGuard, DC_index, (leftGuard+N_zc+2):Nfft].';


% ---- Resource Grid Visualisation ----
resourceGridType = zeros(Nfft, NofdmSyms); % Initialise with zeros (nulls)
resourceGridType(setdiff(activeIndices, DC_index), :) = 2; % Active data
resourceGridType(setdiff(nullIndices, DC_index), :) = 0; % Null 
resourceGridType(DC_index, :) = 1;  % DC carrier

figure();
hasNulls = any(resourceGridType(:) == 0);

if hasNulls
    customColormap = [0.88 0.88 0.88;  % Gray for Guard 
                      0.95 0.55 0.55;  % Red for DC
                      0.30 0.75 0.60]; % Green for ZC Sequence
    imagesc(resourceGridType); 
    colormap(customColormap);
    c = colorbar;
    c.Ticks = [0.33 1 1.67];
    c.TickLabels = {'Guard', 'DC', 'Sequence'};
else
   
    customColormap = [0.95 0.55 0.55;  % Red for DC
                      0.30 0.75 0.60]; % Green for ZC Sequence
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