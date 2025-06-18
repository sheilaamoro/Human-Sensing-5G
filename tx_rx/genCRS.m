%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — CRS Generation Function
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: genCRS.m
%
% Description :
% This function generates the LTE Cell-specific Reference Signal (CRS) defined in
% 3GPP TS 36.211, Section 6.10.1. It creates a pseudo-random QPSK sequence based on 
% the cell ID, OFDM symbol number, slot number, and frame configuration. The output 
% is a complex sequence used for channel estimation and synchronization in LTE systems.
%
% s_crs = genCRS(N_crs, NcellID, Ns, Nl, cp_type, frame_type3)
%
% Input parameters:
%    N_crs: Number of QPSK symbols to generate
%   NcellID: Cell ID used to seed the pseudorandom generator
%   Ns: Slot number within the radio frame
%   Nl: OFDM symbol index within the slot
%   cp_type: Type of cyclic prefix ('normal' or 'extended')
%   frame_type3: Boolean flag indicating type 3 frame structure (true if DRS CRS)
%
% Output:
%    s_crs: Complex vector of QPSK symbols representing the CRS sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s_crs = genCRS(N_crs, NcellID, Ns, Nl, cp_type, frame_type3)

    % 1. Configuration parameters
    Nc = 1600;                       % Initialisation offset
    Nbits = 2 * N_crs + Nc;          % Useful bits + offset
    Ngen = Nbits + 31;              

    % 2. Determine CP type
    if strcmpi(cp_type, 'normal')
        N_CP = 1;
    else
        N_CP = 0;
    end

    % 3. Compute ns depending on frame structure
    if frame_type3
        ns_dash = floor(Ns/10) + mod(Ns, 2);
    else
        ns_dash = Ns;
    end

    % 4. Compute initialisation value c_init according to TS 36.211 Section 6.10.1.1
    c_init = 2^10 * (7 * (ns_dash + 1) * (2 * NcellID + 1) + 2 * NcellID + N_CP);

    % 5. Initialise x1 and x2 sequences
    x1 = [1, zeros(1,30)];
    x2 = de2bi(c_init, 31, 'left-msb');

    % 6. Generate x1 and x2 sequences
    x1_seq = zeros(Ngen,1);
    x2_seq = zeros(Ngen,1);
    x1_seq(1:31) = x1;
    x2_seq(1:31) = x2;

    for n = 1:(Ngen-31)
        x1_seq(n+31) = mod(x1_seq(n+3) + x1_seq(n), 2);
        x2_seq(n+31) = mod(x2_seq(n+3) + x2_seq(n+2) + x2_seq(n+1) + x2_seq(n), 2);
    end

    % 7. Generate c(n)
    c = mod(x1_seq(Nc+1:Nc+2*N_crs) + x2_seq(Nc+1:Nc+2*N_crs), 2);

    % 8. QPSK mapping
    I = 1 - 2*c(1:2:end);
    Q = 1 - 2*c(2:2:end);
    s_crs = (1/sqrt(2)) * (I + 1j*Q);
end
