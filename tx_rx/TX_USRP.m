
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — Transmitter Script
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: TX_USRP.m

% Description:
% This script generates and transmits an OFDM singal using a selected sequence
% (Zadoff-Chu, CRS, or Golay) through a USRP device. The signal is transmitted
% continuously for use in experiments related to human presence detection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Transmitter Configuration Script for OFDM Signal over USRP
clear; close all; clc

%% ---- OFDM Parameters Definition -----

Nfft = 104;                         % FFT Size
Ncp = 7;                            % Cyclic prefix length
ofdmSampleRate = 25e6;              % Sampling rate (Hz)
NofdmSyms = 100;                    % Number of OFDM symbols in the frame
scs = ofdmSampleRate/(1e3*Nfft);    % Subcarrier spacing (kHz)

% Choose the sequence length depending on the type
N_zc = 103;                        % Zadoff-Chu sequence length (active data)
%N_crs = 103;                      % CRS sequence length (active data)
%N_golay = 2^(floor(log2(Nfft)));  % Golay complementary sequences length (active data)

%% ---- Transmit Signal Generation ----

% Uncomment one of the following lines depending on the desired sequence
[txSignal, dataGrid, nullIndices] = generacion_txSignalZC(Nfft,Ncp,ofdmSampleRate,N_zc,NofdmSyms);
%[txSignal, dataGrid, nullIndices] = generacion_txSignalCRS(Nfft,Ncp,ofdmSampleRate,N_crs,NofdmSyms);
%[txSignal, dataGrid, nullIndices] = generacion_txSignalGolayAlterno(Nfft, Ncp, ofdmSampleRate, N_golay, NofdmSyms);

%% ---- USRP Transmitter Setup ----

tx = comm.SDRuTransmitter('192.168.1.150');

% Transmitter parameters
tx.CenterFrequency = 3.75e9;    % Center frequency (Hz)
tx.Gain = 30;                   % Transmit gain (dB)
tx.ClockSource = 'External';    % Use external clock source
tx.PPSSource = 'External';      % Use external PPS signal source
tx.InterpolationFactor = 4;     % Interpolation factor

%% ---- Continuous Transmission Loop ----
while true
   tx(txSignal);
end

