
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: Human detection by 5G radio signals — Receiver Script
% Author: Sheila Moro Robledo
% Institution: Polytechnic School of Engineering of Gijón, University of Oviedo
% Degree: Bachelor's Degree in Telecommunications Engineering
% Date: June, 2025
% File: RX_USRP.m
%
% Description:
% This script receives and stores IQ samples using a USRP device for a
% series of repeated measurements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause(7);       % Initial pause before starting the first measurement
total = 5;
for rep = 1:total
    disp(['Starting measurement ' num2str(rep) ' of ' num2str(total)]);

    % ----- USRP parameters -----
    freq       = 3.75e9;                    % Center frequency (Hz)
    ganancia   = 35;                        % Receiver gain (dB)
    decim      = 4;                         % Decimaton factor
    frameLen   = (Nfft + Ncp)*NofdmSyms;    % Samples per OFDM frame
    maxSamples = 100e6;                    % Max number of samples to capture

    % ---- USRP receiver inicialization ----
    rx = comm.SDRuReceiver('192.168.1.151');
    rx.CenterFrequency    = freq;
    rx.Gain               = ganancia;
    rx.ClockSource        = 'External';
    rx.PPSSource          = 'External';
    rx.DecimationFactor   = decim;
    rx.SamplesPerFrame    = frameLen;
    rx.OutputDataType     = "double";
    fs = rx.MasterClockRate / rx.DecimationFactor;

    % ---- Preallocate storage vector ----
    IQ = complex(zeros(maxSamples, 1)); 
    counter = 0;
    % ---- Main reception loop ----
    disp('Receiving and storing IQ samples... Press Ctrl+C to interrupt');

    try
        while counter < maxSamples
            [data, valid] = rx();
            if valid
                n = length(data);
                remaining = maxSamples - counter;
                if n <= remaining
                    IQ(counter+1 : counter+n) = data;
                    counter = counter + n;
                else
                    IQ(counter+1 : end) = data(1:remaining);
                    counter = maxSamples;
                end
            end
        end
    catch
        disp('Capture interrupted by user');
        IQ = IQ(1:counter);  % Truncate vector if user interrupted
    end

    % ---- Save IQ data ----
    release(rx);
    nombreArchivo = ['IQ_240GOLAY_3PERSONA' datestr(now,'yyyymmdd_HHMMSS') '.mat'];
    save(nombreArchivo, 'IQ', '-v7.3');  
    disp(['Datos guardados en: ' nombreArchivo]);

    if rep < total
        disp('Waiting 5 seconds before the next measurement ....');
        pause(5);
    end
end

disp('All measurements completed.');
