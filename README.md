# Human-Sensing-5G

This repository contains MATLAB code developed for the detection of human presence leveraging CSI and additional signal processing techniques within 5G OFDM framework. The project was carried out as part of the undergraduate thesis **Human Detection by 5G Radio Signals**, for the **Bachelor's Degree in Telecommunication Technologies and Services Engineering** at the **University of Oviedo**.

## Project Overview

This repository provides MATLAB scripts and functions for:
- Generation of sensing sequences and OFDM-modulated signals.
- Signal transmission and reception using USRP hardware (USRP N210), including data acquisition and storage.
- Preprocessing routines for received signals.
- Feature extraction based on CSI and advanced analytical methods for human presence detection.

## System Requirements

- MATLAB R2022b or later, with the following toolboxes: Signal Processing Toolbox, Statistics and Machine Learning Toolbox, Communications Toolbox, NI USRP Radio Support (Communications Toolbox), LTE Toolbox, WLAN Toolbox.
- USRP N210 hardware, configured according to project specifications, is required for executing transmission and reception scripts.

## Repository Structure

- `tx_rx/` — Scripts dedicated to signal transmission, reception and data logging.
- `processing/` — Signal processing algorithms, feature extraction procedures, auxiliary or helper functions. Some functions in `tx_rx/` may be also be required.
- `README.md` — Project documentation.
- `LICENSE` — Software license information.
  
> Please, note that the recorded datasets used during experimentation are not included in this repository due to size constraints.

## Usage instructions

1. Configure key parameters (e.g., carrier frequency, transmission gain) within the scripts located in the `tx_rx/` directory.
2. Launch the transmission and reception processes using the provided USRP control scripts to capture IQ datasets. Two independent MATLAB sessions are requited: one for transmission and one for reception. First, run the transmission script on both ends, then execute the reception script in the corresponding session.
3. Use the scripts in the `processing/` directory to process the acquired data and extract feautres. Each `.mat` file should be sorted into one of four directories, corresponding to the experimental scenarios: `empty`, `1subject`, `2subject`, `3subject`.

## Author

**Sheila Moro Robledo**.

Bachelor's Degree in Telecommunication Technologies and Servicies Engineering, University of Oviedo, 2025.
