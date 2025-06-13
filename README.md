# Human-Sensing-5G

This repository contains MATLAB code developed for the detection of human presence leveraging CSI and additional signal processing techniques within 5G OFDM framework. The project was carried out as part of the undergraduate thesis **Detección de personas mediante señales de radio 5G**, for the **Bachelor's Degree in Telecommunication Technologies and Services Engineering** at the **University of Oviedo**.

## Project Overview

This repository contains scripts for:
- Generation of sensing sequences and OFDM signals.
- Transmission and reception of signals using USRP hardware (USRP N210), alongside data acquisition and storage.
- Preprocessing routines for received signals.
- Feature extraction based on CSI and advanced signal analysis methodologies.

## System Requirements

- MATLAB R2022b or later.
- Signal Processing Toolbox.
- Statistics and Machine Learning Toolbox.
- Communications Toolbox.
- NI USRP Radio Support (Communications Toolbox).
- LTE Toolbox.
- WLAN Toolbox.

## Repository Structure

- `tx_rx/` — Scripts dedicated to signal transmission, reception and data logging.
- `processing/` — Signal processing algorithms, feature extraction procedures, auxiliary or helper functions.
- `README.md` — Projct documentation.
- `LICENSE` — Software license information.
  
> Please, note that the recorded datasets utilised during experimentation are not included in this repository due to size constraints.

## Author

**Sheila Moro Robledo**.

Bachelor's Degree in Telecommunication Technologies and Servicies Engineering, University of Oviedo, 2025.
