%% This code is modified from CLFM_fish_Reconstruction_CPU.m file
%Reference: Zhang, Z., Bai, L., Cong, L. et al. Imaging volumetric dynamics at high speed in mouse and zebrafish brain with confocal light field microscopy. Nat Biotechnol 39, 74–83 (2021).

%% Volumetric reconstruction for any scanning VoV or VoIs.
%Author: Jiazhen Zhai, Xinyue Chen, Yixin Chen, Chi Liu, Junhao Deng, Hao Tang, Weisheng Zhang, Junhao Liang, Chenshuang Zhang, Liangcai Cao, Peifu Tang, Licheng Zhang, Yu Mu, Junle Qu, and Lingjie Kong
%Related to our manuscript: Scanning wide-field tomography for high-speed, mesoscale-volumetric imaging of biodynamics in vivo
% SYSTEM REQUIREMENTS
% Memory: 128 GB RAM
% Graphics: 24 GB RAM
% MATLAB: R2021a

%% This code includes several functional modules: 
%1. Computational optical sectionning by SIM
%2. Pixel-wise alignment
%3. Dual-focus reconstruction
%4. 3D reconstruction by RL deconvolution

%% The demo_data includes:
% 9 raw images of sub-VoV1 to sub-VoV9 in mouse brain cortex in vivo.

Run 'demo_reconstruction.m' to perform volumetric reconstruction for the whole VoV.