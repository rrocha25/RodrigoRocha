clear; clc; close all;

%% 1) Parameters
fs = 44100;
c  = 1500;

dur_s = 80e-3;
t = (0:round(dur_s*fs)-1)' / fs;

SNR_dB = 10;

rng(42);

% Optional bandpass filter
apply_bandpass = false;

% Geometry
fmax_geom  = 850;
lambda_min = c / fmax_geom;
L          = 0.5 * lambda_min;

pos_h = pyramid_array_positions(L);

% Simulated source
az_real = 35;
el_real = 24;
r_real  = 10;

pos_f = sph2cart_deg(r_real, az_real, el_real);

% Output directory for figures
outputDir = './figures';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Standard 2:1 figure dimensions for ALL figures
fig_width  = 1000;
fig_height = 500;

%% 2) Physical maxLag
dmax          = max_pair_distance(pos_h);
maxLagSamples = ceil((dmax/c)*fs) + 2;
fprintf('\nmaxLagSamples = %d samples  %.3f ms\n', ...
        maxLagSamples, 1000*maxLagSamples/fs);

%% 3) Simple 2-tone signal
f1   = 400;
f2   = 800;
A1   = 1.0;
A2   = 0.7;
phi1 = 0;
phi2 = 0;

s = A1*sin(2*pi*f1*t + phi1) + A2*sin(2*pi*f2*t + phi2);

%% --- Ajustar nível de SPL do sinal 2-tons ---
p_ref = 1e-6;   % 1 µPa

%% --- Spectrum plot (Sec. 3) ---
N     = length(s);
S_fft = fft(s);
S_fft = S_fft(1:floor(N/2)+1);
f_fft = (0:floor(N/2))' * fs / N;


%% 4) Simulate reception at 4 hydrophones with relative delays
N = length(s);
dist = zeros(4,1);
for k = 1:4
    dist(k) = norm(pos_f - pos_h(k,:));
end
tau_rel = (dist - dist(1))/c;

guard = maxLagSamples + 4;
s_pad = [zeros(guard,1); s; zeros(guard,1)];

x = zeros(N,4);
for k = 1:4
    sk  = apply_delay_seconds(s_pad, fs, tau_rel(k));
    seg = sk(guard+1:guard+N);
    seg = add_awgn_snr(seg, SNR_dB);
    x(:,k) = seg;
end
%% 5) Plot SPL por Frequência (todos os hidrofones)

NFFT  = 2^nextpow2(N);          % tamanho da FFT (potência de 2)
f_vec = (0:NFFT/2-1) * (fs/NFFT); % vetor de frequências (Hz)

p_ref = 1e-6;                   % referência padrão em underwater acoustics (1 µPa)

figure('Name','SPL por Frequência','NumberTitle','off');
hold on;

colors = lines(4);
labels = cell(4,1);

for k = 1:4
    X_k   = fft(x(:,k), NFFT);          % FFT do canal k
    X_k   = X_k(1:NFFT/2);              % metade positiva

    % Amplitude RMS por bin ? SPL em dB re p_ref
    amp   = abs(X_k) / NFFT;            % normalização
    amp   = amp * sqrt(2);              % pico ? RMS (exceto DC)
    amp(1)= abs(X_k(1)) / NFFT;        % DC: sem fator ?2

    spl_k = 20 * log10(amp / p_ref + eps);

    plot(f_vec/1e3, spl_k, 'Color', colors(k,:), 'LineWidth', 1.2);
    labels{k} = sprintf('Hidrofone %d  (\\tau_{rel} = %.3f ms)', ...
                         k, tau_rel(k)*1e3);
end

xlabel('Frequência (kHz)');
ylabel('SPL  (dB re 1 µPa)');
title('SPL por Frequência – Sinais Recebidos nos Hidrofones');
legend(labels, 'Location','best');
grid on;
xlim([0, fs/2/1e3]);
hold off;

%% Optional bandpass filtering
if apply_bandpass
    fprintf('\nApplying bandpass filter (300-1000 Hz)...\n');
    [b,a] = butter(4, [300, 1000]/(fs/2), 'bandpass');
    x_processed = zeros(size(x));
    for k = 1:4
        x_processed(:,k) = filtfilt(b,a,x(:,k));
    end
else
    fprintf('\nNo bandpass filter.\n');
    x_processed = x;
end

%% Plot: 4 array channels - FORMATO 2:1
t_ref_ms = 38.05;

for k = 1:4
    fig = figure('Position', [100 100 1250 600], 'Color', 'w');

    plot(1e3*t, x(:,k), 'k', 'LineWidth', 1.0); 
    grid on; 
    hold on

    xlabel('Time [ms]', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel(sprintf('Amplitude x_%d', k), 'FontSize', 20, 'FontWeight', 'bold');

    yl = ylim;
    t_del_ms = t_ref_ms + 1e3*tau_rel(k);
    h_ref = plot([t_ref_ms t_ref_ms], yl, 'b--', 'LineWidth', 1.2);
    h_exp = plot([t_del_ms t_del_ms], yl, 'r-',  'LineWidth', 1.2);
    ylim(yl);
    xlim([30 50]);

    legend([h_ref h_exp], {sprintf('Ref: %.2f ms', t_ref_ms), ...
                           sprintf('Exp: %.2f ms', t_del_ms)}, ...
          'Location', 'northeast', 'Interpreter', 'none', 'FontSize', 20);

    set(gca, 'FontSize', 20);
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    hold off

    print(fullfile(outputDir, sprintf('channel_%d.eps', k)), '-depsc')
    %close(fig)
end

%% 5) Estimate TDOAs and plot correlations - FORMATO 2:1
pares = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
Np = size(pares,1);
u_real = [cosd(el_real)*cosd(az_real) cosd(el_real)*sind(az_real) sind(el_real)];

tdoa_phat     = zeros(Np,1);
tdoa_gcc      = zeros(Np,1);
tdoa_th_plane = zeros(Np,1);

for p = 1:Np
    i = pares(p,1);
    j = pares(p,2);

    dp = pos_h(j,:) - pos_h(i,:);
    tdoa_th_plane(p) = -dot(dp, u_real)/c;

    %% GCC-PHAT - FORMATO 2:1
    [rP, lagsP] = gcc_family_local(x_processed(:,j), x_processed(:,i), fs, maxLagSamples, 'phat', 0.1);
    [rPpk, kP] = max(abs(rP));          
    tdoa_phat(p) = lagsP(kP);

    fig = figure('Position', [100 100 1200 600], 'Color', 'w');
    plot(1e3*lagsP, rP, 'k', 'LineWidth', 1.05); 
    grid on; 
    hold on

    yl = ylim;
    h_hat = plot(1e3*[tdoa_phat(p) tdoa_phat(p)], yl, 'm-',  'LineWidth', 1.2);
    h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)], yl, 'b--', 'LineWidth', 1.1);
    h_pk  = plot(1e3*tdoa_phat(p), rP(kP), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
    hold off

    xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
    xlabel('Delay [ms]', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel(sprintf('Corr. (%d,%d)', i, j), 'FontSize', 20, 'FontWeight', 'bold');

    legend([h_hat h_the h_pk], {'Est.', 'Theor.', 'Peak'}, ...
        'Location', 'best', 'Interpreter', 'none', 'FontSize', 20);

    set(gca, 'FontSize', 20);
    set(gca, 'Position', [0.08 0.18 0.90 0.75])

    print(fullfile(outputDir, sprintf('phat_pair_%d%d.eps', i, j)), '-depsc')
    %close(fig)

    %% GCC - FORMATO 2:1
    [rG, lagsG] = gcc_family_local(x_processed(:,j), x_processed(:,i), fs, maxLagSamples, 'gcc',  0.1);
    [rGpk, kG] = max(abs(rG));
    tdoa_gcc(p) = lagsG(kG);

    fig = figure('Position', [100 100 1400 700], 'Color', 'w');
    plot(1e3*lagsG, rG, 'k', 'LineWidth', 1.05); 
    grid on; 
    hold on

    yl = ylim;
    h_hat = plot(1e3*[tdoa_gcc(p) tdoa_gcc(p)], yl, 'm-',  'LineWidth', 1.2);
    h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)], yl, 'b--', 'LineWidth', 1.1);
    h_pk  = plot(1e3*tdoa_gcc(p), rG(kG), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
    hold off

    xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
    xlabel('Delay [ms]', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel(sprintf('Corr. (%d,%d)', i, j), 'FontSize', 20, 'FontWeight', 'bold');

    legend([h_hat h_the h_pk], {'Est.', 'Theor.', 'Peak'}, ...
        'Location', 'best', 'Interpreter', 'none', 'FontSize', 20);

    set(gca, 'FontSize', 20);
    set(gca, 'Position', [0.08 0.18 0.90 0.75])

    print(fullfile(outputDir, sprintf('gcc_pair_%d%d.eps', i, j)), '-depsc')
    %close(fig)
end

fprintf('\nTDOAs (ms) - GCC-PHAT:\n'); disp(1e3*tdoa_phat);
fprintf('TDOAs (ms) - GCC:\n'); disp(1e3*tdoa_gcc);
fprintf('Error (PHAT - theoretical) in ms:\n'); disp(1e3*(tdoa_phat - tdoa_th_plane));
fprintf('Error (GCC  - theoretical) in ms:\n'); disp(1e3*(tdoa_gcc  - tdoa_th_plane));

%% 6) Cost grid for both methods
az_vec = 0:1:359;
el_vec = 0:1:90;

J_phat = doa_cost_grid(pos_h, pares, tdoa_phat, c, az_vec, el_vec);
J_gcc  = doa_cost_grid(pos_h, pares, tdoa_gcc,  c, az_vec, el_vec);

[az_hat_phat, el_hat_phat] = peak_min_cost(J_phat, az_vec, el_vec);
[az_hat_gcc,  el_hat_gcc ] = peak_min_cost(J_gcc,  az_vec, el_vec);

az_real_w = mod(az_real, 360);

%% 7) Plot cost maps - FORMATO 2:1
clim_common = [min([J_phat(:); J_gcc(:)]), max([J_phat(:); J_gcc(:)])];

ms_real = 12; lw_real = 2.2;
ms_est  = 11; lw_est  = 2.4;

% GCC-PHAT cost map
fig = figure('Position', [100 100 1200 600], 'Color', 'w');
imagesc(az_vec, el_vec, J_phat);
set(gca, 'YDir', 'normal');
caxis(clim_common);
colormap(flipud(jet(256)));
cb = colorbar; 
ylabel(cb, 'Cost J', 'FontSize', 20, 'FontWeight', 'bold');
hold on

h_true = plot(az_real_w, el_real, 'o', ...
    'MarkerSize', ms_real, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', lw_real);

plot(az_hat_phat, el_hat_phat, '+', ...
    'Color', 'k', 'MarkerSize', ms_est+2, 'LineWidth', lw_est+1.0);
h_est = plot(az_hat_phat, el_hat_phat, '+', ...
    'Color', 'r', 'MarkerSize', ms_est, 'LineWidth', lw_est);

hold off
xlabel('Azimuth [deg]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Elevation [deg]', 'FontSize', 20, 'FontWeight', 'bold');
legend([h_true h_est], {'True', 'Est.'}, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 20);
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.15 0.81 0.78])  % margem direita ~6%

print(fullfile(outputDir, 'cost_map_phat.eps'), '-depsc')
%close(fig)

% GCC cost map
fig = figure('Position', [100 100 1200 600], 'Color', 'w');
imagesc(az_vec, el_vec, J_gcc);
set(gca, 'YDir', 'normal');
caxis(clim_common);
colormap(flipud(jet(256)));
cb = colorbar; 
ylabel(cb, 'Cost J', 'FontSize', 20, 'FontWeight', 'bold');
hold on

h_true = plot(az_real_w, el_real, 'o', ...
    'MarkerSize', ms_real, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', lw_real);

plot(az_hat_gcc, el_hat_gcc, '+', ...
    'Color', 'k', 'MarkerSize', ms_est+2, 'LineWidth', lw_est+1.0);
h_est = plot(az_hat_gcc, el_hat_gcc, '+', ...
    'Color', 'g', 'MarkerSize', ms_est, 'LineWidth', lw_est);

hold off
xlabel('Azimuth [deg]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Elevation [deg]', 'FontSize', 20, 'FontWeight', 'bold');
legend([h_true h_est], {'True', 'Est.'}, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 20);
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.15 0.81 0.78])  % margem direita ~6%


print(fullfile(outputDir, 'cost_map_gcc.eps'), '-depsc')
%close(fig)

%% Numerical results
u_true = [cosd(el_real)*cosd(az_real_w), cosd(el_real)*sind(az_real_w), sind(el_real)];
u_phat = [cosd(el_hat_phat)*cosd(az_hat_phat), cosd(el_hat_phat)*sind(az_hat_phat), sind(el_hat_phat)];
u_gcc  = [cosd(el_hat_gcc )*cosd(az_hat_gcc ), cosd(el_hat_gcc )*sind(az_hat_gcc ), sind(el_hat_gcc )];

ang_error_phat = acosd(dot(u_true, u_phat));
ang_error_gcc  = acosd(dot(u_true, u_gcc ));

fprintf('\n=== DOA results ===\n');
fprintf('True: az = %.2f deg, el = %.2f deg\n', az_real_w, el_real);
fprintf('Estimated (PHAT): az = %.2f deg, el = %.2f deg, angular error = %.2f deg\n', ...
    az_hat_phat, el_hat_phat, ang_error_phat);
fprintf('Estimated (GCC):  az = %.2f deg, el = %.2f deg, angular error = %.2f deg\n', ...
    az_hat_gcc, el_hat_gcc, ang_error_gcc);

fprintf('\nAll figures saved to: %s\n', outputDir);

%% Comparando Sinal antes e depois da mascara
gcc_phat_mask_effect(x_processed(:,j), x_processed(:,i), fs, maxLagSamples, 0.1);

%% Functions
function pos = pyramid_array_positions(L)
    pos = zeros(4,3);
    pos(1,:) = [0,   0,              0];
    pos(2,:) = [L,   0,              0];
    pos(3,:) = [L/2, L*sqrt(3)/2,    0];
    cx = L/2;
    cy = L*sqrt(3)/6;
    H_tet = L * sqrt(2/3);
    pos(4,:) = [cx, cy, H_tet];
end

function pos = sph2cart_deg(r, az_deg, el_deg)
    pos = [r*cosd(el_deg)*cosd(az_deg), ...
           r*cosd(el_deg)*sind(az_deg), ...
           r*sind(el_deg)];
end

function dmax = max_pair_distance(P)
    dmax = 0;
    for i = 1:size(P,1)
        for j = i+1:size(P,1)
            dmax = max(dmax, norm(P(i,:)-P(j,:)));
        end
    end
end

function y = apply_delay_seconds(x, fs, tau)
    N = length(x);
    n = (0:N-1)';
    y = interp1(n, x, n - tau*fs, 'linear', 0);
end

function y = add_awgn_snr(x, SNR_dB)
    px = mean(x.^2);
    pn = px / (10^(SNR_dB/10));
    y = x + sqrt(pn)*randn(size(x));
end

function J = doa_cost_grid(pos_h, pares, tdoa_med, c, az_vec, el_vec)
    Np = size(pares,1);
    J = zeros(length(el_vec), length(az_vec));

    for iaz = 1:length(az_vec)
        az = az_vec(iaz);
        for iel = 1:length(el_vec)
            el = el_vec(iel);
            u = [cosd(el)*cosd(az), cosd(el)*sind(az), sind(el)];
            tdoa_th = zeros(Np,1);

            for p = 1:Np
                i = pares(p,1);
                j = pares(p,2);
                dp = pos_h(j,:) - pos_h(i,:);
                tdoa_th(p) = -dot(dp, u)/c;
            end

            e = tdoa_med - tdoa_th;
            J(iel, iaz) = sum(e.^2);
        end
    end
end

function [az_hat, el_hat] = peak_min_cost(J, az_vec, el_vec)
    [~, idx] = min(J(:));
    [iel, iaz] = ind2sub(size(J), idx);
    az_hat = az_vec(iaz);
    el_hat = el_vec(iel);
end

function [r, lags_s] = gcc_family_local(x, xref, fs, maxLagSamples, mode, thrFrac)
    x    = x(:)    - mean(x);
    xref = xref(:) - mean(xref);

    N     = length(x);
    Ncorr = 2*N - 1;
    Nfft  = 2^nextpow2(Ncorr);

    X  = fft(x,    Nfft);
    XR = fft(xref, Nfft);

    R = X .* conj(XR);
    mag = abs(R);

    thr  = thrFrac * max(mag);
    mask = (mag >= thr);

    switch lower(mode)
        case 'gcc'
            W = mag ./ (mag);
            Rw = R .* mask .* W;

        case 'phat'
            Rw = (R ./ (mag)) .* mask;

        otherwise
            error('Unknown mode. Use ''gcc'' or ''phat''.');
    end

    cc = ifft(Rw, 'symmetric');
    cc = fftshift(cc);

    lags = (-Nfft/2 : Nfft/2-1).';
    idx  = (lags >= -maxLagSamples) & (lags <= maxLagSamples);

    r      = cc(idx);
    lags_s = lags(idx) / fs;
end

function gcc_phat_mask_effect(x, xref, fs, maxLagSamples, thrFrac, outputDir)
    if nargin < 6
        outputDir = './figures';
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Dimensões padrão 2:1
    fig_width  = 1000;
    fig_height = 500;

    x    = x(:) - mean(x);
    xref = xref(:) - mean(xref);

    N     = length(x);
    Ncorr = 2*N - 1;
    Nfft  = 2^nextpow2(Ncorr);

    X  = fft(x, Nfft);
    XR = fft(xref, Nfft);

    R = X .* conj(XR);
    mag = abs(R);

    thr  = thrFrac * max(mag);
    mask = (mag >= thr);

    mag_mask = mag .* mask;

    W_gcc = mag ./ (mag + eps);
    Rw_gcc_no_mask = R .* W_gcc;
    mag_gcc_no_mask = abs(Rw_gcc_no_mask);
    cc_gcc_no_mask = ifft(Rw_gcc_no_mask, 'symmetric');
    cc_gcc_no_mask = fftshift(cc_gcc_no_mask);

    Rw_gcc_mask = R .* mask .* W_gcc;
    mag_gcc_mask = abs(Rw_gcc_mask);
    cc_gcc_mask = ifft(Rw_gcc_mask, 'symmetric');
    cc_gcc_mask = fftshift(cc_gcc_mask);

    Rw_phat_no_mask = R ./ (mag + eps);
    mag_phat_no_mask = abs(Rw_phat_no_mask);
    cc_phat_no_mask = ifft(Rw_phat_no_mask, 'symmetric');
    cc_phat_no_mask = fftshift(cc_phat_no_mask);

    Rw_phat_mask = (R ./ (mag + eps)) .* mask;
    mag_phat_mask = abs(Rw_phat_mask);
    cc_phat_mask = ifft(Rw_phat_mask, 'symmetric');
    cc_phat_mask = fftshift(cc_phat_mask);

    lags = (-Nfft/2 : Nfft/2-1).';
    idx  = (lags >= -maxLagSamples) & (lags <= maxLagSamples);
    lags_crop = lags(idx) / fs * 1000;

    cc_gcc_no_mask_crop = cc_gcc_no_mask(idx);
    cc_gcc_mask_crop = cc_gcc_mask(idx);

    cc_phat_no_mask_crop = cc_phat_no_mask(idx);
    cc_phat_mask_crop = cc_phat_mask(idx);

    cc_gcc_no_mask_crop = abs(cc_gcc_no_mask_crop) / max(abs(cc_gcc_no_mask_crop));
    cc_gcc_mask_crop = abs(cc_gcc_mask_crop) / max(abs(cc_gcc_mask_crop));

    cc_phat_no_mask_crop = abs(cc_phat_no_mask_crop) / max(abs(cc_phat_no_mask_crop));
    cc_phat_mask_crop = abs(cc_phat_mask_crop) / max(abs(cc_phat_mask_crop));

    freqs = (0:Nfft-1)' * fs / Nfft;
    idx_pos = freqs <= fs/2;
    freqs_pos = freqs(idx_pos);

    color_raw = [0.3 0.3 0.3];
    color_gcc = [0 0.4470 0.7410];
    color_phat = [0.8500 0.3250 0.0980];
    color_mask = [0.6 0.6 0.6];

    %% Figure (a) - Proporção 2:1
    fig = figure('Position', [100 100 1500 750], 'Color', 'w');
    hold on
    h1 = plot(freqs_pos, mag(idx_pos), 'Color', color_raw, 'LineWidth', 2);
    h2 = plot([0 min(2000, fs/2)], [thr thr], 'r--', 'LineWidth', 1.5);
    grid on
    xlabel('Frequency [Hz]', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Magnitude', 'FontSize', 20, 'FontWeight', 'bold')
    xlim([200 1000])
% 	ylim([0 2])
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    legend([h1 h2], {'|R|', 'Threshold'}, 'Location', 'best', 'Interpreter', 'none');
    print(fullfile(outputDir, 'mask_a.eps'), '-depsc')
    %close(fig)

    %% Figure (b) - Proporção 2:1
    fig = figure('Position', [100 100 1500 750], 'Color', 'w');
    hold on
    h1 = plot(freqs_pos, mag_mask(idx_pos), 'Color', color_mask, 'LineWidth', 1.5);
    h2 = plot([0 min(2000, fs/2)], [thr thr], 'r--', 'LineWidth', 1.5);
    grid on
    xlabel('Frequency [Hz]', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Magnitude', 'FontSize', 20, 'FontWeight', 'bold')
    xlim([200 1000])
%     ylim([0 2])
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    legend([h1 h2], {'|R| masked', 'Threshold'}, 'Location', 'best', 'Interpreter', 'none');
    print(fullfile(outputDir, 'mask_b.eps'), '-depsc')
    %close(fig)

    %% Figure (c) - Proporção 2:1
    fig = figure('Position', [100 100 1500 750], 'Color', 'w');
    hold on
    h1 = plot(freqs_pos, mag_gcc_no_mask(idx_pos), 'Color', color_mask, 'LineWidth', 2, 'LineStyle', '-');
    h2 = plot(freqs_pos, mag_gcc_mask(idx_pos), 'Color', color_gcc, 'LineWidth', 2, 'LineStyle', '--');
    grid on
    xlabel('Frequency [Hz]', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Magnitude', 'FontSize', 20, 'FontWeight', 'bold')
    xlim([200 1000])
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    legend([h1 h2], {'No mask', 'With mask'}, 'Location', 'best', 'Interpreter', 'none');
    print(fullfile(outputDir, 'mask_c.eps'), '-depsc')
    %close(fig)

    %% Figure (d) - Proporção 2:1
    fig = figure('Position', [100 100 1500 750], 'Color', 'w');
    hold on
    h1 = plot(freqs_pos, mag_phat_no_mask(idx_pos), 'Color', color_mask, 'LineWidth', 2, 'LineStyle', '-');
    h2 = plot(freqs_pos, mag_phat_mask(idx_pos), 'Color', color_phat, 'LineWidth', 2, 'LineStyle', '--');
    grid on
    xlabel('Frequency [Hz]', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Magnitude', 'FontSize', 20, 'FontWeight', 'bold')
    xlim([200 1000])
%     ylim([0 2])
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    legend([h1 h2], {'No mask', 'With mask'}, 'Location', 'best', 'Interpreter', 'none');
    print(fullfile(outputDir, 'mask_d.eps'), '-depsc')
    %close(fig)

    %% Figure (e) - Proporção 2:1
    fig = figure('Position', [100 100 1200 600], 'Color', 'w');
    hold on
    h1 = plot(lags_crop, abs(cc_gcc_no_mask_crop), 'Color', color_mask, 'LineWidth', 2, 'LineStyle', '-');
    h2 = plot(lags_crop, abs(cc_gcc_mask_crop), 'Color', color_gcc, 'LineWidth', 2, 'LineStyle', '--');
    [peak_gcc_no, idx_gcc_no] = max(abs(cc_gcc_no_mask_crop));
    [peak_gcc_mask, idx_gcc_mask] = max(abs(cc_gcc_mask_crop));
    plot(lags_crop(idx_gcc_no), peak_gcc_no, 'o', 'Color', color_mask, 'MarkerSize', 8, 'LineWidth', 2);
    plot(lags_crop(idx_gcc_mask), peak_gcc_mask, 'o', 'Color', color_gcc, 'MarkerSize', 8, 'LineWidth', 2);
    grid on
    xlabel('Lag [ms]', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Normalized Correlation', 'FontSize', 20, 'FontWeight', 'bold')
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    legend([h1 h2], {'No mask', 'With mask'}, 'Location', 'best', 'Interpreter', 'none');
    print(fullfile(outputDir, 'mask_e.eps'), '-depsc')
    %close(fig)

    %% Figure (f) - Proporção 2:1
    fig = figure('Position', [100 100 1200 600], 'Color', 'w');
    hold on
    h1 = plot(lags_crop, abs(cc_phat_no_mask_crop), 'Color', color_mask, 'LineWidth', 2, 'LineStyle', '-');
    h2 = plot(lags_crop, abs(cc_phat_mask_crop), 'Color', color_phat, 'LineWidth', 2, 'LineStyle', '--');
    [peak_phat_no, idx_phat_no] = max(abs(cc_phat_no_mask_crop));
    [peak_phat_mask, idx_phat_mask] = max(abs(cc_phat_mask_crop));
    plot(lags_crop(idx_phat_no), peak_phat_no, 'o', 'Color', color_mask, 'MarkerSize', 8, 'LineWidth', 2);
    plot(lags_crop(idx_phat_mask), peak_phat_mask, 'o', 'Color', color_phat, 'MarkerSize', 8, 'LineWidth', 2);
    grid on
    xlabel('Lag [ms]', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('Normalized Correlation', 'FontSize', 20, 'FontWeight', 'bold')
    set(gca, 'FontSize', 20)
    set(gca, 'Position', [0.08 0.18 0.90 0.75])
    legend([h1 h2], {'No mask', 'With mask'}, 'Location', 'best', 'Interpreter', 'none');
    print(fullfile(outputDir, 'mask_f.eps'), '-depsc')
    %close(fig)

    %% Summary
    fprintf('\n========== MASK EFFECT SUMMARY ==========\n');
    fprintf('Threshold: %.2f × max(|R|)\n', thrFrac);
    fprintf('Frequency bins retained: %d/%d (%.1f%%)\n\n', sum(mask), length(mask), 100*sum(mask)/length(mask));

    fprintf('--- GCC ---\n');
    fprintf('  WITHOUT mask: peak at %.2f ms (amplitude: %.4f)\n', lags_crop(idx_gcc_no), peak_gcc_no);
    fprintf('  WITH mask:    peak at %.2f ms (amplitude: %.4f)\n', lags_crop(idx_gcc_mask), peak_gcc_mask);
    fprintf('  Enhancement:  %.1f%%\n\n', 100*(peak_gcc_mask/peak_gcc_no - 1));

    fprintf('--- PHAT ---\n');
    fprintf('  WITHOUT mask: peak at %.2f ms (amplitude: %.4f)\n', lags_crop(idx_phat_no), peak_phat_no);
    fprintf('  WITH mask:    peak at %.2f ms (amplitude: %.4f)\n', lags_crop(idx_phat_mask), peak_phat_mask);
    fprintf('  Enhancement:  %.1f%%\n', 100*(peak_phat_mask/peak_phat_no - 1));
    fprintf('============================================\n\n');

    fprintf('Mask figures saved to: %s\n', outputDir);
end

