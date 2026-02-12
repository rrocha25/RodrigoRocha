clear; clc; close all;

%% 1) Parameters
fs = 44100;
c = 1500;

dur_s = 80e-3;
t = (0:round(dur_s*fs)-1)'/fs;

SNR_dB = 10;

rng(42);

% Optional bandpass filter
apply_bandpass = false;  % true/false

% Geometry
fmax_geom = 850;
lambda_min = c / fmax_geom;
L = 0.5 * lambda_min;

pos_h = pyramid_array_positions(L);

% Simulated source
az_real = 35;
el_real = 24;
r_real  = 10;

pos_f = sph2cart_deg(r_real, az_real, el_real);

% Margins and spacings for custom subplot layout
left_margin   = 0.05;
right_margin  = 0.02;
bottom_margin = 0.15;
top_margin    = 0.10;
h_spacing     = 0.05;
v_spacing     = 0.05;

%% 2) Physical maxLag
dmax = max_pair_distance(pos_h);
maxLagSamples = ceil((dmax/c)*fs) + 2;
fprintf('\nmaxLagSamples = %d samples  %.3f ms\n', maxLagSamples, 1000*maxLagSamples/fs);

%% 3) Simple 2-tone signal
f1 = 400;
f2 = 800;
A1 = 1.0;
A2 = 0.7;
phi1 = 0;
phi2 = 0;

s = A1*sin(2*pi*f1*t + phi1) + A2*sin(2*pi*f2*t + phi2);
% s = s .* hann(length(s));
% s = s - mean(s);

% Plot Sinal
% Calcular FFT
N = length(s);
S_fft = fft(s);
S_fft = S_fft(1:floor(N/2)+1);
f_fft = (0:floor(N/2)) * fs / N;

% Calcular magnitude em dB
magnitude_dB = 20*log10(abs(S_fft) + eps);

% Criar figura
figure('Position', [100 100 1400 900])

% Definir margens
left_margin = 0.05;
right_margin = 0.02;
top_margin = 0.05;
bottom_margin = 0.05;

% Calcular dimensões do plot
plot_width  = 1 - left_margin - right_margin;
plot_height = 1 - top_margin - bottom_margin;

% Criar axes com posição customizada
axes('Position', [left_margin, bottom_margin, plot_width, plot_height]);

% Índice de zoom
idx_zoom = f_fft <= 1000;

% Plot
plot(f_fft(idx_zoom), magnitude_dB(idx_zoom), 'b-', 'LineWidth', 1.5)
grid on
xlabel('Frequency [Hz]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 12, 'FontWeight', 'bold')
title('Spectrum of 2-Tone Signal (400 Hz and 800 Hz)', 'FontSize', 14, 'FontWeight', 'bold')
xlim([0 1000])
ylim([-50 100])

%% 4) Simulate reception at 4 hydrophones with relative delays
N = length(s);
dist = zeros(4,1);
for k = 1:4
    dist(k) = norm(pos_f - pos_h(k,:));
end
tau_rel = (dist - dist(1))/c; % s

guard = maxLagSamples + 4;
s_pad = [zeros(guard,1); s; zeros(guard,1)];

x = zeros(N,4);
for k = 1:4
    sk  = apply_delay_seconds(s_pad, fs, tau_rel(k));
    seg = sk(guard+1:guard+N);
    seg = add_awgn_snr(seg, SNR_dB);
    x(:,k) = seg;
end

%% Optional bandpass filtering (before correlations)
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

%% Plot: 4 array channels (time-domain)
t_ref_ms = 37.89;
figure('Position',[120 80 950 700], 'Color','w');

subplot_width  = 1 - left_margin - right_margin;
subplot_height = (1 - top_margin - bottom_margin - (3 * v_spacing)) / 4;

for k = 1:4
    current_bottom = bottom_margin + (4 - k) * (subplot_height + v_spacing);
    ax = subplot('Position', [left_margin, current_bottom, subplot_width, subplot_height]);

    plot(1e3*t, x(:,k), 'k', 'LineWidth', 1.0); grid on; hold on
    ylabel(sprintf('Amplitude $x_{%d}$', k), 'Interpreter','latex');
    title(sprintf('Channel %d true delay relative to channel 1: %.3f ms', k, 1e3*tau_rel(k)), ...
        'Interpreter','latex');

    yl = ylim(ax);
    t_del_ms = t_ref_ms + 1e3*tau_rel(k);
    h_ref = plot([t_ref_ms t_ref_ms], yl, 'b--', 'LineWidth', 1.2);
    h_exp = plot([t_del_ms t_del_ms], yl, 'r-',  'LineWidth', 1.2);
    ylim(ax, yl);
    xlim([30 50]);

    legend([h_ref h_exp], {sprintf('Reference marker: %.2f ms', t_ref_ms), ...
                           sprintf('Expected marker:  %.2f ms', t_del_ms)}, ...
          'Location','best', 'Interpreter','latex');

    if k == 4
        xlabel('Time (ms)','Interpreter','latex');
    else
        set(ax,'XTickLabel',[]);
    end
    hold off
end

annotation('textbox',[0 0.975 1 0.03], ...
  'String','Simulated time-domain signals at the four hydrophones', ...
  'Interpreter','none', 'EdgeColor','none', 'HorizontalAlignment','center');

%% 5) Estimate TDOAs and plot correlations
pares = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
Np = size(pares,1);
u_real = [cosd(el_real)*cosd(az_real) cosd(el_real)*sind(az_real) sind(el_real)];

figure('Position',[120 80 1200 900], 'Color','w');

nRows = Np;
nCols = 2;
subplot_width_corr  = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
subplot_height_corr = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;

tdoa_phat     = zeros(Np,1);
tdoa_gcc      = zeros(Np,1);
tdoa_th_plane = zeros(Np,1);

for p = 1:Np
    i = pares(p,1);
    j = pares(p,2);

    row = p;
    current_bottom = bottom_margin + (nRows - row) * (subplot_height_corr + v_spacing);

    % Plane-wave theoretical TDOA
    dp = pos_h(j,:) - pos_h(i,:);
    tdoa_th_plane(p) = -dot(dp, u_real)/c;

    %% LEFT COLUMN: GCC-PHAT
    col = 1;
    current_left = left_margin + (col - 1) * (subplot_width_corr + h_spacing);
    ax = subplot('Position', [current_left, current_bottom, subplot_width_corr, subplot_height_corr]);

    [rP,lagsP] = gcc_family_local(x_processed(:,j), x_processed(:,i), fs, maxLagSamples, 'phat', 0.1, 0.05);
    [rPpk, kP] = max(abs(rP));          
    tdoa_phat(p) = lagsP(kP);


    plot(1e3*lagsP, rP, 'k', 'LineWidth', 1.05); grid on; hold on
    yl = ylim;
    h_hat = plot(1e3*[tdoa_phat(p) tdoa_phat(p)], yl, 'm-',  'LineWidth', 1.2);
    h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)], yl, 'b--', 'LineWidth', 1.1);
    h_pk  = plot(1e3*tdoa_phat(p), rP(kP), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
    hold off

    xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
    if p == 1, title('GCC-PHAT', 'Interpreter','none'); end
    ylabel(sprintf('pair (%d,%d)', i, j), 'Interpreter','none');

    if p == Np
        xlabel('Delay (ms)', 'Interpreter','none');
    else
        set(ax,'XTickLabel',[]);
    end

    text(0.01, 0.84, sprintf('tau_hat=%.3f ms | tau_th=%.3f ms | |r|max=%.3f', ...
         1e3*tdoa_phat(p), 1e3*tdoa_th_plane(p), rPpk), ...
         'Units','normalized','Interpreter','none','FontSize',9);

    if p == 1
        legend([h_hat h_the h_pk], {'tau_hat','tau_th','peak'}, ...
            'Location','best', 'Interpreter','none');
    end

    %% RIGHT COLUMN: GCC
    col = 2;
    current_left = left_margin + (col - 1) * (subplot_width_corr + h_spacing);
    ax = subplot('Position', [current_left, current_bottom, subplot_width_corr, subplot_height_corr]);

    [rG,lagsG] = gcc_family_local(x_processed(:,j), x_processed(:,i), fs, maxLagSamples, 'gcc',  0.1, 0.05);
    [rGpk, kG] = max(abs(rG));          % <-- defines rGpk
    tdoa_gcc(p) = lagsG(kG);


    plot(1e3*lagsG, rG, 'k', 'LineWidth', 1.05); grid on; hold on
    yl = ylim;
    h_hat = plot(1e3*[tdoa_gcc(p) tdoa_gcc(p)], yl, 'm-',  'LineWidth', 1.2);
    h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)], yl, 'b--', 'LineWidth', 1.1);
    h_pk  = plot(1e3*tdoa_gcc(p), rG(kG), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
    hold off

    xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
    if p == 1, title('GCC', 'Interpreter','none'); end
    set(ax,'YTickLabel',[]);

    if p == Np
        xlabel('Delay (ms)', 'Interpreter','none');
    else
        set(ax,'XTickLabel',[]);
    end

    text(0.01, 0.84, sprintf('tau_hat=%.3f ms | tau_th=%.3f ms | |rho|max=%.3f', ...
         1e3*tdoa_gcc(p), 1e3*tdoa_th_plane(p), rGpk), ...
         'Units','normalized','Interpreter','none','FontSize',9);

    if p == 1
        legend([h_hat h_the h_pk], {'tau_hat','tau_th','peak'}, ...
            'Location','best', 'Interpreter','none');
    end
end

annotation('textbox',[0 0.975 1 0.02], ...
  'String','Each row = a hydrophone pair. Left: GCC-PHAT. Right: GCC. Blue line = plane-wave theoretical delay.', ...
  'Interpreter','none','EdgeColor','none','HorizontalAlignment','center');

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

%% 7) Plot cost maps using the predefined margins (no logic change)
figure('Position', [100, 100, 1200, 800], 'Color', 'w');

% common color axis for both maps (for comparison)
clim_common = [min([J_phat(:); J_gcc(:)]), max([J_phat(:); J_gcc(:)])];

% layout using margins/spacings
nCols = 2;
nRows = 1;
ax_w = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
ax_h = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
ax_bottom = bottom_margin;

ax_left_1 = left_margin;
ax_left_2 = left_margin + ax_w + h_spacing;

% marker style (black outline effect via double plotting)
ms_real = 12; lw_real = 2.2;
ms_est  = 11; lw_est  = 2.4;

% ---- Left axes: GCC-PHAT cost map ----
ax1 = axes('Position',[ax_left_1 ax_bottom ax_w ax_h]); 
imagesc(ax1, az_vec, el_vec, J_phat);
set(ax1, 'YDir', 'normal');
caxis(ax1, clim_common);
colormap(ax1, flipud(jet(256)));  % red = lower cost (visual only)
cb1 = colorbar(ax1); ylabel(cb1,'Cost J (lower is better)','Interpreter','none');
hold(ax1,'on');

% True DOA: white circle with black edge
h_true = plot(ax1, az_real_w, el_real, 'o', ...
    'MarkerSize', ms_real, 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth', lw_real);

% Estimated (PHAT): red plus with black outline (two layers)
plot(ax1, az_hat_phat, el_hat_phat, '+', ...
    'Color','k', 'MarkerSize', ms_est+2, 'LineWidth', lw_est+1.0);
h_est = plot(ax1, az_hat_phat, el_hat_phat, '+', ...
    'Color','r', 'MarkerSize', ms_est, 'LineWidth', lw_est);

hold(ax1,'off');
title(ax1,'GCC-PHAT cost map','Interpreter','none');
xlabel(ax1,'Azimuth (deg)','Interpreter','none');
ylabel(ax1,'Elevation (deg)','Interpreter','none');
legend(ax1, [h_true h_est], {'True DOA','Estimated DOA'}, 'Location','best', 'Interpreter','none');

% ---- Right axes: GCC cost map ----
ax2 = axes('Position',[ax_left_2 ax_bottom ax_w ax_h])
imagesc(ax2, az_vec, el_vec, J_gcc);
set(ax2, 'YDir', 'normal');
caxis(ax2, clim_common);
colormap(ax2, flipud(jet(256)));  % red = lower cost (visual only)
cb2 = colorbar(ax2); ylabel(cb2,'Cost J (lower is better)','Interpreter','none');
hold(ax2,'on');

h_true2 = plot(ax2, az_real_w, el_real, 'o', ...
    'MarkerSize', ms_real, 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth', lw_real);

% Estimated (GCC): green plus with black outline (two layers)
plot(ax2, az_hat_gcc, el_hat_gcc, '+', ...
    'Color','k', 'MarkerSize', ms_est+2, 'LineWidth', lw_est+1.0);
h_est2 = plot(ax2, az_hat_gcc, el_hat_gcc, '+', ...
    'Color','g', 'MarkerSize', ms_est, 'LineWidth', lw_est);

hold(ax2,'off');
title(ax2,'GCC cost map','Interpreter','none');
xlabel(ax2,'Azimuth (deg)','Interpreter','none');
ylabel(ax2,'Elevation (deg)','Interpreter','none');
legend(ax2, [h_true2 h_est2], {'True DOA','Estimated DOA'}, 'Location','best', 'Interpreter','none');

% Global annotation title (avoid sgtitle/subtitle)
if apply_bandpass
    filter_status = 'Bandpass: ON (300-1000 Hz)';
else
    filter_status = 'Bandpass: OFF';
end

annotation('textbox',[0 0.965 1 0.03], ...
    'String', sprintf('DOA cost maps (red = lower cost). %s', filter_status), ...
    'EdgeColor','none', 'HorizontalAlignment','center', 'Interpreter','none', ...
    'FontSize', 13, 'FontWeight','bold');

%% Numerical results (3D angular error)
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


%% ========================================================================
% Functions

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

function [r, lags_s] = gcc_family_local(x, xref, fs, maxLagSamples, mode, thrFrac, alphaFrac)
%GCC_FAMILY_LOCAL  Unified GCC / GCC-PHAT with identical processing pipeline.
%
% Inputs:
%   x, xref        : signals (vectors)
%   fs             : sampling rate (Hz)
%   maxLagSamples  : maximum lag (samples)
%   mode           : 'gcc' or 'phat'
%   thrFrac        : threshold fraction for spectral masking (e.g., 0.1)
%   alphaFrac      : regularization fraction for alpha (e.g., 0.05)
%
% Outputs:
%   r              : correlation values for lags within ±maxLagSamples
%   lags_s         : corresponding lags in seconds

    % 1) DC removal
    x    = x(:)    - mean(x);
    xref = xref(:) - mean(xref);

    % 2) FFT sizing for linear correlation
    N     = length(x);
    Ncorr = 2*N - 1;
    Nfft  = 2^nextpow2(Ncorr);

    % 3) FFTs
    X  = fft(x,    Nfft);
    XR = fft(xref, Nfft);

    % 4) Cross-spectrum
    R = X .* conj(XR);
    mag = abs(R);

    % 5) Same robustness: spectral masking + regularization
    thr  = thrFrac * max(mag);
    mask = (mag >= thr);

    alpha = alphaFrac * median(mag);

    % ------------------- ONLY DIFFERENCE IS THIS WEIGHTING -------------------
    switch lower(mode)
        case 'gcc'
            % "Robust GCC": keep magnitude, just suppress weak/noisy bins smoothly
            W = mag ./ (mag + alpha);          % soft attenuation (0..1)
            Rw = R .* mask .* W;

        case 'phat'
            % "Robust PHAT": normalize by magnitude (phase-based), regularized
            Rw = (R ./ (mag + alpha)) .* mask;

        otherwise
            error('Unknown mode. Use ''gcc'' or ''phat''.');
    end
    % ------------------------------------------------------------------------

    % 6) Back to time-lag domain
    cc = ifft(Rw, 'symmetric');
    cc = fftshift(cc);

    % 7) Lags and physical cropping
    lags = (-Nfft/2 : Nfft/2-1).';
    idx  = (lags >= -maxLagSamples) & (lags <= maxLagSamples);

    r      = cc(idx);
    lags_s = lags(idx) / fs;
end

