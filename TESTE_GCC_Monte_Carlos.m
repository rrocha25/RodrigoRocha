clear; clc; close all;

%% Parâmetros gerais
fs = 44100;              % Hz
c  = 1500;               % m/s
dur_s = 80e-3;           % 80 ms
t = (0:round(dur_s*fs)-1)'/fs;

% Geometria do arranjo
fmax_geom  = 850;        % Hz
lambda_min = c / fmax_geom;
L = 0.5 * lambda_min;
H = L;
pos_h = pyramid_array_positions(L, H);

% Frequências do sinal
f1 = 400;  A1 = 1.0;  phi1 = 0;
f2 = 800;  A2 = 0.7;  phi2 = 0;

% maxLag físico
dmax = max_pair_distance(pos_h);
maxLagSamples = ceil((dmax/c)*fs) + 2;

% Pares de hidrofones
pares = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];

% Distância fixa (far-field)
r_real = 10;             % m

% Grid de busca
az_vec = 0:2:358;   % passo de 2°
el_vec = 0:2:90;    % passo de 2°

% Número de realizações Monte Carlo
N_MC = 50;

% Sinal limpo
s = A1*sin(2*pi*f1*t + phi1) + A2*sin(2*pi*f2*t + phi2);
s = s .* hann(length(s));
s = s - mean(s);
%% ------------------------------------------------------------------------
% Experimento 1: Robustez vs SNR (direção fixa)
% -------------------------------------------------------------------------
SNR_list = -20:2:20;    % dB
N_SNR    = numel(SNR_list);

az_fix = 35;             % graus
el_fix = 24;             % graus
pos_f_fix = sph2cart_deg(r_real, az_fix, el_fix);

erro_snr_mean = zeros(N_SNR,1);
erro_snr_std  = zeros(N_SNR,1);

fprintf('\n=== Experimento 1: Erro vs SNR (az = %.1f, el = %.1f) ===\n', az_fix, el_fix);

for iS = 1:N_SNR
    SNR_dB = SNR_list(iS);
    erros = zeros(N_MC,1);

    for mc = 1:N_MC
        % Simular captação
        x = simulate_array_signals(s, fs, c, pos_h, pos_f_fix, maxLagSamples, SNR_dB);

        % Estimar TDOAs
        tdoa = estimate_tdoas(x, pares, fs, maxLagSamples);

        % Estimar DOA
        J = doa_cost_grid(pos_h, pares, tdoa, c, az_vec, el_vec);
        [az_hat, el_hat] = peak_min_cost(J, az_vec, el_vec);

        % Erro angular
        erros(mc) = angular_error(az_fix, el_fix, az_hat, el_hat);
    end

    erro_snr_mean(iS) = mean(erros);
    erro_snr_std(iS)  = std(erros);

    fprintf('  SNR = %+3d dB: erro médio = %.2f°, desvio padrão = %.2f°\n', ...
        SNR_dB, erro_snr_mean(iS), erro_snr_std(iS));
end

%% Figura: Erro vs SNR (com subplot)
figure('Position',[100 100 700 800]);

% Subplot 1: SNR completo (-20 a 20 dB)
subplot(2,1,1);
errorbar(SNR_list, erro_snr_mean, erro_snr_std, 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Erro angular (graus)', 'Interpreter', 'latex');
title(sprintf('Erro de localização vs SNR - Faixa completa (az = %d graus, el = %d graus, N = %d)', ...
    round(az_fix), round(el_fix), N_MC), 'Interpreter', 'none');
legend('Erro medio ± desvio padrao', 'Location', 'best');
set(gca, 'FontSize', 11);

% Subplot 2: SNR de -10 a 20 dB (zoom)
subplot(2,1,2);
% Filtrar apenas os dados de SNR >= -10
idx_zoom = SNR_list >= -10;
errorbar(SNR_list(idx_zoom), erro_snr_mean(idx_zoom), erro_snr_std(idx_zoom), 'o-', ...
    'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Erro angular (graus)', 'Interpreter', 'latex');
title('Erro de localização vs SNR - Zoom (SNR >= -10 dB)', 'Interpreter', 'none');
legend('Erro medio ± desvio padrao', 'Location', 'best');
set(gca, 'FontSize', 11);
xlim([-10 20]);


%% ------------------------------------------------------------------------
% Experimento 2: Robustez vs direção (SNR fixo)
% -------------------------------------------------------------------------
SNR_dir = 0;             % dB (pode trocar para 10 dB se quiser)
az_list = 0:60:359;      % azimutes
el_list = [0, 30, 60 , 90];  % elevações

N_az = numel(az_list);
N_el = numel(el_list);

erro_dir_mean = zeros(N_el, N_az);
erro_dir_std  = zeros(N_el, N_az);

s = A1*sin(2*pi*f1*t + phi1) + A2*sin(2*pi*f2*t + phi2);
s = s .* hann(length(s));
s = s - mean(s);

fprintf('\n=== Experimento 2: Erro vs direção (SNR = %d dB) ===\n', SNR_dir);
         
for ie = 1:N_el
    for ia = 1:N_az
        az = az_list(ia);
        el = el_list(ie);
        pos_f = sph2cart_deg(r_real, az, el);
        erros = zeros(N_MC,1);

        for mc = 1:N_MC
            x = simulate_array_signals(s, fs, c, pos_h, pos_f, maxLagSamples, SNR_dir);
            tdoa = estimate_tdoas(x, pares, fs, maxLagSamples);
            J = doa_cost_grid(pos_h, pares, tdoa, c, az_vec, el_vec);
            [az_hat, el_hat] = peak_min_cost(J, az_vec, el_vec);

            erros(mc) = angular_error(az, el, az_hat, el_hat);
        end

        erro_dir_mean(ie, ia) = mean(erros);
        erro_dir_std(ie, ia)  = std(erros);

        fprintf('  az = %3d°, el = %2d°: erro médio = %.2f°, desvio padrão = %.2f°\n', ...
            az, el, erro_dir_mean(ie,ia), erro_dir_std(ie,ia));
    end
end


%% ========== FIGURA: Erro vs azimute para todas as elevações ==========
figure('Position',[100 100 1000 900]);

annotation('textbox', [0 0.95 1 0.05], ...
    'String', sprintf('Erro angular medio vs azimute | SNR = %d dB, N = %d', SNR_dir, N_MC), ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 13, 'FontWeight', 'bold');
N_el = length(el_list);

for ie = 1:N_el
    subplot(N_el, 1, ie);

    % Definir cor: vermelho para último subplot, azul para os demais
    if ie == N_el
        cor = 'r';
        ylim_max = 200;
    else
        cor = 'b';
        ylim_max = 15;
    end

    % Plotar erro médio com barras de desvio padrão
    errorbar(az_list, erro_dir_mean(ie, :), erro_dir_std(ie, :), ...
        'o-', 'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', cor, 'Color', cor);

    grid on;
    xlabel('Azimute (graus)');
    ylabel('Erro angular (graus)');
    title(sprintf('Elevacao = %d graus', el_list(ie)), 'FontWeight', 'bold');
    set(gca, 'FontSize', 10);
    xlim([min(az_list), max(az_list)]);
    ylim([0 ylim_max]);
end




%% ----------------- FUNÇÕES AUXILIARES -----------------
function pos = pyramid_array_positions(L, H)
    pos = zeros(4,3);
    pos(1,:) = [0, 0, 0];
    pos(2,:) = [L, 0, 0];
    pos(3,:) = [L/2, L*sqrt(3)/2, 0];
    cx = L/2; cy = L*sqrt(3)/6;
    pos(4,:) = [cx, cy, H];
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

function x = simulate_array_signals(s, fs, c, pos_h, pos_f, maxLagSamples, SNR_dB)
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
        sk = apply_delay_seconds(s_pad, fs, tau_rel(k));
        seg = sk(guard+1:guard+N);
        seg = add_awgn_snr(seg, SNR_dB);
        x(:,k) = seg;
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

function tdoa = estimate_tdoas(x, pares, fs, maxLagSamples)
    Np = size(pares,1);
    tdoa = zeros(Np,1);
    for p = 1:Np
        i = pares(p,1);
        j = pares(p,2);
        [r, lags] = xcorr(x(:,j)-mean(x(:,j)), x(:,i)-mean(x(:,i)), maxLagSamples, 'coeff');
        [~, kpk] = max(abs(r));
        tdoa(p) = lags(kpk)/fs;
    end
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
                i = pares(p,1); j = pares(p,2);
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

function err = angular_error(az_true, el_true, az_est, el_est)
    % Normaliza azimutes para [0, 360)
    az_true = mod(az_true, 360);
    az_est  = mod(az_est, 360);

    % Diferença mínima em azimute (considerando círculo)
    d_az = az_est - az_true;
    d_az = mod(d_az + 180, 360) - 180;  % agora está em [-180, 180]

    % Diferença simples em elevação
    d_el = el_est - el_true;

    % Erro euclidiano no plano (d_az, d_el)
    err = sqrt(d_az.^2 + d_el.^2);
end

