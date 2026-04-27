clear; clc; close all;

%% =========================================================
%% 1) PARÂMETROS GLOBAIS
%% =========================================================

fs   = 44100;   % Taxa de amostragem (Hz)
c    = 1500;    % Velocidade do som na água (m/s)

SNR_dB = -14;   % SNR alvo (dB) — ajuste conforme necessário

% --- Geometria do arranjo ---
fmax_geom  = 1000;
lambda_min = c / fmax_geom;
L          = 0.5 * lambda_min;
pos_h      = pyramid_array_positions(L);

% --- Fonte simulada ---
az_real = 35;   % Azimute real (graus)
el_real = 24;   % Elevação real (graus)
r_real  = 10;   % Distância (m)
pos_f   = sph2cart_deg(r_real, az_real, el_real);

% --- Diretório de saída ---
outputDir = './figures';
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

fig_width  = 1000;
fig_height = 500;

%% =========================================================
%% 2) PARÂMETROS DE VARREDURA
%% =========================================================

% Varredura A: duração de cada segmento
dur_s_list    = [20e-3, 40e-3, 80e-3, 160e-3];
% dur_s_list    = [80e-3];
num_dur_cases = length(dur_s_list);

% Varredura B: realizações de ruído (semente do rng)
num_seg_list0 = 1:1:200;
num_cases0    = length(num_seg_list0);

% Varredura C: número de segmentos por realização
num_seg_list  = [1, 8, 32, 64, 128, 256];
% num_seg_list  = [1];
num_cases     = length(num_seg_list);
max_num_segments = max(num_seg_list);

% Vetor unitário real (fixo para todos os loops)
u_real_vec = [cosd(el_real)*cosd(az_real);
              cosd(el_real)*sind(az_real);
              sind(el_real)];

% Pares de hidrofones (fixo)
pares = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
Np    = size(pares, 1);

%% =========================================================
%% 3) PRÉ-ALOCAÇÃO DAS MATRIZES DE ACÚMULO
%%    Dimensões: (num_cases x num_cases0 x num_dur_cases)
%%    idx_case  = índice do número de segmentos
%%    idx_rng   = índice da realização de ruído
%%    idx_dur   = índice da duração do segmento
%% =========================================================

% DOA mediana por (num_seg, realização, dur_s)
DOA_az_6_median_all = nan(num_cases, num_cases0, num_dur_cases);
DOA_el_6_median_all = nan(num_cases, num_cases0, num_dur_cases);

DOA_az_5_median_all = nan(num_cases, num_cases0, num_dur_cases);
DOA_el_5_median_all = nan(num_cases, num_cases0, num_dur_cases);

DOA_az_4_median_all = nan(num_cases, num_cases0, num_dur_cases);
DOA_el_4_median_all = nan(num_cases, num_cases0, num_dur_cases);

% DOA média por (num_seg, realização, dur_s)
DOA_az_6_mean_all   = nan(num_cases, num_cases0, num_dur_cases);
DOA_el_6_mean_all   = nan(num_cases, num_cases0, num_dur_cases);

DOA_az_5_mean_all   = nan(num_cases, num_cases0, num_dur_cases);
DOA_el_5_mean_all   = nan(num_cases, num_cases0, num_dur_cases);

DOA_az_4_mean_all   = nan(num_cases, num_cases0, num_dur_cases);
DOA_el_4_mean_all   = nan(num_cases, num_cases0, num_dur_cases);

% Erro angular por (num_seg, realização, dur_s)
ang_err_6_median_all = nan(num_cases, num_cases0, num_dur_cases);
ang_err_5_median_all = nan(num_cases, num_cases0, num_dur_cases);
ang_err_4_median_all = nan(num_cases, num_cases0, num_dur_cases);

ang_err_6_mean_all   = nan(num_cases, num_cases0, num_dur_cases);
ang_err_5_mean_all   = nan(num_cases, num_cases0, num_dur_cases);
ang_err_4_mean_all   = nan(num_cases, num_cases0, num_dur_cases);

% DOAs por segmento — cell array para acomodar vetores de tamanhos variáveis
%   az_seg_6_cell{idx_case, idx_rng, idx_dur} = vetor (num_segments x 1)
az_seg_6_cell = cell(num_cases, num_cases0, num_dur_cases);
el_seg_6_cell = cell(num_cases, num_cases0, num_dur_cases);
az_seg_5_cell = cell(num_cases, num_cases0, num_dur_cases);
el_seg_5_cell = cell(num_cases, num_cases0, num_dur_cases);
az_seg_4_cell = cell(num_cases, num_cases0, num_dur_cases);
el_seg_4_cell = cell(num_cases, num_cases0, num_dur_cases);

%% =========================================================
%% 4) LEITURA DO SINAL REAL E DO RUÍDO
%% =========================================================

% ----------------------------------------------------------
% Sinal de interesse (lancha)
% ----------------------------------------------------------
[s_raw, fs_wav] = audioread('H3_lancha2.wav');

% Converte para mono se estéreo
if size(s_raw, 2) > 1
    s_raw = mean(s_raw, 2);
end

% Reamostra para fs do sistema se necessário
if fs_wav ~= fs
    s_raw = resample(s_raw, fs, fs_wav);
    fprintf('Sinal reamostrado de %d Hz para %d Hz\n', fs_wav, fs);
end

s = s_raw(:);   % garante coluna

fprintf('Sinal lido: %d amostras  %.1f ms  fs=%d Hz\n', ...
        length(s), 1e3*length(s)/fs, fs);

% ----------------------------------------------------------
% Ruído de ambiente gravado
% ----------------------------------------------------------
[noise_env_raw, fs_noise] = audioread('Enseada dos Anjos.wav');

% Converte para mono se estéreo
if size(noise_env_raw, 2) > 1
    noise_env_raw = mean(noise_env_raw, 2);
end

% Reamostra para fs do sistema se necessário
if fs_noise ~= fs
    noise_env_raw = resample(noise_env_raw, fs, fs_noise);
    fprintf('Ruido reamostrado de %d Hz para %d Hz\n', fs_noise, fs);
end

noise_env_raw = noise_env_raw(:);

fprintf('Ruido lido: %d amostras  %.1f ms\n', ...
        length(noise_env_raw), 1e3*length(noise_env_raw)/fs);

%% =========================================================
%% 5) maxLag físico
%% =========================================================

dmax          = max_pair_distance(pos_h);
maxLagSamples = ceil((dmax/c)*fs) + 2;

fprintf('\nmaxLagSamples = %d samples  %.3f ms\n', ...
        maxLagSamples, 1000*maxLagSamples/fs);

%% 6) LOOP PRINCIPAL:
%    - sobre durações de segmento (dur_s_list)
%    - sobre realizações de ruído (num_seg_list0)
%    - sobre número de segmentos (num_seg_list)
% =========================================================

for idx_dur = 1:num_dur_cases

    dur_s = dur_s_list(idx_dur);

    N_seg_base = round(dur_s * fs);  % #amostras por segmento

    for idx_rng = 1:num_cases0

        rng(num_seg_list0(idx_rng));   % semente baseada no índice/valor
        fprintf('  Realização de ruído %d (seed = %d) (dur_s = %.1f ms)', idx_rng, num_seg_list0(idx_rng));

        for idx_case = 1:num_cases

            num_segments = num_seg_list(idx_case);
            fprintf('    num_segments = %4d\n', num_segments);

            %% 6.1) Simulação do sinal recebido nos 4 hidrofones

            N_seg   = N_seg_base;
            N_total = N_seg * num_segments;

            x = zeros(N_total, 4);  % [N_total x 4] canais simulados

            % Distâncias / atrasos geométricos (fixos para esta geometria)
            dist = zeros(4,1);
            for k = 1:4
                dist(k) = norm(pos_f - pos_h(k,:));
            end
            tau_rel = (dist - dist(1))/c;  % atrasos relativos ao hidrofone 1

            guard = maxLagSamples + 4;

            for seg_idx = 1:num_segments

                idx_start = (seg_idx-1)*N_seg + 1;
                idx_end   = idx_start + N_seg - 1;

                % --- 6.1.1) Segmento do sinal da lancha ---
                if idx_end > length(s)
                    error('Segmento ultrapassa comprimento de s. Ajuste dur_s ou num_segments.');
                end
                s_seg_raw_current = s(idx_start:idx_end);
                s_seg_raw_current = s_seg_raw_current(:);

                % Normalização por segmento
                max_s_seg = max(abs(s_seg_raw_current));
                if max_s_seg == 0
                    s_seg_normalized = s_seg_raw_current;  % tudo zero
                else
                    s_seg_normalized = s_seg_raw_current / max_s_seg;
                end
                ps_normalized = mean(s_seg_normalized.^2);

                % Aplica atrasos geométricos
                s_seg_pad = [zeros(guard,1); s_seg_normalized; zeros(guard,1)];

                for k = 1:4
                    sk_delayed  = apply_delay_seconds(s_seg_pad, fs, tau_rel(k));
                    seg_delayed_signal = sk_delayed(guard+1:guard+N_seg);

                    % --- 6.1.2) Ruído para este canal/segmento ---
                    if length(noise_env_raw) < N_seg
                        error('Ruído gravado muito curto para formar um segmento.');
                    end

                    i0 = randi(length(noise_env_raw) - N_seg + 1);
                    n_seg_raw_current = noise_env_raw(i0:i0+N_seg-1);
                    n_seg_raw_current = n_seg_raw_current(:);

                    max_n_seg = max(abs(n_seg_raw_current));
                    if max_n_seg == 0
                        n_seg_normalized = n_seg_raw_current;
                    else
                        n_seg_normalized = n_seg_raw_current / max_n_seg;
                    end
                    pn_normalized = mean(n_seg_normalized.^2);

                    % Ajuste de SNR
                    if pn_normalized == 0
                        scale_noise = 0;
                    else
                        scale_noise = sqrt(ps_normalized / (pn_normalized * 10^(SNR_dB/10)));
                    end
                    n_seg_scaled = scale_noise * n_seg_normalized;

                    % Sinal + ruído no canal k
                    x(idx_start:idx_end, k) = seg_delayed_signal + n_seg_scaled;
                end
            end % seg_idx

            %% 6.2) Estimativa de TDOA (GCC-PHAT) por segmento e par

            tdoa_phat = nan(num_segments, Np);  % [seg_idx x par]

            % TDOA teórico para referência (opcional)
            tdoa_th_plane = zeros(Np,1);
            for p = 1:Np
                i = pares(p,1);
                j = pares(p,2);
                dp = pos_h(i,:) - pos_h(j,:);
                tdoa_th_plane(p) = dot(dp, u_real_vec)/c;
            end

            for seg_idx = 1:num_segments

                idx_start = (seg_idx-1)*N_seg + 1;
                idx_end   = idx_start + N_seg - 1;

                x_seg = x(idx_start:idx_end, :);   % [N_seg x 4]

                for p = 1:Np
                    i = pares(p,1);
                    j = pares(p,2);

                    [rP, lagsP] = gcc_family_local(x_seg(:,j), x_seg(:,i), fs, maxLagSamples, 'phat');
                    [~, kP]     = max(abs(rP));
                    tdoa_phat(seg_idx, p) = lagsP(kP);   % em segundos
                end
            end

            %% 6.3) Estimativa de DOA (LS) para 6 / 5 / 4 pares por segmento

            az_seg_6 = nan(num_segments, 1);
            el_seg_6 = nan(num_segments, 1);
            az_seg_5 = nan(num_segments, 1);
            el_seg_5 = nan(num_segments, 1);
            az_seg_4 = nan(num_segments, 1);
            el_seg_4 = nan(num_segments, 1);

            for seg_idx = 1:num_segments

                current_tdoas = tdoa_phat(seg_idx, :).';  % [Np x 1]

                % ---------- 6 pares ----------
                [az6, el6, ~] = calculateDoA_LS(current_tdoas, pares, pos_h, c);
                az_seg_6(seg_idx) = az6;
                el_seg_6(seg_idx) = el6;

                % ---------- 5 pares (remove 1 outlier) ----------
                Np_total     = Np;
                err_sq_out   = zeros(Np_total,1);
                az_sub_5     = zeros(Np_total,1);
                el_sub_5     = zeros(Np_total,1);
                uDoA_sub_5   = zeros(Np_total,3);

                for idx_out = 1:Np_total
                    idx_keep  = setdiff(1:Np_total, idx_out);
                    pairs_sub = pares(idx_keep, :);
                    tdoa_sub  = current_tdoas(idx_keep);
                    tdoa_out  = current_tdoas(idx_out);

                    [az_tmp, el_tmp, u_tmp] = calculateDoA_LS(tdoa_sub, pairs_sub, pos_h, c);

                    az_sub_5(idx_out) = az_tmp;
                    el_sub_5(idx_out) = el_tmp;
                    uDoA_sub_5(idx_out,:) = u_tmp(:).';

                    ii = pares(idx_out,1);
                    jj = pares(idx_out,2);
                    delta_p_ij_over_c = (pos_h(ii,:) - pos_h(jj,:))/c;
                    tau_pred = dot(delta_p_ij_over_c, u_tmp);
                    err_sq_out(idx_out) = (tdoa_out - tau_pred)^2;
                end

                [~, idx_first_outlier] = max(err_sq_out);   % MELHOR “outlier” a descartar
                az_seg_5(seg_idx) = az_sub_5(idx_first_outlier);
                el_seg_5(seg_idx) = el_sub_5(idx_first_outlier);

                % ---------- 4 pares (remove 2 outliers) ----------
                idx_keep_5    = setdiff(1:Np_total, idx_first_outlier);
                pares_5_init  = pares(idx_keep_5, :);
                tdoa_5_init   = current_tdoas(idx_keep_5);
                Np_5_init     = length(idx_keep_5);   % sempre 5

                err_sq_out_2  = zeros(Np_5_init, 1);
                az_sub_4      = zeros(Np_5_init, 1);
                el_sub_4      = zeros(Np_5_init, 1);
                uDoA_sub_4    = zeros(Np_5_init, 3);

                for k2 = 1:Np_5_init

                    idx_keep_4 = setdiff(1:Np_5_init, k2);  % índices dentro de pares_5_init
                    pares_4    = pares_5_init(idx_keep_4, :);
                    tdoa_4     = tdoa_5_init(idx_keep_4);

                    [az4_tmp, el4_tmp, u4_tmp] = calculateDoA_LS(tdoa_4, pares_4, pos_h, c);

                    az_sub_4(k2)   = az4_tmp;
                    el_sub_4(k2)   = el4_tmp;
                    uDoA_sub_4(k2,:) = u4_tmp(:).';

                    % Par removido na segunda etapa
                    idx_out2_original = idx_keep_5(k2);   % índice original (1..Np)
                    tdoa_out2 = current_tdoas(idx_out2_original);

                    ii2 = pares(idx_out2_original,1);
                    jj2 = pares(idx_out2_original,2);
                    delta_p_ij2_over_c = (pos_h(ii2,:) - pos_h(jj2,:))/c;
                    tau_pred2 = dot(delta_p_ij2_over_c, u4_tmp);

                    err_sq_out_2(k2) = (tdoa_out2 - tau_pred2)^2;
                end

                [~, best_k2] = max(err_sq_out_2);
                az_seg_4(seg_idx) = az_sub_4(best_k2);
                el_seg_4(seg_idx) = el_sub_4(best_k2);
            end

            %% 6.4) Guarda DOAs por segmento em cell arrays

            az_seg_6_cell{idx_case, idx_rng, idx_dur} = az_seg_6;
            el_seg_6_cell{idx_case, idx_rng, idx_dur} = el_seg_6;

            az_seg_5_cell{idx_case, idx_rng, idx_dur} = az_seg_5;
            el_seg_5_cell{idx_case, idx_rng, idx_dur} = el_seg_5;

            az_seg_4_cell{idx_case, idx_rng, idx_dur} = az_seg_4;
            el_seg_4_cell{idx_case, idx_rng, idx_dur} = el_seg_4;

            %% 6.5) Calcula média e mediana das DOAs (por segmentos) para cada caso

            % 6 pares
            DOA_az_6_median = median(az_seg_6, 'omitnan');
            DOA_el_6_median = median(el_seg_6, 'omitnan');
            DOA_az_6_mean   = mean( az_seg_6, 'omitnan');
            DOA_el_6_mean   = mean( el_seg_6, 'omitnan');

            DOA_az_6_median_all(idx_case, idx_rng, idx_dur) = DOA_az_6_median;
            DOA_el_6_median_all(idx_case, idx_rng, idx_dur) = DOA_el_6_median;
            DOA_az_6_mean_all(  idx_case, idx_rng, idx_dur) = DOA_az_6_mean;
            DOA_el_6_mean_all(  idx_case, idx_rng, idx_dur) = DOA_el_6_mean;

            % 5 pares
            DOA_az_5_median = median(az_seg_5, 'omitnan');
            DOA_el_5_median = median(el_seg_5, 'omitnan');
            DOA_az_5_mean   = mean( az_seg_5, 'omitnan');
            DOA_el_5_mean   = mean( el_seg_5, 'omitnan');

            DOA_az_5_median_all(idx_case, idx_rng, idx_dur) = DOA_az_5_median;
            DOA_el_5_median_all(idx_case, idx_rng, idx_dur) = DOA_el_5_median;
            DOA_az_5_mean_all(  idx_case, idx_rng, idx_dur) = DOA_az_5_mean;
            DOA_el_5_mean_all(  idx_case, idx_rng, idx_dur) = DOA_el_5_mean;

            % 4 pares
            DOA_az_4_median = median(az_seg_4, 'omitnan');
            DOA_el_4_median = median(el_seg_4, 'omitnan');
            DOA_az_4_mean   = mean( az_seg_4, 'omitnan');
            DOA_el_4_mean   = mean( el_seg_4, 'omitnan');

            DOA_az_4_median_all(idx_case, idx_rng, idx_dur) = DOA_az_4_median;
            DOA_el_4_median_all(idx_case, idx_rng, idx_dur) = DOA_el_4_median;
            DOA_az_4_mean_all(  idx_case, idx_rng, idx_dur) = DOA_az_4_mean;
            DOA_el_4_mean_all(  idx_case, idx_rng, idx_dur) = DOA_el_4_mean;

            %% 6.6) Calcula erro angular entre DOA real e as DOAs média/mediana

            % 6 pares
            u_6_med = azel2unitvec(DOA_az_6_median, DOA_el_6_median);
            u_6_mean = azel2unitvec(DOA_az_6_mean, DOA_el_6_mean);
            ang_err_6_median_all(idx_case, idx_rng, idx_dur) = ang_err_deg(u_real_vec, u_6_med);
            ang_err_6_mean_all(  idx_case, idx_rng, idx_dur) = ang_err_deg(u_real_vec, u_6_mean);

            % 5 pares
            u_5_med = azel2unitvec(DOA_az_5_median, DOA_el_5_median);
            u_5_mean = azel2unitvec(DOA_az_5_mean, DOA_el_5_mean);
            ang_err_5_median_all(idx_case, idx_rng, idx_dur) = ang_err_deg(u_real_vec, u_5_med);
            ang_err_5_mean_all(  idx_case, idx_rng, idx_dur) = ang_err_deg(u_real_vec, u_5_mean);

            % 4 pares
            u_4_med = azel2unitvec(DOA_az_4_median, DOA_el_4_median);
            u_4_mean = azel2unitvec(DOA_az_4_mean, DOA_el_4_mean);
            ang_err_4_median_all(idx_case, idx_rng, idx_dur) = ang_err_deg(u_real_vec, u_4_med);
            ang_err_4_mean_all(  idx_case, idx_rng, idx_dur) = ang_err_deg(u_real_vec, u_4_mean);

        end % idx_case
    end % idx_rng
end % idx_dur

%% =========================================================
%% 7) RESULTADO COMPLETO — todas as combinações de num_seg x dur_s
%% =========================================================

% --- Definição dos métodos ---
num_methods = 6;

method_names = {
    'LS-Mean (6)  ';
    'LS-Mean (5)  ';
    'LS-Mean (4)  ';
    'LS-Median (6)';
    'LS-Median (5)';
    'LS-Median (4)';
};

fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('  RESULTADO COMPLETO — todas as combinacoes seg x ms\n');
fprintf('  az=%.0f deg  el=%.0f deg  SNR=%.0f dB  Realizacoes=%d\n', ...
        az_real, el_real, SNR_dB, num_cases0);
fprintf('%s\n\n', repmat('=', 1, 80));

for id = 1:num_dur_cases
    for is = 1:num_cases

        num_seg_atual = num_seg_list(is);
        dur_ms_atual  = dur_s_list(id) * 1e3;
        update_s      = num_seg_atual * dur_s_list(id);

        fprintf('  --- %d seg @ %.0f ms  (update: %.3f s) ---\n', ...
                num_seg_atual, dur_ms_atual, update_s);

        err_mean6   = squeeze(ang_err_6_mean_all(  is, :, id));
        err_mean5   = squeeze(ang_err_5_mean_all(  is, :, id));
        err_mean4   = squeeze(ang_err_4_mean_all(  is, :, id));
        err_median6 = squeeze(ang_err_6_median_all(is, :, id));
        err_median5 = squeeze(ang_err_5_median_all(is, :, id));
        err_median4 = squeeze(ang_err_4_median_all(is, :, id));

        all_errs = {err_mean6; err_mean5; err_mean4; ...
                    err_median6; err_median5; err_median4};

        fprintf('  %-15s  %10s  %10s\n', 'Method', 'RMSE (deg)', 'Std (deg)');
        fprintf('  %s\n', repmat('-', 1, 40));

        for m = 1:num_methods
            if m == 4
                fprintf('  %s\n', repmat('-', 1, 40));
            end
            e = all_errs{m};
            fprintf('  %-15s  %8.2f      %8.2f\n', ...
                strtrim(method_names{m}), ...
                sqrt(mean(e.^2, 'omitnan')), ...
                std(e, 0, 'omitnan'));
        end
        fprintf('\n');
    end
end

%% 8) SCATTER + ELIPSE DE CONFIANÇA
% =========================================================

target_seg = 1;
target_ms  = 80;
conf_level = 0.15;

idx_seg_target = find(num_seg_list == target_seg, 1);
idx_dur_target = find(abs(dur_s_list - target_ms*1e-3) < 1e-9, 1);

if isempty(idx_seg_target) || isempty(idx_dur_target)
    error('Configuracao %d seg @ %d ms nao encontrada.', target_seg, target_ms);
end

is = idx_seg_target;
id = idx_dur_target;

% ---------------------------------------------------------
% 8.1) Extrai DOAs brutas
% ---------------------------------------------------------
az6 = zeros(num_cases0, 1);  el6 = zeros(num_cases0, 1);
az5 = zeros(num_cases0, 1);  el5 = zeros(num_cases0, 1);
az4 = zeros(num_cases0, 1);  el4 = zeros(num_cases0, 1);

for idx_rng = 1:num_cases0
    az6(idx_rng) = az_seg_6_cell{is, idx_rng, id};
    el6(idx_rng) = el_seg_6_cell{is, idx_rng, id};
    az5(idx_rng) = az_seg_5_cell{is, idx_rng, id};
    el5(idx_rng) = el_seg_5_cell{is, idx_rng, id};
    az4(idx_rng) = az_seg_4_cell{is, idx_rng, id};
    el4(idx_rng) = el_seg_4_cell{is, idx_rng, id};
end

% ---------------------------------------------------------
% 8.2) Erros, centros e elipses
% ---------------------------------------------------------
err6 = [az6 - az_real, el6 - el_real];
err5 = [az5 - az_real, el5 - el_real];
err4 = [az4 - az_real, el4 - el_real];

% Centro = mediana dos erros + DOA real
az6_c = az_real + median(err6(:,1), 'omitnan');
el6_c = el_real + median(err6(:,2), 'omitnan');
az5_c = az_real + median(err5(:,1), 'omitnan');
el5_c = el_real + median(err5(:,2), 'omitnan');
az4_c = az_real + median(err4(:,1), 'omitnan');
el4_c = el_real + median(err4(:,2), 'omitnan');

% Escala chi2 com conf_level
scale = sqrt(chi2inv(conf_level, 2));

% Círculo unitário base
theta       = linspace(0, 2*pi, 300);
unit_circle = [cos(theta); sin(theta)];

% Covariância + Cholesky ? elipse
C6 = cov(err6(all(~isnan(err6),2), :));
C5 = cov(err5(all(~isnan(err5),2), :));
C4 = cov(err4(all(~isnan(err4),2), :));

L6 = chol(C6, 'lower');
L5 = chol(C5, 'lower');
L4 = chol(C4, 'lower');

ell6 = scale * L6 * unit_circle;   % [2 x 300]
ell5 = scale * L5 * unit_circle;
ell4 = scale * L4 * unit_circle;

fprintf('\n[8] (num_seg=%d, dur=%d ms, conf=%.0f%%)\n', ...
        target_seg, target_ms, 100*conf_level);
fprintf('  6 Pairs: centro=(%.2f, %.2f)\n', az6_c, el6_c);
fprintf('  5 Pairs: centro=(%.2f, %.2f)\n', az5_c, el5_c);
fprintf('  4 Pairs: centro=(%.2f, %.2f)\n', az4_c, el4_c);

% ---------------------------------------------------------
% 8.3) Figura
% ---------------------------------------------------------
fig = figure('Position', [100 100 1000 1000], 'Color', 'w');
hold on; grid on;

% Nuvens de pontos
scatter(az6, el6, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.4, ...
        'DisplayName', '6 Pairs');
scatter(az5, el5, 50, 'g', 'filled', 'MarkerFaceAlpha', 0.4, ...
        'DisplayName', '5 Pairs');
scatter(az4, el4, 50, 'r', 'filled', 'MarkerFaceAlpha', 0.4, ...
        'DisplayName', '4 Pairs');

% Elipses de confiança
plot(az6_c + ell6(1,:), el6_c + ell6(2,:), 'b-', 'LineWidth', 3, ...
     'HandleVisibility', 'off');
plot(az5_c + ell5(1,:), el5_c + ell5(2,:), 'g-', 'LineWidth', 3, ...
     'HandleVisibility', 'off');
plot(az4_c + ell4(1,:), el4_c + ell4(2,:), 'r-', 'LineWidth', 3, ...
     'HandleVisibility', 'off');

 % DOA real
plot(az_real, el_real, 'p', 'MarkerSize', 40, ...
     'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 2, ...
     'DisplayName', 'Real DOA');
% Centros
plot(az6_c, el6_c, 'bo', 'MarkerSize', 20,'MarkerFaceColor', 'b','MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 P. Center');
plot(az5_c, el5_c, 'go', 'MarkerSize', 20,'MarkerFaceColor', 'g','MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 P. Center');
plot(az4_c, el4_c, 'ro', 'MarkerSize', 20,'MarkerFaceColor', 'r','MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '4 P. Center');



% ---------------------------------------------------------
% 8.4) Formatação — igual ao plot de referência
% ---------------------------------------------------------
xlabel('Azimuth (degrees)',   'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Elevation (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Limites dinâmicos baseados nos dados
all_az = [az6; az5; az4];
all_el = [el6; el5; el4];

xlim ([ 0 70 ])
ylim([0 50])

h_legend = legend('Location', 'northeast');
set(h_legend, 'FontSize', 30);
set(gca, 'FontSize', 30);
set(gca, 'Position', [0.12 0.15 0.84 0.80]);

hold off;

% ---------------------------------------------------------
% 8.5) Exportação
% ---------------------------------------------------------
fname = sprintf('scatter_ellipse_conf%d_%dseg_%dms', ...
                round(100*conf_level), target_seg, target_ms);
print(fullfile(outputDir, [fname '.eps']), '-depsc');
print(fullfile(outputDir, [fname '.png']), '-dpng', '-r300');

%% 9) PLOT: AZIMUTE VS. ELEVAÇÃO — gradiente por num_seg_list
% =========================================================
target_ms_9 = 80;
idx_dur_9   = find(abs(dur_s_list - target_ms_9*1e-3) < 1e-9, 1);

% ---------------------------------------------------------
% 9.1) Cores base — idênticas ao plot 8
% ---------------------------------------------------------
color_base_6 = [0.00 0.00 1.00];
color_base_5 = [0.00 1.00 0.00];
color_base_4 = [1.00 0.00 0.00];

num_color_steps  = length(num_seg_list);
lightness_factor = 0.05;

colors_gradient_6 = zeros(num_color_steps, 3);
colors_gradient_5 = zeros(num_color_steps, 3);
colors_gradient_4 = zeros(num_color_steps, 3);

for k = 1:num_color_steps
    alpha = (k - 1) / max(num_color_steps - 1, 1);
    colors_gradient_6(k,:) = (1-alpha)*((1-lightness_factor)*[1 1 1] + lightness_factor*color_base_6) + alpha*color_base_6;
    colors_gradient_5(k,:) = (1-alpha)*((1-lightness_factor)*[1 1 1] + lightness_factor*color_base_5) + alpha*color_base_5;
    colors_gradient_4(k,:) = (1-alpha)*((1-lightness_factor)*[1 1 1] + lightness_factor*color_base_4) + alpha*color_base_4;
end

% ---------------------------------------------------------
% 9.2) Figura
% ---------------------------------------------------------
figure('Position', [100 100 1000 2000], 'Color', 'w');
hold on; grid on; box on;
xlim([10 60]); ylim([5 30]);

% Handles de legenda — placeholders com cor base (sempre válidos)
h_real      = plot(NaN, NaN, 'p',  'MarkerSize', 40, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 3,   'DisplayName', 'Real DOA');
% h_est_6     = plot(NaN, NaN, 'o',  'MarkerSize', 8,  'Color', color_base_6);
% h_est_5     = plot(NaN, NaN, 'o',  'MarkerSize', 8,  'Color', color_base_5);
% h_est_4     = plot(NaN, NaN, 'o',  'MarkerSize', 8,  'Color', color_base_4);
% h_med_6     = plot(NaN, NaN, '^',  'MarkerSize', 30, 'MarkerFaceColor', color_base_6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 P. Median');
h_med_5     = plot(NaN, NaN, '^',  'MarkerSize', 30, 'MarkerFaceColor', color_base_5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 P. Median');
h_med_4     = plot(NaN, NaN, '^',  'MarkerSize', 30, 'MarkerFaceColor', color_base_4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '4 P. Median');
% h_mean_6    = plot(NaN, NaN, 'd',  'MarkerSize', 30, 'MarkerFaceColor', color_base_6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 P. Mean');
h_mean_5    = plot(NaN, NaN, 'd',  'MarkerSize', 30, 'MarkerFaceColor', color_base_5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 P. Mean');
h_mean_4    = plot(NaN, NaN, 'd',  'MarkerSize', 30, 'MarkerFaceColor', color_base_4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '4 P. Mean');

% ---------------------------------------------------------
% 9.3) Loop por num_seg_list
% ---------------------------------------------------------
for idx_case = 1:num_color_steps

    current_num_segments = num_seg_list(idx_case);
    c6 = colors_gradient_6(idx_case, :);
    c5 = colors_gradient_5(idx_case, :);
    c4 = colors_gradient_4(idx_case, :);

    % Estimativas individuais — todos os segmentos x todas as realizações
    for idx_rng = 1:num_cases0
        az6_v = az_seg_6_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
        el6_v = el_seg_6_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
        az5_v = az_seg_5_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
        el5_v = el_seg_5_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
        az4_v = az_seg_4_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
        el4_v = el_seg_4_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);

%         plot(az6_v(~isnan(az6_v)), el6_v(~isnan(el6_v)), 'o', 'MarkerSize', 8, 'Color', c6, 'HandleVisibility', 'off');
%         plot(az5_v(~isnan(az5_v)), el5_v(~isnan(el5_v)), 'o', 'MarkerSize', 8, 'Color', c5, 'HandleVisibility', 'off');
%         plot(az4_v(~isnan(az4_v)), el4_v(~isnan(el4_v)), 'o', 'MarkerSize', 8, 'Color', c4, 'HandleVisibility', 'off');
    end

    % Mediana e média agregadas sobre realizações
    az6_med  = median(squeeze(DOA_az_6_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    el6_med  = median(squeeze(DOA_el_6_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    az5_med  = median(squeeze(DOA_az_5_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    el5_med  = median(squeeze(DOA_el_5_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    az4_med  = median(squeeze(DOA_az_4_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    el4_med  = median(squeeze(DOA_el_4_median_all(idx_case,:,idx_dur_9)), 'omitnan');

    az6_mn   = mean(squeeze(DOA_az_6_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    el6_mn   = mean(squeeze(DOA_el_6_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    az5_mn   = mean(squeeze(DOA_az_5_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    el5_mn   = mean(squeeze(DOA_el_5_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    az4_mn   = mean(squeeze(DOA_az_4_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    el4_mn   = mean(squeeze(DOA_el_4_mean_all(idx_case,:,idx_dur_9)), 'omitnan');

%     plot(az6_med, el6_med, '^', 'MarkerSize', 30, 'MarkerFaceColor', c6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az5_med, el5_med, '^', 'MarkerSize', 30, 'MarkerFaceColor', c5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az4_med, el4_med, '^', 'MarkerSize', 30, 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

%     plot(az6_mn,  el6_mn,  'd', 'MarkerSize', 30, 'MarkerFaceColor', c6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az5_mn,  el5_mn,  'd', 'MarkerSize', 30, 'MarkerFaceColor', c5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az4_mn,  el4_mn,  'd', 'MarkerSize', 30, 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
end

% DOA real — por cima de tudo
plot(az_real, el_real, 'p', 'MarkerSize', 40, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'HandleVisibility', 'off');

% ---------------------------------------------------------
% 9.4) Formatação
% ---------------------------------------------------------
xlabel('Azimuth (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Elevation (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');

legend([h_real, h_est_6, h_est_5, h_est_4, h_med_6, h_med_5, h_med_4, h_mean_6, h_mean_5, h_mean_4], ...
       'Location', 'southwest');
set(legend, 'FontSize', 30);
set(gca, 'FontSize', 30);
set(gca, 'Position', [0.14 0.15 0.63 0.78]);

% ---------------------------------------------------------
% 9.5) Legenda de segmentos (patch manual)
% ---------------------------------------------------------
ax_seg = axes('Position', [0.80 0.15 0.10 0.25], 'Visible', 'off');
text(-0.005, 1.1, 'Segments:', 'FontSize', 25, 'FontWeight', 'bold', 'Parent', ax_seg);
for k = 1:num_color_steps
    ypos = 1 - 0.15*k;
    patch([0 0.10 0.10 0], [ypos ypos ypos+0.10 ypos+0.10], ...
          colors_gradient_6(k,:), 'EdgeColor', 'k', 'Parent', ax_seg);
    text(0.15, ypos+0.05, sprintf('%d', num_seg_list(k)), ...
         'FontSize', 20, 'Parent', ax_seg);
end

hold off;

% ---------------------------------------------------------
% 9.6) Exportação
% ---------------------------------------------------------
print(fullfile(outputDir, 'doa_az_el_6_5_4_pairs_segments_color.eps'), '-depsc');
print(fullfile(outputDir, 'doa_az_el_6_5_4_pairs_segments_color.png'), '-dpng', '-r300');

%% =========================================================
%% 9) PLOT: AZIMUTE VS. ELEVAÇÃO — gradiente por num_seg_list
% =========================================================

target_ms_9 = 80;
idx_dur_9   = find(abs(dur_s_list - target_ms_9*1e-3) < 1e-9, 1);

% ---------------------------------------------------------
% 9.1) Cores base — idênticas ao plot 8
% ---------------------------------------------------------
color_base_6 = [0.00 0.00 1.00];
color_base_5 = [0.00 1.00 0.00];
color_base_4 = [1.00 0.00 0.00];

num_color_steps  = length(num_seg_list);
lightness_factor = 0.05;

colors_gradient_6 = zeros(num_color_steps, 3);
colors_gradient_5 = zeros(num_color_steps, 3);
colors_gradient_4 = zeros(num_color_steps, 3);

for k = 1:num_color_steps
    alpha = (k - 1) / max(num_color_steps - 1, 1);
    colors_gradient_6(k,:) = (1-alpha)*((1-lightness_factor)*[1 1 1] + lightness_factor*color_base_6) + alpha*color_base_6;
    colors_gradient_5(k,:) = (1-alpha)*((1-lightness_factor)*[1 1 1] + lightness_factor*color_base_5) + alpha*color_base_5;
    colors_gradient_4(k,:) = (1-alpha)*((1-lightness_factor)*[1 1 1] + lightness_factor*color_base_4) + alpha*color_base_4;
end

% ---------------------------------------------------------
% 9.2) Figura
% ---------------------------------------------------------
figure('Position', [100 100 1000 1000], 'Color', 'w');
hold on; grid on; box on;
xlim([10 60]); ylim([5 30]);

% Placeholders de legenda (criados antes do loop)
h_real   = plot(NaN, NaN, 'p',  'MarkerSize', 40, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 3,   'DisplayName', 'Real DOA');
% h_est_6  = plot(NaN, NaN, 'o',  'MarkerSize', 8,  'Color', color_base_6,  'DisplayName', '6 Pairs');
% h_est_5  = plot(NaN, NaN, 'o',  'MarkerSize', 8,  'Color', color_base_5,  'DisplayName', '5 Pairs');
% h_est_4  = plot(NaN, NaN, 'o',  'MarkerSize', 8,  'Color', color_base_4,  'DisplayName', '4 Pairs');
% h_med_6  = plot(NaN, NaN, '^',  'MarkerSize', 30, 'MarkerFaceColor', color_base_6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 P. Median');
h_med_5  = plot(NaN, NaN, '^',  'MarkerSize', 30, 'MarkerFaceColor', color_base_5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 P. Median');
h_med_4  = plot(NaN, NaN, '^',  'MarkerSize', 30, 'MarkerFaceColor', color_base_4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '4 P. Median');
% h_mean_6 = plot(NaN, NaN, 'd',  'MarkerSize', 30, 'MarkerFaceColor', color_base_6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 P. Mean');
h_mean_5 = plot(NaN, NaN, 'd',  'MarkerSize', 30, 'MarkerFaceColor', color_base_5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 P. Mean');
h_mean_4 = plot(NaN, NaN, 'd',  'MarkerSize', 30, 'MarkerFaceColor', color_base_4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '4 P. Mean');

% ---------------------------------------------------------
% 9.3) PASSO 1 — plota todos os pontos individuais primeiro
% ---------------------------------------------------------
for idx_case = 1:num_color_steps

    current_num_segments = num_seg_list(idx_case);
    % c6 = colors_gradient_6(idx_case, :);
    c5 = colors_gradient_5(idx_case, :);
    c4 = colors_gradient_4(idx_case, :);

    for idx_rng = 1:num_cases0
        % az6_v = az_seg_6_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
        % el6_v = el_seg_6_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
%         az5_v = az_seg_5_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
%         el5_v = el_seg_5_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
%         az4_v = az_seg_4_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
%         el4_v = el_seg_4_cell{idx_case, idx_rng, idx_dur_9}(1:current_num_segments);
% 
%         % plot(az6_v(~isnan(az6_v)), el6_v(~isnan(el6_v)), 'o', 'MarkerSize', 8, 'Color', c6, 'HandleVisibility', 'off');
%         plot(az5_v(~isnan(az5_v)), el5_v(~isnan(el5_v)), 'o', 'MarkerSize', 8, 'Color', c5, 'HandleVisibility', 'off');
%         plot(az4_v(~isnan(az4_v)), el4_v(~isnan(el4_v)), 'o', 'MarkerSize', 8, 'Color', c4, 'HandleVisibility', 'off');
    end
end
% DOA real 
plot(az_real, el_real, 'p', 'MarkerSize', 40, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'HandleVisibility', 'off');

% ---------------------------------------------------------
% 9.4) PASSO 2 — plota medianas e médias por cima dos pontos
% ---------------------------------------------------------
for idx_case = 1:num_color_steps

    % c6 = colors_gradient_6(idx_case, :);
    c5 = colors_gradient_5(idx_case, :);
    c4 = colors_gradient_4(idx_case, :);

    % az6_med  = median(squeeze(DOA_az_6_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    % el6_med  = median(squeeze(DOA_el_6_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    az5_med  = median(squeeze(DOA_az_5_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    el5_med  = median(squeeze(DOA_el_5_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    az4_med  = median(squeeze(DOA_az_4_median_all(idx_case,:,idx_dur_9)), 'omitnan');
    el4_med  = median(squeeze(DOA_el_4_median_all(idx_case,:,idx_dur_9)), 'omitnan');

    % az6_mn   = mean(squeeze(DOA_az_6_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    % el6_mn   = mean(squeeze(DOA_el_6_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    az5_mn   = mean(squeeze(DOA_az_5_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    el5_mn   = mean(squeeze(DOA_el_5_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    az4_mn   = mean(squeeze(DOA_az_4_mean_all(idx_case,:,idx_dur_9)), 'omitnan');
    el4_mn   = mean(squeeze(DOA_el_4_mean_all(idx_case,:,idx_dur_9)), 'omitnan');

    % plot(az6_med, el6_med, '^', 'MarkerSize', 30, 'MarkerFaceColor', c6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az5_med, el5_med, '^', 'MarkerSize', 30, 'MarkerFaceColor', c5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az4_med, el4_med, '^', 'MarkerSize', 30, 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

    % plot(az6_mn, el6_mn, 'd', 'MarkerSize', 30, 'MarkerFaceColor', c6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az5_mn, el5_mn, 'd', 'MarkerSize', 30, 'MarkerFaceColor', c5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    plot(az4_mn, el4_mn, 'd', 'MarkerSize', 30, 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
end

% ---------------------------------------------------------
% 9.5) Formatação
% ---------------------------------------------------------
xlabel('Azimuth (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Elevation (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');

legend([h_real, h_med_5, h_med_4, h_mean_5, h_mean_4], ...
       'Location', 'northeast');
set(legend, 'FontSize', 30);
set(gca, 'FontSize', 30);
set(gca, 'Position', [0.12 0.14 0.85 0.83]);

% ---------------------------------------------------------
% 9.6) Legenda de segmentos (patch manual)
% ---------------------------------------------------------
ax_seg = axes('Position', [0.71 0.15 0.24 0.35], 'Visible', 'off');

% Fundo branco com borda preta
patch([0 1 1 0 0], [0 0 1.3 1.3 0], 'w', ...
      'EdgeColor', 'k', 'LineWidth', 1.5, ...
      'Parent', ax_seg);
  
text(0, 1.15, 'Segments:', 'FontSize', 25, 'FontWeight', 'bold', 'Parent', ax_seg);

bar_h = 0.09;   % altura da barra = mesma altura do texto
for k = 1:num_color_steps
    ypos = 1 - 0.15*k;
    patch([0.02 0.62 0.62 0.02], [ypos ypos ypos+bar_h ypos+bar_h], ...
          colors_gradient_5(k,:), 'EdgeColor', 'k', 'Parent', ax_seg);
    text(0.64, ypos + bar_h/2, sprintf('%d', num_seg_list(k)), ...
         'FontSize', 25, 'Parent', ax_seg, 'VerticalAlignment', 'middle');
end

hold off;

% ---------------------------------------------------------
% 9.7) Exportação
% ---------------------------------------------------------
print(fullfile(outputDir, 'doa_az_el_5_4_pairs_segments_color.eps'), '-depsc');
print(fullfile(outputDir, 'doa_az_el_5_4_pairs_segments_color.png'), '-dpng', '-r300');
%% Functions
function u = azel2unitvec(az_deg, el_deg)
    az_rad = deg2rad(az_deg);
    el_rad = deg2rad(el_deg);
    u = [cos(el_rad).*cos(az_rad);
         cos(el_rad).*sin(az_rad);
         sin(el_rad)];
end

function err_deg = ang_err_deg(u_ref, u_est)
    % u_ref e u_est são vetores unitários 3x1 ou 3xN
    dot_product = dot(u_ref, u_est);
    % Garantir que o produto escalar esteja dentro do domínio de acos [-1, 1]
    dot_product = max(-1, min(1, dot_product));
    err_rad = acos(dot_product);
    err_deg = rad2deg(err_rad);
end

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
    
function [r, lags_s] = gcc_family_local(x, xref, fs, maxLagSamples, mode)
    x    = x(:)    - mean(x);
    xref = xref(:) - mean(xref);

    N     = length(x);
    Ncorr = 2*N - 1;
    Nfft  = 2^nextpow2(Ncorr); % Sua abordagem, que é a mais robusta

    X  = fft(x,    Nfft);
    XR = fft(xref, Nfft);

    % Calcula o espectro de potência cruzada e seu módulo
    Sxy = X .* conj(XR);
    mag_Sxy = abs(Sxy); % Isso é o abs(X1.*X2) do código de referência

    switch lower(mode)
        case 'gcc'
            % ---- GCC convencional ----
            cc = real(ifft(Sxy, Nfft));       % time domain
            cc = fftshift(cc);               % centraliza
            % cc = cc / max(abs(cc) + eps);    % Normalização opcional para GCC

        case 'phat'
            % ---- GCC-PHAT (espelhando o código de referência) ----

            % Denominador PHAT com suavização (std(abs(X1.*X2))/1000+abs(X1.*X2))
            % Note que 'std' de um vetor complexo não é o que queremos aqui.
            % O código de referência usa std(abs(X1.*X2)), que é o std do módulo.
            % Vamos usar o std de 'mag_Sxy' para replicar isso.

            % Para evitar problemas com std de um vetor de zeros ou constante,
            % podemos adicionar um pequeno valor fixo se std for muito pequeno.
            std_mag_Sxy_term = std(mag_Sxy);
            if std_mag_Sxy_term < eps % Se o desvio padrão for quase zero
                std_mag_Sxy_term = 1e-6; % Um pequeno valor fixo
            end

            denom = (std_mag_Sxy_term / 1000) + mag_Sxy;

            % Aplica o filtro PHAT e IFFT
            cc = real(ifft(Sxy ./ denom, Nfft)); % Sxy / denom
            cc = fftshift(cc);                   % centraliza

            % Normalização final para 1, como no código de referência
            cc = cc / max(abs(cc) + eps);        

        otherwise
            error('Modo de GCC não reconhecido. Use ''gcc'' ou ''phat''.');
    end

    % Recorta para o intervalo de lags fisicamente possível
    lags = (-Nfft/2 : Nfft/2-1).';
    idx  = (lags >= -maxLagSamples) & (lags <= maxLagSamples);

    r      = cc(idx);
    lags_s = lags(idx) / fs;
end

function [az_hat, el_hat, uDoA_vector] = calculateDoA_LS(tdoa_values, pairs_indices, microphone_positions, sound_speed)
% calculateDoA_LS Calcula a Direção de Chegada (DoA) usando o método de Mínimos Quadrados (LS).
%
%   [az_hat, el_hat, uDoA_vector] = calculateDoA_LS(tdoa_values, pairs_indices, microphone_positions, sound_speed)
%
%   Esta função implementa o cálculo da DoA (azimute e elevação) usando o
%   método de Mínimos Quadrados (Least Squares - LS), conforme descrito
%   no artigo "Sobre a escolha de sinais em arranjos de microfones
%   estimando DoA com GCC-PhaT".
%
%   Entradas:
%     tdoa_values: Um vetor coluna de TDoAs (Time Difference of Arrival)
%                  em segundos para os pares de microfones selecionados.
%                  Cada elemento corresponde a um par em pairs_indices.
%     pairs_indices: Uma matriz Np x 2, onde Np é o número de pares de
%                    microfones. Cada linha [i, j] representa os índices
%                    dos microfones que formam um par. Os índices devem
%                    corresponder às linhas em microphone_positions.
%     microphone_positions: Uma matriz M x 3, onde M é o número total de
%                           microfones. Cada linha [x, y, z] representa
%                           as coordenadas espaciais de um microfone.
%     sound_speed: A velocidade do som no meio (em metros por segundo).
%
%   Saídas:
%     az_hat: O azimute estimado da DoA em graus.
%     el_hat: A elevação estimada da DoA em graus.
%     uDoA_vector: O vetor unitário da DoA [ax, ay, az].
%
%   Exemplo de Uso:
%     % Supondo que você tenha 7 microfones e 21 pares (C(7,2))
%     % Exemplo de posições de microfones (substitua pelas suas reais)
%     mic_pos = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];
%     
%     % Exemplo de TDoAs para 6 pares (substitua pelos seus valores reais)
%     % Estes TDoAs devem ser calculados previamente (e.g., via GCC-PhaT)
%     tdoas_exemplo = [0.001; 0.0005; -0.0002; 0.0008; -0.0001; 0.0003]; 
%
%     % Exemplo de pares de microfones (substitua pelos seus reais)
%     % Para 6 pares, você precisaria selecionar quais são.
%     % Aqui, um exemplo arbitrário:
%     pares_exemplo = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; 
%
%     c_som = 343; % Velocidade do som
%
%     [azimute, elevacao, vetor_doa] = calculateDoA_LS(tdoas_exemplo, pares_exemplo, mic_pos, c_som);
%     fprintf('Azimute estimado: %.2f graus\n', azimute);
%     fprintf('Elevação estimada: %.2f graus\n', elevacao);
%     fprintf('Vetor DoA: [%.4f, %.4f, %.4f]\n', vetor_doa(1), vetor_doa(2), vetor_doa(3));

    Np = size(pairs_indices, 1); % Número de pares de microfones

    % ----------------------------------------------------------
    % Matriz de geometria (Deltapb)
    % Conforme a equação ?_ij = a_T_?,? * ?p_ij, onde ?p_ij = (p_i - p_j) / v_som
    % ----------------------------------------------------------
    Deltapb = zeros(Np, 3);
    for p = 1:Np
        ii = pairs_indices(p, 1); % Índice do primeiro microfone no par
        jj = pairs_indices(p, 2); % Índice do segundo microfone no par

        % Ajusta para índices base 1 se microphone_positions for indexado de 1 a M
        pos_i = microphone_positions(ii, :);
        pos_j = microphone_positions(jj, :);

        Deltapb(p, :) = (pos_i - pos_j) / sound_speed;
    end

    % ----------------------------------------------------------
    % Cálculo da DoA usando Mínimos Quadrados (LS)
    % A_full = Deltapb' * Deltapb
    % b_phat_full = Deltapb' * tdoa_phat_mean
    % uDoA = A_full \ b_phat_full
    % ----------------------------------------------------------
    A_matrix = Deltapb' * Deltapb;
    b_vector = Deltapb' * tdoa_values; % tdoa_values já é o vetor de TDoAs para os pares atuais

    uDoA_vector = A_matrix \ b_vector;

    % Normalizar o vetor DoA para que seja unitário
    uDoA_vector = uDoA_vector / norm(uDoA_vector);

    % Converter uDoA_vector para azimute e elevação
    % O artigo usa a?,? = [?sin ? cos ? ?sin ? sin ? ?cos ?]T
    % e depois ? = cos?1 az e ? = tan?1 ay/ax.
    % Isso implica que az é a componente Z do vetor unitário.
    % E que a elevação é o ângulo com o plano XY (90 - ângulo rasante).
    % O artigo define ? como "ângulo rasante ou 90? - ângulo de elevação".
    % Se ? é o ângulo rasante, então cos(?) = az.
    % Se ? é 90 - elevação, então elevação = 90 - ?.
    % Vamos usar a convenção mais comum onde elevação é o ângulo com o plano XY.
    %
    % uDoA_vector = [ax; ay; az]
    % azimute (?) = atan2(ay, ax)
    % elevação (?_e) = asin(az)

    az_hat = atan2(uDoA_vector(2), uDoA_vector(1)) * 180/pi;
    el_hat = asin(uDoA_vector(3)) * 180/pi; % Ângulo de elevação em relação ao plano XY

    % O artigo usa ? como "ângulo rasante ou 90? - ângulo de elevação".
    % Se você precisar do "ângulo rasante" (theta_rasante), seria:
    % theta_rasante = acos(uDoA_vector(3)) * 180/pi;
    % E o azimute (phi) seria:
    % phi_rasante = atan2(-uDoA_vector(2), -uDoA_vector(1)) * 180/pi;
    % No entanto, para consistência com a maioria das aplicações e a sua solicitação,
    % estou retornando azimute e elevação padrão (elevação em relação ao plano XY).
end

function plotcorr(seg_idx_to_plot)
% plot_all_correlations_for_segment Plota as funções de correlação cruzada PHAT para TODOS os pares
%                                   de um segmento específico.
%   Assume que as variáveis necessárias (lagsP, rP_all, tdoa_phat, tdoa_th_plane,
%   fs, pares, pos_h, c) já estão disponíveis no workspace base.
%
%   plot_all_correlations_for_segment(seg_idx_to_plot)
%
%   Argumentos:
%     seg_idx_to_plot: Índice do segmento a ser plotado (ex: 1 para o primeiro segmento).

    % Acessa variáveis do workspace base.
    try
        lagsP = evalin('base', 'lagsP');
        rP_all = evalin('base', 'rP_all'); % rP_all deve ser 3D: (num_lags x num_pares x num_segments)
        tdoa_phat_all = evalin('base', 'tdoa_phat'); % tdoa_phat deve ser 2D: (num_segments x num_pares)
        tdoa_th_plane = evalin('base', 'tdoa_th_plane');
        fs = evalin('base', 'fs');
        pares = evalin('base', 'pares');
        % outputDir = evalin('base', 'outputDir'); % Descomente se for usar para salvar
    catch ME
        error('Erro: Variável "%s" não encontrada no workspace base. Certifique-se de que seu script principal foi executado e gerou todas as variáveis necessárias.', ME.identifier(find(ME.identifier == ':', 1, 'last')+1:end));
    end

    num_segments = size(tdoa_phat_all, 1);
    num_pares = size(tdoa_phat_all, 2); % Np

    % Validação básica do índice do segmento
    if seg_idx_to_plot < 1 || seg_idx_to_plot > num_segments
        error('Índice de segmento inválido. Deve estar entre 1 e %d.', num_segments);
    end

    fprintf('Plotando todas as %d correlações para o segmento %d...\n', num_pares, seg_idx_to_plot);

    % Loop para plotar cada par do segmento especificado
    for p_to_plot = 1:num_pares
        % Extrair dados para o segmento e par específicos
        rP = rP_all(:, p_to_plot, seg_idx_to_plot); % Correlação para o par e segmento
        tdoa_phat_val = tdoa_phat_all(seg_idx_to_plot, p_to_plot); % TDOA estimado para o par e segmento
        tdoa_th_val = tdoa_th_plane(p_to_plot); % TDOA teórico para o par

        % Encontrar o pico da correlação para o marcador
        [~, kP] = max(rP);

        % Obter os índices dos microfones do par
        i = pares(p_to_plot, 1);
        j = pares(p_to_plot, 2);

        % Plotagem
        fig = figure('Position', [100 100 1200 600], 'Color', 'w');
        plot(1e3*lagsP, rP, 'k', 'LineWidth', 1.05); grid on; hold on;

        yl = ylim;
        h_hat = plot(1e3*[tdoa_phat_val tdoa_phat_val], yl, 'm-',  'LineWidth', 1.2);
        h_the = plot(1e3*[tdoa_th_val tdoa_th_val],         yl, 'b--', 'LineWidth', 1.1);
        h_pk  = plot(1e3*tdoa_phat_val, rP(kP), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
        hold off;

        % Ajustar limites do eixo X
        maxLagSamples = max(abs(lagsP * fs)); % Calcula o máximo de amostras de atraso
        xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
        ylim([-1 1]);

        xlabel('Delay [ms]', 'FontSize', 20, 'FontWeight', 'bold');
        ylabel(sprintf('Corr. (%d,%d)', i, j), 'FontSize', 20, 'FontWeight', 'bold');

        legend([h_hat h_the h_pk], { ...
            sprintf('Est.   \\tau_{%d%d} = %.3f ms', i, j, 1e3*tdoa_phat_val), ...
            sprintf('Theor. \\tau_{%d%d} = %.3f ms', i, j, 1e3*tdoa_th_val), ...
            'Peak'}, ...
            'Location', 'best', 'Interpreter', 'tex', 'FontSize', 20);

        set(gca, 'FontSize', 20);
        set(gca, 'Position', [0.08 0.18 0.90 0.75])

        Se quiser salvar, descomente a linha abaixo e defina 'outputDir' no seu script principal
        print(fullfile(outputDir, sprintf('phat_pair_%d%d_seg%02d.eps', i, j, seg_idx_to_plot)), '-depsc')
    end
    fprintf('Plotagem concluída para o segmento %d.\n', seg_idx_to_plot);
end
