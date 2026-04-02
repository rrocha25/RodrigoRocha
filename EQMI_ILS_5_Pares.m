% clear; clc; close all;
clear; close all

%% 1) Parameters
fs = 44100;
c  = 1500;

dur_s = 80e-3;
num_segments = 90; % Número de segmentos para o loop
total_dur_s = num_segments * dur_s; % Duração total simulada

% t = (0:round(dur_s*fs)-1)' / fs;

SNR_dB = -14;

% Geometry
fmax_geom  = 1000;
lambda_min = c / fmax_geom;
L          = 0.5 * lambda_min;

pos_h = pyramid_array_positions(L);

% Simulated source
az_real = 0;
el_real = 45;
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


%% Plot 3D da geometria do arranjo piramidal
%{
figure('Position', [100 100 900 700], 'Color', 'w');
hold on; grid on; axis equal;

% Define as margens e espaçamentos para os subplots
% Estes valores podem ser ajustados para o seu layout preferido
left_margin = 0; % Margem esquerda da figura
right_margin = 0; % Margem direita da figura
bottom_margin = 0.08; % Margem inferior da figura (para a colorbar e labels)
top_margin = 0.01; % Margem superior da figura (para o título)
h_spacing = 0.07; % Espaçamento horizontal entre subplots
v_spacing = 0.08; % Espaçamento vertical entre subplots


% Calcula a posição para o único eixo
ax_width_3d = 1 - left_margin - right_margin;
ax_height_3d = 1 - top_margin - bottom_margin;
set(gca, 'Position', [left_margin, bottom_margin, ax_width_3d, ax_height_3d]);


% Calcular altura do tetraedro (apenas para uso interno no plot)
H_calc = L * sqrt(2/3);

% Plotar os 4 hidrofones
scatter3(pos_h(:,1), pos_h(:,2), pos_h(:,3), 150, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Rotular cada hidrofone
for k = 1:4
    text(pos_h(k,1), pos_h(k,2), pos_h(k,3) + 0.05, sprintf('  H%d', k), ...
        'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');
end

% Desenhar as arestas da base (triângulo equilátero)
base_idx = [1 2; 2 3; 3 1];
for i = 1:size(base_idx, 1)
    p1 = pos_h(base_idx(i,1), :);
    p2 = pos_h(base_idx(i,2), :);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'k-', 'LineWidth', 1.5);
end

% Desenhar as arestas da pirâmide (do ápice para a base)
apex_idx = [1 4; 2 4; 3 4];
for i = 1:size(apex_idx, 1)
    p1 = pos_h(apex_idx(i,1), :);
    p2 = pos_h(apex_idx(i,2), :);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], 'b--', 'LineWidth', 1.2);
end

% Adicionar anotações de distâncias importantes
% Lado da base
mid_12 = (pos_h(1,:) + pos_h(2,:))/2;
% text(mid_12(1), mid_12(2), mid_12(3) - 0.05, sprintf('$L = %.3f$ m', L), ...
%     'FontSize', 18, 'Color', 'k', 'Interpreter', 'latex');

% Comprimento de onda mínimo
% text(0.1, 0.1, 0.05, sprintf('$\\lambda_{\\mathrm{min}}/2 = %.3f$ m', lambda_min/2), ...
%     'FontSize', 18, 'Color', 'b', 'Interpreter', 'latex');

% Distância máxima entre pares
text(0.1, 0.1, H_calc + 0.15, sprintf('$d_{\\mathrm{max}} = %.3f$ m', L), ...
    'FontSize', 20, 'Color', [0 0.5 0], 'Interpreter', 'latex', 'FontWeight', 'bold');

% Configurações dos eixos
xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', 20);

% Ajustar visualização
view(45, 30);  % Ângulo de visão 3D
xlim([-0.01 L+0.1]);
ylim([-0.01 L*sqrt(3)/2 + 0.1]);
zlim([0 H_calc + 0.2]);
set(gca, 'FontSize', 20);
box on;

% Adicionar legenda
legend({'Hydrophones', 'Base edges', 'Tetrahedron edges'}, ...
    'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);

% Adicionar seta no topo do eixo Y
y_max = L*sqrt(3)/2 + 0.1;  % mesmo valor do ylim superior

% Seta 3D apontando na direção Y a partir do limite superior
quiver3(0, y_max*2/3, 0, ...   % ponto de origem
        0, y_max/3, 0, ...            % direção e magnitude (dx, dy, dz)
        0, ...                     % AutoScale off (fator = 0)
        'k', ...
        'LineWidth', 2, ...
        'MaxHeadSize', 0.8);
hold off;

% Output directory for figures
outputDir = './figures';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
print(fullfile(outputDir, 'array_MATLAB.png'), '-depsc')
%}

%% 2) Physical maxLag
dmax          = max_pair_distance(pos_h);
maxLagSamples = ceil((dmax/c)*fs) + 2;


%% 3) Sinal real — leitura do arquivo WAV

% ----------------------------------------------------------
% Lê o sinal de interesse
% ----------------------------------------------------------
[s_raw, fs_wav] = audioread('H3_lancha2.wav');

% Converte para mono se estéreo
if size(s_raw, 2) > 1
    s_raw = mean(s_raw, 2);
end

% Reamostra para fs do sistema se necessário
if fs_wav ~= fs
    s_raw = resample(s_raw, fs, fs_wav);
end

s = s_raw(:);   % garante coluna

% ----------------------------------------------------------
% Lê o ruído de ambiente gravado
% ----------------------------------------------------------
[noise_env_raw, fs_noise] = audioread('Enseada dos Anjos.wav');

% Converte para mono se estéreo
if size(noise_env_raw, 2) > 1
    noise_env_raw = mean(noise_env_raw, 2);
end

% Reamostra para fs do sistema se necessário
if fs_noise ~= fs
    noise_env_raw = resample(noise_env_raw, fs, fs_noise);
end

noise_env_raw = noise_env_raw(:);

%% 3.x) Espectro e Espectrograma do Sinal da Lancha e do Ruído (formato 2:1)
%{
% % ----------------------------------------------------------
% % Parâmetros comuns
% % ----------------------------------------------------------
% p_ref = 1e-6;               % referência 1 µPa
% fs_plot = fs;               % mesma Fs do sistema
% 
% % Parâmetros para espectrograma
% window_length = round(0.025 * fs_plot);      % 25 ms
% overlap       = round(0.75 * window_length); % 75% overlap
% nfft_spec     = 2^nextpow2(window_length*4); % FFT para espectrograma
% 
% % ----------------------------------------------------------
% % 3.x.1 – Espectro da lancha (sinal de interesse)
% % ----------------------------------------------------------
% N_s     = length(s);
% Nfft_s  = 2^nextpow2(N_s);
% S_fft   = fft(s, Nfft_s);
% S_fft   = S_fft(1:floor(Nfft_s/2)+1);
% f_s_pos = (0:floor(Nfft_s/2))' * fs_plot / Nfft_s;
% 
% amp_s   = abs(S_fft) / Nfft_s;
% amp_s   = amp_s * sqrt(2);          % pico?RMS (exceto DC)
% amp_s(1)= abs(S_fft(1)) / Nfft_s;   % DC sem fator sqrt(2)
% 
% SPL_s   = 20*log10(amp_s/p_ref + eps);
% 
% idx_zoom_s = f_s_pos <= 1000;       % zoom até 1 kHz
% 
% fig = figure('Position', [100 100 1500 750], 'Color', 'w');
% plot(f_s_pos(idx_zoom_s), SPL_s(idx_zoom_s), ...
%      'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
% grid on
% xlabel('Frequency [Hz]', 'FontSize', 40, 'FontWeight', 'bold');
% ylabel('SPL [dB re 1 \muPa]', 'FontSize', 40, 'FontWeight', 'bold');
% xlim([0 1000]);
% ylim([0 100]);  % ajuste se necessário
% set(gca, 'FontSize', 40);
% set(gca, 'Position', [0.15 0.25 0.80 0.72]);
% 
% print(fullfile(outputDir, 'boat_spectrum.eps'), '-depsc');
% % close(fig);
% 
% % ----------------------------------------------------------
% % 3.x.2 – Espectro do ruído ambiente
% % ----------------------------------------------------------
% N_n      = length(noise_env_raw);
% Nfft_n   = 2^nextpow2(N_n);
% N_fft    = fft(noise_env_raw, Nfft_n);
% N_fft    = N_fft(1:floor(Nfft_n/2)+1);
% f_n_pos  = (0:floor(Nfft_n/2))' * fs_plot / Nfft_n;
% 
% amp_n    = abs(N_fft) / Nfft_n;
% amp_n    = amp_n * sqrt(2);
% amp_n(1) = abs(N_fft(1)) / Nfft_n;
% 
% SPL_n    = 20*log10(amp_n/p_ref + eps);
% 
% idx_zoom_n = f_n_pos <= 1000;
% 
% fig = figure('Position', [100 100 1500 750], 'Color', 'w');
% plot(f_n_pos(idx_zoom_n), SPL_n(idx_zoom_n), ...
%      'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
% grid on
% xlabel('Frequency [Hz]', 'FontSize', 40, 'FontWeight', 'bold');
% ylabel('SPL [dB re 1 \muPa]', 'FontSize', 40, 'FontWeight', 'bold');
% xlim([0 1000]);
% ylim([0 100]);  % ajuste se necessário
% set(gca, 'FontSize', 40);
% set(gca, 'Position', [0.15 0.25 0.80 0.72]);
% 
% print(fullfile(outputDir, 'noise_spectrum.eps'), '-depsc');
% % close(fig);
% 
% % ----------------------------------------------------------
% % 3.x.3 – Espectrograma da lancha
% % ----------------------------------------------------------
% fig = figure('Position', [100 100 1500 750], 'Color', 'w');
% 
% spectrogram(s, hamming(window_length), overlap, nfft_spec, fs_plot, 'yaxis');
% ylim([0 1]);  % até 1 kHz
% 
% clim_max = max(get(gca,'CLim'));
% caxis([clim_max - 60, clim_max]);
% colormap(hot);
% 
% cb = colorbar;
% ylabel(cb, 'Power/Frequency [dB/Hz]', 'FontSize', 40, 'FontWeight', 'bold');
% 
% xlabel('Time [s]',        'FontSize', 40, 'FontWeight', 'bold');
% ylabel('Frequency [kHz]', 'FontSize', 40, 'FontWeight', 'bold');
% set(gca, 'FontSize', 40);
% set(gca, 'Position', [0.11 0.21 0.69 0.73]);
% 
% print(fullfile(outputDir, 'boat_spectrogram.eps'), '-depsc');
% % close(fig);
% 
% % ----------------------------------------------------------
% % 3.x.4 – Espectrograma do ruído ambiente
% % ----------------------------------------------------------
% fig = figure('Position', [100 100 1500 750], 'Color', 'w');
% 
% spectrogram(noise_env_raw, hamming(window_length), overlap, nfft_spec, fs_plot, 'yaxis');
% ylim([0 1]);  % até 1 kHz
% 
% clim_max = max(get(gca,'CLim'));
% caxis([clim_max - 60, clim_max]);
% colormap(hot);
% 
% cb = colorbar;
% ylabel(cb, 'Power/Frequency [dB/Hz]', 'FontSize', 40, 'FontWeight', 'bold');
% 
% xlabel('Time [s]',        'FontSize', 40, 'FontWeight', 'bold');
% ylabel('Frequency [kHz]', 'FontSize', 40, 'FontWeight', 'bold');
% set(gca, 'FontSize', 40);
% set(gca, 'Position', [0.11 0.21 0.69 0.73]);
% 
% print(fullfile(outputDir, 'noise_spectrogram.eps'), '-depsc');
% 
% fprintf('Boat and noise spectrum/spectrogram figures saved to: %s\n', outputDir);

%}

%% Varredura diferentes ruidos
num_seg_list0 = 1:1:300;    % por exemplo, de 25 até 200
num_cases0    = length(num_seg_list0);

for idx_rng = 1:num_cases0
    rng(2*idx_rng);

    %% Varredura no número de segmentos
    num_seg_list = 1:1:1;    % por exemplo, de 25 até 200
    num_cases    = length(num_seg_list);

    err_med_6   = zeros(num_cases,1); % erro angular médio (ou mediano) 6 pares
    err_med_5   = zeros(num_cases,1); % erro angular médio (ou mediano) 5 pares

    for idx_case = 1:num_cases
        num_segments = num_seg_list(idx_case);
        %% 4) Simulate reception at 4 hydrophones with relative delays

        N_seg = round(dur_s * fs);           % amostras por segmento
        N_total = N_seg * num_segments;      % total usado para TDOA/LOO

        x = zeros(N_total, 4);

        % Distâncias / atrasos geométricos (fixos)
        dist = zeros(4,1);
        for k = 1:4
            dist(k) = norm(pos_f - pos_h(k,:));
        end
        tau_rel = (dist - dist(1))/c;

        maxLagSamples = ceil((dmax/c)*fs) + 2;
        guard = maxLagSamples + 4;

        for seg_idx = 1:num_segments

            % Índices deste segmento na linha do tempo global
            idx_start = (seg_idx-1)*N_seg + 1;
            idx_end   = idx_start + N_seg - 1;

            % 1) Segmento de sinal da lancha (80 ms)
            s_seg_raw_current = s(idx_start:idx_end);
            s_seg_raw_current = s_seg_raw_current(:);

            % --- NORMALIZAÇÃO DO SINAL POR TRECHO ---
            max_s_seg = max(abs(s_seg_raw_current));
            s_seg_normalized = s_seg_raw_current / max_s_seg;

            % Calcula a potência do sinal normalizado UMA VEZ por segmento
            ps_normalized = mean(s_seg_normalized.^2);

            % 4) Aplicar atrasos ao SINAL PURO E NORMALIZADO para cada hidrofone
            s_seg_pad = [zeros(guard,1); s_seg_normalized; zeros(guard,1)];

            for k = 1:4
                % Aplica o atraso ao sinal puro e normalizado para o hidrofone 'k'
                sk_delayed  = apply_delay_seconds(s_seg_pad, fs, tau_rel(k));
                seg_delayed_signal = sk_delayed(guard+1:guard+N_seg);

                % --- 2) Segmento de ruído aleatório (80 ms, posição aleatória no ruído gravado) ---
                % A seleção do trecho de ruído agora é feita INDIVIDUALMENTE para cada canal 'k'
                i0 = randi(length(noise_env_raw) - N_seg + 1);
                n_seg_raw_current = noise_env_raw(i0:i0+N_seg-1);
                n_seg_raw_current = n_seg_raw_current(:);

                % --- NORMALIZAÇÃO DO RUÍDO POR TRECHO E POR CANAL ---
                max_n_seg = max(abs(n_seg_raw_current));
                n_seg_normalized = n_seg_raw_current / max_n_seg;

                % --- 3) Ajustar SNR para ESTE segmento e ESTE CANAL ---
                % A potência do ruído é calculada a partir do ruído JÁ NORMALIZADO para este canal
                pn_normalized = mean(n_seg_normalized.^2);

                % Calcula o fator de escala para o ruído para atingir o SNR desejado
                % (ps_normalized é o mesmo para todos os canais neste segmento)
                scale_noise = sqrt(ps_normalized / (pn_normalized * 10^(SNR_dB/10)));

                % Aplica o fator de escala ao ruído normalizado deste canal
                n_seg_scaled = scale_noise * n_seg_normalized;

                % Adiciona o ruído escalado específico deste canal ao sinal atrasado
                x(idx_start:idx_end, k) = seg_delayed_signal + n_seg_scaled;
            end
        end

        %% 5) Estimate TDOAs in segments and (optionally) plot correlations

        pares = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
        Np    = size(pares,1);

        u_real = [cosd(el_real)*cosd(az_real), ...
                  cosd(el_real)*sind(az_real), ...
                  sind(el_real)];

        % TDOA teórico (plano de onda) para cada par (é o mesmo para todos segmentos)
        tdoa_th_plane = zeros(Np,1);
        for p = 1:Np
            i = pares(p,1);
            j = pares(p,2);
            dp = pos_h(j,:) - pos_h(i,:);
            tdoa_th_plane(p) = -dot(dp, u_real)/c;
        end

        % Parâmetros de segmentação
        N_seg = round(dur_s * fs);        % amostras por segmento (80 ms)
        total_samples = size(x,1);        % total de amostras em cada canal

        % Matriz para guardar TDOA estimado por segmento e por par
        % tdoa_phat(seg_idx, p)
        tdoa_phat = nan(num_segments, Np);

        for seg_idx = 1:num_segments

            % índices de amostras deste segmento
            idx_start = (seg_idx-1)*N_seg + 1;
            idx_end   = idx_start + N_seg - 1;

            x_seg = x(idx_start:idx_end, :);   % [N_seg x 4]

            for p = 1:Np
                i = pares(p,1);
                j = pares(p,2);

                % GCC-PHAT neste segmento
                [rP, lagsP] = gcc_family_local(x_seg(:,j), x_seg(:,i), fs, maxLagSamples, 'phat');
                [~, kP]     = max(abs(rP));
                tdoa_phat(seg_idx, p) = lagsP(kP);  
%{
               %--- plotar a correlação de um segmento específico:
%                 if seg_idx == 1  % por exemplo, só o primeiro segmento
%                     fig = figure('Position', [100 100 1200 600], 'Color', 'w');
%                     plot(1e3*lagsP, rP, 'k', 'LineWidth', 1.05); grid on; hold on;
%     
%                     yl = ylim;
%                     h_hat = plot(1e3*[tdoa_phat(seg_idx,p) tdoa_phat(seg_idx,p)], yl, 'm-',  'LineWidth', 1.2);
%                     h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)],         yl, 'b--', 'LineWidth', 1.1);
%                     h_pk  = plot(1e3*tdoa_phat(seg_idx,p), rP(kP), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
%                     hold off;
%     
%                     xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
%                     ylim([-1 1]);
%                     xlabel('Delay [ms]', 'FontSize', 20, 'FontWeight', 'bold');
%                     ylabel(sprintf('Corr. (%d,%d)', i, j), 'FontSize', 20, 'FontWeight', 'bold');
%     
%                     legend([h_hat h_the h_pk], { ...
%                         sprintf('Est.   \\tau_{%d%d} = %.3f ms', i, j, 1e3*tdoa_phat(seg_idx,p)), ...
%                         sprintf('Theor. \\tau_{%d%d} = %.3f ms', i, j, 1e3*tdoa_th_plane(p)), ...
%                         'Peak'}, ...
%                         'Location', 'best', 'Interpreter', 'tex', 'FontSize', 20);
%     
%                     set(gca, 'FontSize', 20);
%                     set(gca, 'Position', [0.08 0.18 0.90 0.75])
%                     % print(fullfile(outputDir, sprintf('phat_pair_%d%d_seg%02d.eps', i, j, seg_idx)), '-depsc')
%                 end
%}     
            end

        end

        %% 6) Least Squares DOA Estimation — PHAT

        % ----------------------------------------------------------
        % Matriz de geometria
        % ----------------------------------------------------------
        Deltapb = zeros(Np, 3);
        for p = 1:Np
            ii = pares(p,1);
            jj = pares(p,2);
            Deltapb(p,:) = (pos_h(ii,:) - pos_h(jj,:)) / c;
        end

        % ----------------------------------------------------------
        %% 6.1) LS com TODAS as 6 correlacoes
        % ----------------------------------------------------------
        % Pré-alocar arrays para armazenar as DOAs por segmento para cada método
        az_per_segment_6_pairs = zeros(num_segments, 1);el_per_segment_6_pairs = zeros(num_segments, 1);
        az_per_segment_5_pairs = zeros(num_segments, 1);el_per_segment_5_pairs = zeros(num_segments, 1);
        az_per_segment_4_pairs = zeros(num_segments, 1);el_per_segment_4_pairs = zeros(num_segments, 1);

        for seg_idx = 1:num_segments
            current_tdoas = tdoa_phat(seg_idx, :)'; % TDOAs para este segmento [Np x 1]

            % ==========================================================
            % 6.1) LS com TODAS as 6 correlações (por segmento)
            % ==========================================================
            [az_seg_6, el_seg_6, ~] = calculateDoA_LS(current_tdoas, pares, pos_h, c);
            az_per_segment_6_pairs(seg_idx) = az_seg_6;
            el_per_segment_6_pairs(seg_idx) = el_seg_6;

            % ==========================================================
            % 6.2) Detecção de outlier via subconjuntos de 5 pares (EQMI adaptado) (por segmento)
            %      Para cada segmento, identificamos o outlier e usamos os 5 pares restantes.
            % ==========================================================
            Np_total = size(pares, 1);
            err_sq_out_seg  = zeros(Np_total, 1);
            az_sub_seg      = zeros(Np_total, 1);
            el_sub_seg      = zeros(Np_total, 1);
            uDoA_sub_seg    = zeros(Np_total, 3);

            for idx_out = 1:Np_total
                idx_keep = setdiff(1:Np_total, idx_out); % Índices dos pares que FICAM (5 pares)

                pairs_sub_seg      = pares(idx_keep, :);
                tdoa_sub_seg       = current_tdoas(idx_keep);
                tdoa_left_out_seg  = current_tdoas(idx_out);

                [az_hat_sub_seg, el_hat_sub_seg, uDoA_hat_sub_seg] = calculateDoA_LS(tdoa_sub_seg, pairs_sub_seg, pos_h, c);

                az_sub_seg(idx_out)   = az_hat_sub_seg;
                el_sub_seg(idx_out)   = el_hat_sub_seg;
                uDoA_sub_seg(idx_out, :) = uDoA_hat_sub_seg(:).';

                ii = pares(idx_out, 1);
                jj = pares(idx_out, 2);
                delta_p_ij_over_c = (pos_h(ii, :) - pos_h(jj, :)) / c;
                tau_pred_seg = dot(delta_p_ij_over_c, uDoA_hat_sub_seg);
                err_sq_out_seg(idx_out) = (tdoa_left_out_seg - tau_pred_seg)^2;
            end
            all_err_sq_out_per_pair(:, seg_idx) = err_sq_out_seg; % Erro para ESTE segmento

            % Encontra o par que, ao ser removido, minimiza o erro quadrático dos pares restantes
            [~, max_eqmi_idx] = max(err_sq_out_seg);

            az_per_segment_5_pairs(seg_idx) = az_sub_seg(max_eqmi_idx);
            el_per_segment_5_pairs(seg_idx) = el_sub_seg(max_eqmi_idx);
            idx_outlier_seg(seg_idx)        = max_eqmi_idx;

            % ==========================================================
            % 6.3) Detecção de outlier via subconjuntos de 4 pares (EQMI adaptado) (por segmento)
            %      Remove o par identificado na etapa 6.2 e mais um dos 5 pares restantes.
            % ==========================================================

            % O primeiro outlier já foi identificado na etapa 6.2 para este segmento
            first_outlier_idx = idx_outlier_seg(seg_idx); 

            % Começamos com os 5 pares que NÃO foram o primeiro outlier
            % 'pares_5_init' são os pares que sobraram após remover o primeiro outlier
            % 'tdoa_5_init' são os TDOAs correspondentes
            idx_pares_5_init = setdiff(1:Np_total, first_outlier_idx);
            pares_5_init     = pares(idx_pares_5_init, :);
            tdoa_5_init      = current_tdoas(idx_pares_5_init); % TDOAs dos 5 pares restantes

            % Np_5_init é o número de pares que restaram (sempre 5)
            Np_5_init = size(pares_5_init, 1); 

            % Inicializa arrays para armazenar resultados para cada possível segundo outlier
            % Agora, iteramos sobre os 5 pares restantes (pares_5_init)
            err_sq_out_2_pairs_seg = zeros(Np_5_init, 1); % Erro para o segundo outlier
            az_sub_4_pairs_seg     = zeros(Np_5_init, 1); % DOA calculada com 4 pares
            el_sub_4_pairs_seg     = zeros(Np_5_init, 1); % DOA calculada com 4 pares
            uDoA_sub_4_pairs_seg   = zeros(Np_5_init, 3); % Vetor uDoA com 4 pares

            % Itera sobre cada um dos 5 pares restantes, removendo um por vez para formar 4 pares
            for k_second_outlier = 1:Np_5_init % k_second_outlier agora é o índice dentro de pares_5_init (de 1 a 5)

                % idx_keep_4_pairs são os índices DENTRO de pares_5_init que FICAM (4 pares)
                idx_keep_4_pairs = setdiff(1:Np_5_init, k_second_outlier); 

                % O par que está sendo removido NESTA iteração (o segundo outlier)
                second_outlier_pair_idx_in_5_init = k_second_outlier; % Índice dentro de pares_5_init
                second_outlier_original_idx = idx_pares_5_init(second_outlier_pair_idx_in_5_init); % Índice original do par

                % Subconjunto de 4 pares e seus TDOAs
                pares_4_pairs = pares_5_init(idx_keep_4_pairs, :);
                tdoa_4_pairs  = tdoa_5_init(idx_keep_4_pairs);

                % Calcula a DoA com os 4 pares restantes
                [az_hat_4_pairs, el_hat_4_pairs, uDoA_hat_4_pairs] = calculateDoA_LS(tdoa_4_pairs, pares_4_pairs, pos_h, c);

                % Armazena os resultados da DoA para este subconjunto de 4 pares
                az_sub_4_pairs_seg(k_second_outlier)   = az_hat_4_pairs;
                el_sub_4_pairs_seg(k_second_outlier)   = el_hat_4_pairs;
                uDoA_sub_4_pairs_seg(k_second_outlier, :) = uDoA_hat_4_pairs;

                % Calcula o TDOA predito para o segundo par removido, usando a DoA dos 4 pares
                tdoa_left_out_second_pair = tdoa_5_init(second_outlier_pair_idx_in_5_init); % TDOA real do segundo par removido

                % Geometria do segundo par removido (usando o índice original)
                ii2 = pares(second_outlier_original_idx, 1);
                jj2 = pares(second_outlier_original_idx, 2);
                delta_p_ij2_over_c = (pos_h(ii2, :) - pos_h(jj2, :)) / c;

                % TDOA predito para o segundo par removido
                tau_pred_second_pair = dot(delta_p_ij2_over_c, uDoA_hat_4_pairs); 

                % Calcula o erro quadrático para o segundo par removido
                err_sq_out_2_pairs_seg(k_second_outlier) = (tdoa_left_out_second_pair - tau_pred_second_pair)^2;
            end

            % AQUI: Identifica qual dos 5 "segundos outliers" resulta no menor erro
            % e seleciona a DOA correspondente calculada com os 4 pares
            [~, best_k_second_outlier] = max(err_sq_out_2_pairs_seg);

            az_per_segment_4_pairs(seg_idx) = az_sub_4_pairs_seg(best_k_second_outlier);
            el_per_segment_4_pairs(seg_idx) = el_sub_4_pairs_seg(best_k_second_outlier);
                                  
        end % Fim do loop de segmentos
   
        DOA_az_6_all(idx_case, idx_rng) = az_per_segment_6_pairs;
        DOA_el_6_all(idx_case, idx_rng) = el_per_segment_6_pairs;
        
        DOA_az_5_all(idx_case, idx_rng) = az_per_segment_5_pairs;
        DOA_el_5_all(idx_case, idx_rng) = el_per_segment_5_pairs;
        
        DOA_az_4_all(idx_case, idx_rng) = az_per_segment_4_pairs;
        DOA_el_4_all(idx_case, idx_rng) = el_per_segment_4_pairs;

        % --- CÁLCULOS NECESSÁRIOS PARA O ERRO ANGULAR (permanecem os mesmos) ---
        % Vetor unitário da DOA real
        u_real_vec = [cosd(el_real)*cosd(az_real); ...
                  cosd(el_real)*sind(az_real); ...
                  sind(el_real)];

        % --- Para 6 Pares ---
        u_6 = azel2unitvec(az_per_segment_6_pairs, el_per_segment_6_pairs);
        ang_err_6 = ang_err_deg(u_real_vec, u_6);
        err_med_6_all(idx_rng) = ang_err_6; % Armazena o erro angular (não o quadrado)

        % --- Para 5 Pares (Outlier por Segmento) ---
        u_5 = azel2unitvec(az_per_segment_5_pairs, el_per_segment_5_pairs);
        ang_err_5 = ang_err_deg(u_real_vec, u_5);
        err_med_5_all(idx_rng) = ang_err_5; % Armazena o erro angular (não o quadrado)

        % --- Para 4 Pares (2 Outliers por Segmento) ---
        u_4 = azel2unitvec(az_per_segment_4_pairs, el_per_segment_4_pairs);
        ang_err_4 = ang_err_deg(u_real_vec, u_4);
        err_med_4_all(idx_rng) = ang_err_4; % Armazena o erro angular (não o quadrado)
        
    end
end

%% CÁLCULOS FINAIS PARA A TABELA (APÓS TODOS OS LOOPS)
% Para 6 Pares
final_rmse_6         = sqrt(mean(err_med_6_all.^2)); % RMSE: Raiz da Média dos Quadrados dos Erros
final_std_error_6    = std(err_med_6_all);           % Desvio Padrão dos Erros

% Para 5 Pares
final_rmse_5         = sqrt(mean(err_med_5_all.^2));
final_std_error_5    = std(err_med_5_all);

% Para 4 Pares
final_rmse_4         = sqrt(mean(err_med_4_all.^2));
final_std_error_4    = std(err_med_4_all);

% --- Exibir os resultados para a tabela ---
fprintf('\n--- Resultados para Tabela de Performance da DOA (SNR = %.1f) ---\n', SNR_dB);
fprintf('Configuração |   RMSE (°)  | Desvio Padrão (°) \n');
fprintf('-----------------------------------------------\n');
fprintf('6 pares      |    %.1f      |       %.1f\n', final_rmse_6, final_std_error_6);
fprintf('5 pares      |    %.1f      |       %.1f\n', final_rmse_5, final_std_error_5);
fprintf('4 pares      |    %.1f      |       %.1f\n', final_rmse_4, final_std_error_4);
fprintf('-----------------------------------------------\n');


%% CÁLCULOS PARA PLOTAGEM (MEDIANAS E MÉDIAS DE AZIMUTE/ELEVAÇÃO)
% Estas variáveis são necessárias para plotar os pontos de mediana/média no gráfico de dispersão
DOA_az_6_all_median = median(DOA_az_6_all);
DOA_el_6_all_median = median(DOA_el_6_all);
DOA_az_6_all_mean   = mean(DOA_az_6_all);
DOA_el_6_all_mean   = mean(DOA_el_6_all);

DOA_az_5_all_median = median(DOA_az_5_all);
DOA_el_5_all_median = median(DOA_el_5_all);
DOA_az_5_all_mean   = mean(DOA_az_5_all);
DOA_el_5_all_mean   = mean(DOA_el_5_all);

DOA_az_4_all_median = median(DOA_az_4_all);
DOA_el_4_all_median = median(DOA_el_4_all);
DOA_az_4_all_mean   = mean(DOA_az_4_all);
DOA_el_4_all_mean   = mean(DOA_el_4_all);
%% Plot do Erro Angular vs. Número de Pares 

figure('Position', [100 100 1000 2500], 'Color', 'w');
hold on;
grid on;
% Ajuste os limites se necessário, ou remova para que o MATLAB ajuste automaticamente
ylim ([10 50])
xlim ([-40 40])

% Caso 6 pares (azul)
plot(DOA_az_6_all, DOA_el_6_all, 'bo', 'DisplayName', '6 Pairs');
% Caso 5 pares (verde)
plot(DOA_az_5_all, DOA_el_5_all, 'go', 'DisplayName', '5 Pairs');
% Caso 4 pares (vermelho)
plot(DOA_az_4_all, DOA_el_4_all, 'ro', 'DisplayName', '4 Pairs');

% Plotar a DOA real (estrela amarela)
plot(az_real, el_real, 'p', 'MarkerSize', 40, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', 'Real DOA');

% Plotar as medianas para cada caso (cruz)
plot(DOA_az_6_all_median, DOA_el_6_all_median, '^', 'MarkerSize', 30, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', '6 P. Median');
plot(DOA_az_5_all_median, DOA_el_5_all_median, '^', 'MarkerSize', 30, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', '5 P. Median');
plot(DOA_az_4_all_median, DOA_el_4_all_median, '^', 'MarkerSize', 30, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', '4 P. Median');
% Plotar as médias para cada caso (asterisco)
plot(DOA_az_6_all_mean, DOA_el_6_all_mean, 'd', 'MarkerSize', 30, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', '6 P. Mean');
plot(DOA_az_5_all_mean, DOA_el_5_all_mean, 'd', 'MarkerSize', 30, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', '5 P. Mean');
plot(DOA_az_4_all_mean, DOA_el_4_all_mean, 'd', 'MarkerSize', 30, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', '4 P. Mean');


xlabel('Azimuth (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Elevation (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Create the legend using the DisplayNames from the plots
h_legend = legend('Location', 'northeast');
set(h_legend, 'FontSize', 30);
set(gca, 'FontSize', 30);
set(gca, 'Position', [0.12 0.15 0.83 0.78]);

hold off;
print(fullfile(outputDir, 'erro_vs_num_segments_6_5_4_pairs_avg_noise.eps'), '-depsc');



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

        % Se quiser salvar, descomente a linha abaixo e defina 'outputDir' no seu script principal
        % print(fullfile(outputDir, sprintf('phat_pair_%d%d_seg%02d.eps', i, j, seg_idx_to_plot)), '-depsc')
    end
    fprintf('Plotagem concluída para o segmento %d.\n', seg_idx_to_plot);
end
