clear; clc; %close all;

%% 1) Parameters
fs = 44100;
c  = 1500;

dur_s = 80e-3;
num_segments = 90; % Número de segmentos para o loop
total_dur_s = num_segments * dur_s; % Duração total simulada

% t = (0:round(dur_s*fs)-1)' / fs;

SNR_dB = -12;

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
figure('Position', [100 100 900 900], 'Color', 'w');
hold on; grid on; axis equal;

% Margens (apenas 1 eixo na figura)
left_margin   = 0;
right_margin  = 0;
bottom_margin = 0.08;
top_margin    = 0.01;

ax_width_3d  = 1 - left_margin - right_margin;
ax_height_3d = 1 - top_margin  - bottom_margin;
set(gca, 'Position', [left_margin, bottom_margin, ax_width_3d, ax_height_3d]);

% Altura do tetraedro (para limites e anotação)
H_calc = L * sqrt(2/3);

% Plotar os 4 hidrofones
h_hydrophones = scatter3(pos_h(:,1), pos_h(:,2), pos_h(:,3), ...
    150, 'filled', ...
    'MarkerFaceColor', [0 0.4470 0.7410], ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);

% Rotular cada hidrofone
for k = 1:4
    text(pos_h(k,1), pos_h(k,2), pos_h(k,3) + 0.05, ...
        sprintf('  H%d', k), ...
        'FontSize', 20, 'FontWeight', 'bold', ...
        'Color', 'k', 'Interpreter', 'latex');
end

% Desenhar TODAS as 6 arestas do tetraedro
edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for i = 1:size(edges, 1)
    p1 = pos_h(edges(i,1), :);
    p2 = pos_h(edges(i,2), :);
    plot3([p1(1) p2(1)], ...
          [p1(2) p2(2)], ...
          [p1(3) p2(3)], ...
          'k-', 'LineWidth', 1.5);
end

% Distância máxima entre pares (anotação)
text(0.1, 0.1, H_calc + 0.15, ...
    sprintf('$d_{\\mathrm{max}} = %.3f$ m', L), ...
    'FontSize', 20, 'Color', [0 0.5 0], ...
    'Interpreter', 'latex', 'FontWeight', 'bold');

% Projeções de H3
H3    = pos_h(3,:);           % [x3, y3, z3]
H3_xy = [H3(1), H3(2), 0];    % projeção no plano XY

% Projeção de H3 no eixo X (a partir do plano XY)
plot3([H3_xy(1) H3_xy(1)], ...
      [H3_xy(2) 0], ...
      [0 0], ...
      'b--', 'LineWidth', 1.2);

% Projeção de H3 no eixo Y (a partir do plano XY)
plot3([H3_xy(1) (L + 0.1)], ...
      [H3_xy(2) H3_xy(2)], ...
      [0 0], ...
      'b--', 'LineWidth', 1.2);

% Projeções de H4
H4    = pos_h(4,:);           % [x4, y4, z4]
H4_xy = [H4(1), H4(2), 0];    % projeção no plano XY

% (a) Projeção de H4 no plano XY – linha pontilhada vertical
plot3([H4(1) H4_xy(1)], ...
      [H4(2) H4_xy(2)], ...
      [H4(3) H4_xy(3)], ...
      'g--', 'LineWidth', 1.2);

% (b) Projeção de H4 no plano YZ
H4_yz = [0, H4(2), H4(3)];
plot3([H4(1) H4_yz(1)], ...
      [H4(2) H4_yz(2)], ...
      [H4(3) H4_yz(3)], ...
      'g--', 'LineWidth', 1.2);

% (c) Do ponto no plano YZ até o eixo Z
H4_z = [0, 0, H4(3)];
plot3([H4_yz(1) H4_z(1)], ...
      [H4_yz(2) H4_z(2)], ...
      [H4_yz(3) H4_z(3)], ...
      'g--', 'LineWidth', 1.2);

% Eixos
xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', 20);

% Ajustar visualização
view(45, 30);
xlim([-0.01, L + 0.1]);
ylim([-0.01, L*sqrt(3)/2 + 0.1]);
zlim([0,     H_calc + 0.2]);
set(gca, 'FontSize', 20);
box on;

% Legenda – apenas hidrofones
legend(h_hydrophones, {'Hydrophones'}, ...
    'Location', 'northeast', ...
    'Interpreter', 'latex', 'FontSize', 20);

hold off;

% Salvar figura
outputDir = './figures';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
print(fullfile(outputDir, 'array_MATLAB.eps'), '-depsc');
%}
%% 2) Physical maxLag
dmax          = max_pair_distance(pos_h);
maxLagSamples = ceil((dmax/c)*fs) + 2;
fprintf('\nmaxLagSamples = %d samples  %.3f ms\n', ...
        maxLagSamples, 1000*maxLagSamples/fs);

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
    fprintf('Sinal reamostrado de %d Hz para %d Hz\n', fs_wav, fs);
end

s = s_raw(:);   % garante coluna

fprintf('Sinal lido: %d amostras  %.1f ms  fs=%d Hz\n', ...
        length(s), 1e3*length(s)/fs, fs);

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
    fprintf('Ruido reamostrado de %d Hz para %d Hz\n', fs_noise, fs);
end

noise_env_raw = noise_env_raw(:);

fprintf('Ruido lido: %d amostras  %.1f ms\n', ...
        length(noise_env_raw), 1e3*length(noise_env_raw)/fs);


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
num_seg_list0 = 1:1:1;    % por exemplo, de 25 até 200
num_cases0    = length(num_seg_list0);

for idx_rng = 1:num_cases0
    rng(idx_rng);

%% Varredura no número de segmentos
num_seg_list = [1 8 32 64 128 256];    % por exemplo, de 25 até 200
num_cases    = length(num_seg_list);
max_num_segments = max(num_seg_list);

err_med_6   = zeros(num_cases,1); % erro angular médio (ou mediano) 6 pares
err_med_5   = zeros(num_cases,1); % erro angular médio (ou mediano) 5 pares

% Use NaN para preencher, assim zeros não afetam a média/mediana
az_per_segment_6_pairs = nan(max_num_segments, num_cases); % (segmentos x casos x realizacoes)
el_per_segment_6_pairs= nan(max_num_segments, num_cases);

az_per_segment_5_pairs = nan(max_num_segments, num_cases);
el_per_segment_5_pairs = nan(max_num_segments, num_cases);

az_per_segment_4_pairs = nan(max_num_segments, num_cases);
el_per_segment_4_pairs = nan(max_num_segments, num_cases);

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
%               --- plotar a correlação de um segmento específico:
    %             if seg_idx == 1  % por exemplo, só o primeiro segmento
    %                 fig = figure('Position', [100 100 1200 600], 'Color', 'w');
    %                 plot(1e3*lagsP, rP, 'k', 'LineWidth', 1.05); grid on; hold on;
    % 
    %                 yl = ylim;
    %                 h_hat = plot(1e3*[tdoa_phat(seg_idx,p) tdoa_phat(seg_idx,p)], yl, 'm-',  'LineWidth', 1.2);
    %                 h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)],         yl, 'b--', 'LineWidth', 1.1);
    %                 h_pk  = plot(1e3*tdoa_phat(seg_idx,p), rP(kP), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
    %                 hold off;
    % 
    %                 xlim(1e3*[-maxLagSamples maxLagSamples]/fs);
    %                 ylim([-1 1]);
    %                 xlabel('Delay [ms]', 'FontSize', 20, 'FontWeight', 'bold');
    %                 ylabel(sprintf('Corr. (%d,%d)', i, j), 'FontSize', 20, 'FontWeight', 'bold');
    % 
    %                 legend([h_hat h_the h_pk], { ...
    %                     sprintf('Est.   \\tau_{%d%d} = %.3f ms', i, j, 1e3*tdoa_phat(seg_idx,p)), ...
    %                     sprintf('Theor. \\tau_{%d%d} = %.3f ms', i, j, 1e3*tdoa_th_plane(p)), ...
    %                     'Peak'}, ...
    %                     'Location', 'best', 'Interpreter', 'tex', 'FontSize', 20);
    % 
    %                 set(gca, 'FontSize', 20);
    %                 set(gca, 'Position', [0.08 0.18 0.90 0.75])
    %                 % print(fullfile(outputDir, sprintf('phat_pair_%d%d_seg%02d.eps', i, j, seg_idx)), '-depsc')
    %             end
    %}
            end

        end
        
        %% Adicionar erro com aleatoriedade (Ruído Uniforme)
        %{ 
        erro_indice_par = 6; 
        % erro_s = 0 * 1e-3;
        % 
        % var_erro_s = 0 * 1e-3; %[-var_erro_s, +var_erro_s]
        % 
        % % rng(42); % Opcional: para reprodutibilidade
        % 
        % random_noise_s = var_erro_s * (2 * rand(num_segments, 1) - 1);
        % tdoa_phat(:, erro_indice_par) = tdoa_phat(:, erro_indice_par) + erro_s + random_noise_s;
        % 
        % fprintf('Erro base de %.3f ms (com ruído Uniforme de +/- %.3f ms) adicionado ao TDOA do par %d em todos os %d segmentos.\n', ...
        %         erro_s*1e3, var_erro_s*1e3, erro_indice_par, num_segments);
        %} 

        %% Plot do espectro do ruído puro 
        %{
        px = mean(s.^2);
        pn = px / (10^(SNR_dB/10));

        rng(123);
        noise_ref = sqrt(pn) * randn(4*N, 1);   % ruído longo para estimativa estável


        % (não regera com rng(42), que é outra realização)
        noise_plot = noise_ref(:) - mean(noise_ref);   % mesmo noise_ref do thr
        M_half     = floor(length(noise_plot) / 2);

        seg1_plot  = noise_plot(1        : M_half) - mean(noise_plot(1:M_half));
        seg2_plot  = noise_plot(M_half+1 : 2*M_half);
        seg2_plot  = seg2_plot - mean(seg2_plot);

        % FFT com o mesmo Nfft do gcc_family_local
        Xn1_plot = fft(seg1_plot, Nfft);
        Xn2_plot = fft(seg2_plot, Nfft);

        R_noise_plot   = Xn1_plot .* conj(Xn2_plot);
        mag_noise_plot = abs(R_noise_plot);

        % Vetor de frequências
        freqs_plot   = (0 : Nfft-1)' * fs / Nfft;
        idx_pos_plot = freqs_plot <= fs/2;
        freqs_pos_plot = freqs_plot(idx_pos_plot);

        % Diagnóstico: verifica se thr está dentro do ylim
        max_mag_plot = max(mag_noise_plot(idx_pos_plot));
        fprintf('max(|R_noise|) plotado : %.4e\n', max_mag_plot);

        % ----------------------------------------------------------
        % FIGURA: Cross-power spectrum do ruído puro
        % ----------------------------------------------------------

        fig = figure('Position', [100 100 2000 800], 'Color', 'w');
        hold on

        h1 = plot(freqs_pos_plot, mag_noise_plot(idx_pos_plot), ...
                  'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);


        grid on
        xlabel('Frequency [Hz]', 'FontSize', 40, 'FontWeight', 'bold')
        ylabel('|R_{noise}|',    'FontSize', 40, 'FontWeight', 'bold')
        xlim([0 1000])

        set(gca, 'FontSize', 40)
        set(gca, 'Position', [0.15 0.25 0.73 0.65])


        hold off
        print(fullfile(outputDir, 'mask_noise_only.eps'), '-depsc')

        %% -------------------------------------------------------
        %  FIGURA: Comparação sinal+ruído vs ruído puro
        % -------------------------------------------------------

        x1 = x(:,1) - mean(x(:,1));
        x2 = x(:,2) - mean(x(:,2));
        X1 = fft(x1, Nfft);
        X2 = fft(x2, Nfft);
        R_sig   = X1 .* conj(X2);
        mag_sig = abs(R_sig);

        % Limites do plot definidos antes das linhas horizontais
        x_lim = [0 1000];
        y_top = 1.2 * max(mag_sig(idx_pos_plot));

        fig2 = figure('Position', [100 100 1500 750], 'Color', 'w');
        hold on

        % |R| sinal + ruído (par 1-2)
        h1 = plot(freqs_pos_plot, mag_sig(idx_pos_plot), ...
                  'Color', [0 0.4470 0.7410], 'LineWidth', 2);

        % |R_noise| ruído puro — usa mag_noise_plot já calculado
        h2 = plot(freqs_pos_plot, mag_noise_plot(idx_pos_plot), ...
                  'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

        % Threshold único (calculado na banda 0-1000 Hz)

        grid on
        xlabel('Frequency [Hz]', 'FontSize', 40, 'FontWeight', 'bold')
        ylabel('|R|',            'FontSize', 40, 'FontWeight', 'bold')
        xlim(x_lim)

        set(gca, 'FontSize', 40)
        set(gca, 'Position', [0.10 0.25 0.85 0.65])


        hold off
        print(fullfile(outputDir, 'mask_signal_vs_noise.eps'), '-depsc')

        %}

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

        for seg_idx = 1:num_segments
            current_tdoas = tdoa_phat(seg_idx, :)'; % TDOAs para este segmento [Np x 1]

            % ==========================================================
            % 6.1) LS com TODAS as 6 correlações (por segmento)
            % ==========================================================
            [az_seg_6, el_seg_6, ~] = calculateDoA_LS(current_tdoas, pares, pos_h, c);
            az_per_segment_6_pairs(seg_idx, idx_case) = az_seg_6;
            el_per_segment_6_pairs(seg_idx, idx_case) = el_seg_6;

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

            az_per_segment_5_pairs(seg_idx, idx_case) = az_sub_seg(max_eqmi_idx);
            el_per_segment_5_pairs(seg_idx, idx_case) = el_sub_seg(max_eqmi_idx);
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

            az_per_segment_4_pairs(seg_idx, idx_case) = az_sub_4_pairs_seg(best_k_second_outlier);
            el_per_segment_4_pairs(seg_idx, idx_case) = el_sub_4_pairs_seg(best_k_second_outlier);
            
        end % Fim do loop de segmentos
        
      
        % ==========================================================
        % Mediana Padrão das DOAs
        % ==========================================================
        DOA_az_6_median = nanmedian(az_per_segment_6_pairs, 1);
        DOA_el_6_median = nanmedian(el_per_segment_6_pairs, 1);
        
        DOA_az_5_median = nanmedian(az_per_segment_5_pairs, 1);
        DOA_el_5_median = nanmedian(el_per_segment_5_pairs, 1);
        
        DOA_az_4_median = nanmedian(az_per_segment_4_pairs, 1);
        DOA_el_4_median = nanmedian(el_per_segment_4_pairs, 1);
        
        % ==========================================================
        % Média Padrão das DOAs
        % ==========================================================
        
        DOA_az_6_mean = nanmean(az_per_segment_6_pairs, 1);
        DOA_el_6_mean = nanmean(el_per_segment_6_pairs, 1);
        
        DOA_az_5_mean = nanmean(az_per_segment_5_pairs, 1);
        DOA_el_5_mean = nanmean(el_per_segment_5_pairs, 1);
        
        DOA_az_4_mean = nanmean(az_per_segment_4_pairs, 1);
        DOA_el_4_mean = nanmean(el_per_segment_4_pairs, 1);
        
        % Erro e STD  
        
        % --- CÁLCULOS NECESSÁRIOS PARA O ERRO ANGULAR (permanecem os mesmos) ---
        % Vetor unitário da DOA real
        u_real_vec = [cosd(el_real)*cosd(az_real); ...
                  cosd(el_real)*sind(az_real); ...
                  sind(el_real)];
        
        % --- Para 6 Pares ---
        u_6_mean = azel2unitvec(DOA_az_6_mean, DOA_el_6_mean);
        ang_err_6_mean = ang_err_deg(u_real_vec, u_6_mean);
        
        u_6_median = azel2unitvec(DOA_az_6_median, DOA_el_6_median);
        ang_err_6_median = ang_err_deg(u_real_vec, u_6_median);
        
        % --- Para 5 Pares ---
        u_5_mean = azel2unitvec(DOA_az_5_mean, DOA_el_5_mean);
        ang_err_5_mean = ang_err_deg(u_real_vec, u_5_mean);
        
        u_5_median = azel2unitvec(DOA_az_5_median, DOA_el_5_median);
        ang_err_5_median = ang_err_deg(u_real_vec, u_5_median);
        
        % --- Para 4 Pares ---
        u_4_mean = azel2unitvec(DOA_az_4_mean, DOA_el_4_mean);
        ang_err_4_mean = ang_err_deg(u_real_vec, u_4_mean);
        
        u_4_median = azel2unitvec(DOA_az_4_median, DOA_el_4_median);
        ang_err_4_median = ang_err_deg(u_real_vec, u_4_median);
%         
        
    end

end

%% =========================================================================
% PLOT: AZIMUTE VS. ELEVAÇÃO (DOA ESTIMATES, REAL, MEDIAN, MEAN)
% COM VARIAÇÃO DE COR POR NÚMERO DE SEGMENTOS (APENAS 6 E 5 PARES)
% =========================================================================

figure('Position', [100 100 1000 2000], 'Color', 'w');
hold on;
grid on;
box on;

xlim ([-40 40])
ylim ([10 50])

% Definindo cores base para cada número de pares
color_base_6 = [0 0 1];     % Azul para 6 pares
color_base_5 = [1 0 0]; % Verde para 5 pares

% Gerar gradientes de cores para cada base
num_color_steps = length(num_seg_list); % Número de variações de segmentos
lightness_factor = 0.05; % Quão clara a cor mais clara será (0=cor base, 1=branco)

colors_gradient_6 = zeros(num_color_steps, 3);
colors_gradient_5 = zeros(num_color_steps, 3);
% colors_gradient_4 = zeros(num_color_steps, 3); (REMOVIDO)

for k = 1:num_color_steps
    alpha = (k - 1) / (num_color_steps - 1); % alpha de 0 (mais claro) a 1 (cor base)

    start_color_6 = (1 - lightness_factor) * [1 1 1] + lightness_factor * color_base_6;
    start_color_5 = (1 - lightness_factor) * [1 1 1] + lightness_factor * color_base_5;
    % start_color_4 = (1 - lightness_factor) * [1 1 1] + lightness_factor * color_base_4; (REMOVIDO)

    colors_gradient_6(k,:) = (1 - alpha) * start_color_6 + alpha * color_base_6;
    colors_gradient_5(k,:) = (1 - alpha) * start_color_5 + alpha * color_base_5;
    % colors_gradient_4(k,:) = (1 - alpha) * start_color_4 + alpha * color_base_4; (REMOVIDO)
end


% --- 2. Plotar todas as estimativas individuais por segmento ---
h_est_6_legend = []; h_est_5_legend = []; % h_est_4_legend = []; (REMOVIDO)

% --- 3. Plotar as medianas e médias para cada caso (para cada num_segments) ---
% Coletar handles para a legenda
h_med_6_legend = []; h_med_5_legend = []; % h_med_4_legend = []; (REMOVIDO)
h_mean_6_legend = []; h_mean_5_legend = []; % h_mean_4_legend = []; (REMOVIDO)


for idx_case = 1:length(num_seg_list)
    current_num_segments = num_seg_list(idx_case);
    current_color_6 = colors_gradient_6(idx_case, :);
    current_color_5 = colors_gradient_5(idx_case, :);
    % current_color_4 = colors_gradient_4(idx_case, :); (REMOVIDO)

    % --- Plotar Estimativas Individuais ---
    % Para 6 Pares
    az_6_valid = az_per_segment_6_pairs(1:current_num_segments, idx_case);
    el_6_valid = el_per_segment_6_pairs(1:current_num_segments, idx_case);
    valid_idx_6 = ~isnan(az_6_valid) & ~isnan(el_6_valid);
    if any(valid_idx_6)
        h_plot = plot(az_6_valid(valid_idx_6), el_6_valid(valid_idx_6), 'o', 'MarkerSize', 8, 'Color', current_color_6, 'HandleVisibility', 'off');
        if isempty(h_est_6_legend) % Guarda o handle do primeiro plot para a legenda
            h_est_6_legend = h_plot;
        end
    end

    % Para 5 Pares
    az_5_valid = az_per_segment_5_pairs(1:current_num_segments, idx_case);
    el_5_valid = el_per_segment_5_pairs(1:current_num_segments, idx_case);
    valid_idx_5 = ~isnan(az_5_valid) & ~isnan(el_5_valid);
    if any(valid_idx_5)
        h_plot = plot(az_5_valid(valid_idx_5), el_5_valid(valid_idx_5), 'o', 'MarkerSize', 8, 'Color', current_color_5, 'HandleVisibility', 'off');
        if isempty(h_est_5_legend)
            h_est_5_legend = h_plot;
        end
    end

    % Para 4 Pares (REMOVIDO)
    % az_4_valid = az_per_segment_4_pairs(1:current_num_segments, idx_case);
    % el_4_valid = el_per_segment_4_pairs(1:current_num_segments, idx_case);
    % valid_idx_4 = ~isnan(az_4_valid) & ~isnan(el_4_valid);
    % if any(valid_idx_4)
    %     h_plot = plot(az_4_valid(valid_idx_4), el_4_valid(valid_idx_4), 'o', 'MarkerSize', 8, 'Color', current_color_4, 'HandleVisibility', 'off');
    %     if isempty(h_est_4_legend)
    %         h_est_4_legend = h_plot;
    %     end
    % end
    
    % --- 1. Plotar a DOA real (estrela amarela) ---
    h_real = plot(az_real, el_real, 'p', 'MarkerSize', 40, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 3, 'DisplayName', 'Real DOA');


    % --- Plotar Mediana para o num_segments atual ---
    % Para 6 Pares
    if ~isnan(DOA_az_6_median(idx_case)) && ~isnan(DOA_el_6_median(idx_case))
        h_plot = plot(DOA_az_6_median(idx_case), DOA_el_6_median(idx_case), '^', 'MarkerSize', 30, 'MarkerFaceColor', current_color_6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
        if isempty(h_med_6_legend)
            h_med_6_legend = h_plot;
        end
    end
    % Para 5 Pares
    if ~isnan(DOA_az_5_median(idx_case)) && ~isnan(DOA_el_5_median(idx_case))
        h_plot = plot(DOA_az_5_median(idx_case), DOA_el_5_median(idx_case), '^', 'MarkerSize', 30, 'MarkerFaceColor', current_color_5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
        if isempty(h_med_5_legend)
            h_med_5_legend = h_plot;
        end
    end
    % Para 4 Pares (REMOVIDO)
    % if ~isnan(DOA_az_4_median(idx_case)) && ~isnan(DOA_el_4_median(idx_case))
    %     h_plot = plot(DOA_az_4_median(idx_case), DOA_el_4_median(idx_case), '^', 'MarkerSize', 15, 'MarkerFaceColor', current_color_4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    %     if isempty(h_med_4_legend)
    %         h_med_4_legend = h_plot;
    %     end
    % end

    % --- Plotar Média para o num_segments atual ---
    % Para 6 Pares
    if ~isnan(DOA_az_6_mean(idx_case)) && ~isnan(DOA_el_6_mean(idx_case))
        h_plot = plot(DOA_az_6_mean(idx_case), DOA_el_6_mean(idx_case), 'd', 'MarkerSize', 30, 'MarkerFaceColor', current_color_6, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
        if isempty(h_mean_6_legend)
            h_mean_6_legend = h_plot;
        end
    end
    % Para 5 Pares
    if ~isnan(DOA_az_5_mean(idx_case)) && ~isnan(DOA_el_5_mean(idx_case))
        h_plot = plot(DOA_az_5_mean(idx_case), DOA_el_5_mean(idx_case), 'd', 'MarkerSize', 30, 'MarkerFaceColor', current_color_5, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
        if isempty(h_mean_5_legend)
            h_mean_5_legend = h_plot;
        end
    end
    % Para 4 Pares (REMOVIDO)
    % if ~isnan(DOA_az_4_mean(idx_case)) && ~isnan(DOA_el_4_mean(idx_case))
    %     h_plot = plot(DOA_az_4_mean(idx_case), DOA_el_4_mean(idx_case), 'd', 'MarkerSize', 15, 'MarkerFaceColor', current_color_4, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    %     if isempty(h_mean_4_legend)
    %         h_mean_4_legend = h_plot;
    %     end
    % end
end

% Definir DisplayName e HandleVisibility para os handles de legenda coletados
% Estimativas Individuais
if ~isempty(h_est_6_legend)
    set(h_est_6_legend, 'DisplayName', '6 Pairs', 'HandleVisibility', 'on');
    % MUDANÇA AQUI: Definir a cor do marcador na legenda para a cor mais escura (último passo do gradiente)
    set(h_est_6_legend, 'Color', colors_gradient_6(end,:)); % <--- MUDANÇA
else
    h_est_6_legend = plot(NaN, NaN, 'o', 'MarkerSize', 30, 'Color', colors_gradient_6(end,:), 'DisplayName', '6 Pairs (Segment Estimates)'); % Placeholder com cor mais escura
end
if ~isempty(h_est_5_legend)
    set(h_est_5_legend, 'DisplayName', '5 Pairs', 'HandleVisibility', 'on');
    % MUDANÇA AQUI: Definir a cor do marcador na legenda para a cor mais escura (último passo do gradiente)
    set(h_est_5_legend, 'Color', colors_gradient_5(end,:)); % <--- MUDANÇA
else
    h_est_5_legend = plot(NaN, NaN, 'o', 'MarkerSize', 30, 'Color', colors_gradient_5(end,:), 'DisplayName', '5 Pairs (Segment Estimates)'); % Placeholder com cor mais escura
end

% Mediana
if ~isempty(h_med_6_legend)
    set(h_med_6_legend, 'DisplayName', '6 P. Median', 'HandleVisibility', 'on');
    % MUDANÇA AQUI: Definir a cor do marcador na legenda para a cor mais escura (último passo do gradiente)
    set(h_med_6_legend, 'MarkerFaceColor', colors_gradient_6(end,:)); % <--- MUDANÇA
else
    h_med_6_legend = plot(NaN, NaN, '^', 'MarkerSize', 30, 'MarkerFaceColor', colors_gradient_6(end,:), 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 Pairs (Median)'); % Placeholder com cor mais escura
end
if ~isempty(h_med_5_legend)
    set(h_med_5_legend, 'DisplayName', '5 P. Median', 'HandleVisibility', 'on');
    % MUDANÇA AQUI: Definir a cor do marcador na legenda para a cor mais escura (último passo do gradiente)
    set(h_med_5_legend, 'MarkerFaceColor', colors_gradient_5(end,:)); % <--- MUDANÇA
else
    h_med_5_legend = plot(NaN, NaN, '^', 'MarkerSize', 30, 'MarkerFaceColor', colors_gradient_5(end,:), 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 Pairs (Median)'); % Placeholder com cor mais escura
end

% Média
if ~isempty(h_mean_6_legend)
    set(h_mean_6_legend, 'DisplayName', '6 P. Mean', 'HandleVisibility', 'on');
    % MUDANÇA AQUI: Definir a cor do marcador na legenda para a cor mais escura (último passo do gradiente)
    set(h_mean_6_legend, 'MarkerFaceColor', colors_gradient_6(end,:)); % <--- MUDANÇA
else
    h_mean_6_legend = plot(NaN, NaN, 'd', 'MarkerSize', 30, 'MarkerFaceColor', colors_gradient_6(end,:), 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '6 Pairs (Mean)'); % Placeholder com cor mais escura
end
if ~isempty(h_mean_5_legend)
    set(h_mean_5_legend, 'DisplayName', '5 P. Mean', 'HandleVisibility', 'on');
    % MUDANÇA AQUI: Definir a cor do marcador na legenda para a cor mais escura (último passo do gradiente)
    set(h_mean_5_legend, 'MarkerFaceColor', colors_gradient_5(end,:)); % <--- MUDANÇA
else
    h_mean_5_legend = plot(NaN, NaN, 'd', 'MarkerSize', 30, 'MarkerFaceColor', colors_gradient_5(end,:), 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', '5 Pairs (Mean)'); % Placeholder com cor mais escura
end


xlabel('Azimuth (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('Elevation (degrees)', 'FontSize', 40, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Criar a legenda com todos os elementos
h_legend = legend([h_real, h_est_6_legend, h_est_5_legend, ... % Removido h_est_4_legend
                  h_med_6_legend, h_med_5_legend, ... % Removido h_med_4_legend
                  h_mean_6_legend, h_mean_5_legend], ... % Removido h_mean_4_legend
                  'Location', 'southwest');
set(h_legend, 'FontSize', 30);
set(gca, 'FontSize', 30);
set(gca, 'Position', [0.14 0.15 0.83 0.78]);

% Adicionar a legenda de cores para os segmentos (manual)
ax_seg_legend = axes('Position', [0.75 0.15 0.1 0.25], 'Visible', 'off'); % Posição ajustável
text(-0.005, 1.1, 'Segments:', 'FontSize', 25, 'FontWeight', 'bold', 'Parent', ax_seg_legend);

for k = 1:num_color_steps
    current_seg_val = num_seg_list(k);
    current_color = colors_gradient_6(k,:); % Usamos o gradiente de 6 pares como exemplo

    patch([0 0.1 0.1 0], [1 - 0.15*k, 1 - 0.15*k, 1 - 0.15*k + 0.1, 1 - 0.15*k + 0.1], current_color, ...
          'EdgeColor', 'k', 'Parent', ax_seg_legend);
    text(0.15, 1 - 0.15*k + 0.05, sprintf('%d', current_seg_val), 'FontSize', 20, 'Parent', ax_seg_legend);
end

hold off;

print(fullfile(outputDir, 'doa_az_el_plot_6_5_pairs_segments_color.eps'), '-depsc');


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
function u_vec = azel2unitvec(az_deg, el_deg)
    % Converte azimute e elevação em graus para um vetor unitário 3D.
    % az_deg e el_deg podem ser escalares ou vetores.
    u_vec(1,:) = cosd(el_deg) .* cosd(az_deg);
    u_vec(2,:) = cosd(el_deg) .* sind(az_deg);
    u_vec(3,:) = sind(el_deg);
end

function ang_err = ang_err_deg(u_real, u_est)
    % Calcula o erro angular em graus entre dois vetores unitários.
    % u_real: vetor unitário real (3x1)
    % u_est: vetor(es) unitário(s) estimado(s) (3xN)
    % Retorna um escalar ou vetor de erros angulares.

    % Garante que u_real seja um vetor coluna
    u_real = u_real(:); 

    % Calcula o produto escalar e limita para evitar erros de floating point
    dot_product = u_real' * u_est;
    dot_product = max(-1, min(1, dot_product));

    ang_err = acosd(dot_product);
end