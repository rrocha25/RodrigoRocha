    % =========================================================================
%% SÍNTESE DE ÁUDIO A PARTIR DO ESPECTRO DO CSV (NÍVEL DE SPL EXATO)
% =========================================================================
clear; clc; close all;

% -------------------------------------------------------------------------
%% 1. IMPORTAR CSV
% -------------------------------------------------------------------------
dados = importdata('plot-data.csv');
freq = dados.data(:,1);
SPL = dados.data(:,2);

[freq, idx] = unique(freq,'stable'); %unique: Remove valores duplicados do vetor freq
SPL = SPL(idx);

% -------------------------------------------------------------------------
%% 2. SPL ? Pressão (Pa)
% -------------------------------------------------------------------------
p_ref = 1e-6; %Pressão de referência padrão em água = 1 µPa
P = p_ref * 10.^(SPL/20); %Vetor P com amplitudes de pressão em Pa para cada frequência

% -------------------------------------------------------------------------
%% 3. GRID UNIFORME PARA IFFT
% -------------------------------------------------------------------------
fs = 44100;
dur = 60;
N = fs * dur;
f_ifft = linspace(0, fs/2, floor(N/2)+1)'; 
% A IFFT requer frequências igualmente espaçadas

P_interp = interp1(freq, P, f_ifft, 'linear', 'extrap');
P_interp(P_interp<0) = 0; 
% Atribui zero: Pressão não pode ser negativa

% -------------------------------------------------------------------------
%% 4. CONSTRUÇÃO DO ESPECTRO COMPLEXO
% -------------------------------------------------------------------------
fase = 2*pi*rand(size(P_interp));
%Cria fases aleatórias para cada componente de frequência (em radianos)
% Fases aleatórias uniformes geram ruído realista com espectro controlado
% Garante propriedades estatísticas adequadas para simular ruído ambiente

S = P_interp .* exp(1i*fase); 
% Vetor complexo com: Magnitude: P_interp (pressão em Pa), Fase: fase (aleatória)

S_full = [S ; conj(S(end-1:-1:2))];
% Inverte o vetor S, excluindo primeiro (DC) e último (Nyquist) elementos
% Calcula o conjugado complexo (inverte o sinal da parte imaginária)
% Concatena verticalmente

% -------------------------------------------------------------------------
%% 5. IFFT ? SINAL TEMPORAL
% -------------------------------------------------------------------------
x = real(ifft(S_full,'symmetric'));
x = x(1:N); % Garante exatamente 60 segundos de áudio

% -------------------------------------------------------------------------
%% 6. AJUSTE ABSOLUTO: ÁUDIO FICA COM O MESMO SPL DO CSV
% -------------------------------------------------------------------------
energia_espectro = mean(P_interp.^2);  % Potência acústica média por unidade de área 

energia_audio = mean(x.^2);
% Energia média do sinal temporal em Pa²

ganho = sqrt(energia_espectro / energia_audio);
% Raiz da Razão entre energias (adimensional)para converter de energia (potência) para amplitude

x = x * ganho;
% Ajusta o nível absoluto para corresponder aos valores SPL do CSV

% -------------------------------------------------------------------------
%% 7. SALVAR ÁUDIO (SINAL PURO)
% -------------------------------------------------------------------------
audiowrite('auv_chen_exato.wav', x, fs);

% -------------------------------------------------------------------------
%% 8. CALCULAR ESPECTRO DO SINAL PURO
% -------------------------------------------------------------------------
Y = fft(x);
f_fft_all = (0:N-1)' * fs / N;
% Criação do Vetor de Frequências Completo

f_fft_pos = f_fft_all(1:floor(N/2)+1);
% Extração das Frequências Positivas

PSD_sinal = (abs(Y(1:floor(N/2)+1)).^2) / N;
% Cálculo da Densidade Espectral de Potência (PSD)

PSD_sinal(2:end-1) = 2 * PSD_sinal(2:end-1);
% Correção para Espectro Unilateral
% (2*) porque estamos descartando as frequências negativas
% Exceções:
%   DC (índice 1): Não multiplica (não tem par negativo)
%   Nyquist (end): Não multiplica (é único, não tem par)
% Resultado: PSD unilateral com energia total preservada

SPL_sinal = 10 * log10(PSD_sinal / (p_ref^2) + eps);
%  Conversão PSD ? SPL

% -------------------------------------------------------------------------
%% 9. CARREGAR RUÍDO DE FUNDO
% -------------------------------------------------------------------------
[ruido, fs_ruido] = audioread('Enseada dos Anjos.wav');

if size(ruido, 2) > 1
    ruido = ruido(:, 1);
end
% Conversão para Mono

if fs_ruido ~= fs
    ruido = resample(ruido, fs, fs_ruido);
end
% Reamostragem pois a taxa de amostragem do ruído difere de 44100 Hz

if length(ruido) < N
    ruido = repmat(ruido, ceil(N/length(ruido)), 1);
end
%  Extensão do Ruído para 60s 

ruido = ruido(1:N);
% Garante: Exatamente N amostras (60 segundos)

% -------------------------------------------------------------------------
%% 10. AJUSTAR SNR
% -------------------------------------------------------------------------
SNR_desejada_dB = 20;

potencia_sinal = mean(x.^2);
potencia_ruido = mean(ruido.^2);
ganho_ruido = sqrt(potencia_sinal / (potencia_ruido * 10^(SNR_desejada_dB/10)));
%  Calcula ganho necessário para escalar o ruído até atingir a SNR desejada

ruido = ruido * ganho_ruido;

SNR_final = 10*log10(potencia_sinal / mean(ruido.^2));
% Verifica a SNR final obtida

% -------------------------------------------------------------------------
%% 11. COMBINAR SINAL + RUÍDO
% -------------------------------------------------------------------------
sinal_com_ruido = x + ruido;
audiowrite('auv_chen_com_ruido.wav', sinal_com_ruido, fs);

% -------------------------------------------------------------------------
%% 12. CALCULAR ESPECTRO DO RUÍDO
% -------------------------------------------------------------------------
Y_ruido = fft(ruido);
PSD_ruido = (abs(Y_ruido(1:floor(N/2)+1)).^2) / N;
PSD_ruido(2:end-1) = 2 * PSD_ruido(2:end-1);
SPL_ruido = 10 * log10(PSD_ruido / (p_ref^2) + eps);

% -------------------------------------------------------------------------
%% 13. CALCULAR ESPECTRO DO SINAL COM RUÍDO
% -------------------------------------------------------------------------
Y_total = fft(sinal_com_ruido);
PSD_total = (abs(Y_total(1:floor(N/2)+1)).^2) / N;
PSD_total(2:end-1) = 2 * PSD_total(2:end-1);
SPL_total = 10 * log10(PSD_total / (p_ref^2) + eps);

% -------------------------------------------------------------------------
%% 14. PLOTAGEM COMPARATIVA
% -------------------------------------------------------------------------
figure('Position', [100 100 1400 900])
idx_zoom = f_fft_pos <= 1000;

subplot(3, 1, 1)
plot(f_fft_pos(idx_zoom), SPL_sinal(idx_zoom), 'b-', 'LineWidth', 1.5)
grid on
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 12, 'FontWeight', 'bold')
title('Espectro do Sinal Puro (AUV)', 'FontSize', 14, 'FontWeight', 'bold')
xlim([0 1000])
ylim([0 100])
set(gca, 'XTickLabel', [])

subplot(3, 1, 2)
plot(f_fft_pos(idx_zoom), SPL_ruido(idx_zoom), 'r-', 'LineWidth', 1.5)
grid on
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 12, 'FontWeight', 'bold')
title('Espectro do Ruído de Fundo (Enseada dos Anjos)', 'FontSize', 14, 'FontWeight', 'bold')
xlim([0 1000])
ylim([0 100])
set(gca, 'XTickLabel', [])

subplot(3, 1, 3)
plot(f_fft_pos(idx_zoom), SPL_total(idx_zoom), 'k-', 'LineWidth', 1.5)
grid on
xlabel('Frequência [Hz]', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 12, 'FontWeight', 'bold')
title(sprintf('Espectro do Sinal com Ruído (SNR = %.1f dB)', SNR_final), 'FontSize', 14, 'FontWeight', 'bold')
xlim([0 1000])
ylim([0 100])

%% 15 - Espectro + Espectrograma do Sinal AUV e Ruído
figure('Position', [100 100 1600 1000], 'Color', 'w')

% Definir margens e gaps
left_margin   = 0.06;
right_margin  = 0.02;
top_margin    = 0.05;
bottom_margin = 0.06;
vertical_gap  = 0.06;
horizontal_gap = 0.08;

% Dimensões dos subplots
plot_width  = (1 - left_margin - right_margin - horizontal_gap) / 2;
plot_height = (1 - top_margin - bottom_margin - vertical_gap) / 2;

% Índice de zoom para espectro (até 1 kHz)
idx_zoom = f_fft_pos <= 1000;

% Parâmetros para espectrograma (baseado na literatura consultada)
window_length = round(0.025 * fs);  % 25 ms (compromisso tempo-frequência)
overlap = round(0.75 * window_length);  % 75% overlap (padrão para boa resolução)
nfft = 2^nextpow2(window_length * 4);  % Zero-padding para melhor resolução espectral

% ========== COLUNA ESQUERDA: SINAL AUV ==========

% --- Subplot 1: Espectro do Sinal AUV ---
ax1 = axes('Position', [left_margin, ...
                        bottom_margin + plot_height + vertical_gap, ...
                        plot_width, ...
                        plot_height]);
plot(f_fft_pos(idx_zoom), SPL_sinal(idx_zoom), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
grid on
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 11, 'FontWeight', 'bold')
title('AUV-like Signal Spectrum', 'FontSize', 13, 'FontWeight', 'bold')
xlim([0 1000])
ylim([0 100])
set(gca, 'XTickLabel', [], 'FontSize', 10)

% --- Subplot 3: Espectrograma do Sinal AUV ---
ax3 = axes('Position', [left_margin, ...
                        bottom_margin, ...
                        plot_width, ...
                        plot_height]);
spectrogram(x, hamming(window_length), overlap, nfft, fs, 'yaxis');
ylim([0 1])  % Limita até 1 kHz
% Ajuste dinâmico de contraste (últimos 60 dB)
clim_max = max(get(gca,'CLim'));
caxis([clim_max - 60, clim_max])
% Colormap vermelho intenso -> amarelo claro
colormap(ax3, hot)
colorbar('FontSize', 10)
xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Frequency [kHz]', 'FontSize', 11, 'FontWeight', 'bold')
title('AUV-like Signal Spectrogram', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

% ========== COLUNA DIREITA: RUÍDO ==========

% --- Subplot 2: Espectro do Ruído ---
ax2 = axes('Position', [left_margin + plot_width + horizontal_gap, ...
                        bottom_margin + plot_height + vertical_gap, ...
                        plot_width, ...
                        plot_height]);
plot(f_fft_pos(idx_zoom), SPL_ruido(idx_zoom), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
grid on
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 11, 'FontWeight', 'bold')
title('Background Noise Spectrum (Enseada dos Anjos)', 'FontSize', 13, 'FontWeight', 'bold')
xlim([0 1000])
ylim([0 100])
set(gca, 'XTickLabel', [], 'FontSize', 10)

% --- Subplot 4: Espectrograma do Ruído ---
ax4 = axes('Position', [left_margin + plot_width + horizontal_gap, ...
                        bottom_margin, ...
                        plot_width, ...
                        plot_height]);
spectrogram(ruido, hamming(window_length), overlap, nfft, fs, 'yaxis');
ylim([0 1])  % Limita até 1 kHz
% Ajuste dinâmico de contraste
clim_max = max(get(gca,'CLim'));
caxis([clim_max - 60, clim_max])
% Colormap vermelho intenso -> amarelo claro
colormap(ax4, hot)
colorbar('FontSize', 10)
xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Frequency [kHz]', 'FontSize', 11, 'FontWeight', 'bold')
title('Noise Spectrogram (Enseada dos Anjos)', 'FontSize', 13, 'FontWeight', 'bold')
set(gca, 'FontSize', 10)

% Link eixos X dos espectros para zoom/pan sincronizado
linkaxes([ax1, ax2], 'x');


% =========================================================================
%% 16. ESPECTROGRAMAS COM CONTRASTE AUTOMÁTICO
% =========================================================================
janela = hamming(4096);
overlap = round(0.75*length(janela));
nfft_spec = 8192;

figure('Position',[100 100 1400 900])

subplot(3, 1, 1)
spectrogram(x , janela , overlap , nfft_spec , fs , 'yaxis');
ylim([0 1]);
title('Espectrograma do sinal puro');
ylabel('Frequência em kHz');
xlabel('Tempo em s');
colorbar;

subplot(3, 1, 2)
spectrogram(ruido , janela , overlap , nfft_spec , fs , 'yaxis');
ylim([0 1]);
title('Espectrograma do ruído de fundo');
ylabel('Frequência em kHz');
xlabel('Tempo em s');
colorbar;

subplot(3, 1, 3)
spectrogram(sinal_com_ruido , janela , overlap , nfft_spec , fs , 'yaxis');
ylim([0 1]);
title('Espectrograma do sinal somado ao ruído');
ylabel('Frequência em kHz');
xlabel('Tempo em s');
colorbar;

% =========================================================================
%% 17. REPRESENTAÇÃO TEMPORAL DOS TRÊS SINAIS
% =========================================================================
t = (0:N-1)'/fs;

figure('Position',[100 100 1400 900])

subplot(3,1,1)
plot(t, x, 'b-', 'LineWidth', 0.4)
grid on
xlim([0 dur])
ylabel('Amplitude em Pa', 'FontSize', 12, 'FontWeight', 'bold')
title('Sinal puro sintetizado no domínio do tempo', 'FontSize', 14, 'FontWeight', 'bold')
set(gca, 'XTickLabel', [])

subplot(3,1,2)
plot(t, ruido, 'r-', 'LineWidth', 0.4)
grid on
xlim([0 dur])
ylabel('Amplitude em Pa', 'FontSize', 12, 'FontWeight', 'bold')
title('Ruído ambiental no domínio do tempo', 'FontSize', 14, 'FontWeight', 'bold')
set(gca, 'XTickLabel', [])

subplot(3,1,3)
plot(t, sinal_com_ruido, 'k-', 'LineWidth', 0.4)
grid on
xlim([0 dur])
xlabel('Tempo em s', 'FontSize', 12, 'FontWeight', 'bold')
ylabel('Amplitude em Pa', 'FontSize', 12, 'FontWeight', 'bold')
title('Sinal combinado para a razão sinal-ruído prescrita', 'FontSize', 14, 'FontWeight', 'bold')


