    % =========================================================================
%% SÍNTESE DE ÁUDIO A PARTIR DO ESPECTRO DO CSV (NÍVEL DE SPL EXATO)
% =========================================================================
clear; clc; close all;

% Output directory for figures
outputDir = './figures';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
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


%% 15 - Espectro + Espectrograma do Sinal AUV e Ruído (figuras separadas 2:1)

% Parâmetros para espectrograma
window_length = round(0.025 * fs);      % 25 ms
overlap       = round(0.75 * window_length);  % 75% overlap
nfft          = 2^nextpow2(window_length * 4);

% Índice de zoom para espectro (até 1 kHz)
idx_zoom = f_fft_pos <= 1000;

% Dimensão comum 2:1 (largura ~ 2 * altura)
fig_width  = 1000;
fig_height = 500;

% (a) Espectro do Sinal AUV
fig = figure('Position', [100 100 1200 600], 'Color', 'w');

plot(f_fft_pos(idx_zoom), SPL_sinal(idx_zoom), ...
     'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
grid on
xlabel('Frequency [Hz]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 20, 'FontWeight', 'bold');
xlim([0 1000]);
ylim([0 100]);
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.18 0.88 0.77]);

print(fullfile(outputDir, 'auv_spectrum.eps'), '-depsc');
%close(fig);
%%
% (b) Espectro do Ruído
fig = figure('Position', [100 100 1200 600], 'Color', 'w');

plot(f_fft_pos(idx_zoom), SPL_ruido(idx_zoom), ...
     'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
grid on
xlabel('Frequency [Hz]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 20, 'FontWeight', 'bold');
xlim([0 1000]);
ylim([0 100]);
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.18 0.88 0.77]);

print(fullfile(outputDir, 'noise_spectrum.eps'), '-depsc');
%close(fig);
%%
% (c) Espectrograma do Sinal AUV
fig = figure('Position', [100 100 1200 600], 'Color', 'w');

spectrogram(x, hamming(window_length), overlap, nfft, fs, 'yaxis');
ylim([0 1]);  % até 1 kHz
clim_max = max(get(gca,'CLim'));
caxis([clim_max - 60, clim_max]);
colormap(hot);

cb = colorbar;
ylabel(cb, 'Power/Frequency [dB/Hz]', 'FontSize', 20, 'FontWeight', 'bold');

xlabel('Time [s]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.15 0.78 0.82]);

print(fullfile(outputDir, 'auv_spectrogram.eps'), '-depsc');
%close(fig);
%%
% (d) Espectrograma do Ruído
fig = figure('Position', [100 100 1200 600], 'Color', 'w');

spectrogram(ruido, hamming(window_length), overlap, nfft, fs, 'yaxis');
ylim([0 1]);
clim_max = max(get(gca,'CLim'));
caxis([clim_max - 60, clim_max]);
colormap(hot);

cb = colorbar;
ylabel(cb, 'Power/Frequency [dB/Hz]', 'FontSize', 20, 'FontWeight', 'bold');

xlabel('Time [s]', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.15 0.78 0.82]);

print(fullfile(outputDir, 'noise_spectrogram.eps'), '-depsc');
%close(fig);

fprintf('Spectrum and spectrogram figures (2:1) saved to: %s\n', outputDir);

%% 3) Simple 2-tone signal
f1   = 400;
f2   = 800;
A1   = 1.0;
A2   = 0.7;
phi1 = 0;
phi2 = 0;
t = (0:round(dur*fs)-1)'/fs;

s = A1*sin(2*pi*f1*t + phi1) + A2*sin(2*pi*f2*t + phi2);

%% --- Spectrum plot (Sec. 3) ---
N     = length(s);
S_fft = fft(s);
S_fft = S_fft(1:floor(N/2)+1);
f_fft = (0:floor(N/2))' * fs / N;

% Calcula PSD e SPL antes do ajuste
PSD_sinal          = (abs(S_fft).^2) / N;
PSD_sinal(2:end-1) = 2 * PSD_sinal(2:end-1);
SPL_sinal          = 10 * log10(PSD_sinal / (p_ref^2) + eps);

% Ajusta ganho para o PICO ESPECTRAL ficar em SPL_pico_desejado
SPL_pico_desejado = 95;                          % dB re 1 µPa
SPL_pico_atual    = max(SPL_sinal);
ganho_dB          = SPL_pico_desejado - SPL_pico_atual;
ganho_lin         = 10^(ganho_dB / 20);

% Aplica ganho no domínio do tempo e recalcula espectro
s                  = s * ganho_lin;
S_fft              = S_fft * ganho_lin;          % evita refazer FFT
PSD_sinal          = (abs(S_fft).^2) / N;
PSD_sinal(2:end-1) = 2 * PSD_sinal(2:end-1);
SPL_sinal          = 10 * log10(PSD_sinal / (p_ref^2) + eps);

fprintf('Pico espectral apos ajuste: %.1f dB re 1 uPa\n', max(SPL_sinal));

% Plot
fig      = figure('Position', [100 100 1200 600], 'Color', 'w');
idx_zoom = f_fft <= 1000;

plot(f_fft(idx_zoom), SPL_sinal(idx_zoom), 'b-', 'LineWidth', 1.5)
grid on
xlabel('Frequency [Hz]',      'FontSize', 20, 'FontWeight', 'bold')
ylabel('SPL [dB re 1 \muPa]', 'FontSize', 20, 'FontWeight', 'bold')
xlim([0 1000])
ylim([0 100])
set(gca, 'FontSize', 20)
set(gca, 'Position', [0.08 0.18 0.88 0.77])
print(fullfile(outputDir, 'signal_spectrum.eps'), '-depsc')
% close(fig)


%% --- Spectrogram (Sec. 3) ---

fig = figure('Position', [100 100 1200 600], 'Color', 'w');

spectrogram(s, hamming(window_length), overlap, nfft, fs, 'yaxis');
ylim([0 1]);

clim_max = max(get(gca, 'CLim'));
caxis([clim_max - 60, clim_max]);
colormap(hot);

cb = colorbar;
ylabel(cb, 'Power/Frequency [dB/Hz]', 'FontSize', 20, 'FontWeight', 'bold');

xlabel('Time [s]',        'FontSize', 20, 'FontWeight', 'bold');
ylabel('Frequency [kHz]', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 20);
set(gca, 'Position', [0.08 0.15 0.78 0.82]);  % same as Sec. 15

print(fullfile(outputDir, 'signal_spectrogram.eps'), '-depsc');
% close(fig);
% =========================================================================