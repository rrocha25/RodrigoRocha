clear; clc; close all;

%% 1) Parametros
fs = 44100;
c = 1500;

dur_s = 80e-3;
t = (0:round(dur_s*fs)-1)'/fs;

SNR_dB = 0;

% Geometria
fmax_geom = 850;
% Define a frequência máxima considerada para a geometria do arranjo.

lambda_min = c / fmax_geom;
% Usado para obter o espaçamento máximo permitido entre hidrofones para evitar ambiguidades espaciais.
% Picos do sinal em 400 e 800Hz

L = 0.5 * lambda_min;
% Largura limite de lambida/2
H = L;

pos_h = pyramid_array_positions(L, H);

% Fonte simulada
az_real = 35;
el_real = 24;
r_real = 10;

pos_f = sph2cart_deg(r_real, az_real, el_real);
% Converte as coordenadas esféricas para coordenadas cartesianas

%% 2) maxLag fisico
dmax = max_pair_distance(pos_h);
% Calcula a maior distância entre quaisquer dois hidrofones do arranjo(função)

maxLagSamples = ceil((dmax/c)*fs) + 2;
% Define o atraso maximo fisicamente possivel entre dois hidrofones em amostras
% Limitando o lag da correlacao cruzada ao intervalo plausivel e reduzindo picos espurios.

fprintf('\nmaxLagSamples = %d amostras  %.3f ms\n', maxLagSamples, 1000*maxLagSamples/fs);

%% 3) Sinal simples com duas frequencias
f1 = 400;
f2 = 800;

A1 = 1.0;
A2 = 0.7;

phi1 = 0;
phi2 = 0;

s = A1*sin(2*pi*f1*t + phi1) + A2*sin(2*pi*f2*t + phi2);
% Cria um sinal do tipo CW com 400 e 800Hz

s = s .* hann(length(s));
% Aplica uma janela de Hann ao sinal

s = s - mean(s);
% Remove o médio do sinal, elimininando a componente de DC 
% Melhora a correlação e evita enviesamento do pico

%% Sinal chirp usado no primeiro teste

% s = chirp(t, 300, t(end), 1200, 'linear');
% % Cria um sinal do tipo chirp 300 a 1200Hz
% 
% s = s .* hann(length(s));
% % Aplica uma janela de Hann ao sinal
% 
% s = s - mean(s);
% % Remove o médio do sinal, elimininando a componente de DC 
% % Melhora a correlação e evita enviesamento do pico

%% 4) Simular captacao nos 4 hidrofones com atrasos relativos
N = length(s);
dist = zeros(4,1);
for k = 1:4
    dist(k) = norm(pos_f - pos_h(k,:));
    % Calcula a distância euclidiana entre a fonte e cada hidrofone
end
tau_rel = (dist - dist(1))/c; % s
% Constrói atrasos relativos tomando o hidrofone 1 como referência

guard = maxLagSamples + 4;
% Define uma margem de segurança em amostras para evitar problemas de borda ao aplicar atrasos

s_pad = [zeros(guard,1); s; zeros(guard,1)];
% Cria uma versão do sinal com preenchimento de zeros antes e depois

x = zeros(N,4);
for k = 1:4
    sk = apply_delay_seconds(s_pad, fs, tau_rel(k));
    % Aplica atraso temporal ao sinal
    
    seg = sk(guard+1:guard+N);
    % Recorta do sinal atrasado exatamente N amostras sem os zeros
    
    seg = add_awgn_snr(seg, SNR_dB);
    % Adiciona ruído branco gaussiano aditivo ao canal k com SNR definido

    x(:,k) = seg;
    % Armazena o sinal de cada hidrofone
end

%% Plot dos 4 sinais do arranjo apos aplicar os delays
t_ref_ms = 37.89;

figure('Position',[120 80 950 700]);

for k = 1:4
    ax = subplot(4,1,k);
    plot(1e3*t, x(:,k), 'k', 'LineWidth', 1.0); grid on; hold on
    ylabel(sprintf('$x_{%d}$',k),'Interpreter','latex');
    title(sprintf('Canal %d atraso real relativo ao canal 1: %.3f ms', k, 1e3*tau_rel(k)), 'Interpreter','latex');

    yl = ylim(ax);
    t_del_ms = t_ref_ms + 1e3*tau_rel(k);

    h_ref = plot([t_ref_ms t_ref_ms], yl, 'b--', 'LineWidth', 1.2);
    h_exp = plot([t_del_ms t_del_ms], yl, 'r-',  'LineWidth', 1.2);
    ylim(ax, yl);

    legend([h_ref h_exp], {sprintf('Pico referencia: %.2f ms', t_ref_ms), sprintf('Pico estimado: %.2f ms', t_del_ms)}, 'Location','best', 'Interpreter','latex');

    if k == 4
        xlabel('Tempo em ms','Interpreter','latex');
    else
        set(ax,'XTickLabel',[]);
    end

    hold off
end

%% 5) Estimar TDOAs com xcorr com maxLag fisico e plotar a correlacao por par
pares = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
Np = size(pares,1);

tdoa = zeros(Np,1);
tdoa_th_plane = zeros(Np,1);

u_real = [cosd(el_real)*cosd(az_real) cosd(el_real)*sind(az_real) sind(el_real)];

figure('Position',[140 120 1100 650]);

for p = 1:Np
    i = pares(p,1);
    j = pares(p,2);

    [r,lags] = xcorr(x(:,j)-mean(x(:,j)), x(:,i)-mean(x(:,i)), maxLagSamples, 'coeff');
    [rpk,kpk] = max(abs(r));
    tdoa(p) = lags(kpk)/fs;

    dp = pos_h(j,:) - pos_h(i,:);
    tdoa_th_plane(p) = -dot(dp, u_real)/c;

    subplot(3,2,p);
    plot(1e3*lags/fs, r, 'k', 'LineWidth', 1.1); grid on; hold on

    yl = ylim;
    h_hat = plot(1e3*[tdoa(p) tdoa(p)], yl, 'm-', 'LineWidth', 1.2);
    h_the = plot(1e3*[tdoa_th_plane(p) tdoa_th_plane(p)], yl, 'b--', 'LineWidth', 1.1);
    h_pk  = plot(1e3*tdoa(p), r(kpk), 'mo', 'MarkerSize', 6, 'LineWidth', 1.2);
    ylim(yl);

    xlim([-0.6 0.2]);
    xlabel('Atraso em ms','Interpreter','latex');
    ylabel(sprintf('$\\rho_{%d%d}$', j, i),'Interpreter','latex');
    title(sprintf('Par %d %d', i, j), 'Interpreter','latex');

    legend([h_hat h_the h_pk], {sprintf('$\\hat{\\tau}$ %.3f ms',1e3*tdoa(p)), sprintf('$\\tau$ teorico %.3f ms',1e3*tdoa_th_plane(p)), sprintf('$|\\rho|_{max}$ %.3f',rpk)}, 'Location','best','Interpreter','latex');

    hold off
end

annotation('textbox',[0 0.97 1 0.03],'String','Correlacao cruzada normalizada por par de hidrofone modelo de onda plana','Interpreter','latex','EdgeColor','none','HorizontalAlignment','center');

fprintf('\nTDOAs medidos em ms:\n'); disp(1e3*tdoa);
fprintf('\nErro estimacao menos teorico em ms:\n'); disp(1e3*(tdoa - tdoa_th_plane));


%% 6) Estimar DOA por varredura em 360 graus
az_vec = 0:1:359;
el_vec = 0:1:90;

J = doa_cost_grid(pos_h, pares, tdoa, c, az_vec, el_vec);
[az_hat, el_hat] = peak_min_cost(J, az_vec, el_vec);

az_real_w = mod(az_real, 360);
az_hat_w = mod(az_hat, 360);

fprintf('\n=== DOA ===\n');
fprintf('Real az = %.2f el = %.2f\n', az_real_w, el_real);
fprintf('Est  az = %.2f el = %.2f\n', az_hat_w, el_hat);

%% 7) Plot do mapa de custo em 4 setores de 90 graus
setores = [0 90; 90 180; 180 270; 270 360];
clim = [min(J(:)) max(J(:))];

figure('Position',[100 100 1100 650]);

for s = 1:4
    a0 = setores(s,1);
    a1 = setores(s,2);

    idx = find(az_vec >= a0 & az_vec < a1);
    az_s = az_vec(idx);
    J_s = J(:, idx);

    ax = subplot(2,2,s);
    imagesc(az_s, el_vec, J_s);
    set(ax,'YDir','normal');
    set(ax,'CLim',clim);
    grid on
    hold on

    plot(az_real_w, el_real, 'k+', 'MarkerSize', 12, 'LineWidth', 2);
    plot(az_hat_w,  el_hat,  'wo', 'MarkerSize', 8,  'LineWidth', 2);

    xlim([a0 a1]);
    ylim([min(el_vec) max(el_vec)]);

    xlabel('Azimute em graus','Interpreter','latex');
    ylabel('Elevacao em graus','Interpreter','latex');
    title(sprintf('Matriz do erro quadrático total: setor %.0f a %.0f', a0, a1));
    
    hE = plot(az_hat_w,  el_hat,  'o', 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 1 1]);
    hR = plot(az_real_w, el_real, '+', 'Color',[0.85 0 0], 'MarkerSize', 14, 'LineWidth', 3);

    lg = legend([hR hE], {sprintf('Real az %.2f el %.2f', az_real_w, el_real), sprintf('Estimado az %.2f el %.2f', az_hat_w, el_hat)}, 'Location','best', 'Interpreter','latex');
    set(lg,'Color','w','EdgeColor','k');
    hold off
end

cmap = flipud(jet(256));
colormap(gcf, cmap);
cb = colorbar;
cb.Label.String = 'Erro quadratico';
cb.Label.Interpreter = 'latex';

%% Figura 2 mapa de custo completo 0 a 360 em escala log
az_vec_360 = [az_vec 360];
J_360 = [J J(:,1)];

Jlog_360 = log10(J_360 + eps);
clim_log = [min(Jlog_360(:)) max(Jlog_360(:))];

figure('Position',[120 120 980 520]);
imagesc(az_vec_360, el_vec, Jlog_360);
set(gca,'YDir','normal');
set(gca,'CLim',clim_log);
colormap(cmap);
cb = colorbar;
cb.Label.String = 'Erro quadratico';
cb.Label.Interpreter = 'latex';grid on
hold on

plot(az_real_w, el_real, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
plot(az_hat_w,  el_hat,  'wo', 'MarkerSize', 8,  'LineWidth', 2);

xlim([0 360]);
ylim([min(el_vec) max(el_vec)]);

xlabel('Azimute em graus','Interpreter','latex');
ylabel('Elevacao em graus','Interpreter','latex');
title('Matriz logaritima do erro quadratico total J de 0 a 360');

hE = plot(az_hat_w,  el_hat,  'o', 'MarkerSize', 10, 'LineWidth', 2.5, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 1 1]);
hR = plot(az_real_w, el_real, '+', 'Color',[0.85 0 0], 'MarkerSize', 14, 'LineWidth', 3);

lg = legend([hR hE], {sprintf('Real az %.2f el %.2f', az_real_w, el_real), sprintf('Estimado az %.2f el %.2f', az_hat_w, el_hat)}, 'Location','best', 'Interpreter','latex');
set(lg,'Color','w','EdgeColor','k');
hold off



%% ========================================================================
%  Funcoes

function pos = pyramid_array_positions(L, H)
% Retorna as coordenadas 3D de um arranjo piramidal com 4 hidrofones, dado o lado L e a altura H.

    pos = zeros(4,3);
% Prealoca a matriz 4x3 que armazena as posicoes cartesianas de cada hidrofone.

    pos(1,:) = [0, 0, 0];
% Define o hidrofone 1 na origem do sistema de coordenadas.

    pos(2,:) = [L, 0, 0];
% Define o hidrofone 2 no eixo x, a uma distancia L do hidrofone 1.

    pos(3,:) = [L/2, L*sqrt(3)/2, 0];
% Define o hidrofone 3 no plano z igual a zero, formando triangulo equilatero com lado L.

    cx = L/2; cy = L*sqrt(3)/6;
% Calcula as coordenadas do centroide do triangulo equilatero da base.

    pos(4,:) = [cx, cy, H];
% Define o hidrofone 4 no apice da piramide, acima do centroide, com altura H.

end

function pos = sph2cart_deg(r, az_deg, el_deg)
% Converte coordenadas esfericas em graus para coordenadas cartesianas 3D.

    pos = [r*cosd(el_deg)*cosd(az_deg), r*cosd(el_deg)*sind(az_deg), r*sind(el_deg)];
% Calcula x, y e z a partir do raio, azimute e elevacao, usando funcoes trigonometricas em graus.

end

function dmax = max_pair_distance(P)
% Calcula a maior distancia entre quaisquer dois pontos de um conjunto de posicoes 3D.

    dmax = 0;
% Inicializa a maior distancia com zero para permitir atualizacao por maximo.

    for i = 1:size(P,1)
% Percorre cada ponto como primeiro elemento do par.

        for j = i+1:size(P,1)
% Percorre os pontos seguintes para formar pares unicos sem repeticao.

            dmax = max(dmax, norm(P(i,:)-P(j,:)));
% Atualiza dmax com a maior norma euclidiana encontrada entre os pares avaliados.

        end
% Finaliza o loop interno sobre o segundo indice do par.

    end
% Finaliza o loop externo sobre o primeiro indice do par.

end

function y = apply_delay_seconds(x, fs, tau)
% Aplica um atraso temporal tau em segundos a um sinal discreto por interpolacao linear.

    N = length(x);
% Define o numero de amostras do sinal de entrada.

    n = (0:N-1)';
% Cria o vetor de indices de amostra em coluna para uso na interpolacao.

    y = interp1(n, x, n - tau*fs, 'linear', 0);
% Interpola o sinal em indices deslocados por tau vezes fs, preenchendo fora do suporte com zero.

end

function y = add_awgn_snr(x, SNR_dB)
% Adiciona ruido branco gaussiano ao sinal para atingir uma SNR especificada em dB.

    px = mean(x.^2);
% Estima a potencia media do sinal como a media do quadrado das amostras.

    pn = px / (10^(SNR_dB/10));
% Calcula a potencia de ruido necessaria a partir da SNR desejada.

    y = x + sqrt(pn)*randn(size(x));
% Soma ruido gaussiano de variancia pn ao sinal de entrada.

end

function J = doa_cost_grid(pos_h, pares, tdoa_med, c, az_vec, el_vec)
% Calcula o mapa de custo J em uma grade de azimute e elevacao a partir de TDOAs medidos.

    Np = size(pares,1);
% Armazena o numero de pares de hidrofones considerados na estimativa.

    J = zeros(length(el_vec), length(az_vec));
% Prealoca a matriz de custo com linhas para elevacao e colunas para azimute.

    for iaz = 1:length(az_vec)
% Varre todos os valores de azimute da grade.

        az = az_vec(iaz);
% Seleciona o azimute atual da varredura.

        for iel = 1:length(el_vec)
% Varre todos os valores de elevacao da grade.

            el = el_vec(iel);
% Seleciona a elevacao atual da varredura.

            u = [cosd(el)*cosd(az), cosd(el)*sind(az), sind(el)];
% Calcula o vetor unitario de direcao no sistema cartesiano para o par az, el.

            tdoa_th = zeros(Np,1);
% Prealoca o vetor de TDOAs teoricos para todos os pares na direcao considerada.

            for p = 1:Np
% Percorre cada par de hidrofones para calcular o TDOA teorico correspondente.

                i = pares(p,1); j = pares(p,2);
% Extrai os indices dos dois hidrofones que compoem o par p.

                dp = pos_h(j,:) - pos_h(i,:);
% Calcula o vetor separacao entre os hidrofones do par, do sensor i para o sensor j.

                tdoa_th(p) = -dot(dp, u)/c;
% Calcula o TDOA teorico por projecao da linha de base na direcao u e divisao por c.

            end
% Finaliza o loop sobre pares para os TDOAs teoricos.

            e = tdoa_med - tdoa_th;
% Calcula o vetor de erros entre TDOAs medidos e teoricos para a direcao atual.

            J(iel, iaz) = sum(e.^2);
% Define o custo como soma do erro quadraticos sobre todos os pares.

        end
    end
end

function [az_hat, el_hat] = peak_min_cost(J, az_vec, el_vec)
% Encontra o par azimute e elevacao que minimiza o mapa de custo J.

    [~, idx] = min(J(:));
% Localiza o indice linear do menor valor de custo ao vetorizar a matriz J.

    [iel, iaz] = ind2sub(size(J), idx);
% Converte o indice linear para indices de linha e coluna correspondentes em J.

    az_hat = az_vec(iaz);
% Mapeia a coluna do minimo para o valor de azimute na grade.

    el_hat = el_vec(iel);
% Mapeia a linha do minimo para o valor de elevacao na grade.

end

