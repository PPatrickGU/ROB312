%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           TP Filtre de Kalman partie 2
%           Nicolas Merlinge (ONERA, TP ENSTA ROB312)
%           Version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear % effacer toutes les variables
close all % fermer toutes les fenetres
clc % effacer la ligne de commande
rng(123456) % imposer la graine de generation de nombres pseudo-aleatoire pour la repetabilite

% Parametres initiaux et de simulation
X_reel = [0, 0, 0]' + sqrt(diag([10, 10, 10*pi/180].^2))*randn(3,1); % etat vrai (inconnu du filtre)    
d = size(X_reel,1); % dimension de l'etat
dt = 0.1; % pas de temps
Ri = diag([10,1*pi/180].^2); % matrice de covariance du bruit de mesure pour chaque amer
amers(1).P = [10;00]; % position de l'amer 1 (point de repere)
amers(2).P = [0;50]; % position de l'amer 2 (point de repere)
amers(3).P = [10;20]; % position de l'amer 2 (point de repere)
% amers(4).P = [10;5]; % position de l'amer 2 (point de repere)
% amers(5).P = [0;20]; % position de l'amer 2 (point de repere)
% amers(6).P = [0;0]; % position de l'amer 2 (point de repere)
dm = size(amers,2)*size(Ri,1); % dimension du vecteur de mesures
R = [];
for i = 1:size(amers,2)
    R = blkdiag(R,Ri);
end

% Initialisation des variables de stockage des donnees
tk=1;
t_sim(tk) = 0;
X_reel_sim(:,tk) = X_reel;

%% Boucle de simulation physique de la trajectoire
T = 10; % duree (s)
for tk = 2:(T/dt)
    % commande (definit la trajectoire)
    V = 10; % Vitesse
    omega = cos(tk*2*pi/100); % vitesse angulaire
    
    % simulation de l'etat vrai (attention, inconnu du filtre)
    X_reel = [X_reel(1) + V*dt*cos(X_reel(3)); X_reel(2) + V*dt*sin(X_reel(3)); X_reel(3) + dt*omega]; % propagation de l'etat reel
    X_reel(3) = mod(X_reel(3), 2*pi);
    
    % enregistrement des variables (pour plot)
    X_reel_sim(:,tk) = X_reel;
end

%% Boucle de simulation physique des mesures
dt_mesure = 1;
for tk = 2:(T/dt)
%     if mod(tk, dt_mesure/dt) == 0 
    % Récupération de l'etat vrai
        X_reel = X_reel_sim(:,tk);
%     end
        % generation de la mesure reelle
    Y = [];
    for i = 1:size(amers,2)
        P_amer = amers(i).P;
        Delta = P_amer(1:2)-X_reel(1:2); % attention, Delta reel inconnu du filtre
        Yi = [norm(Delta); atan2(Delta(2),Delta(1))-X_reel(3)];
        Yi(2) = mod(Yi(2),2*pi);
        Y = [Y;Yi];
    end
    noiseY = randn(dm,1);
    Y = Y + sqrt(R)*noiseY;
    
    % enregistrement des variables (pour plot)
    Y_sim(:,tk) = Y;
    noiseY_sim(:,tk) = noiseY;
end

%% Boucle du filtre de Kalman

% Initialisation du filtre
P_hat = diag([10, 10, 10*pi/180].^2); % matrice de covariance initiale
X_hat = [0, 0, 0]'; % estimation initial
Qf = diag([0.1,0.1,0.001]); % matrice de covariance de bruit de dynamique
% Qf = 0.1 * diag([0.1,0.1,0.001]); % matrice de covariance de bruit de dynamique
% Rfi = diag([3,1*pi/180].^2); % matrice de covariance du bruit de mesure pour chaque amer
Rfi = diag([10,1*pi/180].^2); % matrice de covariance du bruit de mesure pour chaque amer
Rf = []; % matrice de covariance du bruit de mesure pour l'ensemble des amers
for i = 1:size(amers,2)
    Rf = blkdiag(Rf,Rfi);
end

% Initialisation des variables de stockage des donnees
tk=1;
K = zeros(d,dm);
K_sim(:,:,tk) = K;
inno = zeros(dm,1);
inno_sim(:,tk) = inno;
P_sim(:,:,tk) = P_hat;
Pdiag_sim(:,:,tk) = diag(P_hat);
X_hat_sim(:,tk) = X_hat;

for tk = 2:(T/dt)
    % commande (definit la trajectoire)
    V = 10; % Vitesse
    omega = cos(tk*2*pi/100); % vitesse angulaire

    % prediction (?completer: variables X_hat et P_hat)
    F = eye(3) ; % matrice de dynamique (a completer)   
    F(1,3) = -V*dt*cos(X_hat(3,1));
    F(2,3) = V*dt*cos(X_hat(3,1));
    
    X_hat = [X_hat(1) + V*dt*cos(X_hat(3)); X_hat(2) + V*dt*sin(X_hat(3)); X_hat(3) + dt*omega]; % ¨¦tat estim¨¦ pr¨¦dit
    X_hat(3) = (mod(X_hat(3), 2*pi));  % modulo sur l'angle theta (3e composante de X_hat)
    
    P_hat = F*P_hat*F.' + Qf; % covariance prédite
    
    % Récupération de la mesure reelle
    Y = Y_sim(:,tk);
    
    % validit?de la mesure reelle (gestion de la frequence)
    dt_mesure = 0.1; % (s)
    if mod(tk*dt, dt_mesure) == 0
        is_measurementValid = true;
    else
        is_measurementValid = false;
    end
    
    % correction (a completer: variables K, P_hat, inno, X_hat)
    H = [];
    if is_measurementValid
        Y_hat = [];
        for i = 1:size(amers,2)
            P_amer = amers(i).P;
            Delta_hat = P_amer(1:2) - X_hat(1:2);
            
            Y_hati = [norm(Delta_hat); atan2(Delta_hat(2),Delta_hat(1)) - X_hat(3)];
            Y_hati(2) = mod(Y_hati(2),2*pi); % modulo de l'angle phi
            Y_hat = [Y_hat; Y_hati];
            
            Hi = [-Delta_hat(1)/norm(Delta_hat), -Delta_hat(2)/norm(Delta_hat), 0;...
                   Delta_hat(2)/norm(Delta_hat)^2, -Delta_hat(1)/norm(Delta_hat)^2, -1];
            H = [H; Hi];
        end
        
        % matrice d'observation:
        K = P_hat*H.'*inv(Rf + H*P_hat*H.');
        P_hat = (eye(d) - K*H)*P_hat; 
        inno = Y - Y_hat;
%         inno(2) = (mod(inno(2), 2*pi));
        X_hat = X_hat + K*inno; 
%         X_hat(3) = (mod(X_hat(3), 2*pi));
    end
    
    % enregistrement des variables (pour plot)
    K_sim(:,:,tk) = K;
    inno_sim(:,tk) = inno;
    P_sim(:,:,tk) = P_hat;
    Pdiag_sim(:,:,tk) = diag(P_hat);
    X_hat_sim(:,tk) = X_hat;
    t_sim(tk) = tk;
    
    % plot instantan?
    figure(2)
    clf
    hold on
    plot(X_reel_sim(1,1:tk), X_reel_sim(2,1:tk), '-b')
    plot(X_hat_sim(1,1:tk), X_hat_sim(2,1:tk), '-r')
    scatter(X_reel_sim(1,tk), X_reel_sim(2,tk), 'b')
    scatter(X_hat_sim(1,tk), X_hat_sim(2,tk), '.r')
    plotcov(X_hat_sim(:,tk),3^2*P_sim(:,:,tk),'r')
    for i = 1:size(amers,2)
        P_amer = amers(i).P;
        scatter(P_amer(1), P_amer(2), 'ok')
        Delta = P_amer(1:2)-X_reel_sim(1:2,tk);
        plot([X_reel_sim(1,tk), X_reel_sim(1,tk) + Delta(1)], [X_reel_sim(2,tk), X_reel_sim(2,tk) + Delta(2)], 'g');
    end
    grid on
    axis equal
    legend('Position vraie', 'position estimee','Location','best')
    drawnow
end

figure(1)
labels = {'x (m)','y (m)','\theta (rad)'};
for i = 1:d
    subplot(3,1,i)
    hold on
    fill([t_sim flip(t_sim,2)],[X_hat_sim(i,:) - 3*sqrt(Pdiag_sim(i,:)), X_hat_sim(i,end:-1:1) + 3*sqrt(Pdiag_sim(i,end:-1:1))], [7 7 7]/8);
    plot(t_sim, X_reel_sim(i,:), 'b')
    plot(t_sim, X_hat_sim(i,:), 'r')
    grid on
    xlabel('time (s)')
    ylabel(labels(i))
end
legend('uncertainty (3\sigma)', 'actual state', 'estimated state')


