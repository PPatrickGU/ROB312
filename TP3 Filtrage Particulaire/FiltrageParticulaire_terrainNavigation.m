%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           TP Filtre Particulaire
%           Nicolas Merlinge (ONERA, TP ENSTA ROB312)
%           Version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear % effacer toutes les variables
close all % fermer toutes les fenetres
clc % effacer la ligne de commande
rng(123457) % imposer la graine de generation de nombres pseudo-aleatoire pour la reetabilite


% Parametres initiaux et de simulation
erreurInitiale = sqrt(diag([5000, 5000, 100, 20*pi/180].^2))*randn(4,1);
X_reel = [230000, 90000, 1000, 150*pi/180]' + erreurInitiale; % etat vrai initial (inconnu du filtre)    
d = size(X_reel,1); % dimension de l'etat
dt = 1; % pas de temps
R = 20.^2; % matrice de covariance du bruit de mesure reelle (inconnu du filtre)
dm = size(R,1); % dimension du vecteur de mesures

% Chargement des donnees
load carte.mat
params.pasx_reel = pasx_reel;
params.pasy_reel = pasy_reel;
params.nrow_h = size(h_MNT,1);
params.dxh_MNT = diff(h_MNT,1,2)/pasx_reel;
params.dyh_MNT = diff(h_MNT)/pasy_reel;
params.h_MNT = h_MNT;
params.x_MNT = x_MNT;
params.y_MNT = y_MNT;

% Initialisation des variables de stockage des donnees
tk=1;
t_sim(tk) = 0;
Y_sim(:,tk) = 0;
X_reel_sim(:,tk) = X_reel;

%% Boucle de simulation physique de la trajectoire
T = 80; % duree (s)
for tk = 2:(T/dt)
    % commande (definit la trajectoire)
    V = 300; % Vitesse (connue)
    omega = -0.01; % vitesse angulaire (rad/s) (connue)
    
    % simulation de l'etat vrai (attention, inconnu du filtre)
    t = dt*(tk-1); % temps courant
    X_reel = [X_reel(1) + V*dt*cos(X_reel(4)); X_reel(2) + V*dt*sin(X_reel(4)); X_reel(3); X_reel(4) + dt*omega]; % propagation de l'etat reel (a completer)
    X_reel(4) = mod(X_reel(4), 2*pi); % modulo 2pi sur le cap

    X_reel_sim(:,tk) = X_reel;
    t_sim(tk) = t;
end

%% Boucle de simulation physique des mesures
T = 80; % duree (s)
for tk = 2:(T/dt)
    % Reup�ration de l'etat
    X_reel = X_reel_sim(:,tk);
    
    % generation de la mesure reelle    
    Y = X_reel(3,:) - lectureCarte(X_reel, params) + sqrt(R)*randn(dm,1);
    
    Y_sim(:,tk) = Y;
end

%% Boucle du filtre particulaire

% Parametres du filtre
N = 3000; % nombre de particules
P_hat = diag([5000, 5000, 100, 20*pi/180].^2); % matrice de covariance initiale
X_hat = [230000, 90000, 1000, 150*pi/180]'; % estime initial (x, y, z, theta)
Qf = diag([3, 3, 0.6, 0.001*180/pi].^2); % matrice de covariance de bruit de dynamique
Rf = 20.^2; % covariance du bruit de mesure du filtre
threshold_resampling = 1; % seuil de reehantillonnage (theta_eff)

% Initialisation des variables de stockage des donnees
tk=1;
P_sim(:,:,tk) = P_hat;
Pdiag_sim(:,tk) = diag(P_hat);
X_hat_sim(:,tk) = X_hat;

% Initialisation du filtre
Xp = X_hat*ones(1,N) + sqrt(P_hat)*randn(d,N); % Tirage des particules autours de X_hat initial (a completer)
wp = 1/N*ones(1,N); % poids initiaux (a completer)

T = 80; % duree (s)
for tk = 2:(T/dt)
    % commande (definit la trajectoire)
    V = 300; % Vitesse (connue)
    omega = -0.01; % vitesse angulaire (rad/s) (connue)
    
    % Reup�ration de la mesure reelle
    Y = Y_sim(:,tk);
    
    % prediction (a completer: variables Xp, X_hat et P_hat)
    Xp = [Xp(1, :) + V*dt*cos(Xp(4,:)); Xp(2,:) + V*dt*sin(Xp(4,:)); Xp(3,:); Xp(4,:) + dt*omega]; % particules predites
    Xp = Xp + sqrtm(Qf)*randn(d,N);
    
    X_hat =  Xp*wp';  % etat estim predit
    X_hat(4,:) = mod(X_hat(4,:), 2*pi); % modulo 2pi sur le cap
    P_hat = ((Xp-X_hat(:,ones(1,N))).*wp(ones(d,1),:))*(Xp-X_hat(:,ones(1,N)))'; % matrice de covariance predite
    
    
    % validite de la mesure reelle (a completer pour la gestion des fr�quences et des trous de mesures)
    is_measurementValid = true;
%    Q7
%     is_measurementValid = true;
%     if tk > 50/dt && tk < 75/dt
%         is_measurementValid = false;
%     end

% Q8
%     is_measurementValid = false;
%     if mod(tk, 10) == 0
%         is_measurementValid = true;
%     end
     
    % correction (a completer)
    if is_measurementValid
        % definition de la mesure predite (a completer)%         
        
        Y_hat = Xp(3,:) - lectureCarte(Xp, params) + sqrt(Rf) * randn(1,N);
%         Y_hat = Xp(3,:) - lectureCarte(Xp, params);
        
        % correction des poids des particules (A completer)
        inno = Y - Y_hat;
        likelihood = exp(-inno.^2/(2*Rf));
        wp = wp.*likelihood; % correction des poids
        wp = wp/sum(wp); % normalisation des poids (la somme doit etre �gale ?1)
        
    end
    
    % Reehantillonnage (critere de seuil a completer, puis coder un
    % autre algorithme de reehantillonnage de votre choix pour la
    % derniere question, en substitution du fichier select.p)
    criterionResampling = 1/(sum(wp.^2)); % a completer
    if criterionResampling < N*threshold_resampling
%         Xp = Xp(:,select(wp)); % selection des nouvelles particules selon l'algorithme de reehantillonnage multinomial
        Xp = Xp(:,stratified(wp)); % selection des nouvelles particules selon l'algorithme de reehantillonnage multinomial
        wp = 1/N*ones(1,N); % reinitialisation des poids (a completer)
    end
    
    % enregistrement des variables (pour plot)
    P_sim(:,:,tk) = P_hat;
    Pdiag_sim(:,tk) = diag(P_hat);
    X_hat_sim(:,tk) = X_hat;
    
    % plot instantanee(ne pas hesiter ?passer <is_temporalPlot> ?false pour gagner du temps d'exeution)
    is_temporalPlot = true;
    if is_temporalPlot
        figure(2)
        clf
        hold on
        imagesc(params.x_MNT*params.pasx_reel/1000, params.y_MNT*params.pasx_reel/1000, h_MNT)
        xlabel('km'); ylabel('km'); 
        title(['Erreur position: ', num2str(norm(X_hat(1:3) - X_reel(1:3))), ' m'])
        grid
        colorbar
        hold on
        plot(X_reel_sim(1,1:tk)./1000, X_reel_sim(2,1:tk)./1000,'.k')
        plot(X_hat_sim(1,1:tk)./1000, X_hat_sim(2,1:tk)/1000,'.r')
        scatter(Xp(1,:)./1000, Xp(2,:)./1000, '.y')
        scatter(X_reel_sim(1,tk)./1000, X_reel_sim(2,tk)./1000, '.k')
        scatter(X_hat_sim(1,tk)./1000, X_hat_sim(2,tk)./1000, '.r')
        grid on
        ylim([50, 150])
        xlim([170, 250])
        legend('Position vraie', 'position estimee','particules')
        drawnow
        


    end

end

% Plot des resultats
figure(1)
labels = {'x (m)','y (m)','z (m)','\theta (rad)'};
for i = 1:d
    subplot(4,1,i)
    hold on
    fill([t_sim flip(t_sim,2)],[X_hat_sim(i,:) - 3*sqrt(Pdiag_sim(i,:)), X_hat_sim(i,end:-1:1) + 3*sqrt(Pdiag_sim(i,end:-1:1))], [7 7 7]/8);
    plot(t_sim, X_reel_sim(i,:), 'b')
    plot(t_sim, X_hat_sim(i,:), 'r')
    grid on
    xlabel('time (s)')
    ylabel(labels(i))
end
legend('uncertainty (3\sigma)', 'actual state', 'estimated state')

% figure(3)
% histogram(wp)
