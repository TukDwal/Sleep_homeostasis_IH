function [error,P] = fit_S_dual_stages(params, hypno, t, obs_t, obs_P)

if nargin<=3
   obs_t = [1];
   obs_P = [0];
end
% Version of the model where decay and build-up is constant, but with a
% stage-dependent factor 
% H(t+1) = f_rise * [UA - (UA - H(t)) * exp(-dt / tau_r)] + 
%          f_decay * [(H(t) - LA) * exp(-dt / tau_d) + LA]
% where f_rise and f_decay are scalars that take up different values in
% every sleep stage, they are equal to 0 at lowest and the sum of both
% should be 1 (so only f_rise can be defined, other is 1-f_rise)

S0 = params(1);
tau_d_nrem = params(2);
tau_r_nrem = params(3);
tau_d_rem = params(4);
tau_r_rem = params(5);
tau_d_w = params(6);
tau_r_w = params(7);
LA = params(8);
UA = params(9);

dt = t(2)-t(1);
% Initialisation des variables
P = zeros(length(hypno), 1); % Vecteur pour stocker la pression de sommeil modélisée
P(1) = S0; % Pression initiale

% Boucle temporelle pour simuler la pression de sommeil à chaque transition
for i = 2:length(P)
    
    if hypno(i) == 1  % Si on est en veille
        % Croissance exponentielle de la pression de sommeil
        P(i) =  P(i-1) - ((P(i-1) - (UA - (UA - P(i-1)) * exp(-dt / tau_r_w)))+...
            (P(i-1) - ((P(i-1) - LA) * exp(-dt / tau_d_w) + LA)));
        
    elseif hypno(i) == 2
        
        % Croissance exponentielle de la pression de sommeil
        P(i) = P(i-1) - ((P(i-1) - (UA - (UA - P(i-1)) * exp(-dt / tau_r_rem)))+...
            (P(i-1) - ((P(i-1) - LA) * exp(-dt / tau_d_rem) + LA)));
        
    elseif hypno(i) > 2
        
        % Croissance exponentielle de la pression de sommeil
        P(i) = P(i-1) - ((P(i-1) - (UA - (UA - P(i-1)) * exp(-dt / tau_r_nrem)))+...
            (P(i-1) - ((P(i-1) - LA) * exp(-dt / tau_d_nrem) + LA)));
        
    end
end

% Calcul de l'erreur (somme des carrés des différences entre les observations et le modèle)
weights = obs_P;

[~,idmin] = min(abs(t-obs_t));
error = sqrt( sum( ((P(idmin)' - obs_P).^2).*weights,'omitnan')/sum(weights));
end