% Number of cells
L = 7;
% Normalized cell radius
Radius_cell = 1;
% Pathloss coefficient
beta_loss = 3.7;
% Minimum radius where UEs are distributed
R_min = 2/3*Radius_cell;
% Maximum radius where UEs are distributed
R_max = 2/3*Radius_cell;
% Rician factor
kappa_rice = 4;
% Number of users per cell
numOfUEs = 10;
% Reset the rand generator
rand('state',0); randn('state',0); %#ok<RAND>

for ii=1:L    
    % UEs positions
    user_distance{ii} = sqrt((R_max^2-R_min^2)*rand(numOfUEs,1)+R_min^2); %#ok<*SAGROW>
    % Rician factor
    kappa{ii} = kappa_rice*ones(numOfUEs,1);
    % Angle of arrivals
    AoA{ii} = -pi + 2*pi*rand(numOfUEs,1);
end
% Training SNR in dB
rho_training = 6;
% Training SNR in linear
rho_training_lin = 10^(rho_training/10);
% SNR in dB
rho = 10;
% SNR in linear
rho_lin = 10^(rho/10);
% Transmit power
Power = 1;
% Noise variance in linear
noise_lin = Power*rho_lin;


