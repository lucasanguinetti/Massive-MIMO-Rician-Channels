% This Matlab script can be used to reproduce the spectral efficiecny of a
% multicell Massive MIMO system with uncorellated Rician fading channels as
% a function of the Rician factor. In particular, it can be used to
% regenerate Fig. 5 of the following work, which is currently under review 
% for publication in TCOM:
%
% L. Sanguinetti, A. Kammoun, M. Debbah "Theoretical Perfomance Limits of 
% Massive MIMO with Uncorrelated Rician Fading Channels" under review.
%
%
% This is version 1.0 (Last edited: 2018-5-22)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% work as described above.


%Empty workspace and close figures
clear; close all
% Add path of functions
addpath ./functions_DL

% Number of iterations
num_iter = 10000;
% parpool(40);

% Load system parameters
define

% Number of antennas per BSs;
numOfAnt = 100;

% Scaling coefficient for RZF
lambda = numOfUEs/(numOfAnt*rho_lin);

% Values of Rician factors to simulate
kappa_tab = [0:2:20,40,60,80,100];

% Loop over the values of Rician factors
for index_kappa=1:length(kappa_tab)
    
    % Load the current value of Rician factor
    kappa_rice = kappa_tab(index_kappa);
    
    % Overwrite the kappa values given in define.m
    for ii=1:L
        kappa{ii} = kappa_rice*ones(numOfUEs,1);
    end
    
    % Generate coordinates of BS of interest
    coor0 = [0,0]; 
    % Generate coordinates of interfering BSs 
    coor1 = 2*[cos(pi/6) sin(pi/6)]; coor2 = 2*[cos(pi/6) -sin(pi/6)]; coor3 = 2*[0,-1];
    coor4 = 2*[-cos(pi/6),-sin(pi/6)]; coor5 = 2*[-cos(pi/6) sin(pi/6)]; coor6 = [0,2];
    % Save coordinates of BSs
    coor_cells={coor0,coor1,coor2,coor3,coor4,coor5,coor6};
    
    % Prepare the steering vector of ULA
    ULA = (0:numOfAnt-1)'*pi;
    
    % Loop over the cells
    for ii=1:L
        
        % Generate positions of UEs
        coor_users{ii}=[coor_cells{ii}(1)*ones(numOfUEs,1) coor_cells{ii}(2)*ones(numOfUEs,1)]+[user_distance{ii}.*cos(AoA{ii}) user_distance{ii}.*sin(AoA{ii})];
        
        % Compute distances between the generated users and all the other cells
        % distance{ii}(:,jj) accounts for the between UEs in cell ii and cell jj
        
        somme_pathloss{ii}=0;
        for jj=1:L
            distance{ii}(:,jj)= sqrt(sum([ coor_users{ii}-ones(numOfUEs,1)*[coor_cells{jj}(1),coor_cells{jj}(2)]].^2,2));
            pathloss{ii}(:,jj)= 1./(distance{ii}(:,jj).^beta_loss);
            if (ii==jj)
                pathloss{ii}(:,jj)=pathloss{ii}(:,jj)./(1+kappa{ii});
                A{ii}=exp(-1i*ULA*sin(AoA{ii}'));
                H_bar{ii}=(ones(numOfAnt,1)*(sqrt(pathloss{ii}(:,ii).*kappa{ii}))').*A{ii};
            end
            somme_pathloss{ii}=somme_pathloss{ii}+ pathloss{ii}(:,jj);
        end
        for jj=1:L
            Phi{ii}(:,jj)=pathloss{ii}(:,ii).*pathloss{ii}(:,jj)./(1/rho_training_lin+somme_pathloss{ii});
        end
    end
    
    % Loop over the channel realizations
    parfor iter=1:num_iter
        
        % Generate results for RZF and MRT
        [RZF_results,MRT_results] = system_random(H_bar,Phi,pathloss,somme_pathloss,lambda,Power,rho_training_lin);
        
        % Save RZF results
        energy_signal_cell{iter} = RZF_results.energy_signal_cell;
        pilot_contamination{iter} = RZF_results.pilot_contamination;
        intra_cell_inter_cell_interference{iter} = RZF_results.intra_cell_inter_cell_interference;
        first_term_interference_emp{iter} = RZF_results.first_term_interference_emp;
        second_term_interference_emp{iter}=RZF_results.second_term_interference_emp;
        third_term_interference_emp{iter}= RZF_results.third_term_interference_emp;
        
        % Save MRT results
        energy_signal_cell_MRT{iter}=MRT_results.energy_signal_cell;
        pilot_contamination_MRT{iter}=MRT_results.pilot_contamination;
        intra_cell_inter_cell_interference_MRT{iter}=MRT_results.intra_cell_inter_cell_interference ;
        intra_cell_interference_MRT{iter}=MRT_results.intra_cell_interference_MRT;
        
    end
    
    for l=1:L
        
        energy_signal_cell_moy{l}=0;
        pilot_contamination_moy{l}=0;
        intra_cell_inter_cell_interference_moy{l}=0;
        first_term_interference_emp_moy{l}=0;
        second_term_interference_emp_moy{l}=0;
        third_term_interference_emp_moy{l}=0;
        
        energy_signal_cell_MRT_moy{l}=0;
        pilot_contamination_MRT_moy{l}=0;
        intra_cell_inter_cell_interference_MRT_moy{l}=0;
        intra_cell_interference_MRT_moy{l} = 0;
        
    end
    
    for l=1:L
        for iter=1:num_iter
            
            energy_signal_cell_moy{l}=energy_signal_cell_moy{l}+ energy_signal_cell{iter}{l};
            pilot_contamination_moy{l}= pilot_contamination_moy{l} +pilot_contamination{iter}{l};
            intra_cell_inter_cell_interference_moy{l}= intra_cell_inter_cell_interference_moy{l}+ intra_cell_inter_cell_interference{iter}{l};
            first_term_interference_emp_moy{l}=first_term_interference_emp{iter}{l}+  first_term_interference_emp_moy{l};
            second_term_interference_emp_moy{l}= second_term_interference_emp_moy{l} + second_term_interference_emp{iter}{l};
            third_term_interference_emp_moy{l}= third_term_interference_emp_moy{l} + third_term_interference_emp{iter}{l};
            
            energy_signal_cell_MRT_moy{l}= energy_signal_cell_MRT_moy{l}+ energy_signal_cell_MRT{iter}{l};
            pilot_contamination_MRT_moy{l}= pilot_contamination_MRT_moy{l}+ pilot_contamination_MRT{iter}{l};
            intra_cell_inter_cell_interference_MRT_moy{l}=intra_cell_inter_cell_interference_MRT_moy{l} + intra_cell_inter_cell_interference_MRT{iter}{l};
            intra_cell_interference_MRT_moy{l} = intra_cell_interference_MRT_moy{l}+ intra_cell_interference_MRT{iter}{l};
            
        end
    end
    
    % RZF
    energy_signal_cell_moy = cellfun(@(x) x./(num_iter),energy_signal_cell_moy,'UniformOutput',false);
    pilot_contamination_moy = cellfun(@(x) x./(num_iter),pilot_contamination_moy,'UniformOutput',false);
    intra_cell_inter_cell_interference_moy = cellfun(@(x) x./(num_iter),intra_cell_inter_cell_interference_moy,'UniformOutput',false);
    
    % MRT
    energy_signal_cell_MRT_moy = cellfun(@(x) x./(num_iter),energy_signal_cell_MRT_moy,'UniformOutput',false);
    pilot_contamination_MRT_moy = cellfun(@(x) x./(num_iter),pilot_contamination_MRT_moy,'UniformOutput',false);
    intra_cell_inter_cell_interference_MRT_moy = cellfun(@(x) x./(num_iter),intra_cell_inter_cell_interference_MRT_moy,'UniformOutput',false);
    intra_cell_interference_MRT_moy = cellfun(@(x) x./(num_iter),intra_cell_interference_MRT_moy,'UniformOutput',false);
    
    for ii=1:L
        total_interference{ii}=intra_cell_inter_cell_interference_moy{ii} + pilot_contamination_moy{ii};
    end
    
    snr_signal_cell_moy = cellfun(@(x,y) x./(y+1/(numOfAnt*rho_lin)),energy_signal_cell_moy,total_interference,'UniformOutput',false);
    
    
    for ii=1:L
        snr_signal_cell_moy_tab(index_kappa,ii) = mean(snr_signal_cell_moy{ii});
    end
    for ii=1:L
        total_interference_MRT{ii} = intra_cell_inter_cell_interference_MRT_moy{ii} + pilot_contamination_MRT_moy{ii};
    end
    snr_signal_cell_MRT_moy=cellfun(@(x,y) x./(y+1/(numOfAnt*rho_lin)),energy_signal_cell_MRT_moy,total_interference_MRT,'UniformOutput',false);
    for ii=1:L
        snr_signal_cell_MRT_moy_tab(index_kappa,ii)=mean(snr_signal_cell_MRT_moy{ii});
    end
    [snr_signal_cell_RZF_th,snr_signal_cell_MRT_th]= results_mrt_rzf_th(H_bar,Phi,pathloss,somme_pathloss,lambda,Power,rho_training_lin,rho_lin);
    for ii=1:L
        snr_signal_cell_MRT_th_tab(index_kappa,ii)=mean(snr_signal_cell_MRT_th{ii});
        snr_signal_cell_RZF_th_tab(index_kappa,ii)=mean(snr_signal_cell_RZF_th{ii});
    end
    
end



figure(1); hold on; box on;

plot(kappa_tab,(mean(log2(1+real(snr_signal_cell_MRT_th_tab)),2)),'k-.','LineWidth',2);
plot(kappa_tab,(mean(log2(real(snr_signal_cell_RZF_th_tab)+1),2)),'k:','LineWidth',2);
plot(kappa_tab,(mean(log2(1+real(snr_signal_cell_MRT_moy_tab)),2)),'kd','LineWidth',1);
plot(kappa_tab,(mean(log2(real(snr_signal_cell_moy_tab)+1),2)),'ko','LineWidth',1);


legend('MRC Approx','S-MMSE Approx','MRC Sim','S-MMSE Sim')

xlabel('Rician factor');
ylabel('Average sum SE [bit/s/Hz/UE]');
ylim([3 7]);
