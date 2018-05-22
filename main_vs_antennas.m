% This Matlab script can be used to reproduce the spectral efficiecny of a
% multicell Massive MIMO system with uncorellated Rician fading channels as
% a function of BS antennas. In particular, it can be used to
% regenerate Fig. 4 of the following work, which is currently under review
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
num_iter = 1000;
% parpool(40);

% Load system parameters
define


% Values of BS antennas to simulate
numOfAnt_tab = [20:20:100,150,200,300,400];

for index_tab=1:length(numOfAnt_tab)
    
    % Load the current value of BS antennas
    numOfAnt=numOfAnt_tab(index_tab);
    
    % Scaling coefficient for RZF
    lambda = numOfUEs/(numOfAnt*rho_lin);
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
                % H_bar{ii}=zeros(numOfAnt,numOfUEs);
            end
            somme_pathloss{ii}=somme_pathloss{ii}+ pathloss{ii}(:,jj);
        end
        
        for jj=1:L
            Phi{ii}(:,jj)=pathloss{ii}(:,ii).*pathloss{ii}(:,jj)./(1/rho_training_lin+somme_pathloss{ii});
        end
        
    end
    
    parfor iter=1:num_iter
        
        %[energy_signal_cell{iter},pilot_contamination{iter},intra_cell_inter_cell_interference{iter},first_term_interference_emp{iter},second_term_interference_emp{iter},third_term_interference_emp{iter}]=system_random(H_bar,Phi,pathloss,somme_pathloss,lambda,Power,rho_training_lin);
        [RZF_results,MRT_results]= system_random(H_bar,Phi,pathloss,somme_pathloss,lambda,Power,rho_training_lin);
        energy_signal_cell{iter}=RZF_results.energy_signal_cell;
        pilot_contamination{iter}=RZF_results.pilot_contamination;
        intra_cell_inter_cell_interference{iter}= RZF_results.intra_cell_inter_cell_interference;
        first_term_interference_emp{iter}=RZF_results.first_term_interference_emp;
        second_term_interference_emp{iter}=RZF_results.second_term_interference_emp;
        third_term_interference_emp{iter}= RZF_results.third_term_interference_emp;
        energy_signal_cell_MRT{iter}=MRT_results.energy_signal_cell;
        pilot_contamination_MRT{iter}=MRT_results.pilot_contamination;
        intra_cell_inter_cell_interference_MRT{iter}=MRT_results.intra_cell_inter_cell_interference ;
        
        
        %[energy_signal_cell{iter},pilot_contamination{iter},intra_cell_inter_cell_interference{iter},first_term_interference_emp{iter},second_term_interference_emp{iter},third_term_interference_emp{iter}]=system_random(H_bar,Phi,pathloss,somme_pathloss,lambda,Power,rho_training_lin);
        
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
            
        end
    end
    
    
    energy_signal_cell_moy=cellfun(@(x) x./(num_iter),energy_signal_cell_moy,'UniformOutput',false);
    pilot_contamination_moy=cellfun(@(x) x./(num_iter),pilot_contamination_moy,'UniformOutput',false);
    intra_cell_inter_cell_interference_moy=cellfun(@(x) x./(num_iter),intra_cell_inter_cell_interference_moy,'UniformOutput',false);
    
    
    energy_signal_cell_MRT_moy=cellfun(@(x) x./(num_iter),energy_signal_cell_MRT_moy,'UniformOutput',false);
    pilot_contamination_MRT_moy= cellfun(@(x) x./(num_iter),pilot_contamination_MRT_moy,'UniformOutput',false);
    intra_cell_inter_cell_interference_MRT_moy=cellfun(@(x) x./(num_iter),intra_cell_inter_cell_interference_MRT_moy,'UniformOutput',false);
    
    
    
    first_term_interference_emp_moy=cellfun(@(x) x./(num_iter),first_term_interference_emp_moy,'UniformOutput',false);
    second_term_interference_emp_moy=cellfun(@(x) x./(num_iter),second_term_interference_emp_moy,'UniformOutput',false);
    third_term_interference_emp_moy=cellfun(@(x) x./(num_iter),third_term_interference_emp_moy,'UniformOutput',false);
    
    
    %%%%%% Equivalent deterministe: Checking of each term %%%%%%%%%%%%%%%%%%%
    
    delta_tilde_tot=zeros(L,1);
    delta_tot=zeros(L,1);
    gamma_d_tab=zeros(L,1);
    gamma_tilde_d_tab=zeros(L,1);
    F=zeros(L,1);
    Delta=zeros(L,1);
    %ubar=zeros(numOfUEs,L);
    
    
    for ii=1:L
        [delta_tilde,delta,T,T_tilde,Psi,Psi_tilde]=deterministic_matrix(H_bar{ii}/sqrt(numOfAnt),diag(Phi{ii}(:,ii)),lambda);
        
        delta_tilde_tot(ii)=delta_tilde;
        delta_tot(ii)=delta;
        T_tab{ii}=T;
        T_tilde_tab{ii}=T_tilde;
        gamma_d_tab(ii)=1/numOfAnt*trace(T^2);
        D_tilde=diag(Phi{ii}(:,ii));
        gamma_tilde_d_tab(ii)=1/numOfAnt*trace(diag(Phi{ii}(:,ii))*T_tilde*diag(Phi{ii}(:,ii))*T_tilde);
        F(ii)=(1/(numOfAnt^2))*trace(T^2*H_bar{ii}*diag((Phi{ii}(:,ii))./((1+delta*Phi{ii}(:,ii)).^2))*H_bar{ii}');
        F1(ii)=(1/(numOfAnt^2))*trace(T^2*H_bar{ii}*D_tilde*inv(eye(numOfUEs)+delta*D_tilde)^2*H_bar{ii}');
        Delta(ii)=(1-F(ii))^2-lambda^2*gamma_d_tab(ii)*gamma_tilde_d_tab(ii);
        nu_d(ii)=(1/Delta(ii))*(1/numOfAnt)*trace(T^2);
        Psi_bar(ii)=(numOfUEs/numOfAnt) /(delta-lambda*nu_d(ii)); %%%% Equivalent deterministe of betaii/N ***********Checking DONE **********
        
        %%%% Slow implementation
        %for k=1:numOfUEs
        %   ubar(k,ii)= delta_tot(ii)*Phi{ii}(k,ii) + 1/(lambda * T_tilde(k,k)) * (H_bar{ii}(:,k))'*T*H_bar{ii}(:,k) /(numOfAnt*(1+delta*Phi{ii}(k,ii)));
        %end
        
        %%%% Fast implementation
        ubar{ii}= delta_tot(ii) * Phi{ii}(:,ii) + (diag(H_bar{ii}'*T*H_bar{ii}).*(1./(lambda*diag(T_tilde))))./(numOfAnt*(1+delta*Phi{ii}(:,ii)));
        omega{ii}=delta_tot(ii) * pathloss{ii}(:,ii) + diag(H_bar{ii}'*T*H_bar{ii}).*(1./(lambda*diag(T_tilde)))./(numOfAnt*(1+delta*Phi{ii}(:,ii)));
        
        
        energy_signal_cell_th{ii}=Psi_bar(ii)*ubar{ii}.^2./((1+ubar{ii}).^2); %%%% Checking of the signal term %%%%%**********Checking DONE **********
        
    end
    
    for l=1:L
        
        zeta(:,l)= (1-F(l))/Delta(l) * (diag(H_bar{l}'*T_tab{l}^2*H_bar{l}).*(1./(lambda*diag(T_tilde_tab{l})).^2))./(numOfAnt*(1+delta_tot(l)*Phi{l}(:,l)).^2);
        
        for k=1:numOfUEs
            for ii=1:numOfUEs
                if (ii~=k)
                    zeta(k,l) =zeta(k,l) + gamma_d_tab(l)/Delta(l)*Phi{l}(ii,l)/(1+Phi{l}(ii,l)*delta_tot(l))^2 *(1/(numOfAnt^2)) * (H_bar{l}(:,k))'*T_tab{l}*H_bar{l}(:,ii)*(H_bar{l}(:,ii))'*T_tab{l}*H_bar{l}(:,k)/(lambda^2*(T_tilde_tab{l}(k,k))^2*(1+delta_tot(l)*Phi{l}(k,l))^2);
                end
            end
        end
        
    end
    %%%%%%%%% Checking of pilot contamination: Some difference due to the fact that the values are small, but I think results are correct  %%%%%%%%%%%
    for l=1:L
        pilot_contamination_th{l}=0;
        for lcomp=1:L
            if (lcomp~=l)
                pilot_contamination_th{l} = pilot_contamination_th{l} + Psi_bar(lcomp)*delta_tot(lcomp)^2*Phi{lcomp}(:,l).^2./((1+ubar{lcomp}).^2);
            end
        end
    end
    
    
    %%%%% Checking first term interference  %%%%%%%%%%%%%% %%%%%%%%%%%%% First term checked: Good results %%%%%%%%%%%%%%
    for l=1:L
        first_term_inter_intra_interference{l}=0;
        for lcomp=1:L
            if (lcomp~=l)
                first_term_inter_intra_interference{l}= first_term_inter_intra_interference{l} + Psi_bar(lcomp)*(pathloss{lcomp}(:,l)*delta_tot(lcomp)-delta_tot(lcomp)^2*Phi{lcomp}(:,l).^2./((1+ubar{lcomp})));
            else
                first_term_inter_intra_interference{l}= first_term_inter_intra_interference{l}+ Psi_bar(l)*(omega{l}-ubar{l}.^2./(1+ubar{l}));
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for l=1:L
        second_term_inter_intra_interference{l}=0;
        for lcomp=1:L
            if (lcomp~=l)
                second_term_inter_intra_interference{l} = second_term_inter_intra_interference{l} +  Psi_bar(lcomp) * (Phi{lcomp}(:,l)*delta_tot(lcomp)./(1+ubar{lcomp})).^2;
            else
                second_term_inter_intra_interference{l}=  second_term_inter_intra_interference{l} +  Psi_bar(lcomp)*(ubar{lcomp}./(1+ubar{lcomp})).^2;
            end
        end
    end
    
    for l=1:L
        third_term_inter_intra_interference{l}=0;
        for lcomp=1:L
            if (lcomp~=l)
                mu=pathloss{lcomp}(:,l) * nu_d(lcomp) - Phi{lcomp}(:,l).^2.*delta_tot(lcomp).*(2*nu_d(lcomp)*(1+ubar{lcomp})-delta_tot(lcomp)*(Phi{lcomp}(:,lcomp)*nu_d(lcomp)+zeta(:,lcomp)))./((1+ubar{lcomp}).^2);
                third_term_inter_intra_interference{l}= third_term_inter_intra_interference{l} + Psi_bar(lcomp) *mu;
            else
                mu=pathloss{l}(:,l)*nu_d(l)+zeta(:,l)- ubar{l}.*(Phi{l}(:,l)*nu_d(l)+zeta(:,l)) .* (2+ubar{l})./((1+ubar{l}).^2);
                third_term_inter_intra_interference{l}= third_term_inter_intra_interference{l}+ Psi_bar(l)*mu;
            end
        end
        
    end
    
    for ii=1:L
        total_interference_th{ii}= first_term_inter_intra_interference{ii} - lambda* third_term_inter_intra_interference{ii}-second_term_inter_intra_interference{ii} +pilot_contamination_th{ii};
        intra_cell_inter_cell_interference_th{ii} = first_term_inter_intra_interference{ii} - lambda* third_term_inter_intra_interference{ii} -second_term_inter_intra_interference{ii};
        
    end
    for ii=1:L
        total_interference{ii}=intra_cell_inter_cell_interference_moy{ii} + pilot_contamination_moy{ii};
    end
    snr_signal_cell_moy=cellfun(@(x,y) x./(y+1/(numOfAnt*rho_lin)),energy_signal_cell_moy,total_interference,'UniformOutput',false);
    
    snr_signal_cell_th=cellfun(@(x,y) x./(y+1/(numOfAnt*rho_lin)),energy_signal_cell_th,total_interference_th,'UniformOutput',false);
    
    for ii=1:L
        snr_signal_cell_moy_tab(index_tab,ii)=mean(snr_signal_cell_moy{ii});
        snr_signal_cell_th_tab(index_tab,ii)=mean(snr_signal_cell_th{ii});
    end
    
    
    
    %%%%%%%%%%%%%%Theoretical values MRT %%%%%%%%%%%%%HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    for ii=1:L
        Theta_bar(ii)=1/(1/numOfUEs*(sum(Phi{ii}(:,ii)) + 1/numOfAnt * trace((H_bar{ii})'*H_bar{ii})));
    end
    
    for l=1:L
        energy_signal_cell_MRT_th{l}=Theta_bar(l)*(Phi{l}(:,l)+1/numOfAnt*diag((H_bar{l})'*H_bar{l})).^2;
        pilot_contamination_MRT_th{l}=0;
        for lcomp=1:L
            if (lcomp~=l)
                pilot_contamination_MRT_th{l}=pilot_contamination_MRT_th{l}+Theta_bar(lcomp)*Phi{lcomp}(:,l).^2;
            end
        end
    end
    
    for l=1:L
        first_term_interference_MRT{l}=0;
        for lcomp=1:L
            first_term_interference_MRT{l}=first_term_interference_MRT{l} + Theta_bar(lcomp)/numOfAnt*(pathloss{lcomp}(:,l))*(sum(Phi{lcomp}(:,lcomp)) + 1/numOfAnt*trace((H_bar{lcomp})'*H_bar{lcomp}));
        end
    end
    
    for l=1:L
        second_term_interference_MRT{l}= 1/(numOfAnt^2)*Theta_bar(l)*(sum(Phi{l}(:,l))*diag((H_bar{l})'*H_bar{l})-diag((H_bar{l})'*H_bar{l}).*Phi{l}(:,l)+diag((H_bar{l})'*H_bar{l}*(H_bar{l})'*H_bar{l})-diag((H_bar{l})'*H_bar{l}).^2);
    end
    
    
    for l=1:L
        interference_MRT_th{l}= first_term_interference_MRT{l} + second_term_interference_MRT{l};
        total_interference_MRT_th{l}= interference_MRT_th{l} + pilot_contamination_MRT_th{l};
    end
    
    
    for ii=1:L
        total_interference_MRT{ii}=intra_cell_inter_cell_interference_MRT_moy{ii} + pilot_contamination_MRT_moy{ii};
    end
    snr_signal_cell_MRT_moy=cellfun(@(x,y) x./(y+1/(numOfAnt*rho_lin)),energy_signal_cell_MRT_moy,total_interference_MRT,'UniformOutput',false);
    
    snr_signal_cell_MRT_th=cellfun(@(x,y) x./(y+1/(numOfAnt*rho_lin)),energy_signal_cell_MRT_th,total_interference_MRT_th,'UniformOutput',false);
    
    [total_interference_th_bis,energy_signal_cell_th_bis,Psi_bar_bis,snr_signal_cell_RZF_th_bis]= results_mrt_rzf_th_new_formula(H_bar,Phi,pathloss,somme_pathloss,lambda,Power,rho_training_lin,rho_lin);
    
    for ii=1:L
        snr_signal_cell_MRT_moy_tab(index_tab,ii)=mean(snr_signal_cell_MRT_moy{ii});
        snr_signal_cell_MRT_th_tab(index_tab,ii)=mean(snr_signal_cell_MRT_th{ii});
    end
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); hold on; box on;
plot(numOfAnt_tab,mean(log2(1+real(snr_signal_cell_MRT_th_tab)),2),'k-.','LineWidth',2);
plot(numOfAnt_tab,mean(log2(real(snr_signal_cell_th_tab)+1),2),'k:','LineWidth',2);
plot(numOfAnt_tab,mean(log2(1+real(snr_signal_cell_MRT_moy_tab)),2),'kd','LineWidth',1);
plot(numOfAnt_tab,mean(log2(1+real(snr_signal_cell_moy_tab)),2),'ko','LineWidth',1);

legend('MRC Approx','S-MMSE Approx','MRC Sim','S-MMSE Sim')

xlabel('Number of antennas (N)');
ylabel('Average sum SE [bit/s/Hz/UE]');
ylim([0 9]);

