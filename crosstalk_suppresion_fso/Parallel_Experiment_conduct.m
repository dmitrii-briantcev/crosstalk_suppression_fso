%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmitrii Briantcev, 2020
%
% Conduct params.iter independent system runs, each over 
% params.time_lim time steps
%
% Output is written in a .mat file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Parallel_Experiment_conduct (SNR_db, r0, poolnum)

tic

% Outer initialization
params_main = params_setup(r0, SNR_db);
tlim_r = params_main.time_lim-params_main.delay;
CoeffArray_1 = cell(1, params_main.iter);
CoeffArray_2 = cell(1, params_main.iter);
Fidelity_1 = cell(1, params_main.iter);
Efficiency_1 = cell(1, params_main.iter);
Fidelity_2 = cell(1, params_main.iter);
Efficiency_2 = cell(1, params_main.iter);
E1 = cell(1, params_main.iter);
E2 = cell(1, params_main.iter);
TF_1 = cell(1, params_main.iter);
TF_2 = cell(1, params_main.iter);
phz = cell(1, params_main.iter);
seed = cell(1, params_main.iter);
Precoding_1 = cell(1, params_main.iter);
Precoding_2 = cell(1, params_main.iter);
C1 = cell(1, params_main.iter);
C2 = cell(1, params_main.iter);
E = cell(1, params_main.iter);
Temp_vec = cell(1, params_main.iter);

parfor iter = 1:params_main.iter 
    
    disp(iter)
    params(iter) = params_main;

    % Initialization for a run
    CoeffArray_1{iter} = zeros(length(params(iter).m_alph));
    CoeffArray_2{iter} = zeros(length(params(iter).m_alph));
    Fidelity_1{iter} = zeros(length(params(iter).m_alph));
    Efficiency_1{iter} = zeros(1,tlim_r);
    Fidelity_2{iter} = zeros(length(params(iter).m_alph));
    Efficiency_2{iter} = zeros(1,tlim_r);
    	
    % Random spectrum initialization
	seed{iter} = iter*100 + poolnum;
    [C1{iter}, C2{iter}] = gen_spec(params(iter),seed{iter});
    
    % Random starting time shift initialization
    shift = randi([0,1000],1,1);
	  
    for t = 1:tlim_r
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % First transciever
        
        Precoding_1{iter} = eye(length(params(iter).m_alph));
        
        for j=1:length(params(iter).m_alph)
            
            E{iter} = zeros(size(params(iter).x,2),size(params(iter).y,2));
            
            % Pre-encoding
            for k=1:length(params(iter).m_alph)
                E{iter} = E{iter} + OptMode (params(iter).m_alph(k), params(iter).n_alph(k), params(iter))*Precoding_1{iter}(j,k);
            end
            
            % Normalization & SNR
            E{iter} = E{iter}/norm(E{iter}, 'fro').*sqrt(params(iter).SNR);
            
            % Propagation
            for i = 1:params(iter).z_div
                E{iter} = propagate(E{iter},params(iter),params(iter).dz);
                phz{iter} = Phz(params(iter), t+shift, C1{iter}{i}, C2{iter}{i});
                E{iter} = E{iter}.*exp(1i*phz{iter});
            end
            
            % Decomposition
            CoeffArray_1{iter}(j,:) = Decompose (E{iter}, params(iter));
			
			% Saving Coeffs without Errors
            Coeff_mem_1{iter}(:,:,t) = CoeffArray_1{iter};
            
            % Introduce responsivity
            CoeffArray_1{iter}(j,:) = CoeffArray_1{iter}(j,:)*sqrt(params(iter).responsivity);
            
            % Introduce noise and phase estimation error
            for k = 1:length(params(iter).m_alph)
                CoeffArray_1{iter}(j,k) = CoeffArray_1{iter}(j,k)*exp(1i*normrnd(0,params(iter).phase_error));
                CoeffArray_1{iter}(j,k) = CoeffArray_1{iter}(j,k) + (normrnd(0,1)+1i*normrnd(0,1))/sqrt(2);
            end
            
            % Get mode fidelity for a time frame
            Fidelity_1{iter}(j,:) = abs(CoeffArray_1{iter}(j,:)).^2;
              
        end
        
        % Get overall- and mode- efficiency
        Efficiency_1{iter}(t) = trace(Fidelity_1{iter})/sum(sum(Fidelity_1{iter}));
        for k = 1:length(params(iter).m_alph)
            E1{iter}(k,t) = Fidelity_1{iter}(k,k)/sum(Fidelity_1{iter}(k,:));
			Temp_vec{iter} = Fidelity_1{iter}(k,:);
            Temp_vec{iter}(k) = 0;
            TF_1{iter}(k,t) = Fidelity_1{iter}(k,k)/(Fidelity_1{iter}(k,k)+max(Temp_vec{iter}));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Second transciever
        
        Precoding_2{iter} = pinv(CoeffArray_1{iter});
     
        for j=1:length(params(iter).m_alph)
            
            % Pre-encoding
            E{iter} = zeros(size(params(iter).x,2),size(params(iter).y,2));
            for k=1:length(params(iter).m_alph)
                E{iter} = E{iter} + OptMode (params(iter).m_alph(k), params(iter).n_alph(k), params(iter))*Precoding_2{iter}(j,k);
            end
            
            % Normalization & SNR
            E{iter} = E{iter}/norm(E{iter}, 'fro').*sqrt(params(iter).SNR);
            
            % Propagation
            for i = 1:params(iter).z_div
                E{iter} = propagate(E{iter},params(iter),params(iter).dz);
                phz{iter} = Phz(params(iter), t+shift+params(iter).delay, C1{iter}{i}, C2{iter}{i});
                E{iter} = E{iter}.*exp(1i*phz{iter});
            end
            
            % Decomposition
            CoeffArray_2{iter}(j,:) = Decompose (E{iter}, params(iter));
			
			% Saving Coeffs without Errors
            Coeff_mem_2{iter}(:,:,t) = CoeffArray_2{iter};
            
            % Introduce responsivity
            CoeffArray_2{iter}(j,:) = CoeffArray_2{iter}(j,:)*sqrt(params(iter).responsivity);
            
            % Introduce noise and phase estimation error
            for k = 1:length(params(iter).m_alph)
                CoeffArray_2{iter}(j,k) = CoeffArray_2{iter}(j,k)*exp(1i*normrnd(0,params(iter).phase_error));
                CoeffArray_2{iter}(j,k) = CoeffArray_2{iter}(j,k) + (normrnd(0,1)+1i*normrnd(0,1))/sqrt(2);
            end
            
            % Get fidelity for a time frame
            Fidelity_2{iter}(j,:) = abs(CoeffArray_2{iter}(j,:)).^2;
              
        end
        
        % Get overall- and mode- efficiency, transmit fidelity for the
        % system
        Efficiency_2{iter}(t) = trace(Fidelity_2{iter})/sum(sum(Fidelity_2{iter}));
        for k = 1:length(params(iter).m_alph)
            E2{iter}(k,t) = Fidelity_2{iter}(k,k)/sum(Fidelity_2{iter}(k,:));
            Temp_vec{iter} = Fidelity_2{iter}(k,:);
            Temp_vec{iter}(k) = 0;
            TF_2{iter}(k,t) = Fidelity_2{iter}(k,k)/(Fidelity_2{iter}(k,k)+max(Temp_vec{iter}));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end

% Save the Workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = strcat(num2str(params_main.SNR_dB),'-',num2str(params_main.r0),'-',num2str(poolnum),'.mat');
save(name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc
















































