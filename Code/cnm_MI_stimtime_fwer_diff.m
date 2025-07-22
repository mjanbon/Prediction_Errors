function [MI_diff,MI_A, MI_B, sigMask] = cnm_MI_stimtime_fwer_diff(data_A, data_B, stim, nPerm)
    % Dimensions <trials x samples>
    nSamples = size(data_A', 2);
    
    % Number of stimuli
    Ns = length(unique(stim));
    
    % Perform copnorm
    cdata_A = copnorm(data_A');
    cdata_B = copnorm(data_B');
    
    % Compute MI
    MI_A = nan(nSamples, 1);
    for tIdx = 1 : nSamples
        MI_A(tIdx) = mi_model_gd(cdata_A(:, tIdx), stim, Ns);
    end

     % Compute MI
    MI_B = nan(nSamples, 1);
    for tIdx = 1 : nSamples
        MI_B(tIdx) = mi_model_gd(cdata_B(:, tIdx), stim, Ns);
    end

    MI_diff = MI_B - MI_A;
    
 if nPerm > 0
        MI_perm = zeros(nSamples, nPerm);   
        parfor permIdx = 1 : nPerm
            permIdx;
            
            % Randomize trials
            randorder = randperm(size(data_A', 1));
            for tIdx = 1 : nSamples
                MI_perm_A(tIdx, permIdx) = mi_model_gd(cdata_A(:, tIdx), stim(randorder), Ns);
                MI_perm_B(tIdx, permIdx) = mi_model_gd(cdata_B(:, tIdx), stim(randorder), Ns);
            end
        end
 end
    MI_perm = MI_perm_B - MI_perm_A

    % Take the max across the permutation dimension (1st dimension)
    Iperm_max = max(MI_perm, [], 1);  % Result is now 1 x timepoints x timepoints
    Iperm_min = min(MI_perm, [], 1);  % Result is now 1 x timepoints x timepoints

    %Flatten Iperm_max to a vector to include all timepoint pairs
    Iperm_max_vector = Iperm_max(:);  % Converts to a vector of all max values
    Iperm_min_vector = Iperm_min(:);  % Converts to a vector of all max values

    %Compute the 95th percentile for FWER thresholding
    fwer_thresh_max = prctile(Iperm_max_vector, 97.5);
    fwer_thresh_min = prctile(Iperm_min_vector, 2.5);
    
    
    %Initialize a significance mask for the observed CoI matrix
    sigMask = zeros(size(MI_diff));  % Assuming CoI_obs is your observed 1 x timepoints x timepoints matrix

    %Apply the threshold to create the significance mask
    sigMask(MI_diff > fwer_thresh_max) = 1;
    sigMask(MI_diff < fwer_thresh_min) = -1;
    