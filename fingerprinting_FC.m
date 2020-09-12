function [Accuracy,Accuracy_chance] = fingerprinting_FC(FC_all)
% Fingerprinting analysis based on the functional connectivity matrices from a specific preprocessing pipeline

% FC_all = [nScans x nSubjects x nEdges] matrix, consisting of the upper triangular functional connectivity values for a specific scan and subject.
% nScans = number of scans being evaluated. For instance, in our paper we had 4 resting-state scans from the HCP data (Rest1_LR, Rest1_RL, Rest2_LR, Rest2_RL)
% nSubjects = number of subjects
% nEdges = number of upper triangular funcional connectivity values, it will depend on the atlas used

[nScans, nSubjects, nEdges] = size(FC_all);
FC_flat = reshape(FC_all(:,1:nSubjects,:),nScans*nSubjects,nEdges);   % scan1 subj1, scan2 subj1, scan3 subj1, scan4 subj1, scan1 subj2 ...
Accuracy = zeros(nScans,nScans);
Accuracy_chance = zeros(nScans,nScans);

for scan_i = 1:nScans
    FC_i = FC_flat(scan_i:nScans:end,:);   
    scans = 1:nScans;
    scans(scan_i) = [];
    for scan_j = scans          
        FC_j = FC_flat(scan_j:nScans:end,:);
        C_ij = corr(FC_i',FC_j');
        [~,ID] = max(C_ij,[],2);
        Accuracy(scan_i,scan_j) = (length(find(ID-(1:nSubjects)'==0))/nSubjects)*100;
        Accuracy_chance(scan_i,scan_j) = (length(find(ID-randi(nSubjects,nSubjects,1)==0))/nSubjects)*100;           % randperm(nSubj,nSubj), without replacement; randi(nSubj,nSubj,1), permutation with replacement
    end
end

end

