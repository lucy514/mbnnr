function KV=similarity_disease(interaction)
%SIMILARITY_DISEASE Calculate Gaussian similarity matrix for diseases
%   KV = SIMILARITY_DISEASE(interaction) calculates the Gaussian kernel
%   similarity matrix between diseases based on their interaction patterns
%
%   Inputs:
%       interaction - interaction matrix where columns represent diseases
%
%   Outputs:
%       KV - disease-disease similarity matrix

[~,nm]=size(interaction);
for i=1:nm
    sv(i)=norm(interaction(:,i))^2;  % calculate gamma for Gaussian kernel calculation
end
gamam=nm/sum(sv');
for i=1:nm
    for j=1:nm
        kv(i,j)=exp(-gamam*(norm(interaction(:,i)-interaction(:,j)))^2);  % calculate Gaussian kernel similarity
    end
end 
KV=kv;  % the integrated similarity between diseases

% Note: This function calculates disease-disease similarity using Gaussian kernel
% based on interaction patterns to avoid data leakage in cross-validation
