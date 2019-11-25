function sCl = kmeansoptrnd(E,N,NCl)

[Np,Ne] = size(E);

%% Randomization
ERnd = zeros(Np,Ne);
for i = 1:Ne
    ERnd(:,i) = E(randperm(Np),i);
end

%% Covariace matrix
M = zeros(Ne,Ne);
parfor i = 1:Ne
    for j = 1:Ne
        M(i,j) = covnorm(ERnd(:,i),ERnd(:,j),0);
    end
end
M(isnan(M)) = 0;


%% k-means
parfor k = 1:N
    IDX = kmeans(M,NCl);
    s = silh(M,IDX);
    IDX0(k,:) = IDX;
    S(k) = median(s);
end
% keep best silhouette
[~,ClOK] = max(S);
IDX = IDX0(ClOK,:);
s = silh(M,IDX);
sCl = zeros(1,NCl);
for i = 1:NCl
    sCl(i) = median(s(IDX==i));
end
sCl = max(sCl);