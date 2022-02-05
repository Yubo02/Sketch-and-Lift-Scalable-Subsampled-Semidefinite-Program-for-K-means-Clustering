 [data.lev, levhdr]=fca_readfcs('Levine_32dim.fcs'); % Extract the data
 data.lev1=data.lev(1:72463,5:36); % Select the data we needed
 idx_ori=data.lev(1:72463,40);  % Select the correponding ground truth classifications 
 X0=data.lev1';

 n=72463; % Number of data points
 K=14;  % Number of groups
 p=32;  % Dimensions
 gama = 0.02; % Set the sub-sampling parameter
 r=1; % Number of rounds for WSL

% K-means++ method
tic
idx_hat_0=kmeansplus(X0,K);
t1=toc;


% Set the initialization for r-th round WSL
ind_old=idx_hat_0';
idd=ind_old;

for o=1:r

tic    
sumsub=zeros(1,K);
sumsub_v=ind_old;

for k=1:K
sumsub(k)=sum(ind_old==k);
sumsub_v(ind_old==k)=sumsub(k);
end


columns_new=(rand(1,n) <(gama*n/K)./sumsub_v' );



X_hat = X0(:,columns_new); %% New matrix with dimension p*q


X_2=kmeans_sdp_1(X_hat, K);

[U_2,S_2,V_2]=svd(X_2);

X_rel_hat= U_2(:,1:K);

idx_hat = kmeans(X_rel_hat,K);  % Get the cluster





X_hat_1=X_hat';
C_hat=zeros(p,K);
for cc=1:K
linearIndices = find(idx_hat==cc);
inter=mean(X_hat_1(linearIndices,:)); % Select the rows of tranposed X_hat.
C_hat(:,cc)=inter';

end

% Assign Xi to nearesr centroid of X_hat


idx_hat_3 = zeros(n,1);

for j=1:n

fmv=zeros(1,K);
for i=1:K
fmv(1,i)=norm(X0(:,j)-C_hat(:,i)); % Every point compared with centers
end
[mv,mp]=min(fmv);

idx_hat_3(j)=mp; % Assigned to the position of center

end

idd = [idd idx_hat_3]; % index of groups after kmeans_sdp & sketching


ind_old = idx_hat_3;
idd=[idd ind_old];
t2 = toc;

end


% Find the approximation error rate
error_1=zeros(1,3000);
error_2=zeros(1,3000);
for i=1:3000
Per=randsample(K,K); 

idx_loop0=ones(1,sum(idx_ori==1)).*Per(1);
for j=2:K
    idx_loop=ones(1,sum(idx_ori==j)).*Per(j);
    idx_loop0=[idx_loop0 idx_loop];    
end

error_1(i)=sum(idx_loop0~=idx_hat_0)/n;
error_2(i)=sum(idx_loop0'~=idx_hat_3)/n;
end
error_11=min(error_1); % Error rate for K-means++
error_12=min(error_2); % Error rate for WSL
t1; % Time cost for K-means++
t2; % Time cost for WSL




