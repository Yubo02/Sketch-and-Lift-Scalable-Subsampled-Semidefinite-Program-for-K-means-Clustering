
n=2000; % # of data points
K=4;  % # of groups
pch=[200 500 2000 4000 6000]; % dim of a points
sigma=1; % sd of the error

% iterate 10 times
AA=[0;0;0];

for ll=1:length(pch)
p=pch(ll); 
gama = 0.1;
l = 1.2;    % choose based on simple case
delt=sqrt((1+sqrt(1+K*p/log(n)/gama/n))*(4*log(n)))*l; % how faraway of the centers


mu_base = eye(p).*delt/sqrt(2); %Construct mu1,..,muK from rows of mu_base
Sigma1 = eye(p).*sigma;
part_preset=[n/K/2 n/K/2 n/K*1.5 n/K*1.5]; % Preset partition

tic
%Construct X matrix
X1=mvnrnd(mu_base(:,1),Sigma1,part_preset(1))';
for i=2:K
X2=mvnrnd(mu_base(:,i),Sigma1,part_preset(i))';
X1= [X1 X2];
end

X=X1;


ind_old=kmeans(X',K);

idd=ind_old;

for o=1:4
    
sumsub=[0 0 0 0];
sumsub_v=ind_old;

for k=1:K
sumsub(k)=sum(ind_old==k);
sumsub_v(ind_old==k)=sumsub(k);
end


columns_new=(rand(1,n) <(gama*n/K)./sumsub_v' );



X_hat = X(:,columns_new); %% New matrix with dimension p*q


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
fmv(1,i)=norm(X(:,j)-C_hat(:,i)); % Every point compared with centers
end
[mv,mp]=min(fmv);

idx_hat_3(j)=mp; % Assigned to the position of center

end

idd = [idd idx_hat_3]; % index of groups after kmeans_sdp & sketching


ind_old = idx_hat_3;


end

timeElapsed_5 = toc;








% Determine the error
order = 1:K;
Per = perms(order);
error_1=zeros(1,size(Per,1));
for i=1:size(Per,1)
    
 % index of groups for origin centers idx_ori=[ones(1,10) 2*ones(1,4) 3*ones(1,6)]'
idx_ori=ones(1,part_preset(1)).*Per(i,1);
for j=2:K
    idx_loop=ones(1,part_preset(j)).*Per(i,j);
    idx_ori=[idx_ori idx_loop];    
end
error_4(i)=sum(idx_ori'~=idd(:,5))/n;


end

error_ks_sk4=min(error_4);
Error= error_ks_sk4;
Time=timeElapsed_5;
BB=[Error;Time;p];


AA=[AA BB];
end

err5_p=[AA(1,2:size(AA,2)); AA(3,2:size(AA,2))]; % Error rate for different p
tim5_p=[AA(2,2:size(AA,2)); AA(3,2:size(AA,2))]; % Time cost for different p