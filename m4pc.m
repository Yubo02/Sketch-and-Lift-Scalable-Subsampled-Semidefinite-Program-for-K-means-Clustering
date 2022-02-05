
n=2000; % # of data points
K=4;  % # of groups
pch=[200 500 2000 4000 6000]; % dim of a points
sigma=1; % sd of the error

% initiate the matrix for storing errors and time cost with respect to p
AA=[0;0;0];

for ll=1:length(pch)
p=pch(ll); 
gama = 0.1;
l = 1.2;    % choose based on simple case
delt=sqrt((1+sqrt(1+K*p/log(n)/gama/n))*(4*log(n)))*l; % how faraway of the centers


mu_base = eye(p).*delt/sqrt(2); %Construct mu1,..,muK from rows of mu_base
Sigma1 = eye(p).*sigma;
part_preset=[n/K n/K n/K n/K]; % Preset partition

%Construct X matrix
X1=mvnrnd(mu_base(:,1),Sigma1,part_preset(1))';
for i=2:K
X2=mvnrnd(mu_base(:,i),Sigma1,part_preset(i))';
X1= [X1 X2];
end

X=X1;




tic
% Second do kmeans sdp for shrinked X_hat's

q =gama*n;  % Random select q data points

n_q=n/q;


C_whole=zeros(p,1);
for ac=1:n_q
columns = randperm(size(X,2),q);
X_hat = X(:,columns); %% New matrix with dimension p*q


X_2=kmeans_sdp_1(X_hat, K);

[U_2,S_2,V_2]=svd(X_2);

X_rel_hat= U_2(:,1:K);

idx_hat = kmeans(X_rel_hat,K);  % Get the cluster

% Get the centers

X_hat_1=X_hat';

C_hat=zeros(p,K);
for cc=1:K
linearIndices = find(idx_hat==cc);
inter=mean(X_hat_1(linearIndices,:)); % Select the rows of tranposed X_hat.
C_hat(:,cc)=inter';
end

C_whole=[C_whole C_hat];
X(:,columns)=[];
end

C_whole_tr=C_whole(:,2:size(C_whole,2));

C_final_cent=optim(C_whole_tr,K,n_q);  % Get the centers by opmatch
          
% Assign Xi to nearesr centroid of X_hat's


idx_hat_f = zeros(n,1);

for j=1:n

fmv=zeros(1,K);
for i=1:K
fmv(1,i)=norm(X1(:,j)-C_final_cent(:,i)); % Every point compared with centers
end
[mv,mp]=min(fmv);

idx_hat_f(j)=mp; % Assigned to the position of center

end

idx_hat_f; % index of groups after kmeans_sdp & sketching
timeElapsed_4 = toc;


% Determine the error
order = 1:K;
Per = perms(order);
error_1=zeros(1,size(Per,1));
for i=1:size(Per,1);
    
 % index of groups for origin centers idx_ori=[ones(1,10) 2*ones(1,4) 3*ones(1,6)]'
idx_ori=ones(1,part_preset(1)).*Per(i,1);
for j=2:K
    idx_loop=ones(1,part_preset(j)).*Per(i,j);
    idx_ori=[idx_ori idx_loop];    
end

error_1(i)=sum(idx_ori'~=idx_hat_f)/n;

end
error_ks_sk1=min(error_1);


Error=error_ks_sk1;

Time=timeElapsed_4;

BB=[Error;Time;p];


AA=[AA BB];
end

err4_p=[AA(1,2:size(AA,2)); AA(3,2:size(AA,2))]; % Error rate for different p
tim4_p=[AA(2,2:size(AA,2)); AA(3,2:size(AA,2))]; % Time cost for different p