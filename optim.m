
% Given matrix of Centers, get the final centers


function X=optim(C_whole_tr,K,n_q)


Sum_ce=C_whole_tr(:,1:K);

for l=1:n_q-1
Cost_m=zeros(K,K); % Write down the cost matrix
for j=1:K
for i=1:K
    
    Cost_m(i,j)=norm(C_whole_tr(:,i+l*K)-C_whole_tr(:,j)); 
end            
end              
  Matc = matchpairs(Cost_m,10000000,'min');

for k=1:K
  Sum_ce(:,k)=Sum_ce(:,k)+C_whole_tr(:,l*K+Matc(k,1)); % Get the sum of each matched centers
end

end

X=Sum_ce/n_q;
end



