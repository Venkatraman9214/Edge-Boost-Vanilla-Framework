function C_allocation=Cachingpolicy(StaOfNode,alpha,beta,NumOfNode,lambda)

min_Cost=min_ODE_Cost(StaOfNode,alpha,beta,NumOfNode,lambda);

for ii=1:numel(StaOfNode(1,:))
    %calculate the current number of replicas
    C_current(ii)=sum(find(StaOfNode(:,ii)==3));
    %calculate the minimum caching cost
    C_minC(ii)=min_Cost(ii);
    %calculate the caching gap
    C_Gap(ii)= C_minC(ii)-C_current(ii);
end

%calculate the remaining caching space
C_c=numel(StaOfNode(find(StaOfNode==3)));
%Let caching space equal to 10;
C_total=numel(StaOfNode(:,1))*10;
C_remain=C_total-C_c;

%Caching allocation policy
for ii=1:numel(StaOfNode(1,:))
    %find the caching gap larger than 0
    N=sum((C_Gap(find(C_Gap>0))));
    %calculate the caching allocation portion
    if(C_Gap(ii)>0)
        C_allocation(ii)=C_remain*C_Gap(ii)/N;
    else 
        C_allocation(ii)=0;
    end

end

 
