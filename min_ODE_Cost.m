function opt=min_ODE_Cost(StaOfNode,alpha,beta,NumOfNode,lambda)
%extract the number of nodes in different state, ODE not yet, using a rand
%to replace,assuming 10 samples of different initial state
interval=10;
%ODE=rand(100,5,interval);

%-------------------------------
%draft for obtaining the ODE function arrey
%require the ODE function output as a [M,5] matrix, where M indicates the
%time point, 5 indicates the state

for jj=1:numel(StaOfNode(1,:))
    
    NumOfOu = (NumOfNode - numel(StaOfNode(:,jj)))/NumOfNode;
    NumOfIn = sum(StaOfNode(:,jj)==0)/NumOfNode;
    NumOfRe = sum(StaOfNode(:,jj)==1)/NumOfNode;
    NumOfSa = sum(StaOfNode(:,jj)==3)/NumOfNode;
    NumOfCa = sum(StaOfNode(:,jj)==2)/NumOfNode;
    %derive the initial state of ODE
    INI=[NumOfIn,NumOfRe,NumOfCa,NumOfSa,NumOfOu];
    %
    %ODE=zeros(100,5,10);
    for ii=1:interval
        INI_C=INI;
        INI_C(1)=INI_C(1)-NumOfIn*ii/interval;
        INI_C(3)=INI_C(3)+NumOfIn*ii/interval;
        ODE_R=ODE_ICME(INI_C,lambda(jj));
        %find the maximum value of requesting status
        sta_max=find(ODE_R(:,2)==(max(ODE_R(:,2))));
        ODE(jj,:,ii)=ODE_R(sta_max,:);
    end
end

%-------------------------------
%find the minimum Cost
cost=alpha*ODE(:,2,:)+beta*ODE(:,3,:);
cost_function=reshape(cost,[numel(StaOfNode(1,:)),10]);
for ii=1:numel(ODE(:,1,1))
    [opt_cost(ii),opt_num(ii)]=min(cost_function(ii,:));
    opt(ii)=ODE(ii,3,opt_num(ii));
end
