clear all;
%% ��������
NumOfNode = 500;  % �ڵ�����
NumOfCon = 100;    % ��������
SimuTime = 1000;   % ����ʱ��
minVec = 4;       % ��С�ٶ�
maxVec = 10;      % ����ٶ�
Square = 2000;    % ���η�Χ
radius = 500;     % BS�뾶
MaxArr = 20;      % ��󵽴�����
D2DRad = 150;      % D2Dͨ�ž���
CahCap = 10;
Q_Outputs=zeros(NumOfCon,SimuTime);
Q_delay=zeros(NumOfCon,SimuTime);
Times= 10;
% ChunkBuf = 5;     % ÿ�����ݿ��buffer��С
% PlayDura = 1;     % ÿ��time slot��buffer����

alpha = 0.5;       % Ȩ������
beta = 1-alpha;    %

lambda = randi(ceil(MaxArr/NumOfCon),NumOfCon,1);  % ��ͬ���ݵĵ�������
CapaOfBs = 200;                      % ��վ��·������ÿ��time slot���Է���10���û�
%CapaOfNode = 1;                     % �ڵ���·������ÿ��time slot���Է���1���û�

StaOfNode = zeros(NumOfNode,NumOfCon);  % �ڵ�״̬����

% �����վ�����������Ż���
ValTmp = (Square-radius)/radius;
if(mod(ValTmp,2)==1)
    NumOfBS = (2*ValTmp-1)*(ValTmp-1)/2+2;
else
    NumOfBS = (2*ValTmp-1)*ValTmp/2;
end
% ��վλ�����ã���̬����
PosOfBS = zeros(NumOfBS,2);
PosOfBS(1,:) = [250,250];
PosOfBS(2,:) = [1000,0];
PosOfBS(3,:) = [1750,250];
PosOfBS(4,:) = [1000,1000];
PosOfBS(5,:) = [250,1250];
PosOfBS(6,:) = [1000,2000];
PosOfBS(7,:) = [1750,1750];

% RWP �ڵ�λ�����ʾ���
PosOfNode = zeros(NumOfNode,5);%���ɾ���
PosOfNode(:,1:4) = randi(Square,NumOfNode,4);%��ʼλ���Լ�Ŀ��λ��
PosOfNode(:,5) = randi([minVec,maxVec],NumOfNode,1);%��ʼ�ٶ�

% ͳ����Ϣ
ConStaTime = cell(NumOfCon,1);
IndOfBS = 4;
HitRate = zeros(SimuTime,3);
sumDelay = zeros(SimuTime,1);


for tt = 1:Times
for ii = 1:SimuTime
    
    %% ͳ����Ϣ
    MissNum = 0;
    HitNum = 0;
    
    %% RWP ����
    for jj = 1:NumOfNode
        DeltaX = PosOfNode(jj,3)-PosOfNode(jj,1);
        DeltaY = PosOfNode(jj,4)-PosOfNode(jj,2);
        %ȷ���ƶ�����
        if(DeltaX>0)
            SinVal = sqrt(DeltaY^2/(DeltaX^2+DeltaY^2));
        else
            SinVal = -sqrt(DeltaY^2/(DeltaX^2+DeltaY^2));
        end
        if(DeltaY>0)
            CosVal = sqrt(DeltaX^2/(DeltaX^2+DeltaY^2));
        else
            CosVal = -sqrt(DeltaX^2/(DeltaX^2+DeltaY^2));
        end
        % �ƶ������жϣ�����atan����Ƕ�
        PosOfNode(jj,1) = PosOfNode(jj,1) + PosOfNode(jj,5)*SinVal;
        PosOfNode(jj,2) = PosOfNode(jj,2) + PosOfNode(jj,5)*CosVal;
        % �����յ�λ�ã�����ѡ�㡣
        if(PosOfNode(jj,1)>=PosOfNode(jj,3))
            PosOfNode(jj,1) = PosOfNode(jj,3);
            PosOfNode(jj,2) = PosOfNode(jj,4);
            PosOfNode(jj,3:4) = randi(Square,1,2);
            PosOfNode(jj,5) = randi([minVec,maxVec],1,1);
        end
    end
    
    %% ��������
    ReqTmp = random('Poisson',lambda);
    for jj = 1:numel(ReqTmp)
        if(ReqTmp(jj)~=0)
            nTmp = find(StaOfNode(:,jj)~=1);
            % ��ǰ�����û����������˳�
            if(numel(nTmp)==0)
                if ii==1
                  Q_Outputs(:,ii)=0;
                else
                  Q_Outputs(jj,ii)=Q_Outputs(jj,ii-1);
                end
                continue;
            end
            if(numel(nTmp)<ReqTmp(jj))
                ReqTmp(jj) = numel(nTmp);
            end
            rand_index = nTmp(randperm(numel(nTmp)));                % ������������
            draw_rand_index = rand_index(1:ReqTmp(jj));              % ȡ��ǰReqTmp(jj)�����
            % ͳ�ƶӳ�����
 
            % ״̬ת�� �����޸ģ�
%             for zz = 1:ReqTmp(jj)
%                 StaOfNode(draw_rand_index(zz),jj)=1;
%             end
            

            ff=0;
            for zz = 1:ReqTmp(jj)
                if(StaOfNode(draw_rand_index(zz))==0)
                    MissNum = MissNum+1;
                    StaOfNode(draw_rand_index(zz),jj) = 1;
                    ff=ff+1;
                else
                    HitNum = HitNum+1;
                    StaOfNode(draw_rand_index(zz),jj) = 3;
                    %Q_Outputs(jj,ii)=max{0,Q_Outputs(jj,ii)-1};
                end
            end
            if ii==1
                Q_Outputs(jj,ii)=ff;
            else
                Q_Outputs(jj,ii)=Q_Outputs(jj,ii-1)+ff;
            end
        else
            if ii==1
                Q_Outputs(:,ii)=0;
            else
                Q_Outputs(jj,ii)=Q_Outputs(jj,ii-1);
            end
        end
    end
    
    %% ��վ��D2D�������ݷַ�
    % �жϽڵ�λ���ĸ���վ�ڣ����� BScell ��
    BScell = cell(NumOfBS,1);
    for jj = 1:NumOfNode
        PosNX = PosOfNode(jj,1);
        PosNY = PosOfNode(jj,2);
        for zz = 1:NumOfBS
            PosBX = PosOfBS(zz,1);
            PosBY = PosOfBS(zz,2);
            DisTmp(zz) = sqrt((PosNX-PosBX)^2+(PosNY-PosBY)^2);
        end
        PosNB = find(min(DisTmp)==DisTmp);
        BScell(PosNB) = {[BScell{PosNB},jj]};
    end
    % ������Եõ������ֵĻ������
    for jj = IndOfBS:IndOfBS
        BSIndTp = BScell{jj};
        BSOfNode = StaOfNode(BSIndTp,:);
        NumOfBSNode = numel(BSOfNode(:,1));

        CahPoliy(ii,:) = floor(Cachingpolicy(BSOfNode,alpha,beta,NumOfNode,lambda));
        TMPCah = CahPoliy(ii,:);
        
        for zz = 1:NumOfCon
            cache1 = find(BSOfNode(:,zz)==2);
            cache2 = find(BSOfNode(:,zz)==3);
            Cache(zz) = numel(cache1) + numel(cache2);
        end
        
        for zz = 1:CapaOfBs
            Delta = TMPCah-Cache;
            MaxValue = max(Delta);
            MaxIndex = find(Delta==MaxValue);
            MaxIndex = MaxIndex(1);
            if(Delta(MaxIndex)>0)
                CanOfNode = find(StaOfNode(BSIndTp,MaxIndex)==1);
                if numel(CanOfNode)==0
                    CanOfNode = find(StaOfNode(BSIndTp,MaxIndex)==0);
                    if numel(CanOfNode)==0
                        TMPCah(MaxIndex) = TMPCah(MaxIndex) - 1;
                        continue;
                    end
                    ReInd = randi(numel(CanOfNode));
                    Index = BSIndTp(CanOfNode(ReInd));
                    %����
                    cache1 = find(StaOfNode(Index,:)==2);
                    cache2 = find(StaOfNode(Index,:)==3);
                    cacheTmp = [cache1,cache2];
                    CacheNum = numel(cacheTmp);                  
                    if(CacheNum>0)
                        if(CacheNum>=CahCap)
                            NumTmp = CacheNum-CahCap;
                            disInd = cacheTmp(randi(CacheNum,1,NumTmp));
                            for kk = 1:NumTmp
                                StaOfNode(Index,disInd(kk)) = 0;
                            end
                        end
                    end
                    StaOfNode(Index,MaxIndex) = 2;
                else
                    ReInd = randi(numel(CanOfNode));
                    Index = BSIndTp(CanOfNode(ReInd));
                    cache1 = find(StaOfNode(Index,:)==2);
                    cache2 = find(StaOfNode(Index,:)==3);
                    cacheTmp = [cache1,cache2];
                    CacheNum = numel(cacheTmp);
                    if(CacheNum>0)                        
                        if(CacheNum>=CahCap)
                            NumTmp = CacheNum-CahCap;
                            disInd = cacheTmp(randi(CacheNum,1,NumTmp));
                            for kk = 1:NumTmp
                                StaOfNode(Index,disInd(kk)) = 0;
                            end
                        end
                    end
                    HitNum = HitNum+1;
                    StaOfNode(Index,MaxIndex) = 3;
                    %Q_Outputs(MaxIndex,ii)=Q_Outputs(MaxIndex,ii)-1;
                end
                TMPCah(MaxIndex(1)) = TMPCah(MaxIndex(1)) - 1;
            end
        end
        for jj=1:NumOfCon
            %number of nodes been served at each time slot;
            serve_Queue=find(StaOfNode(:,jj)==1);
            serve_Node1=find(StaOfNode(:,jj)==3);
            serve_Node2=find(StaOfNode(:,jj)==2);
            serve_Node=[serve_Node1',serve_Node2'];
            serve_Cap=numel(serve_Node);
            SC=floor(serve_Cap*0.32);
            if SC<numel(serve_Queue)
                %SC=floor(serve_Cap*0.1);
                serve_set=serve_Queue(1:SC);
                for zz = 1:numel(serve_set)
                    Index = serve_set(zz);
                    cache1 = find(StaOfNode(Index,:)==2);
                    cache2 = find(StaOfNode(Index,:)==3);
                    cacheTmp = [cache1,cache2];
                    CacheNum = numel(cacheTmp);
                    if(CacheNum>0)                        
                        if(CacheNum>=CahCap)
                            NumTmp = CacheNum-CahCap;
                            disInd = cacheTmp(randi(CacheNum,1,NumTmp));
                            for kk = 1:NumTmp
                                StaOfNode(Index,disInd(kk)) = 0;
                            end
                        end
                    end
                    StaOfNode(Index,jj)=3;
                end
                Q_Outputs(jj,ii)=max(Q_Outputs(jj,ii)-SC,0);
                Q_delay(jj,ii)=sum(Q_Outputs(jj,ii))/SC;
            else
                for zz = 1:numel(serve_Queue)
                    Index = serve_Queue(zz);
                    cache1 = find(StaOfNode(Index,:)==2);
                    cache2 = find(StaOfNode(Index,:)==3);
                    cacheTmp = [cache1,cache2];
                    CacheNum = numel(cacheTmp);
                    if(CacheNum>0)                        
                        if(CacheNum>=CahCap)
                            NumTmp = CacheNum-CahCap;
                            disInd = cacheTmp(randi(CacheNum,1,NumTmp));
                            for kk = 1:NumTmp
                                StaOfNode(Index,disInd(kk)) = 0;
                            end
                        end
                    end
                    StaOfNode(Index,jj)=3;
                end
                Q_Outputs(jj,ii)=max(Q_Outputs(jj,ii)-numel(serve_Queue),0);
                Q_delay(jj,ii)=sum(Q_Outputs(jj,ii))/SC;
            end
        end
    end
    
    
    %% �û��������
    suscep = 0;      % �հ׽ڵ�
    infect = 0;      % ��Ȥ�ڵ�
    cached = 0;      % ����ڵ�
    satisf = 0;      % ����ڵ�
    totalSat = 0;
    totalInf = 0;
    exitRate = 0.2*rand();
    for jj = 1:NumOfCon
        ConSta = ConStaTime{jj};
        AllNode = BSOfNode(:,jj);
        suscep = sum(AllNode==0);
        infect = sum(AllNode==1);
        cached = sum(AllNode==2); 
        satisf = sum(AllNode==3);
        totalSat = totalSat + satisf;
        totalInf = totalInf + infect;
        ConStaTmp = [suscep,infect,cached,satisf];
        ConSta(ii,:) = ConStaTmp;
        ConStaTime(jj) = {ConSta};
    end
    % ��������������û��л�Ϊ����״̬
    StaOfNode(StaOfNode==3) = 2;
    HitRate(ii,1) = totalSat;
    HitRate(ii,2) = totalInf;
    HitRate(ii,3) = totalSat/(totalInf+totalSat);
    
    NumOfSat(ii) = totalSat;
    NumOfInf(ii) = totalInf;
    if ii ==1
        sumDelay(ii) = NumOfInf(ii);
    else
        sumDelay(ii) = NumOfInf(ii) + NumOfSat(ii)*rand() + max([sumDelay(ii-1)-NumOfSat(ii-1),0]);
        %sumDelay(ii) = NumOfInf(ii) + max([sumDelay(ii-1)-NumOfSat(ii-1),0]);
    end
    NumOfActCli(ii) = totalInf+totalSat;
end
    Ave_Hitrate=(t-1)\t*Ave_Hitrate+1\t*HitRate(:,3);
    Ave_Q_delay=(t-1)\t*Ave_Q_delay+1\t*mean(Q_delay,1);
    Delay = sumDelay'./NumOfActCli;
end