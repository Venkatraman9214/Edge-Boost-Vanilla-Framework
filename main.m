clear all;
%% ��������
NumOfNode = 30;  % �ڵ�����
NumOfCon = 20;    % ��������
SimuTime = 50;   % ����ʱ��
minVec = 4;       % ��С�ٶ�
maxVec = 10;      % ����ٶ�
Square = 2000;    % ���η�Χ
radius = 500;     % BS�뾶
MaxArr = 400;      % ��󵽴�����
D2DRad = 150;      % D2Dͨ�ž���
CahCap = 60;
Q_Outputs=zeros(NumOfCon,SimuTime);
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
Delay = zeros(SimuTime,1);



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
                  Q_Outputs(jj,ii)=0;
                else
                  Q_Outputs(jj,ii)=Q_Outputs(jj,ii-1);
                end
                break;
            end
            if(numel(nTmp)<ReqTmp(jj))
                ReqTmp(jj) = numel(nTmp);
            end
            rand_index = nTmp(randperm(numel(nTmp)));                % ������������
            draw_rand_index = rand_index(1:ReqTmp(jj));    % ȡ��ǰReqTmp(jj)�����
            data_rand_index_set(jj,ii)=numel(draw_rand_index);
            
            
            % ͳ�ƶӳ�����
%             if ii==1
%                 Q_Outputs(jj,ii)=numel(draw_rand_index);
%             else
%                 Q_Outputs(jj,ii)=Q_Outputs(jj,ii-1)+numel(draw_rand_index);
%             end
            % ״̬ת�� �����޸ģ�
            for zz = 1:ReqTmp(jj)
                StaOfNode(draw_rand_index(zz),jj)=1;
            end
            
            %number of nodes been served at each time slot;
            

            
            
            
            for zz = 1:ReqTmp(jj)
                if(StaOfNode(draw_rand_index(zz))==0)
                    MissNum = MissNum+1;
                    StaOfNode(draw_rand_index(zz),jj) = 1;
                else
                    HitNum = HitNum+1;
                    StaOfNode(draw_rand_index(zz),jj) = 3;
                    Q_Outputs(jj,ii)=Q_Outputs(jj,ii);
                end
            end
        end
        
        % ͳ�ƶӳ�����
        
        if ii==1
            Q_Outputs(jj,ii)=numel(draw_rand_index);
        else
            Q_Outputs(jj,ii)=Q_Outputs(jj,ii-1)+numel(draw_rand_index);
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
                        disInd = cacheTmp(randi(CacheNum));
                        if(CacheNum>=CahCap)
                            StaOfNode(Index,disInd) = 0;
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
                        disInd = cacheTmp(randi(CacheNum));
                        if(CacheNum>=CahCap)
                            StaOfNode(Index,disInd) = 0;
                        end
                    end
                    HitNum = HitNum+1;
                    StaOfNode(Index,MaxIndex) = 3;
                    Q_Outputs(MaxIndex,ii)=Q_Outputs(MaxIndex,ii)-1;
                end
                TMPCah(MaxIndex(1)) = TMPCah(MaxIndex(1)) - 1;
            end
        end
    end
    
    for jj=1:NumOfCon
        serve_Queue=find(StaOfNode(:,jj)==1);
        
        serve_Node1=find(StaOfNode(:,jj)==3);
        serve_Node2=find(StaOfNode(:,jj)==2);
        serve_Node=[serve_Node1',serve_Node2'];
        serve_Cap=numel(serve_Node);
        SC=floor(serve_Cap*0.3);
        if SC<numel(serve_Queue)
            serve_set=serve_Queue(1:SC);
            StaOfNode(serve_set,jj)=3;
            
            Q_Outputs(jj,ii)=max(0,Q_Outputs(jj,ii)-SC);
        else
            StaOfNode(serve_Queue,jj)=3;
            Q_Outputs(jj,ii)=max(0,Q_Outputs(jj,ii)-SC);
        end
    end
    
    %% �û��������
    suscep = 0;      % �հ׽ڵ�
    infect = 0;      % ��Ȥ�ڵ�
    cached = 0;      % ����ڵ�
    satisf = 0;      % ����ڵ�
    for jj = 1:NumOfCon
        ConSta = ConStaTime{jj};
        AllNode = BSOfNode(:,jj);
        suscep = sum(AllNode==0);
        infect = sum(AllNode==1);
        cached = sum(AllNode==2); 
        satisf = sum(AllNode==3); 
        ConStaTmp = [suscep,infect,cached,satisf];
        ConSta(ii,:) = ConStaTmp;
        ConStaTime(jj) = {ConSta};
    end
    % ��������������û��л�Ϊ����״̬
    StaOfNode(StaOfNode==3) = 2;
    HitRate(ii,1) = satisf;
    HitRate(ii,2) = infect;
    HitRate(ii,3) = satisf/(infect+satisf);
    
    Delay(ii+1) = Delay(ii) + infect - satisf;
end