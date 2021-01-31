function dU=ICME(t,U,flag,beta,psi)
% beta 流行度
% psi 缓存策略
% beta=0.005;
% psi=0.002;
dU=zeros(5,1);
dU(1)=0.09*U(5)-(beta+psi)*U(1)+0.008*U(3)-0.02*U(1);
dU(2)=beta*U(1)-(0.05*(U(3)+U(4)))*U(2)-0.02*U(2);
dU(3)=psi*U(1)-(0.008+beta+0.02)*U(3);
dU(4)=(0.05*(U(3)+U(4)))*U(2)+beta*U(3)-0.02*U(4);
dU(5)=0.02*(1-U(5))-U(5)*0.09;