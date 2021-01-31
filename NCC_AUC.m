% NCC_AUC() -  AUC optimization model
%
% Usage:
%              >> NCC_AUC(S,group,lamda);
%              >> [feature_index] = NCC_AUC(S,group,lamda);
% Inputs:
%       S   = Expression of genes selected by the network propagation model
%       group   = Labels of genes selected by the network propagation model
%       lamda   = parameter used in the AUC optimization model
% Outputs:
%       feature_index   = Statistic scores of genes
function [feature_index] = NCC_AUC(S,group,lamda)
    n=size(S,2);
    Sp=S(group,:);
    Sn=S(~group,:);
    m1=size(Sp,1);
    m2=size(Sn,1);
    A1=[];
    mean_p=mean(Sp);
    mean_n=mean(Sn);
    %% Linear programming
    C=1/m1/m2;
    f=[lamda*ones(n,1);C*ones(m1*m2,1);zeros(m1*m2,1)];
    for i=1:size(Sn,1)
        A1 = [repmat(Sn(i,:).*(mean_p-mean_n),size(Sp,1),1);A1];
    end
    A2=repmat(Sp*diag(mean_p-mean_n),size(Sn,1),1);
    Aeq=[A2-A1 speye(m1*m2) -speye(m1*m2)];
    beq=ones(m1*m2,1);
    x=linprog(f,[],[],Aeq,beq,zeros(1,2*m1*m2+n),[]);
    auc=sum((1-x(n+1:m1*m2+n)+x(m1*m2+n+1:2*m1*m2+n))>0)/m1/m2;

    %% coefficient
    c1=x(1:n).*(mean(Sp)-mean(Sn))';
    theda=x(1:n);
    score=S*c1;
    %% identify panel of  features
    feature_index=find(abs(theda)>1e-5);
end
