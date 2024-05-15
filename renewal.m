classdef renewal<handle%Construct a renewal system
    properties(GetAccess=public,SetAccess=public)
        P%the probability vector
        n%number of systems
        m%size of P
        ss%record the value of r.v. S[k]
        S%record the number of S[k]
        u%decision, a (m,n,2)-matrix including cost and time
        L%record every step's the average cost(using P to caculate)
        all_u%record all decisions
    end
    methods(Access=public)
        function obj=renewal(P,n)
            if sum(P)~=1
                error("the sum of the probability vector is not 1")
            else
            obj.P=P;%row vector
            obj.m=size(P,2);
            obj.ss=[];
            obj.S=zeros(1,obj.m);
            obj.n=n;
            obj.u=2*ones(1,2);%represents the decision when S[k]=2;(The first column is time, the second is cost)
            %According to the specific problem in the paper, when n=2, we only need u is 2-dimention actually;
            obj.L=[];
            obj.all_u=[];
            end
            if obj.n==1
            obj.u=2*ones(1,2);
            end
        end
        
        function g=forward(obj)%return the observation of gradient for iteration
            s = randsample(1:obj.m, 1, true, obj.P);%generate S[k]
            obj.ss=s;
            obj.S(s)=obj.S(s)+1;%update S
            S_pro=obj.S/sum(obj.S);%caculate the estimator of probability vector
            g=zeros(size(obj.u));
            if obj.n==1
            ave=sum([[1 1];obj.u].*S_pro',1);%Caculate ave time and cost using estimator
            g(2)=S_pro(2)/ave(1);%the gradient of cost
            g(1)=-g(2)*ave(2)/ave(1);%the gradient of time
            elseif obj.n==2
                ave1=sum([[1 1];obj.u].*S_pro',1);
                ave2=sum([obj.u;[1 1]].*S_pro',1);%Caculate ave time and cost of two systems using estimator
                g(2)=S_pro(2)/ave1(1)+S_pro(1)/ave2(1);
                g(1)=-S_pro(2)*ave1(2)/(ave1(1)^2)-S_pro(1)*ave2(2)/(ave2(1)^2);
            end
        end
        
        function obj=update_L(obj)
            
            if obj.n==1
            loss_and_time=sum([[1 1];obj.u].*obj.P',1);%Caculate ave cost and time using true P
            loss_time=loss_and_time(2)/loss_and_time(1);
              if obj.ss==1
                  obj.all_u=[obj.all_u,[1,1]'];
              elseif obj.ss==2
                  obj.all_u=[obj.all_u,obj.u'];
              end
            else
            loss_and_time1=sum([[1 1];obj.u].*obj.P',1);
            loss_and_time2=sum([obj.u;[1 1]].*obj.P',1);
            loss_time=loss_and_time1(2)/loss_and_time1(1)+loss_and_time2(2)/loss_and_time2(1);
              if obj.ss==1
                  obj.all_u=[obj.all_u,[1,1,obj.u]'];
              elseif obj.ss==2
                  obj.all_u=[obj.all_u,[obj.u,1,1]'];
              end
            end
            obj.L=[obj.L,loss_time];%record the current average cost into L        
        end
        function obj=reinitial(obj)
            obj.S=zeros(1,obj.m);            
            obj.u=2*ones(1,2);
            obj.L=[];            
        end
    end
end