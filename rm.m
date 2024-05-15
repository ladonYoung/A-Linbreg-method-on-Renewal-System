classdef rm<handle%the original algorithm(Robbins-Monro) for Renewal system, only for single system.
    properties(GetAccess=public,SetAccess=public)
        r%renewal system
        theta
        stepsize
        itemaxiter
        k
    end
    methods(Access=public)
        function obj=rm(r)
            obj.r=r;
            obj.theta=1;
            obj.itemaxiter=200;
            obj.stepsize=1/3;
            obj.k=1;
        end
        function iterate(obj)
           while(obj.k<=obj.itemaxiter)
               s=randsample(1:obj.r.m, 1, true, obj.r.P);
               if s==1
                   if obj.k==1
                       obj.r.all_u=[1,1]';
                       obj.r.L=[obj.r.L,1];
                   else
                   obj.r.L=[obj.r.L,obj.r.L(size(obj.r.L,2))];%when s=1, u stay the same as the last step
                   obj.r.all_u=[obj.r.all_u,[1,1]'];
                   end
                   obj.theta=max(min(obj.theta+obj.stepsize*(1-obj.theta),1),0.8);%2*sqrt(2)-2 is almost 0.8
                   obj.k=obj.k+1;
                   obj.stepsize=1/(obj.k+2);%stepsize is 1/(k+2)
               elseif s==2
                   obj.r.u=obj.rm_min(obj.theta);%find the min of C[k]-theta[k]T[k]
                   obj.r.all_u=[obj.r.all_u,obj.r.u'];
                   obj.r=obj.r.update_L;
                   obj.theta=max(min(obj.theta+obj.stepsize*(obj.r.u(2)-obj.theta*obj.r.u(1)),1),0.8);
                   obj.k=obj.k+1;
                   obj.stepsize=1/(obj.k+2);
               end
           end
        end
        function obj=reinitial(obj)
            obj.r=obj.r.reinitial;
            obj.itemaxiter=200;
            obj.theta=1;
            obj.k=1;
            obj.stepsize=1/3;
        end
        function setitemaxiter(obj,max)
            obj.itemaxiter=max;            
        end
        function max=getitemaxiter(obj)
            max=obj.itemaxiter;
        end
    end
    methods(Static)
        function u=rm_min(theta)
            f=@(u)u(2)-theta*u(1);
            function [c, ceq] = constraint2(u)
                    x = u(1);
                    y = u(2);
                    c = [-(y - (x - 1)^2 - 1); y - x];
                    ceq = [];
            end            
                lb = [1; 1];  
                ub = [2; 2];                    
                options = optimoptions('fmincon', 'Display', 'notify'); 
                u= fmincon(f, [1.5,1.5], [], [], [], [], lb, ub, @constraint2, options);            
        end
    end
end