classdef linbreg_renew <handle%Linbreg for renewal system

    properties(GetAccess=public,SetAccess=public)
        r%renewal class
        J%regularizer，input=1 or 2 to set the regularizer l1/2-norm 
        tol
        itemaxiter %the maximal iteration
        stepsize
        k%constant step        
    end

    methods(Access=public)
        function obj = linbreg_renew(r,J)
            obj.r=r;
            obj.J=J;            
            obj.tol=1e-3;
            obj.itemaxiter=200;
            obj.stepsize=0.8;%default 0.8 const
            obj.k=1;
        end

        function iterate(obj)
            v=zeros(size(obj.r.u));
            u=obj.r.u;
            while (obj.k <= obj.itemaxiter) %&& ((obj.sens > ...
                %obj.tol) || (obj.k < 3) )
                u_pre=u;                
                g=obj.r.forward;    
                %if obj.k==1||round(log(obj.k)/log(3))==log(obj.k)/log(3)
                uu=u+obj.stepsize*(v-g);
                u=obj.prox(obj.J,  obj.stepsize,  uu);                
                v=v-(u-u_pre+obj.stepsize*g)/obj.stepsize;
                obj.r.u=u;
                %end
                obj.r=obj.r.update_L;%record the average cost, for plots
                %obj.sens=norm(u-u_pre);
                obj.k=obj.k+1;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods(Static)
        function u=prox(J,t,v)
            if J==1
                fun=@(u)0.5 * norm(u - v)^2 + t*norm(u, 1);
            elseif J==2
                fun=@(u)0.5 * norm(u - v)^2+t*norm(u,2);
            end
            function [c, ceq] = constraint1(u)
                    x = u(1);
                    y = u(2);
                    c = [-(y - (x - 1)^2 - 1); y - x];
                    ceq = [];
            end            
                lb = [1; 1];  % lower bound in fmincon
                ub = [2; 2];  % upper bound in fmincon           
                
                options = optimoptions('fmincon', 'Display', 'notify'); % set iteration info
                u= fmincon(fun, [1.5,1.5], [], [], [], [], lb, ub, @constraint1, options);            
            end
        end
        
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      methods(Access=public)
        function settol(obj,tol)
            obj.tol=tol;
        end
        function tol=gettol(obj)
            tol=obj.tol;
        end
        function setitemaxiter(obj,max)
            obj.itemaxiter=max;          
        end
        function max=getitemaxiter(obj)
            max=obj.itemaxiter;
        end
        function setstepsize(obj,step_size)
            obj.stepsize=step_size;
        end
        function step_size=getstepsize(obj)
            step_size=obj.stepsize;
        end
        function obj=reinitial(obj)
            obj.r=obj.r.reinitial;
            obj.tol=1e-3;
            obj.itemaxiter=200;
            obj.stepsize=0.8;%redefault 0.8
            obj.k=1;%the current step            
        end
    end
end