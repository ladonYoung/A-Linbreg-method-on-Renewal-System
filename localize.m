classdef localize<handle%A solving method for target localization
    properties(GetAccess=public,SetAccess=public)
        A%sensor's location:(n+1)*n matrix
        n%dimension
        x_0%the optimal point
        x_sgd
        x_lin
        d%distance,every step add (n+1) values
        d0%the true distance
        k%step
        L_sgd% cost of sgd
        L_lin%cost of linbreg
        itemaxiter
        lin_stepsize   
    end

    methods(Access=public)
        function obj=localize(x0,A0)
            obj.x_0=x0;
            obj.n=size(x0,2);
            obj.A=A0;
            obj.d=zeros(obj.n+1,1);
            v=obj.A-obj.x_0;
            for i=1:obj.n+1
              obj.d0(i)=norm(v(i,:),2);
            end
            obj.d0=obj.d0';
            obj.k=1;
            obj.itemaxiter=200;
            obj.x_sgd=obj.x_0+0.5;
            obj.x_lin=obj.x_0+0.5;
            obj.lin_stepsize=0.1;
          
        end
        function iterate(obj)               
            v=zeros(1,obj.n);
            while obj.k<=obj.itemaxiter
                   g_sgd=zeros(1,obj.n);
                   g_lin=zeros(1,obj.n);
                   w=0.01*randn(obj.n+1,1);
                   d_1=obj.d0+w;
                   obj.d=(obj.k*obj.d+d_1)/(obj.k+1);%avereage distance
                   for j=1:obj.n
                       for i=1:obj.n+1
                       g_sgd(j)=g_sgd(j)+2*(1-obj.d(i)/norm(obj.A(i,:)-obj.x_sgd(obj.k,:),2))*(obj.x_sgd(obj.k,j)-obj.A(i,j));
                       g_lin(j)=g_lin(j)+2*(1-obj.d(i)/norm(obj.A(i,:)-obj.x_lin(obj.k,:),2))*(obj.x_lin(obj.k,j)-obj.A(i,j));
                       end
                   end
                   if obj.n==2
                       sgd_stepsize=1/obj.k;
                   else 
                       sgd_stepsize=0.5/sqrt(obj.k);
                   end
                   obj.x_sgd(obj.k+1,:)=obj.x_sgd(obj.k,:)-sgd_stepsize*g_sgd;                
                    obj.x_lin(obj.k+1,:)=obj.prox(obj.lin_stepsize, obj.x_lin(obj.k,:)+obj.lin_stepsize*(v-g_lin));
                    v=v-(obj.x_lin(obj.k+1,:)-obj.x_lin(obj.k,:)+obj.lin_stepsize*g_lin)/obj.lin_stepsize;
                    obj.L_sgd=[obj.L_sgd,obj.loss(obj.x_sgd(obj.k+1,:))];
                    obj.L_lin=[obj.L_lin,obj.loss(obj.x_lin(obj.k+1,:))];
                    obj.k=obj.k+1;
            end
        end
   
    function y=loss(obj,x)
       y=0;
       for i=1:obj.n+1
           y=y+(norm(obj.A(i,:)-x,2)-obj.d0(i))^2;
       end
    end
    end
 
methods(Static)
   function u=prox(t,v)
            u=v/(t+1);%assume J is l2-norm and let the derivative be 0
   end
   
end

end




