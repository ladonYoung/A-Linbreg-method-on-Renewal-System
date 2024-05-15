function [y,u]=minff(q)%Caculate the min of multiple-systems' loss function(for plotting)
f=@(u)(q*u(2)+1-q)/(q*u(1)+1-q)+((1-q)*u(2)+q)/((1-q)*u(1)+q);
          function [c, ceq] = constraint3(u)
                    x = u(1);
                    y = u(2);
                    c = [-(y - (x - 1)^2 - 1); y - x];
                    ceq = [];
            end
            
                lb = [1; 1];  
                ub = [2; 2];           
                
                options = optimoptions('fmincon', 'Display', 'notify'); 
                u= fmincon(f, [1.5,1.5], [], [], [], [], lb, ub, @constraint3, options);  
                y=f(u);
end