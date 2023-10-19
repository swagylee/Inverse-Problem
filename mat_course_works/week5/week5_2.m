n = 111;

x = linspace(0,1.2,n);

f = zeros(1,length(x)).';

L = zeros(n-1,n);

l = 0;

v_temp = zeros(1,n);
v_temp(1:2) = [-1,1];

for ii = 1:length(x)
   if(x(ii)> 0.25) && (x(ii)<=0.35)
       f(ii) = 2;
   elseif (x(ii) >0.4) && (x(ii) <0.5)
       f(ii) = -1.5;
   elseif (x(ii)> 0.6)&&(x(ii) <0.7)
       f(ii) = 8*(2*x(ii) -1.2);
   elseif (x(ii) >0.75) &&(x(ii)<=1)
       f(ii) = sin(2*pi*x(ii))   ;
   end
end


for ii = 1:n-1
   L(ii,:) = v_temp(1:n);
   v_temp = circshift(v_temp,1); 
end


L2norm = norm(f)^2;

L1norm = sum(abs(f));

LF_l1norm = sum(abs(L*f));

LF_l2norm = norm(L*f)^2;

figure(1)
 clf
 plot(x,f,'k','LineWidth',1.5)
 hold on
 stem(x,f)

