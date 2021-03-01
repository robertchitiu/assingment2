%%  Analytical Series Solution
nx=150;
ny=100;
nx1=nx/2;
vtemp =zeros(nx1*2,ny);
v=zeros(nx1*2,ny);
x= zeros(1,nx1*2);
y = linspace(0,ny,ny);
%Scaling the x direction for proper calculations
for i=1:nx1*2
    if i==1
    x(i)=-nx1;
    else
    x(i)=x(i-1)+1;   
    end
end
for i=1:nx1*2
    for j=1:ny
            for n=1:100
                if (rem(n,2)~=0) %This make the sumation on every odd number
                vtemp(i,j) = vtemp(i,j)+(1/n)*(cosh(n*pi*x(i)/ny)/cosh(n*pi*nx1/ny))*sin(n*pi*j/ny);
                end
            end
               v(i,j)=4/pi*vtemp(i,j);
    end
end
%{
When looking at this solution it is very important to note that it is solved by convergences. This means the actual value intended is never 
actually reached. This takes more computation power and as seen in the figure below is that the end of the graph are curved meaning that the 
exact value needed is not actually correct but the convergence gives a good understanding on how the system works.
%}
figure(1)
surf(v)
title ('Analytical Series Solution Graph');
%%