%% Finite Difference method 
nx=150; %X direction
ny=100; %Y direction
nx1=nx/2;
G = sparse(nx*ny); %initializing sparse matrix
B = zeros(1,nx*ny); %initializing B matrix
cMap = ones(nx,ny); %Conductivity map, all ones no conductivity change
Vmap =zeros(nx,ny); %initalizing Voltage map

%%
%Creating of the G matrix and solving for voltages
for i=1:nx
    for j=1:ny
        n=j+(i-1)*ny;
        %left side
        if i==1 %Left side
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
        elseif i==nx %Right side
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1 %Bottom side
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nyp = j+1+(i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        
        elseif j==ny %Top side
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) +cMap(i,j-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) + cMap(i,j-1))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
    end
end
V=G\B';
%Reverting back to matrix form for plotting
for i=1:nx
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap(i,j) = V(n); 
    end
end
figure(1)
surf(Vmap)
title( 'Electrostatic potential solution');

%% Electrostatic potential saddle
%Creating of the G matrix and solving for voltages Saddle
Vmap1 =zeros(nx,ny); %initalizing Voltage map
for i=1:nx
    for j=1:ny
        n=j+(i-1)*ny;
        %left side
        if i==1 %Left side
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
        elseif i==nx %Right side
            G(n,:) = 0;
            G(n,n) = 1;
             B(n) = 1;
        elseif j==1 %Bottom side
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==ny %Top side
           G(n,:) = 0;
           G(n,n) = 1;
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2;
            rxp = (cMap(i,j) + cMap(i+1,j))/2;
            rym = (cMap(i,j) + cMap(i,j-1))/2;
            ryp = (cMap(i,j) + cMap(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
    end
end
V=G\B';
%Reverting back to matrix form for plotting
for i=1:nx
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap1(i,j) = V(n); 
    end
end
figure(2)
surf(Vmap1)
title( 'Electrostatic potential saddle');

%% Analytical Series Solution
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
            for n=1:300
                if (rem(n,2)~=0) %This make the sumation on every odd number
                vtemp(i,j) = vtemp(i,j)+(1/n)*(cosh(n*pi*x(i)/ny)/cosh(n*pi*nx1/ny))*sin(n*pi*j/ny);
                end
            end
               v(i,j)=4/pi*vtemp(i,j);
    end
end
figure(3)
surf(v)
title ('Analytical Series Solution Graph');
%%
%{
When looking at this solution it is very important to note that it is solved by convergences. This means the actual value intended is never 
actually reached. This takes more computation power and as seen in the figure below is that the end of the graph are curved meaning that the 
exact value needed is not actually correct but the convergence gives a good understanding on how the system works.

This method of solving the electrostatic problem is very easy to adjust the size and provides great base for more complex systems. The 
solution provide is a finite answer which is different than the convergence method. This could include complex conductivity maps that is 
explored in part 2 of this assignment. The major problem is that this system only provides the solution for rectangular shapes and different 
space of space with curve are not able to solved as accurately. The stair case effect is used in order to solve curved areas but this provide 
a certain aspect of error and affects the accuracy. 
%}
%% The difference
figure(4)
surf(abs(Vmap1-v))
title('Analytical VS FD')

%%


















