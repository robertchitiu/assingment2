%% Finite Difference Method with Boxes
%The problem with this program is that functions were not used this makes
%things messy and complicated so in the future in the next assingments this
%program will be changed making all the componets into fuctions in order to
%call them and used them more efficently 
nx=150;
ny=100;
G = sparse(nx*ny);
B = zeros(1,nx*ny);
%% Changing conductivity in regions changes
cMap = ones(nx,ny);
for i=50:100 
    for j=1:ny
        if j>=1 && j<=40
            cMap(i,j)=10^-2;
        elseif j>=60 && j<=100
            cMap(i,j) = 10^-2;
        end
    end
end
cMap1 = ones(nx,ny);
for i=10:40 
    for j=1:ny
        if j>=1 && j<=50
            cMap1(i,j)=10^-2;
        elseif j>=60 && j<=100
            cMap1(i,j) = 10^-2;
        end
    end
end
cMap2 = ones(nx,ny);
for i=20:100 
    for j=1:ny
        if j>=1 && j<=40
            cMap2(i,j)=10^-2;
        elseif j>=60 && j<=100
            cMap2(i,j) = 10^-2;
        end
    end
end
cMap3 = ones(nx,ny);
for i=20:100 
    for j=1:ny
        if j>=1 && j<=20
            cMap3(i,j)=10^-2;
        elseif j>=40 && j<=100
            cMap3(i,j) = 10^-2;
        end
    end
end
%% Making different G matrix for area changes
for i=1:nx %Creating G Matrix
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
Vmap =zeros(nx,ny);
for i=1:nx %%Converting V to matrix to plot
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap(i,j) = V(n); 
    end
end
for i=1:nx %Electric Field calculation
    for j=1:ny
        if i==1
            Ey(i,j) = (Vmap(i+1,j) - Vmap(i,j));
        elseif i==nx
            Ey(i,j) = (Vmap(i,j) - Vmap(i-1,j));
        else
            Ey(i,j) = (Vmap(i+1,j) - Vmap(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex(i,j) = (Vmap(i,j+1) - Vmap(i,j));
        elseif j==ny
            Ex(i,j) = (Vmap(i,j) - Vmap(i,j-1));
        else
            Ex(i,j) = (Vmap(i,j+1) - Vmap(i,j-1)) * 0.5;
        end 
        
    end
end
Ex = -Ex;
Ey = -Ey;
FlowX = cMap .* Ex;
FlowY = cMap .* Ey;
figure(2)
surf(Vmap)
title('Potential Map')
%new g matrix size
nx1=375;
ny1=250;
G1 = sparse(nx1*ny1);
B1 = zeros(1,nx1*ny1);
cMap11 = ones(nx1,ny1);
for i=1:nx1 %Creating G Matrix
    for j=1:ny1
        n=j+(i-1)*ny1;
        %left side
        if i==1 %Left side
            G1(n,:) = 0;
            G1(n,n) = 1;
            B1(n) = 1;
        elseif i==nx1 %Right side
            G1(n,:) = 0;
            G1(n,n) = 1;
        elseif j==1 %Bottom side
            nxm = j+(i-2)*ny1;
            nxp = j+i*ny1;
            nyp = j+1+(i-1)*ny1;
            
            rxm = (cMap11(i,j) + cMap11(i-1,j))/2;
            rxp = (cMap11(i,j) + cMap11(i+1,j))/2;
            ryp = (cMap11(i,j) + cMap11(i,j+1))/2;
            
            G1(n,n) = -(rxm+rxp+ryp);
            G1(n,nxm) = rxm;
            G1(n,nxp) = rxp;
            G1(n,nyp) = ryp;
        
        elseif j==ny1 %Top side
            nxm = j+(i-2)*ny1;
            nxp = j+i*ny1;
            nym = j-1+(i-1)*ny1;
            
            rxm = (cMap11(i,j) + cMap11(i-1,j))/2;
            rxp = (cMap11(i,j) + cMap11(i+1,j))/2;
            rym = (cMap11(i,j) +cMap11(i,j-1))/2;
            
            G1(n,n) = -(rxm+rxp+rym);
            G1(n,nxm) = rxm;
            G1(n,nxp) = rxp;
            G1(n,nym) = rym;
        else
            nxm = j + (i-2)*ny1;
            nxp = j + i*ny1;
            nym = j-1 + (i-1)*ny1;
            nyp = j+1 + (i-1)*ny1;
            
            rxm = (cMap11(i,j) + cMap11(i-1,j))/2;
            rxp = (cMap11(i,j) + cMap11(i+1,j))/2;
            rym = (cMap11(i,j) + cMap11(i,j-1))/2;
            ryp = (cMap11(i,j) + cMap11(i,j+1))/2;
            
            G1(n,n) = -(rxm+rxp+rym+ryp);
            G1(n,nxm) = rxm;
            G1(n,nxp) = rxp;
            G1(n,nym) = rym;
            G1(n,nyp) = ryp;
            
        end
    end
end
V1=G1\B1';
Vmap1 = zeros(nx1,ny1);
for i=1:nx1 %%Converting V to matrix to plot
    for j=1:ny1
        n = j+(i-1)*ny1;
        Vmap1(i,j) = V1(n); 
    end
end
Ex1 = zeros(nx1,ny1);
Ey1 = zeros(nx1,ny1);
for i=1:nx1 %Electric Field calculation
    for j=1:ny1
        if i==1
            Ey1(i,j) = (Vmap1(i+1,j) - Vmap1(i,j));
        elseif i==nx1
            Ey1(i,j) = (Vmap1(i,j) - Vmap1(i-1,j));
        else
            Ey1(i,j) = (Vmap1(i+1,j) - Vmap1(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex1(i,j) = (Vmap1(i,j+1) - Vmap1(i,j));
        elseif j==ny1
            Ex1(i,j) = (Vmap1(i,j) - Vmap1(i,j-1));
        else
            Ex1(i,j) = (Vmap1(i,j+1) - Vmap1(i,j-1)) * 0.5;
        end 
        
    end
end
Ex1 = -Ex1;
Ey1 = -Ey1;
FlowX1 = cMap11 .* Ex1;
FlowY1 = cMap11 .* Ey1;
%new g matrix size
nx2=75;
ny2=50;
G2 = sparse(nx2*ny2);
B2 = zeros(1,nx2*ny2);
cMap12 = ones(nx2,ny2);
for i=1:nx2 %Creating G Matrix
    for j=1:ny2
        n=j+(i-1)*ny2;
        %left side
        if i==1 %Left side
            G2(n,:) = 0;
            G2(n,n) = 1;
            B2(n) = 1;
        elseif i==nx2 %Right side
            G2(n,:) = 0;
            G2(n,n) = 1;
        elseif j==1 %Bottom side
            nxm = j+(i-2)*ny2;
            nxp = j+i*ny2;
            nyp = j+1+(i-1)*ny2;
            
            rxm = (cMap12(i,j) + cMap12(i-1,j))/2;
            rxp = (cMap12(i,j) + cMap12(i+1,j))/2;
            ryp = (cMap12(i,j) + cMap12(i,j+1))/2;
            
            G2(n,n) = -(rxm+rxp+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nyp) = ryp;
        
        elseif j==ny2 %Top side
            nxm = j+(i-2)*ny2;
            nxp = j+i*ny2;
            nym = j-1+(i-1)*ny2;
            
            rxm = (cMap12(i,j) + cMap12(i-1,j))/2;
            rxp = (cMap12(i,j) + cMap12(i+1,j))/2;
            rym = (cMap12(i,j) +cMap12(i,j-1))/2;
            
            G2(n,n) = -(rxm+rxp+rym);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
        else
            nxm = j + (i-2)*ny2;
            nxp = j + i*ny2;
            nym = j-1 + (i-1)*ny2;
            nyp = j+1 + (i-1)*ny2;
            
            rxm = (cMap12(i,j) + cMap12(i-1,j))/2;
            rxp = (cMap12(i,j) + cMap12(i+1,j))/2;
            rym = (cMap12(i,j) + cMap12(i,j-1))/2;
            ryp = (cMap12(i,j) + cMap12(i,j+1))/2;
            
            G2(n,n) = -(rxm+rxp+rym+ryp);
            G2(n,nxm) = rxm;
            G2(n,nxp) = rxp;
            G2(n,nym) = rym;
            G2(n,nyp) = ryp;
        end
    end
end
V2=G2\B2';
Vmap2 =zeros(nx2,ny2);
for i=1:nx2 %%Converting V to matrix to plot
    for j=1:ny2
        n = j+(i-1)*ny2;
        Vmap2(i,j) = V2(n); 
    end
end
Ex2 = zeros(nx2,ny2);
Ey2 = zeros(nx2,ny2);
for i=1:nx2 %Electric Field calculation
    for j=1:ny2
        if i==1
            Ey2(i,j) = (Vmap2(i+1,j) - Vmap2(i,j));
        elseif i==nx2
            Ey2(i,j) = (Vmap2(i,j) - Vmap2(i-1,j));
        else
            Ey2(i,j) = (Vmap2(i+1,j) - Vmap2(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex2(i,j) = (Vmap2(i,j+1) - Vmap2(i,j));
        elseif j==ny2
            Ex2(i,j) = (Vmap2(i,j) - Vmap2(i,j-1));
        else
            Ex2(i,j) = (Vmap2(i,j+1) - Vmap2(i,j-1)) * 0.5;
        end 
        
    end
end
Ex2 = -Ex2;
Ey2 = -Ey2;
FlowX2 = cMap12 .* Ex2;
FlowY2 = cMap12 .* Ey2;
%new g matrix size
nx3=300;
ny3=200;
G3 = sparse(nx3*ny3);
B3 = zeros(1,nx3*ny3);
cMap13 = ones(nx3,ny3);
for i=1:nx3 %Creating G Matrix
    for j=1:ny3
        n=j+(i-1)*ny3;
        %left side
        if i==1 %Left side
            G3(n,:) = 0;
            G3(n,n) = 1;
            B3(n) = 1;
        elseif i==nx3 %Right side
            G3(n,:) = 0;
            G3(n,n) = 1;
        elseif j==1 %Bottom side
            nxm = j+(i-2)*ny3;
            nxp = j+i*ny3;
            nyp = j+1+(i-1)*ny3;
            
            rxm = (cMap13(i,j) + cMap13(i-1,j))/2;
            rxp = (cMap13(i,j) + cMap13(i+1,j))/2;
            ryp = (cMap13(i,j) + cMap13(i,j+1))/2;
            
            G3(n,n) = -(rxm+rxp+ryp);
            G3(n,nxm) = rxm;
            G3(n,nxp) = rxp;
            G3(n,nyp) = ryp;
        
        elseif j==ny3 %Top side
            nxm = j+(i-2)*ny3;
            nxp = j+i*ny3;
            nym = j-1+(i-1)*ny3;
            
            rxm = (cMap13(i,j) + cMap13(i-1,j))/2;
            rxp = (cMap13(i,j) + cMap13(i+1,j))/2;
            rym = (cMap13(i,j) +cMap13(i,j-1))/2;
            
            G3(n,n) = -(rxm+rxp+rym);
            G3(n,nxm) = rxm;
            G3(n,nxp) = rxp;
            G3(n,nym) = rym;
        else
            nxm = j + (i-2)*ny3;
            nxp = j + i*ny3;
            nym = j-1 + (i-1)*ny3;
            nyp = j+1 + (i-1)*ny3;
            
            rxm = (cMap13(i,j) + cMap13(i-1,j))/2;
            rxp = (cMap13(i,j) + cMap13(i+1,j))/2;
            rym = (cMap13(i,j) + cMap13(i,j-1))/2;
            ryp = (cMap13(i,j) + cMap13(i,j+1))/2;
            
            G3(n,n) = -(rxm+rxp+rym+ryp);
            G3(n,nxm) = rxm;
            G3(n,nxp) = rxp;
            G3(n,nym) = rym;
            G3(n,nyp) = ryp;
            
        end
    end
end
V3=G3\B3';
Vmap3 =zeros(nx3,ny3);
for i=1:nx2 %%Converting V to matrix to plot
    for j=1:ny3
        n = j+(i-1)*ny2;
        Vmap3(i,j) = V3(n); 
    end
end
Ex3 = zeros(nx3,ny3);
Ey3 = zeros(nx3,ny3);
for i=1:nx3 %Electric Field calculation
    for j=1:ny3
        if i==1
            Ey3(i,j) = (Vmap3(i+1,j) - Vmap3(i,j));
        elseif i==nx3
            Ey3(i,j) = (Vmap3(i,j) - Vmap3(i-1,j));
        else
            Ey3(i,j) = (Vmap3(i+1,j) - Vmap3(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex3(i,j) = (Vmap3(i,j+1) - Vmap3(i,j));
        elseif j==ny3
            Ex3(i,j) = (Vmap3(i,j) - Vmap3(i,j-1));
        else
            Ex3(i,j) = (Vmap3(i,j+1) - Vmap3(i,j-1)) * 0.5;
        end 
        
    end
end
Ex3 = -Ex3;
Ey3 = -Ey3;
FlowX3 = cMap13 .* Ex3;
FlowY3 = cMap13 .* Ey3;

%% various bottle neck effects
for i=1:nx %Creating G Matrix
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
            
            rxm = (cMap1(i,j) + cMap1(i-1,j))/2;
            rxp = (cMap1(i,j) + cMap1(i+1,j))/2;
            ryp = (cMap1(i,j) + cMap1(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        
        elseif j==ny %Top side
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (cMap1(i,j) + cMap1(i-1,j))/2;
            rxp = (cMap1(i,j) + cMap1(i+1,j))/2;
            rym = (cMap1(i,j) +cMap1(i,j-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap1(i,j) + cMap1(i-1,j))/2;
            rxp = (cMap1(i,j) + cMap1(i+1,j))/2;
            rym = (cMap1(i,j) + cMap1(i,j-1))/2;
            ryp = (cMap1(i,j) + cMap1(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
    end
end
V=G\B';
Vmap11 =zeros(nx,ny);
for i=1:nx %%Converting V to matrix to plot
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap11(i,j) = V(n); 
    end
end
for i=1:nx %Electric Field calculation
    for j=1:ny
        if i==1
            Ey11(i,j) = (Vmap11(i+1,j) - Vmap11(i,j));
        elseif i==nx
            Ey11(i,j) = (Vmap11(i,j) - Vmap11(i-1,j));
        else
            Ey11(i,j) = (Vmap11(i+1,j) - Vmap11(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex11(i,j) = (Vmap11(i,j+1) - Vmap11(i,j));
        elseif j==ny
            Ex11(i,j) = (Vmap11(i,j) - Vmap11(i,j-1));
        else
            Ex11(i,j) = (Vmap11(i,j+1) - Vmap11(i,j-1)) * 0.5;
        end 
        
    end
end
Ex11 = -Ex11;
Ey11 = -Ey11;
FlowX11 = cMap1 .* Ex11;
FlowY11 = cMap1 .* Ey11;
for i=1:nx %Creating G Matrix
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
            
            rxm = (cMap2(i,j) + cMap2(i-1,j))/2;
            rxp = (cMap2(i,j) + cMap2(i+1,j))/2;
            ryp = (cMap2(i,j) + cMap2(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        
        elseif j==ny %Top side
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (cMap2(i,j) + cMap2(i-1,j))/2;
            rxp = (cMap2(i,j) + cMap2(i+1,j))/2;
            rym = (cMap2(i,j) +cMap2(i,j-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap2(i,j) + cMap2(i-1,j))/2;
            rxp = (cMap2(i,j) + cMap2(i+1,j))/2;
            rym = (cMap2(i,j) + cMap2(i,j-1))/2;
            ryp = (cMap2(i,j) + cMap2(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
    end
end
V=G\B';
Vmap12 =zeros(nx,ny);
for i=1:nx %%Converting V to matrix to plot
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap12(i,j) = V(n); 
    end
end
for i=1:nx %Electric Field calculation
    for j=1:ny
        if i==1
            Ey12(i,j) = (Vmap12(i+1,j) - Vmap12(i,j));
        elseif i==nx
            Ey12(i,j) = (Vmap12(i,j) - Vmap12(i-1,j));
        else
            Ey12(i,j) = (Vmap12(i+1,j) - Vmap12(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex12(i,j) = (Vmap12(i,j+1) - Vmap12(i,j));
        elseif j==ny
            Ex12(i,j) = (Vmap12(i,j) - Vmap12(i,j-1));
        else
            Ex12(i,j) = (Vmap12(i,j+1) - Vmap12(i,j-1)) * 0.5;
        end 
        
    end
end
Ex12 = -Ex12;
Ey12 = -Ey12;
FlowX12 = cMap2 .* Ex12;
FlowY12 = cMap2 .* Ey12;
for i=1:nx %Creating G Matrix
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
            
            rxm = (cMap3(i,j) + cMap3(i-1,j))/2;
            rxp = (cMap3(i,j) + cMap3(i+1,j))/2;
            ryp = (cMap3(i,j) + cMap3(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        
        elseif j==ny %Top side
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (cMap3(i,j) + cMap3(i-1,j))/2;
            rxp = (cMap3(i,j) + cMap3(i+1,j))/2;
            rym = (cMap3(i,j) +cMap3(i,j-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap3(i,j) + cMap3(i-1,j))/2;
            rxp = (cMap3(i,j) + cMap3(i+1,j))/2;
            rym = (cMap3(i,j) + cMap3(i,j-1))/2;
            ryp = (cMap3(i,j) + cMap3(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
    end
end
V=G\B';
Vmap =zeros(nx,ny);
for i=1:nx %%Converting V to matrix to plot
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap13(i,j) = V(n); 
    end
end
for i=1:nx %Electric Field calculation
    for j=1:ny
        if i==1
            Ey13(i,j) = (Vmap13(i+1,j) - Vmap13(i,j));
        elseif i==nx
            Ey13(i,j) = (Vmap13(i,j) - Vmap13(i-1,j));
        else
            Ey13(i,j) = (Vmap13(i+1,j) - Vmap13(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex13(i,j) = (Vmap13(i,j+1) - Vmap13(i,j));
        elseif j==ny
            Ex13(i,j) = (Vmap13(i,j) - Vmap13(i,j-1));
        else
            Ex13(i,j) = (Vmap13(i,j+1) - Vmap13(i,j-1)) * 0.5;
        end 
        
    end
end
Ex13 = -Ex13;
Ey13 = -Ey13;
FlowX13 = cMap3 .* Ex13;
FlowY13 = cMap3 .* Ey13;
%% Conductivity changes
cMap111 = ones(nx,ny);
for i=50:100 
    for j=1:ny
        if j>=1 && j<=40
            cMap111(i,j)=0.05;
        elseif j>=60 && j<=100
            cMap111(i,j) = 0.05;
        end
    end
end
for i=1:nx %Creating G Matrix
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
            
            rxm = (cMap111(i,j) + cMap111(i-1,j))/2;
            rxp = (cMap111(i,j) + cMap111(i+1,j))/2;
            ryp = (cMap111(i,j) + cMap111(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
        
        elseif j==ny %Top side
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            
            rxm = (cMap111(i,j) + cMap111(i-1,j))/2;
            rxp = (cMap111(i,j) + cMap111(i+1,j))/2;
            rym = (cMap111(i,j) +cMap111(i,j-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap111(i,j) + cMap111(i-1,j))/2;
            rxp = (cMap111(i,j) + cMap111(i+1,j))/2;
            rym = (cMap111(i,j) + cMap111(i,j-1))/2;
            ryp = (cMap111(i,j) + cMap111(i,j+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            
        end
    end
end
V111=G\B';
Vmap111 =zeros(nx,ny);
for i=1:nx %%Converting V to matrix to plot
    for j=1:ny
        n = j+(i-1)*ny;
        Vmap111(i,j) = V111(n); 
    end
end
for i=1:nx %Electric Field calculation
    for j=1:ny
        if i==1
            Ey111(i,j) = (Vmap111(i+1,j) - Vmap111(i,j));
        elseif i==nx
            Ey111(i,j) = (Vmap111(i,j) - Vmap111(i-1,j));
        else
            Ey111(i,j) = (Vmap111(i+1,j) - Vmap111(i-1,j)) * 0.5;
        end 
        
        if j==1
            Ex111(i,j) = (Vmap111(i,j+1) - Vmap111(i,j));
        elseif j==ny
            Ex111(i,j) = (Vmap111(i,j) - Vmap111(i,j-1));
        else
            Ex111(i,j) = (Vmap111(i,j+1) - Vmap111(i,j-1)) * 0.5;
        end 
        
    end
end
Ex111 = -Ex111;
Ey111 = -Ey111;
FlowX111 = cMap111 .* Ex111;
FlowY111 = cMap111 .* Ey111;

%% Plotting graphs
figure (1)
surf(cMap)
title('Conductivity Map used for solution')
figure (3)
quiver(Ex,Ey)
title ('Electric Field')
figure (4)
quiver(FlowX,FlowY)
title ('Current Flow')

%{
This plots muliple size of area and looks at the current. The current is
always pointing upwards and the scales change the density seen. Since
subplot 2 has the biggest area it has the most current componets and that
is why it look so dense comparative to subplot 2 and 3
%}

figure (5)
subplot(2,2,1)
quiver(FlowX,FlowY)
title ('Current Flow used for calculation')
subplot(2,2,2)
quiver(FlowX1,FlowY1,0.1)
title ('Current Flow')
subplot(2,2,3)
quiver(FlowX2,FlowY2,0.1)
title ('Current Flow')
subplot(2,2,4)
quiver(FlowX2,FlowY2,0.1)
title ('Current Flow')
%{
These plots show the different types of bottle neck and the current
reaction to those. It is noticable that the current concentrates in the
conductive regions on the area.
%}
figure (6)
subplot (2,2,1)
quiver(FlowX,FlowY)
title('CurrentBottle Neck used for solution')
subplot(2,2,2)
quiver(FlowX11,FlowY11)
title('Current Different Bottle Neck')
subplot(2,2,3)
quiver(FlowX12,FlowY12)
title('Current Different Bottle Neck')
subplot(2,2,4)
quiver(FlowX13,FlowY13)
title('Current Different Bottle Neck')
%{
These 2 graphs are almost identical execpt that subplot 2 have a high
conductivity in the insulating region. This mean more current can flow
through and it is seen to be true. This is expects that little current
flows throught the insulating region but some still does so the current
arrows in subplot 1 and very small and in subplot 2 they are slightly
bigger
%}
figure (7)
subplot (2,1,1)
quiver(FlowX,FlowY)
title('Conductiviy current used for solution')
subplot(2,1,2)
quiver(FlowX111,FlowY111)


