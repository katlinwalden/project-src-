%discretizing the domain
m=2*pi; %domain length
nx=200; %1-total nu1mber of x nodes
ny=400; %1-total nu1mber of y nodes
dx=m/nx; %gap between x nodes
dy=m/ny; %gap between y nodes
x=0:dx:m; %x nodes
y=0:dy:m; %y nodes
%GAUSS-SEIDEL portion
    %making initial arrays for my linear system
F=zeros(nx+1,ny+1); %preallocate
u1=zeros(nx+1,ny+1); %preallocate
for i=1:nx+1
    v=cos(x(i)/4+1); %stride 1 access optimization
    for j=1:ny+1
        F(i,j)=v*sin(y(j)/2); %right hand side matrix w/ stride 1 access optimization
    end
end
G=transpose(F); %saved in storage for optimizing loops later
    %running actual iterations
percenterror=1; %predefine for functionality
s=ones(nx+1,ny+1); %predefine for functionality
while percenterror >=.01
    for j=1:ny+1
        for i=1:nx+1
            if i==1
                u1(i,j)=m*y(j); %boundary condition
            elseif j==1
                u1(i,j)=(m-x(i))^2*x(i); %boundary condition
            elseif j==ny+1
                u1(i,j)=(m-x(i))^2*cos(x(i)/2); %boundary condition
            elseif i==nx+1
                u1(i,j)=.5*(F(i,j)*dx^2*dy^2-u1(i,j-1)-u1(i,j+1)-2*u1(i-1,j));
            else
                u1(i,j)=.5*(F(i,j)*dx^2*dy^2-u1(i,j-1)-u1(i,j+1)-u1(i-1,j)-u1(i+1,j)); %defnining u1 matrix
            end
        end
    end
    percenterror=max(max(abs((u1-s)/u1)));
    s=u1;
end
