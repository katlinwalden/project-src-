%discretizing the domain
m=2*pi; %domain length
nx=200; %1-total nu1mber of x nodes
ny=200; %1-total nu1mber of y nodes
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
l=1; %predefine for functionalityones(nx+1,ny+1)
percenterror=ones(nx+1:ny+1);
s=ones(nx+1:ny+1);
while l >=.00001
    for j=1:ny+1 
        if j==1
            u1(:,j)=(m-x).^2.*x; %boundary condition
        elseif j==ny+1
                u1(:,j)=(m-x).^2.*cos(x/2); %boundary condition
        else
        end
        for i=1:nx+1
            if i==1
                u1(i,j)=m*y(j); %boundary condition
            elseif i==nx+1
                u1(i,j)=.5*(F(i,j)*dx^2*dy^2-u1(i,j-1)-u1(i,j+1)-2*u1(i-1,j)); %Neumann condition with Taylor series substitution f(x-h)=f(x+h)
            else
                u1(i,j)=.5*(F(i,j)*dx^2*dy^2-u1(i,j-1)-u1(i,j+1)-u1(i-1,j)-u1(i+1,j)); %defnining u1 matrix
            end
            percenterror=(u1-s)./u1;
            l=max(max(percenterror));
            s=u1;
        end   
    end
end
