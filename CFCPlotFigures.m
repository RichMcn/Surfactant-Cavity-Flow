function [] ...
    = CFCPlotFigures(x,y,psi,SurfactantProfile,psi_y,psi_x,SurfaceVelocity,StressProfile,Vorticity,SurfacePressure,Height)
%CFCPlotFigures Called by the main functions CavictyFlowCalculator which plots figures 

%Contour plot of the stream function
figure(1)
contours=[0,10^(-2),1.5*10^(-2),2*10^(-2),5*10^(-3),10^(-4),10^(-5),10^(-6),10^(-7),10^(-8),10^(-9),10^(-10),10^(-11),10^(-12)];
contours=[contours, -contours]';
contour(x,y,psi,contours)
xlabel("x");
ylabel("y");


%Plots surfactant concentration profile
figure(2)
plot(x,SurfactantProfile)
xlabel("x")
ylabel("\Gamma")


%Contour plot of the velocity fiels in the x direction
figure(3)
contour(x,y,psi_y)
xlabel("x")
ylabel("y")


%Contour plot of the velocity fiels in the y direction
figure(4)
contour(x,y,-psi_x)
xlabel("x")
ylabel("y")


%Plots the surface velocity
figure(5)
plot(x,SurfaceVelocity)
xlabel("x")
ylabel("u_s")



%Plots surface stress
figure(6)
plot(x,StressProfile(end,:),'b:')
hold on
plot(x(10:end-9),StressProfile(end,10:end-9),'b-')
xlabel("x")
ylabel("\tau")
hold off


%Contour plot of the vorticity field
figure(7)
contour(x,y,Vorticity);
xlabel("x")
ylabel("y")

%Plot of the surface pressure profile
figure(8)
plot(x,SurfacePressure)
xlabel("x")
ylabel("P(x,0)")


%Plot of the leading order height deflection
figure(9)
plot(x,Height)
xlabel("x")
ylabel("$\hat{h}$",'Interpreter','Latex')

end

