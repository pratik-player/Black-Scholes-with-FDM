% an approximation to the differential equation (x^2*sig^2*uxx)/2+rxux-ru=f
N = input("Please input the reciprocal of the step size: ");
sig = input("Please input the constant sigma: ");
r = input("Please input the constant r: ");

h=1/N; %step size
points = (h:h:1-h)'; %domain not including boundary points
main_diag = (-1*N.^2*((points).^2.*sig.^2)-r); %main diagnol, taking ones per diagnol
under_diag = (N.^2*(points(2:end)).^2*sig.^2)/2-(r*(points(2:end).*N)/2); %under
right_diag = (N.^2*(points(1:end-1)).^2*sig.^2)/2+(r*(points(1:end-1).*N)/2); %right
A = diag(main_diag, 0) + diag(under_diag, -1) + diag(right_diag, 1);
b=(points).^2.*sig.^2+r.*points.*(2.*(points)-1)-r.*(points.^2-points); %solution matrix b
disp(A);
disp(b);
x = A \ b; %solve for vector x the approximated function values 
disp(x); 

figure(1);
plot([0;points;1], [0;x;0]);
xlabel('x')
ylabel('u(x)')
title('Solution')

%Error calculation
u_exact = @(x) x.*(x-1); %exact function values at each x point
x_values = linspace(0, 1, N+1); 
u_values_exact = u_exact(x_values); 
x_refined = linspace(0, 1, 1000);
hold on
plot(x_refined, u_exact(x_refined));
hold off
figure(2);
plot(x_values, u_values_exact);
u_values_approx = x.';
err = abs(u_values_exact(2:end-1)-u_values_approx);
disp(vecnorm(err)) %vecnorm approximation
