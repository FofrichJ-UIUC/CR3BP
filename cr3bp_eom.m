function [y,tt] = cr3bp_eom(x,t,mu)
%cr3bp_eom: eom for numerical integrator for a cr3bp dynamical model
%  x is nondimensional
%  t is nondimen
%  mu is nondimensional mass ratio of second body divided by sum of system
%  mass
r1_cube = ((x(1)+mu)^2 + x(2)^2 + x(3)^2)^(3/2);
r2_cube = ((x(1)-1+mu)^2 + x(2)^2 + x(3)^2)^(3/2);

y(1) = x(4);
y(2) = x(5);
y(3) = x(6);
y(4) = 2*y(2) + x(1) - (1-mu)*(x(1)+mu)/r1_cube - mu(x(1) -1 +mu)/r2_cube;
y(5) = -2*x(4) +x(2) - (1-mu)*x(2)/r1_cube -mu*x(2)/r2_cube;
y(6) = -(1-mu)*x(3)/r1_cube - mu*x(3)/r2_cube;


end