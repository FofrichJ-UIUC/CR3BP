clear all
close all

%% Selecting which point to perform analysis on 

point = 5

%% Calculating the location of the equilibrium points
m_earth = 5.972 *10^24;
m_moon = 7.347 * 10^22;
m1 = m_earth;
m2 = m_moon;
mstar = m1 + m2;
lstar = 382500; % https://www.nasa.gov/sites/default/files/files/Distance_to_the_Moon.pdf
G = 6.6743E-20;
tstar = (lstar^3/(G*mstar))^0.5
mu = m2/mstar;



syms x

%L1
eqn = x - (1-mu)/((x+mu)^2) + mu/((x-1+mu)^2) == 0;
sol = double(solve(eqn,x));
for j = 1:length(sol)
    if isreal(sol(j)) ==1
        L1=sol(j);
    end
end

%L2
eqn = x - (1-mu)/((x+mu)^2) - mu/((x-1+mu)^2) == 0;
sol = double(solve(eqn,x));
for j = 1:length(sol)
    if isreal(sol(j)) ==1
        L2=sol(j);
    end
end

%L3
eqn = x + (1-mu)/((x+mu)^2) + mu/((x-1+mu)^2) == 0;
sol = double(solve(eqn,x));
for j = 1:length(sol)
    if isreal(sol(j)) ==1
        L3=sol(j);
    end
end

%L4
L4x = 0.5-mu; L4y = sqrt(3)/2;

%L5
L5x = 0.5-mu; L5y = -sqrt(3)/2;


%% Calculating the stability of the equilbrium points

if point == 1
    y = [L1,0,0];
elseif point ==2
    y = [L2,0,0];
elseif point ==3
    y = [L3,0,0];
elseif point ==4
y = [L4x,L4y,0];
elseif point ==5
y = [L5x,L5y,0];
end


r1_cube = ((y(1)+mu)^2 + y(2)^2 + y(3)^2)^(3/2);
r2_cube = ((y(1)-1+mu)^2 + y(2)^2 + y(3)^2)^(3/2);

r1_fif = ((y(1)+mu)^2 + y(2)^2 + y(3)^2)^(5/2);
r2_fif = ((y(1)-1+mu)^2 + y(2)^2 + y(3)^2)^(5/2);

uzz = -(1-mu)/r1_cube - mu/r2_cube + 3*(1-mu)*(y(3)^2)/r1_fif + (3*mu*y(3)^2)/r2_fif;

oscmodepos = 1i*sqrt(abs(uzz)) %just calculating the positive version of it, but need the negative version of it too
oscmodeneg = -1i*sqrt(abs(uzz)) %just calculating the positive version of it, but need the negative version of it too


uxx = 1 - (1-mu)/r1_cube -mu/r2_cube + 3*(1-mu)*(mu+y(1))^2/r1_fif + 3*mu*(y(1)-1+mu)^2/r2_fif;
uyy = 1 - (1-mu)/r1_cube -mu/r2_cube + 3*(1-mu)*y(2)^2/r1_fif +3*mu*y(2)^2/r2_fif;
uxy = 3*(1-mu)*(y(1)+mu)*y(2)/r1_fif + 3*mu*(y(1)-1+mu)*y(2)/r2_fif;

%formulating the A matrix for the planar motion

A11 = zeros(2,2);
A12 = eye(2);
A21 = [uxx uxy;uxy uyy];
A22 = [0 2;-2 0];

A2d = [A11 A12;A21 A22];

[V,D] = eig(A2d) % columns of V correspond to the eigenvectors of the eigenvalues in the D matrix
