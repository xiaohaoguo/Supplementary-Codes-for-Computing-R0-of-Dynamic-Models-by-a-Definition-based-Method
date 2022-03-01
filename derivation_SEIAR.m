close all; clear all; clc;

%% Define all symbolic variables and parameters involved in SEAIR model
syms kappa beta S N p omega1 omega2 gamma1 gamma2 f E A I d_r;  % define symbolic variables

%% Define infected variables
variables = [E; A; I]; % all infected cells
dim = numel(variables);  
F = sym(zeros(dim,1));  % initialize the vector of newly infections
V = sym(zeros(dim,1));  % initialize the vector of transitions (exluding newly infections)

% Define vector F and V
% Note that the original ODEs can be written as  dx/dt = F - V
F(1) = beta*S*(I+kappa*A);
V(1) = p*omega2*E + (1-p)*omega1*E + d_r*E;
V(2) = -p*omega2*E + gamma2*A + d_r*A;
V(3) = -(1-p)*omega1*E + (f + gamma1)*I + d_r*I;



%% construct the next generation matrix and compute its leading eigenvalue.
% compute the jacobian matrix 
JF = sym(zeros(dim)); 
JV = sym(zeros(dim)); 
for i = 1:dim
    for j = 1:dim
        JF(i,j) = diff(F(i),variables(j));
        JV(i,j) = diff(V(i),variables(j));
    end
end

% construct the next-generation matrix Mat = F*V^(-1)
Mat = JF*inv(JV)

% all eigenvalues of the Next-Generation matrix
eigenvalues = eig(Mat)

% maximum real part of eigenvalues 
Reff = eigenvalues(end)

% R0 is obtained at the disease free equalibrium (N, 0, 0, 0, 0, 0, 0)
R0 = subs(Reff,S,N)


pretty(Reff)
latex(Reff)