close all; clear; clc;
format long;

%% Define all symbolic variables and parameters involved in the SIRC model
syms S I R C N;
syms Gamma a p gamma mu beta;

%% Define infected variables
variables = [I C]; % all infected cells
dim = numel(variables);  
F = sym(zeros(dim,1));  % initialize the vector of newly infections
V = sym(zeros(dim,1));  % initialize the vector of transitions (exluding newly infections)

% Define vector F and V
% Note that the original ODEs can be written as  dx/dt = F - V
F(1) = beta*S*(I+C);
V(1) = (a + gamma)*I;
V(2) = -p*gamma*I + (mu + a)*C;




%% construct the next generation matrix and compute its leading eigenvalue.
% compute the jacobian matrices JF, JV via differentiating F, V on variables
JF = sym(zeros(dim)); 
JV = sym(zeros(dim)); 
for i = 1:dim
    for j = 1:dim
        JF(i,j) = diff(F(i),variables(j));
        JV(i,j) = diff(V(i),variables(j));
    end
end

% construct the next-generation matrix Mat = F*V^(-1)
invJV = inv(JV)
Mat = JF*invJV

% all eigenvalues of the Next-Generation matrix
eigenvalues = eig(Mat)

% maximum real part of eigenvalues 
Reff = eigenvalues(end)

% R0 is obtained at the disease free equalibrium (N, 0, 0, 0, 0, 0, 0)
R0 = subs(Reff,S,N)