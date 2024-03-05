function f = ColeHopfTransformation(initialCond, c)
% This function coleHopfTransformation, applies the Cole-Hopf transformation to a given anonymous function
% for a chosen value of the dissipation coefficient. This function receives the initial condition for
% Burgers' equation and converts it into the corresponding initial condition for the Heat equation. The
% given initial condition must be an integrable function i.e. finding a symbolic anti-derivative must be 
% posible.
%
%         initialCond = The chosen initial condition for Burgers' equation given as an anonymous function
%                   c = The diffusion/viscosity constant (positive real number)
%                  
% @author Jarren Ralf

syms v(x) % Declare the symbolic variable of v as a function of x

% The Cole-Hopf tranformation and a condition to solve for the constant that arises after integration
ode = diff(v, x) == -1/2*initialCond(x)/c*v;
cond = v(0) == 1;

% Solve the ODE and output the solution as a function handle (anonymous function)
f = matlabFunction(dsolve(ode, cond));
