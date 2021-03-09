function [x, t, psi, psire, psiim, psimod, prob, v] = ...
sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)
% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
%
% Outputs
%
% x: Vector of x coordinates [nx]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx]
% psire Array of computed psi_re values [nt x nx]
% psiim Array of computed psi_im values [nt x nx]
% psimod Array of computed sqrt(psi psi*) values [nt x nx]
% prob Array of computed running integral values [nt x nx]
% v Array of potential values [nx]

   % Define mesh and derived parameters ...
   nx = 2^level + 1;
   x = linspace(0.0, 1.0, nx);
   dx = x(2) - x(1);
   dt = lambda * dx;
   nt = round(tmax / dt) + 1;
   t = [0 : nt-1] * dt;
   
   
   % Initialize solution, and set initial data ...
   psi = zeros(nt, nx);
    
   if idtype == 0
      psi(1, :) = sin(idpar(1)* pi * x);
   elseif idtype == 1
      x0 = idpar(1);
      delta = idpar(2);
      p = idpar(3);
      psi(1, :) = exp(1i * p .* x) .* exp(-((x - x0) ./ delta) .^ 2);
   else
      fprintf('sch_1d_cn: Invalid idtype=%d\n', idtype);
      return
   end
   
   % Determine Potential Term ---------------------------------------
   v = zeros(1, nx);
   if vtype == 1  % rectangular well
       imin = round(vpar(1) / dx) +1;
       imax = round(vpar(2) / dx) + 1;
       v(imin:imax) = vpar(3);
   elseif vtype == 0
       % do nothing
   else 
       fprintf('sch_1d_cn: Invalid vtype=%d\n', vtype);
   end

   % Initialize storage for sparse matrix and RHS ...
   dl = zeros(nx,1);
   d  = zeros(nx,1);
   du = zeros(nx,1);
   f  = zeros(nx,1);

   % Set up tridiagonal system ...
   dl = 0.5 / dx^2 * ones(nx, 1);
   d  = (1i / dt - 1.0 / dx^2) * ones(nx,1);
   d = d - 0.5 * v.';
   du = dl;
   % Fix up boundary cases ...
   d(1) = 1.0;
   du(2) = 0.0;
   dl(nx-1) = 0.0;
   d(nx) = 1.0;
   % Define sparse matrix ...
   A = spdiags([dl d du], -1:1, nx, nx);

   % Compute solution using CN scheme ...
   for n = 1 : nt-1
      % Define RHS of linear system ...
      f(2:nx-1) = 1i*psi(n, 2:nx-1) / dt - 0.5 * ( ...
         psi(n, 1:nx-2) - 2 * psi(n, 2:nx-1) + psi(n, 3:nx)) / dx^2 ...
         + 0.5 * v(2:nx-1) .* psi(n, 2:nx-1);
      f(1) = 0.0;
      f(nx) = 0.0;
      % Solve system, thus updating approximation to next time 
      % step ...
      psi(n+1, :) = A \ f;

   end
   
   psire = real(psi);
   psiim = imag(psi);
   psimod = psi.*conj(psi);
   prob = zeros(nt, nx);
   for n = 1 : nt
       sum = 0;
        for k = 2 : nx
            next = (psimod(n, k-1)+ psimod(n, k)) * (x(k)-x(k-1));
            sum = sum + next;
            prob(n,k) = 0.5 * sum;
    
        end
   end
   
   psimod = sqrt(psi.*conj(psi));
   
end
   