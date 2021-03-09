function [x, y, t, psi, psire, psiim, psimod, v] = ...
sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)
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
% y: Vector of y coordinates [ny]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx x ny]
% psire Array of computed psi_re values [nt x nx x ny]
% psiim Array of computed psi_im values [nt x nx x ny]
% psimod Array of computed sqrt(psi psi*) values [nt x nx x ny]
% v Array of potential values [nx x ny]

    % construct meshes
    nx = 2^level + 1;
    ny = nx;
    x = linspace(0.0, 1.0, nx);
    y = linspace(0.0, 1.0, ny).';
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    dt = lambda * dx;
    nt = round(tmax/dt) + 1;
    t = [0 : nt-1] * dt;
    
    % construct output arrays
    psi = zeros(nt, nx, ny);
    psire = zeros(nt, nx, ny);
    psiim = zeros(nt, nx, ny);
    psimod = zeros(nt, nx, ny);
    v = zeros(nx, ny);

    % initialize solution, set initial data ---------------------------
    
    if idtype == 0  % exact family
        mx = idpar(1);
        my = idpar(2);
        psi(1,:,:) = sin(mx*pi*x) .* sin(my*pi*y);
    elseif idtype == 1  % boosted gaussian
        x0 = idpar(1);
        y0 = idpar(2);
        delx = idpar(3);
        dely = idpar(4);
        px = idpar(5);
        py = idpar(6);
        psi(1,:,:) = exp(1i .* px .* x -((x - x0).^2 ./ delx^2)) .* ...
            exp(1i .* py .* y - ((y - y0).^2/ dely^2));
       
    else
        fprintf('sch_2d_adi: Invalid idtype=%d\n', idtype);
        return
    end
            
    if vtype == 0  % no potential
        % do nothing
    elseif vtype == 1  % rectangular barrier or well
       iminx = round(vpar(1) / dx) +1;
       imaxx = round(vpar(2) / dx) + 1;
       iminy = round(vpar(3) / dy) +1;
       imaxy = round(vpar(4) / dy) + 1;
       v(iminx:imaxx, iminy:imaxy) = vpar(5);
    elseif vtype == 2  % double slit
        x1 = vpar(1);
        x2 = vpar(2);
        x3 = vpar(3);
        x4 = vpar(4);
        vc = vpar(5);
        it1 = round(x1/ dx) + 1;
        it2 = round(x2/ dx) + 1;
        it3 = round(x3/ dx) + 1;
        it4 = round(x4/ dx) + 1;
        
        jprime = (ny - 1)/4 + 1;
        v(1:it1,jprime:jprime+1) = vc;
        v(it2:it3,jprime:jprime+1) = vc;
        v(it4:nx,jprime:jprime+1) = vc;
    else
        fprintf('sch_2d_adi: Invalid idtype=%d\n', idtype);
        return
    end
    
    % Initialize storage for sparse matrices and RHS ...
    dlhalf = zeros(nx,1);
    dhalf  = zeros(nx,1);
    duhalf = zeros(nx,1);
    fhalf  = zeros(nx,1);
    dl = zeros(ny,1);
    du = zeros(ny,1);
    
    f  = zeros(ny,1);
    
    % set up tridiagonal systems -------------------------------------
    
    dlhalf = - 1i * dt/(2 * dx^2) * ones(nx, 1);
    dhalf  = (1 + (1i*dt)/(dx^2)) * ones(nx, 1);
    duhalf = dlhalf;
    dl = - 1i * dt / (2 * dy^2) * ones(ny, 1);
    du = dl;
    
    
    % Fix up boundary cases ...
    dhalf(1) = 1.0;
    duhalf(2) = 0.0;
    dlhalf(nx-1) = 0.0;
    dhalf(nx) = 1.0;
    
    du(2) = 0.0;
    dl(ny-1) = 0.0;

    
   % Define sparse matrices ...
   Ahalf = spdiags([dlhalf dhalf duhalf], -1:1, nx, nx);


   v = v.';
    
   for n = 1 : nt-1
       psihalf = zeros(nx, ny);
       for j = 2 : ny-1
           delyypsin = (1/(dy^2))*(psi(n,2:nx-1,j+1) - 2 * psi(n,2:nx-1,j) + ...
               psi(n, 2:nx-1, j-1));  % I'm doing this for all x here even tho 
                                % it's only valid for 2:nx-1 but I will
                                % only use the 2:nx-1 values in my
                                % calculations
           
           % expansion of 2nd bracket on rhs of eq 15
           temp = zeros(nx,1);
           fhalf = zeros(nx, 1);
           
           temp(2:nx-1) = psi(n, 2:nx-1, j) - 1i*(dt/2)*v(2:nx-1,j).' .* ...
               psi(n, 2:nx-1, j) + 1i * (dt/2) * delyypsin;
           
           delxxpsin = (1/dx^2) * (temp(3:nx) - 2 * temp(2:nx-1) ...
               + temp(1:nx-2));
           
           fhalf(2:nx-1) = temp(2:nx-1) + 1i * (dt/2) * delxxpsin;
           fhalf(1) = 0.0;
           fhalf(nx) = 0.0;
           psihalf(:,j) = Ahalf \ fhalf; 
       end


       for i = 2 : nx - 1
           f(1) = 0.0;
           f(ny) = 0.0;
           f(2:ny-1) = psihalf(i, 2:ny-1);

           % calc matrix A
           
        
           d = (1 + (1i* (dt/dy^2))) + (1i* (dt/2) .* v(i,:).') .* ones(ny, 1);
           d(1) = 1.0;
           d(ny) = 1.0;
           A = spdiags([dl d du], -1:1, ny, ny);
           psi(n+1,i,:) = A \ f;
       end
   end
   
   psire = real(psi);
   psiim = imag(psi);
   psimod = sqrt(psi.*conj(psi));
    
end