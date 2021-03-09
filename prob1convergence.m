% Convergence test 1 ------------------------------------------------

idtype = 0;
vtype = 0;
idpar = [3];
tmax = 0.25;
lambda = 0.1;
vpar = [0 0 0];
lmin = 6;
lmax = 9;


for l = lmin : lmax
   % Compute the solution  ...
   [x{l}, t{l}, psi{l}, psire{l}, psiim{l}, psimod{l}, prob{l}, v{l}] = ...
       sch_1d_cn(tmax, l, lambda, idtype, idpar, vtype, vpar);
   
   [nt{l}, nx{l}] = size(psi{l});
   
   % If possible, compute exact solution ...
   if idtype == 0
      m = idpar(1);
      psixct{l} = zeros(nt{l}, nx{l});
      psierr = zeros(nt{l}, nx{l});
      errnorm{l} = zeros(nt{1}, 1);
      for n = 1 : nt{l}
         psixct{l}(n,:) = exp(-1i*m^2*pi^2*t{l}(n)) * ...
         sin(m* pi * x{l});
         psierr(n,:) = psi{l}(n,:) - psixct{l}(n,:);
         % error norm
         errnorm{l}(n) = sqrt(sum(sum((psierr .* conj(psierr)))) ...
             / prod(size(psierr)));
      end
   end     

end

% plots of error --------------------------------------------------


clf; hold on; plot(t{6}, errnorm{6}, 'r-.'); plot(t{7}, 4*errnorm{7}, 'g-.');
plot(t{8}, 16*errnorm{8}, 'b-.'); plot(t{9}, 64*errnorm{9}, 'c-.'); 

xlabel('time (t)'); 
ylabel('Error Norm') ;
legend({'level = 6','level = 7', 'level = 8', 'level = 9'});

input("press enter for next graph")

psi{7} = psi{7}(1:2:end, 1:2:end);
psi{8} = psi{8}(1:4:end, 1:4:end);
psi{9} = psi{9}(1:8:end, 1:8:end);



dpsi6 = psi{7} - psi{6};
dpsi7 = psi{8} - psi{7};
dpsi8 = psi{9} - psi{8};

dpsimod6  = zeros(nt{6}, 1); 
dpsimod7  = zeros(nt{6}, 1);
dpsimod8  = zeros(nt{6}, 1);
for i = 1 : nt{6}
    dpsimod6(i) = sqrt(sum(sum((dpsi6(i,:) .* conj(dpsi6(i,:))))) ...
        / prod(size(dpsi6(i,:))));
    dpsimod7(i) = sqrt(sum(sum((dpsi7(i,:) .* conj(dpsi7(i,:))))) ...
        / prod(size(dpsi7(i,:))));
    dpsimod8(i) = sqrt(sum(sum((dpsi8(i,:) .* conj(dpsi8(i,:))))) ...
        / prod(size(dpsi8(i,:))));
end

dpsimod7 = 4 * dpsimod7;
dpsimod8 = 16 * dpsimod8;



% plots of difference ---------------------------------------------
clf; hold on; plot(t{6}, dpsimod6, 'r-.'); plot(t{6}, dpsimod7, 'g-.');
plot(t{6}, dpsimod8, 'b-.'); 
xlabel('time (t)'); 
ylabel('||dpsi||2') ;
legend({'levels = 6, 7','levels = 7, 8', 'levels = 8, 9'});

input("press enter for next graph")
disp(t{6}(5))

% plot of real component of psi ------------------------------------
% get real comp of downsampled psi for every level
psire7 = real(psi{7});
psire8 = real(psi{8});
psire9 = real(psi{9});

clf; hold on; plot(x{6}, psire{6}(5,:), 'r-.');
plot(x{6}, psire7(5,:), 'g-.');
plot(x{6}, psire8(5,:), 'b-.'); 
plot(x{6}, psire9(5,:), 'c-.');
plot(x{6},  real(psixct{6}(5,:)), 'm-.'); 
xlabel('x'); 
ylabel('Real component of Psi') ;
legend({'level = 6','level = 7', 'level = 8', 'level = 9', 'Exact Solution'});

input("press enter for next graph")

% plot of imaginary component of psi --------------------------------
% get imaginary comp of downsampled psi for every level
psiim7 = imag(psi{7});
psiim8 = imag(psi{8});
psiim9 = imag(psi{9});

clf; hold on; plot(x{6}, psiim{6}(5,:), 'r-.');
plot(x{6}, psiim7(5,:), 'g-.');
plot(x{6}, psiim8(5,:), 'b-.'); 
plot(x{6}, psiim9(5,:), 'c-.');
plot(x{6},  imag(psixct{6}(5,:)), 'm-.'); 
xlabel('x'); 
ylabel('Imaginary component of Psi') ;
legend({'level = 6','level = 7', 'level = 8', 'level = 9', 'Exact Solution'});

input("press enter for next graph")

% Convergence test 2 ------------------------------------------------

idtype = 1;
vtype = 0;
idpar = [0.50 0.075 0.0];
tmax = 0.01;
lambda = 0.01;
vpar = [0 0 0];
lmin = 6;
lmax = 9;

for l = lmin : lmax

   % Compute the solution  ...
   [x{l}, t{l}, psi{l}, psire{l}, psiim{l}, psimod{l}, prob{l}, v{l}] = ...
       sch_1d_cn(tmax, l, lambda, idtype, idpar, vtype, vpar);
   
   [nt{l}, nx{l}] = size(psi{l});
   
   % If possible, compute exact solution ...
   if idtype == 0
      psixct{l} = zeros(nt{l}, nx{l});
      for n = 1 : nt{l}
         psixct{l}(n,:) = sin(idpar(1)* pi * x{l});
      end
   end     


end


psi{7} = psi{7}(1:2:end, 1:2:end);
psi{8} = psi{8}(1:4:end, 1:4:end);
psi{9} = psi{9}(1:8:end, 1:8:end);


dpsi6 = psi{7} - psi{6};
dpsi7 = psi{8} - psi{7};
dpsi8 = psi{9} - psi{8};

dpsimod6  = zeros(nt{6}, 1); 
dpsimod7  = zeros(nt{6}, 1);
dpsimod8  = zeros(nt{6}, 1);

for i = 1 : nt{6}
    dpsimod6(i) = sqrt(sum((dpsi6(i,:) .* conj(dpsi6(i,:)))) ...
        / prod(size(dpsi6(i,:))));
    dpsimod7(i) = sqrt(sum((dpsi7(i,:) .* conj(dpsi7(i,:)))) ...
        / prod(size(dpsi7(i,:))));
    dpsimod8(i) = sqrt(sum((dpsi8(i,:) .* conj(dpsi8(i,:)))) ...
        / prod(size(dpsi8(i,:))));

end

dpsimod7 = 4 * dpsimod7;
dpsimod8 = 16 * dpsimod8;

% plots of difference ---------------------------------------------
clf; hold on; plot(t{6}, dpsimod6, 'r-.o'); plot(t{6}, dpsimod7, 'g-.+');
plot(t{6}, dpsimod8, 'b-.*'); 
xlabel('time (t)'); 
ylabel('||dpsi||2') ;
legend({'levels = 6, 7','levels = 7, 8', 'levels = 8, 9'});
input("press enter for next graph")
disp(t{6}(49))


% plot of real component of psi ------------------------------------
% get real comp of downsampled psi for every level
psire7 = real(psi{7});
psire8 = real(psi{8});
psire9 = real(psi{9});

clf; hold on; plot(x{6}, psire{6}(49,:), 'r-.');
plot(x{6}, psire7(49,:), 'g-.');
plot(x{6}, psire8(49,:), 'b-.'); 
plot(x{6}, psire9(49,:), 'c-.');
xlabel('x'); 
ylabel('Real component of Psi') ;
legend({'level = 6','level = 7', 'level = 8', 'level = 9'});

input("press enter for next graph")

% plot of imaginary component of psi --------------------------------
% get imaginary comp of downsampled psi for every level
psiim7 = imag(psi{7});
psiim8 = imag(psi{8});
psiim9 = imag(psi{9});

clf; hold on; plot(x{6}, psiim{6}(49,:), 'r-.');
plot(x{6}, psiim7(49,:), 'g-.');
plot(x{6}, psiim8(49,:), 'b-.'); 
plot(x{6}, psiim9(49,:), 'c-.');
xlabel('x'); 
ylabel('Imaginary component of Psi') ;
legend({'level = 6','level = 7', 'level = 8', 'level = 9'});
