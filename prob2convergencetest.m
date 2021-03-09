idtype = 0;
vtype = 0;
idpar = [2,3];
tmax = 0.05;
lambda = 0.05;
lmin = 6;
lmax = 9;
vpar = [0];

   
[x6, y6, t6, psi6, psire6, psiim6, psimod6, v6] = ...
    sch_2d_adi(tmax, 6, lambda, idtype, idpar, vtype, vpar);
[x7, y7, t7, psi7, psire7, psiim7, psimod7, v7] = ...
    sch_2d_adi(tmax, 7, lambda, idtype, idpar, vtype, vpar);
[x8, y8, t8, psi8, psire8, psiim8, psimod8, v8] = ...
    sch_2d_adi(tmax, 8, lambda, idtype, idpar, vtype, vpar);
[x9, y9, t9, psi9, psire9, psiim9, psimod9, v9] = ...
    sch_2d_adi(tmax, 9, lambda, idtype, idpar, vtype, vpar);

[nt6, nx6, ny6] = size(psi6);
[nt7, nx7, ny7] = size(psi7);
[nt8, nx8, ny8] = size(psi8);
[nt9, nx9, ny9] = size(psi9);


% % plots of error --------------------------------------------------

mx = idpar(1);
my = idpar(2);

psixct6 = zeros(nt6, nx6, ny6);
psierr6 = zeros(nt6, nx6, ny6);
errnorm6 = zeros(nt6, 1);

psixct7 = zeros(nt7, nx7, ny7);
psierr7 = zeros(nt7, nx7, ny7);
errnorm7 = zeros(nt7, 1);

psixct8 = zeros(nt8, nx8, ny8);
psierr8 = zeros(nt8, nx8, ny8);
errnorm8 = zeros(nt8, 1);

psixct9 = zeros(nt9, nx9, ny9);
psierr9 = zeros(nt9, nx9, ny9);
errnorm9 = zeros(nt9, 1);

for n = 1 : nt6
    psixct6(n,:, :) = exp(-1i * (mx^2 + my^2) * pi^2 * t6(n)) ...
        .* sin(mx * pi * x6) .* sin(my * pi * y6);
    psierr6(n,:,:) = psi6(n,:,:) - psixct6(n,:,:);
    % error norm
    errnorm6(n) = sqrt(sum(sum((psierr6(n,:,:) .* conj(psierr6(n,:,:))))) ...
        / prod(size(psierr6(n,:,:))));
    
end


for n = 1 : nt7
    psixct7(n,:, :) = exp(-1i * (mx^2 + my^2) * pi^2 * t7(n)) ...
        .* sin(mx * pi * x7) .* sin(my * pi * y7);
    psierr7(n,:,:) = psi7(n,:,:) - psixct7(n,:,:);
    % error norm
    errnorm7(n) = sqrt(sum(sum((psierr7(n,:,:) .* conj(psierr7(n,:,:))))) ...
        / prod(size(psierr7(n,:,:))));
    
end


for n = 1 : nt8
    psixct8(n,:, :) = exp(-1i * (mx^2 + my^2) * pi^2 * t8(n)) ...
        .* sin(mx * pi * x8) .* sin(my * pi * y8);
    psierr8(n,:,:) = psi8(n,:,:) - psixct8(n,:,:);
    % error norm
    errnorm8(n) = sqrt(sum(sum((psierr8(n,:,:) .* conj(psierr8(n,:,:))))) ...
        / prod(size(psierr8(n,:,:))));
    
end


for n = 1 : nt9
    psixct9(n,:, :) = exp(-1i * (mx^2 + my^2) * pi^2 * t9(n)) ...
        .* sin(mx * pi * x9) .* sin(my * pi * y9);
    psierr9(n,:,:) = psi9(n,:,:) - psixct9(n,:,:);
    % error norm
    errnorm9(n) = sqrt(sum(sum((psierr9(n,:,:) .* conj(psierr9(n,:,:))))) ...
        / prod(size(psierr9(n,:,:))));
    
end

clf; hold on; plot(t6, errnorm6, 'r-.'); plot(t7, 4*errnorm7, 'g-.');
plot(t8, 16*errnorm8, 'b-.'); 
plot(t9, 64*errnorm9, 'c-.'); 
xlabel('time (t)'); 
ylabel('Error Norm') ;
legend({'level = 6','level = 7', 'level = 8', 'level = 9'});


input("press enter for next graph")

% 

psi7 = psi7(1:2:end, 1:2:end, 1:2:end);
psi8 = psi8(1:4:end, 1:4:end, 1:4:end);
psi9 = psi9(1:8:end, 1:8:end, 1:8:end);


dpsi6 = psi7 - psi6;
dpsi7 = psi8 - psi7;
dpsi8 = psi9 - psi8;



dpsimod6  = zeros(nt6, 1); 
dpsimod7  = zeros(nt6, 1);
dpsimod8  = zeros(nt6, 1);
for i = 1 : nt6
    dpsimod6(i) = sqrt(sum(sum((dpsi6(i,:,:) .* conj(dpsi6(i,:,:))))) ...
        / prod(size(dpsi6(i,:,:))));
    dpsimod7(i) = sqrt(sum(sum((dpsi7(i,:,:) .* conj(dpsi7(i,:,:))))) ...
        / prod(size(dpsi7(i,:,:))));
    dpsimod8(i) = sqrt(sum(sum((dpsi8(i,:,:) .* conj(dpsi8(i,:,:))))) ...
        / prod(size(dpsi8(i,:,:))));
end

dpsimod7 = 4 * dpsimod7;
dpsimod8 = 16 * dpsimod8;

clf; hold on; plot(t6, dpsimod6, 'r-.o'); plot(t6, dpsimod7, 'g-.+');
plot(t6, dpsimod8, 'b-.*'); 
xlabel('time (t)'); 
ylabel('||dpsi||2') ;
legend({'levels = 6, 7','levels = 7, 8', 'levels = 8, 9'});