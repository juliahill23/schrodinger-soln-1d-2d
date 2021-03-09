% Experiment 1 : Barrier Survey --------------------------------------
tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 20.0];
vtype = 1;
x1 = 0.8;
x2 = 1.0;

lnv0 = linspace(-2.0, 5.0, 251);
lnFe = zeros(length(lnv0), 1);


for i = 1:length(lnv0)
    vpar = [0.6 0.8 exp(lnv0(i))];
    
    [x, t, psi, psire, psiim, psimod, prob, v] = ...
        sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    
    [nt, nx] = size(psi);
    
    temp_avg = zeros(nx, 1);
    for j = 1:nx
        temp_avg(j) = sum(prob(:,j))/nt;
    end
    
    temp_avg = temp_avg/temp_avg(nx);
    
    dx = x(2) - x(1);    
    x1i = round(x1 / dx) + 1;
    x2i = round(x2 / dx) + 1;
    
    frac = temp_avg(x2i) - temp_avg(x1i);
    lnFe(i) = log(frac/(x2-x1));
    
end

plot(lnv0, lnFe, 'c-')
xlabel('ln(V)'); 
ylabel('ln(Fe)') ;

input("press enter for next graph")

% Experiment 2 : Well Survey -----------------------------------------

tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 0.0];
vtype = 1;
x1 = 0.6;
x2 = 0.8;


lnv0 = linspace(2, 10, 251);
lnFe = zeros(length(lnv0), 1);


for i = 1:length(lnFe)
    vpar = [0.6 0.8 -exp(lnv0(i))];
    
    [x, t, psi, psire, psiim, psimod, prob, v] = ...
        sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    
    [nt, nx] = size(psi);
    
    temp_avg = zeros(nx, 1);
    for j = 1:nx
        temp_avg(j) = sum(prob(:,j))/nt;
    end
    
    temp_avg = temp_avg/temp_avg(nx);
    
    dx = x(2) - x(1);    
    x1i = round(x1 / dx) + 1;
    x2i = round(x2 / dx) + 1;
    
    frac = temp_avg(x2i) - temp_avg(x1i);
    lnFe(i) = log(frac/(x2-x1));
    
end

plot(lnv0, lnFe, 'c-')
xlabel('ln(V)'); 
ylabel('ln(Fe)') ;