% Pere Camps
% TFM - Airfoil Optimisation with Genetic Algorithm
% All rights reserved

close all; clear;

PUASE_TIME = 0.1;

results = textread('Results.txt');
gen = results(:,1);
cl  = results(:,2);
cd  = results(:,3);
eff = results(:,4);
fit = results(:,5);
cm  = results(:,6);

eff_(1) = eff(1);
MAX     = 0;

for i=2:length(eff)
    if(eff(i)>MAX)
        MAX = eff(i);
    end
    if(eff(i)<MAX)
        eff_(i) = eff_(i-1);
    else
        eff_(i) = eff(i);
    end
end

best_gen_aux = find(fit==min(fit));
best_gen     = best_gen_aux(1)-1;

%% LOOP PLOT
n_max      = max(gen);
formatSpec = '%f %f';
sizeA      = [2 199];
aux = 0;
figure (01)
 xlim([0 1]); ylim([-0.5 0.5]);
 grid on; grid minor
 title('Best Arifoil of Generation');
 hold on
for i=0:n_max
   
    name   = sprintf('airfoil_gen%d.dat',i);
    fileID = fopen(name,'r');
    A  = fscanf(fileID,formatSpec,sizeA).';
    x  = A(100:end,1);       
    zu = flipud(A(1:100,2));
    zl = A(100:end,2);
    
    NamePlot = sprintf('Generation %d',i);
    text(0.4,0.3,NamePlot)
    
    if(fit(i+1)== min(fit) && aux == 0)
        plot(x,zu,'b',x,zl,'b');
        aux = 1;
    elseif(aux == 1) 
        plot(x,zu,'r',x,zl,'r');
    else
        plot(x,zu,'k',x,zl,'k');
    end
    
    pause(PUASE_TIME);
    fclose(fileID);
    
    cla reset;
    xlim([0 1]); ylim([-0.5 0.5]);
    grid on; grid minor;
    title('Best Arifoil of Generation');
    hold on  
end
hold off

%% FITS
n2 = sprintf('Best Generation: %d',best_gen);

figure(02)
hold on
grid on; grid minor;
xlabel('Generation');
yyaxis left
    grid minor;
    ylim([0 max(cl)+0.2]);
    text(max(gen)/2-1,max(cl)-0.2,n2)
    plot(gen,cl);
    ylabel('C_l');
yyaxis right
    grid minor;
    plot(gen,cd);
    %ylim([0.009 0.02]);
    ylabel('C_d');
hold off

figure(03)
hold on
grid on; grid minor;
xlabel('Generation');ylabel('Fitness function');
title('Fitness Evolution');
plot(gen,fit);
hold off

figure(04)
grid on; grid minor;
hold on
title('Efficiency Evolution');
xlabel('Generation'); ylabel('Efficiency');
plot(gen,eff,'b--',gen,eff_,'b');
hold off

figure(05)
hold on
grid on; grid minor;
xlabel('Generation');ylabel('C_m');
title('C_m Evolution');
plot(gen,cm);
hold off

%% Best airfoil

name = sprintf('airfoil_gen%d.dat',best_gen);
fileID = fopen(name,'r');
A = fscanf(fileID,formatSpec,sizeA).';

x  = A(100:end,1);
zu = flipud(A(1:100,2));
zl = A(100:end,2);

  
figure (06)
hold on
title('BEST AIRFOIL');
xlim([0 1]); ylim([-0.5 0.5]);
grid on; grid minor;
NamePlot = sprintf('Generation %d',best_gen);
text(0.4,0.3,NamePlot)
plot(x,zu,'b',x,zl,'b');
plot(x,(zu+zl)/2,'g'); 
hold off
