%% Espacio de fases con tantos vehículos como indique FSV
n        = 1; %póngase aquí el diezmado temporal (1 es nada)
lane     = 2; %póngase aquí el carril que se quiere (interesante ver el 2 por los aumentos repentinos de densidad)
fsv      = 1; %el FSV de siempre para adelgazar la representación 

tiempos  = 0:n*param.dt:param.tmax;
I        = ismember(stv_gt(:,2),tiempos) & stv_gt(:,9)==lane;
stv_gt_n = stv_gt(I,:);
[stv_fsv, ids, Nsv] = sensordata10(stv_gt_n,fsv,0,7);

% Diagrama fundamental microscópico o espacio de fases
h = figure;
axes1 = axes('Parent',h,...
             'FontWeight','bold',...
             'ZDir','reverse',...
             'FontSize',14);
         
%trayectorias completas en gris
hold on 
for i = 1:1:Nsv
    I = stv_fsv(:,7) == ids(i); 
    plot(1000./(stv_fsv(I,5)+param.l),3600.*stv_fsv(I,3)./(stv_fsv(I,5)+param.l),'color',[0.5 0.5 0.5]);
end
axis([0 150 0 2500]);
hold off

%puntos en color, que es la velocidad
h=figure;
axes1 = axes('Parent',h,...
             'FontWeight','bold',...
             'ZDir','reverse',...
             'FontSize',14);
scatter(1000./(stv_fsv(:,5)+param.l),3600.*stv_fsv(:,3)./(stv_fsv(:,5)+param.l),1,stv_fsv(:,3),'filled');
axis([0 150 0 2500]);
%% Película con campo FSV y phase portrait 
% esto sí que mola
h=figure;

a1 = subplot(2,2,1,...
     'FontWeight','bold',...
     'FontSize',14);
axis([0 param.tmax 0 param.smax]);
set(a1,'CLim',[0 33]);
box on
hold on
title('STV field');
xlabel('time(s)');
ylabel('space(m)');

a2 = subplot(2,2,2,...
   'FontWeight','bold',...
   'FontSize',14);

a3 = subplot(2,2,[3 4],...
   'FontWeight','bold',...
   'FontSize',14);

tiempos = 150:param.tmax;
tau = 2; %cuantos puntos quieres representar a la vez para intuir la trayectoria en el espacio de fases
for t = 1:numel(tiempos)
    I1 = stv_fsv(:,2) >= tiempos(t) & stv_fsv(:,2) <= tiempos(t) & stv_fsv(:,8) == 1;
    I2 = stv_fsv(:,2) >= tiempos(t) & stv_fsv(:,2) <= tiempos(t) & stv_fsv(:,8) == 2;
    II1 = stv_fsv(:,2) >= (tiempos(t) - tau) & stv_fsv(:,2) <= (tiempos(t)+tau) & stv_fsv(:,8) == 1; 
    II2 = stv_fsv(:,2) >= (tiempos(t) - tau) & stv_fsv(:,2) <= (tiempos(t)+tau) & stv_fsv(:,8) == 2; 

    axes(a1);
    scatter(stv_fsv(I1|I2,2),stv_fsv(I1|I2,1),5,stv_fsv(I1|I2,3),'filled');
    
    axes(a2);
    plot(1000./(stv_fsv(II1,5)+param.l),3600.*stv_fsv(II1,3)./(stv_fsv(II1,5)+param.l),'k.');
    hold on
    plot(1000./(stv_fsv(II2,5)+param.l),3600.*stv_fsv(II2,3)./(stv_fsv(II2,5)+param.l),'ro');
    hold off
    
    axis([0 150 0 2500]);
    title('flow-density phase space');
    xlabel('density (vh/km)');
    ylabel('flow (vh/h)');
    
    axes(a3);
    plot(stv_fsv(I1,1),ones(length(stv_fsv(I1,1)),1),'k.');
    hold on
    plot(stv_fsv(I2,1),ones(length(stv_fsv(I2,1)),1),'r.');
    hold off
    axis([0 1000 0.9 1.1]);
    drawnow;
 end