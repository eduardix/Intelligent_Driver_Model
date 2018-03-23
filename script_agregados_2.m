% Este script calcula diagramas fundamentales en regiones R(t) y R(t,s) con
% distintos parámetros. Los parámetros son el número de lugares de medida,
% el tamaño de los lugares de medida y el tiempo de promedio.
% Adicionalmente, la medida en región R(t) tiene una velocidad mínima de
% activación de la espira.

%% Cálculo de densidad y flujo en una región R(t) 
%diferentes Tmp en una sola gráfica

smax        = param.smax;
tmax        = param.tmax;

dx_espiras  = 1000;
smin        = 400;
anchura     = 14; %un coche
vmin        = 1; %m/s

x_espiras   = smin:dx_espiras:smax;

Tmp         = [30 60 90 120];
color       = [0.7 0.7 0.7; 0.5 0.5 0.5; 0.25 0.25 0.25; 0 0 0];

h = figure;
axes1 = axes('Parent',h,...
    'YTickLabel',{'0', '250', '500', '750', '1000', '1250', '1500', '1750', '2000', '2250', '2500'},...
    'YTick',[0 250 500 750 1000 1250 1500 1750 2000 2250 2500],...
    'XTickLabel',{'0','25','50','75','100','125','150'},...
    'XTick',[0 25 50 75 100 125 150],...
    'FontWeight','bold',...
    'FontSize',14);
axis([0 150 0 2500]);
hold on

for t=1:4
    
    tmin        = Tmp(t);
    t_espiras   = Tmp(t):Tmp(t):tmax-Tmp(t);
    
    N           = zeros(numel(x_espiras),numel(t_espiras)-1);
    Vs          = N;
    Vt          = N;
    varVs       = N;
    sumvs       = N;
    sumvt       = N;
    rho1        = N;

    for i = 1:numel(x_espiras)

        I = stv_gt(:,1) >= (x_espiras(i) - anchura) & stv_gt(:,1) <= (x_espiras(i) + anchura);
        temp = stv_gt(I,:);

        for j = 1:numel(t_espiras)-1
            I = temp(:,2) >= t_espiras(j) & temp(:,2) < t_espiras(j+1);
            a = temp(I,:);

            %Umbral de detección
            a = a(a(:,3)>=vmin,:);

            %Contar sólo una vez
            [~,i1] = unique(a(:,7));
            a = a(i1,:);

            %Cuenta
            N(i,j) = size(a,1);

            %Velocidades
            Vs(i,j) = harmmean(a(:,3)); %media harmonic
            Vt(i,j) = mean(a(:,3));     %media aritmética
            varVs(i,j) = (1./(N(i,j)-1)).*sum((a(:,3)-Vs(i,j)).^2);

            sumvs(i,j) = sum(1./a(:,3));
            sumvt(i,j) = sum(a(:,3));

            %rho(i,j) = 1./t_agregado*V1(i,j); 
            %I1 = a(:,1) <= x_espiras(i);
            %I2 = a(:,1) > x_espiras(i);
            %ids = unique(a(I1,7));
            %I = ismember(ids,a(I2,7));
            %V(i,j) = mean(a(I2,3));
            %N(i,j) = sum(I);
            %O(i,j) = N(i,j)/numel(I);
        end
    end
    Q1 = 3600.*N./Tmp(t); %El flujo medido
    rho1 = 1000.*(1/Tmp(t)).*sumvs; %Densidad calculada a partir de Vs
    %rho = Q./(3.6.*Vs);     %Densidad a partir de Relación Fundamental
    
    scatter(rho1(:),Q1(:),5,'markerfacecolor',color(t,:),'markeredgecolor',color(t,:));
    %scatter(rho1(:),Q1(:),30,'d');
    xlabel('density (vh/km)','fontsize',14,'FontWeight','bold');
    ylabel('flow (vh/h)','fontsize',14,'FontWeight','bold');
    box on
end

clear a i1 I temp i j

%% Cálculo de densidad y flujo en una región R(t,s)

% Submuestreo del campo para 
tiempos = 0:1:param.tmax;
I       = ismember(stv_gt(:,2),tiempos);
stv_gt_n = stv_gt(I,:);

smax        = param.smax;
tmax        = param.tmax;

smin        = 400;

mp = [2 5 7 10];        % Número de coches máximo por segmento tamaño K
prom = [30 60 90 120];  % Promedio temporal
color = [0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2];
% Las regiones R(t,s) están separadas G y tienen unas dimensiones de Tmp x K
D = 1000;
for t = 1:length(prom)
    
    Tmp = prom(t);
    tmin = Tmp;
     
    h = figure;
    axes1 = axes('Parent',h,...
        'YTickLabel',{'0', '250', '500', '750', '1000', '1250', '1500', '1750', '2000', '2250', '2500'},...
        'YTick',[0 250 500 750 1000 1250 1500 1750 2000 2250 2500],...
        'XTickLabel',{'0','25','50','75','100','125','150'},...
        'XTick',[0 25 50 75 100 125 150],...
        'FontWeight','bold',...
        'FontSize',14);
    axis([0 150 0 2500]);
    hold on
    box on
    xlabel('density (vh/km)','fontsize',14,'FontWeight','bold');
    ylabel('flow (vh/h)','fontsize',14,'FontWeight','bold');
    
    for  m = 1:length(mp)
        
        max_pac = mp(m);
        K = max_pac*(param.l+param.s0);  % K en función del número de coches máximo
        G = D - K;

        %Calcula parámetros de la rejilla
        franjas = floor((smax-smin)/K);
        smax    = franjas*K+smin;
        NT      = floor((tmax-tmin)/Tmp);    
        tmax    = NT*Tmp+tmin;
        T       = (tmin:tmax)';
        S       = (smin:smax)';
        Tc      = (tmin:Tmp:tmax)';
        Sc      = (smin:D:smax)';

        disp('------------------------------');
        disp(['Número de regiones R(t,:): ' num2str(numel(Sc))]);
        disp(['Número de regiones R(:,s): ' num2str(NT)])
        disp(['K:                         ' num2str(K)])
        disp(['Tmp:                       ' num2str(Tmp)]);
        disp(['Empaquetamiento máximo:    ' num2str(max_pac)]);
        disp(['Tmax:                      ' num2str(max(T(:)))]);
        disp(['Smax:                      ' num2str(max(S(:)))]);

        Rts = zeros((length(Tc)-1)*length(Sc),4);
        kk = 1;
        for i = 1:length(Sc)
        I = stv_gt_n(:,1) >= Sc(i) & stv_gt_n(:,1) <= Sc(i)+K;
        temp1 = stv_gt_n(I,:);
            for j = 1:length(Tc)-1
                I = temp1(:,2) >= Tc(j) & temp1(:,2) < Tc(j+1);
                temp2   = temp1(I,:);

                %%%%% CALCULOS %%%%%
                tiempos = Tc(j):Tc(j+1)-1;
                rho     = 1000/(Tmp*K)*sum(histc(temp2(:,2),tiempos));
                q1      = zeros(Tmp,1);

                for k = 1:Tmp
                    q1(k) = 1/K*sum(temp2(temp2(:,2) == tiempos(k),3));  
                end

                q = (3600/Tmp)*sum(q1);

                Rts(kk,:) = [i j rho q];
                kk     = kk + 1;

            end
        end
        clear i j k kk q rho I temp1 temp2 tiempos q1

        % Diagrama fundamental formado por algunas trayectorias individuales
        % Diezmado a 1s N=4

        %n        = 1;
        %tiempos  = 0:n*param.dt:param.tmax;
        %I        = ismember(stv_gt(:,2),tiempos);
        %stv_gt_n = stv_gt(I,:);
        %I        = stv_gt_n(:,1) > smin;

        %[stv_fsv, ids, Nsv] = sensordata10(stv_gt_n(I,:),0.001,0,7);

        %rho_micro = 1000./(stv_fsv(:,5)+param.l);

        %for i = 1:Nsv
        %    I = stv_fsv(:,7) == ids(i); 
        %    plot(1000./(stv_fsv(I,5)+param.l),3600.*stv_fsv(I,3)./(stv_fsv(I,5)+param.l),'color',[0.5 0.5 0.5]);
        %end

        %scatter(rho1(:),Q1(:),15,3.6.*Vs(:),'filled');
        %scatter(Rts(:,3),Rts(:,4),5,Rts(:,4)./Rts(:,3),'filled');
        scatter(Rts(:,3),Rts(:,4),5,'o','markeredgecolor',color(m,:));
        
        %title(['K = ' num2str(K) 'm. Tmp = ' num2str(Tmp) 's. G = ' num2str(G) 'm.']);
        
    end
    hold off
end