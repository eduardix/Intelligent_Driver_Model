% Este script sirve para probar una configuración de parámetros en medidas
% R(t) y R(t,s).

%% Cálculo de densidad y flujo en una región R(t)

%¿Qué carril?
stv_estudio = stv_lane2;

smax        = param.smax;
tmax        = param.tmax;
Tmp         = 60;
tmin        = Tmp;
dx_espiras  = 200;
smin        = 0;
anchura     = 7; %un coche
vmin        = 0.7; %m/s

x_espiras   = smin:dx_espiras:smax;
t_espiras   = Tmp:Tmp:tmax-Tmp;

N           = zeros(numel(x_espiras),numel(t_espiras)-1);
Vs          = N;
Vt          = N;
varVs       = N;
sumvs       = N;
sumvt       = N;
rho1        = N;


for i = 1:numel(x_espiras)
    
    I = stv_estudio(:,1) >= (x_espiras(i) - anchura) & stv_estudio(:,1) <= (x_espiras(i) + anchura);
    temp = stv_estudio(I,:);
    
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

Q1 = 3600.*N./Tmp; %El flujo medido
rho1 = 1000.*(1/Tmp).*sumvs; %Densidad calculada a partir de Vs (armónica)
%rho = Q./(3.6.*Vs);     %Densidad a partir de Relación Fundamental

figure, scatter(rho1(:),Q1(:),5,3.6.*Vs(:),'filled');
axis([0 145 0 2500]);
box on
clear a i1 I temp i j

%% Cálculo de densidad y flujo en una región R(t,s)

% Diezmado a 1s -> N=4
n = 4;
tiempos = 0:4*param.dt:param.tmax;
I       = ismember(stv_estudio(:,2),tiempos);
stv_gt_n = stv_estudio(I,:);

% Las regiones R(t,s) están separadas G y tienen unas dimensiones de Tmp x K

max_pac = 3;               % Número de coches máximo por segmento tamaño K
K = max_pac*(param.l+param.s0);  % K en función del número de coches máximo
G = dx_espiras;                        % Separación en función de K

Tmp = 60;                        % Promedio temporal   

%Calcula parámetros de la rejilla
franjas = floor((smax-smin)/K);
smax    = franjas*K+smin;
NT      = floor((tmax-tmin)/Tmp);    
tmax    = NT*Tmp+tmin;
T       = (tmin:tmax)';
S       = (smin:smax)';
Tc      = (tmin:Tmp:tmax)';
Sc      = (smin:K+G:smax)';

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
% Diezmado a 1s -> N=4

n        = 4;
tiempos  = 0:n*param.dt:param.tmax;
I        = ismember(stv_estudio(:,2),tiempos);
stv_gt_n = stv_estudio(I,:);
I        = stv_gt_n(:,1) > smin;

[stv_fsv, ids, Nsv] = sensordata10(stv_gt_n(I,:),0.005,0,7);

rho_micro = 1000./(stv_fsv(:,5)+param.l);

figure
axis([0 150 0 2500]);
hold on

for i = 1:Nsv
    I = stv_fsv(:,7) == ids(i); 
    plot(1000./(stv_fsv(I,5)+param.l),3600.*stv_fsv(I,3)./(stv_fsv(I,5)+param.l),'color',[0.5 0.5 0.5]);
end

%scatter(rho1(:),Q1(:),15,3.6.*Vs(:),'filled');
scatter(Rts(:,3),Rts(:,4),5,Rts(:,4)./Rts(:,3),'filled');
hold off
box on
axis([0 150 0 2500]);
title(['K = ' num2str(K) 'm. Tmp = ' num2str(Tmp) 's. G = ' num2str(G) 'm.']);
xlabel('density rho vh/km');
ylabel('flow Q vh/h');