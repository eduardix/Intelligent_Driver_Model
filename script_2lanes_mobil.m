clear all
close all
%dt, tiempo máximo, longitud carretera
dt      = 0.25;         %(s) 
tmax    = 1000;         %(s)
smax    = 1000;
%Parámetros IDM (malos conductores)

%T       = 1.2;          %(s)       %timehead demasiado pequeño
%a       = 0.8;          %(m/s^2)   %aceleración demasiado pequeña
%b       = 2.5;          %(m/s^2)   %decelación brusca
%ab      = 2*sqrt(a*b);  %(m/s^2)
%s0      = 2;            %(m)
%l       = 5;            %(m)
%v0      = 120/3.6;      %(m/s)

%Parámetros IDM (buenos conductores)
s0      = 2;            %(m)        
T       = 1.6;          %(s)       %timehead recomendado internacional
a       = 1;            %(m/s^2)   %aceleración suave pero suficiente
b       = 1.65;         %(m/s^2)   %deceleración suave
ab      = 2*sqrt(a*b);  %(m/s^2)
l       = 5;            %(m)
v0      = 120/3.6;      %(m/s)

%Aleatoriedad velocidad deseada
desv    = 0/3.6;         %(m/s)

%Parámetros entrada vehículos lanes:
lambda1 = 2; %1 coche por lambda s lane 1 
lambda2 = 5; %1 coche por lambda s lane 2

%Introducción de vehículos para cumplir condición de contorno:
llegadas1 = zeros(1,numel(0:dt:tmax));
llegadas2 = zeros(1,numel(0:dt:tmax));

vin = @(va,s) min(v0,real((sqrt(ab.^2.*T.^2-2.*ab.*T.*va+va.^2-4.*ab.*(s0-s))-ab.*T+va)./2));
da  = @(va) (s0+va*T)/sqrt(1-(va/v0)^4);

%Generación velocidades deseadas aleatorias
V0    = v0+desv.*randn(1,ceil(2/lambda1*tmax)+ceil(2/lambda2*tmax)); %esto asegura que la vmax=v0
I     = V0 <= v0-0.7*desv | V0 > v0+0.7*desv;
V0(I) = v0; %Y que la velocidad mínima no va a ser menor que mu-n*sigma

%parámetros MOBIL
global p b_mobil ath;

p = 0.5;
b_mobil = -4;
ath = 0.2;

%Variables de salida y control CARRIL 1:
%Ax = [S T V a Ds Dv id s0 T a ab v0 l source_lane actual_lane];

%carril 1
id = 1;
A1 = [10000 0 0 0 0 0 -1 0 0 0 0 0 0 1 1];
A1 = [0 0 V0(id) 0 A1(1,1)-l 0 id s0 T a ab v0 l 1 1; A1];

%carril 2
id = 2;
A2 = [10000 0 0 0 0 0 -2 0 0 0 0 0 0 2 2];
A2 = [0 0 V0(id) 0 A2(1,1)-l 0 id s0 T a ab v0 l 2 2; A2];
id = id + 1;

numt = numel(0:dt:tmax);
nums = ceil(smax/(s0+l));
stv_gt = zeros(2*numt*nums,9);

j = 1;
k1 = 1;
k2 = 1;

tic
for t = 0:dt:tmax
   disp(t);
   %CARRIL 1 
   %Llegadas Poissonianas, parámetro lambda1
   t0 = 0;
   arrivaltimes1=[];
   while t0<=dt %
       t0 = t0 - log(rand)*lambda1; %tiempo a la proxima llegada
       arrivaltimes1 = [arrivaltimes1 t0];
   end
   umbral1 = sum(arrivaltimes1<=dt);
   %Introduzco el coche, si cabe
   if umbral1 >= 1
       if A1(1,1) > s0+l
            %A = [S T V a Ds Dv id s0 T a ab v0 l source_lane actual_lane];
            v_temp = vin(A1(1,3),A1(1,1)-A1(1,13));
            %v_temp = vin(A(1,3),0);
            %v_temp = A(1,3);
            A1 = [0 t v_temp 0 A1(1,1)-A1(1,13) v_temp-A1(1,3) id s0 T a ab V0(id) l 1 1; A1];
            %A1 = [0 t v_temp 0 A(1,1)-A(1,13) 0 id s0 T a ab V0(id) l 1 1; A];
            id = id + 1;
            llegadas1(k1) = 1;
       end
   end
   k1 = k1+1;
   
   %CARRIL 2:
   %Llegadas Poissonianas, parámetro lambda2
   t0 = 0;
   arrivaltimes2=[];
   while t0<=dt
       t0 = t0 - log(rand)*lambda2; %tiempo a la proxima llegada
       arrivaltimes2 = [arrivaltimes2 t0];
   end
   umbral2 = sum(arrivaltimes2<=dt);
   %Introduzco el coche, si cabe
   if umbral2 >= 1
       if A2(1,1) > s0+l
            %A = [S T V a Ds Dv id s0 T a ab v0 l source_lane actual_lane];
            v_temp = vin(A2(1,3),A2(1,1)-A2(1,13));
            %v_temp = vin(A(1,3),0);
            %v_temp = A(1,3);
            A2 = [0 t v_temp 0 A2(1,1)-A2(1,13) v_temp-A2(1,3) id s0 T a ab V0(id) l 2 2; A2];
            %A2 = [0 t v_temp 0 A(1,1)-A(1,13) 0 id s0 T a ab V0(id) l 2 2; A];
            id = id + 1;
            llegadas2(k2) = 1;
       end
   end
   k2 = k2+1;
   
   %Paso de simulación IDM + MOBIL
   [A1,A2] = MIXED_engine(A1,A2,dt);
   
   %Elimina coches de la carretera cuando alcanzan el final
   I = A1(:,1) > smax & A1(:,7) > 0;
   A1(I,:) = [];
   
   %Elimina coches de la carretera cuando alcanzan el final
   I = A2(:,1) > smax & A2(:,7) > 0;
   A2(I,:) = [];
   
   %Guarda en el campo
   stv_gt(j:(j-1)+size(A1,1),:) = A1(:,[1:7 14 15]);
   j = j+size(A1,1);
   
   stv_gt(j:(j-1)+size(A2,1),:) = A2(:,[1:7 14 15]);
   j = j+size(A2,1);
end

%Limpieza campo STV
stv_gt = stv_gt(stv_gt(:,7) > 0,:); %quita dummies

param.s0    = s0;
param.T     = T;
param.a     = a;
param.b     = b;
param.ab    = ab;
param.v0    = v0;
param.desv  = desv;
param.l     = l;
param.V0    = V0;

param.dt = dt;
param.tmax = tmax; %en s.
param.t0 = t0;
param.smax = smax;
param.lambda1 = lambda1;
param.llegadas1 = llegadas1;
param.lambda2 = lambda2;
param.llegadas2 = llegadas2;

clear s0 Tt Ss a b ab v0 desv l V0 dt tmax smax T I1 I2 I3 I4 In dIn t0 vin n v_temp
clear umbral lambda A I id k t j numt nums sih plot_rt m Tin  Tf lih arrivaltimes llegadas

toc
%% DIBUJA DISTINTOS CARRILES EN EL MISMO PLOT

%carril1
figure
hold on
colores = ['-k';'-r'];
%pinta:
%A = [S T V a Ds Dv id source_lane actual_lane];
I = stv_gt(:,9) == 1;
stv_lane1 = stv_gt(I,:);
I = stv_gt(:,9) == 2;
stv_lane2 = stv_gt(I,:);
ids1 = unique(stv_lane1(:,7));
ids2 = unique(stv_lane2(:,7));
for i = 1:numel(ids1)
    I = stv_lane1(:,7) == ids1(i);
    plot(stv_lane1(I,2), stv_lane1(I,1), colores(unique(stv_lane1(I,8)),:))
end
axis([0 param.tmax 0 param.smax]);
xlabel('time (min)','fontsize',14);
ylabel('space (km)','fontsize',14);
box on
title('Lane 1');

%carril2
figure
hold on
for i = 1:numel(ids2)
    I = stv_lane2(:,7) == ids2(i);
    plot(stv_lane2(I,2), stv_lane2(I,1), colores(unique(stv_lane2(I,8)),:))
end
axis([0 param.tmax 0 param.smax]);
xlabel('time (min)','fontsize',14);
ylabel('space (km)','fontsize',14);
box on
title('Lane 2');
%% Película
figure
for t = 5:param.dt:param.tmax
    I = stv_gt(:,9)==1 & stv_gt(:,2)==t;
    II = stv_gt(:,9)==2 & stv_gt(:,2)==t;
    
    scatter(stv_gt(I,1),ones(length(stv_gt(I,1)),1),20,stv_gt(I,8),'filled');
    hold on
    scatter(stv_gt(II,1),1.3.*ones(length(stv_gt(II,1)),1),20,stv_gt(II,8),'filled');
    hold off
    axis([0 1000 0.9 1.4]);
    drawnow;
    disp(t)
end