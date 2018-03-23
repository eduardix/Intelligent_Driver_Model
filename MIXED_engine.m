function [Ao,An,cambios] = MIXED_engine(A1,A2,dt)
    global p ath b_mobil gap_d gap_t;
    %A = [S T V a Ds Dv id s0 T a  ab v0 l  source_lane actual_lane];
    %     1 2 3 4 5  6  7  8  9 10 11 12 13          14          15];

    %PARAMETROS: Luego irán en las matrices y podrán cambiar
    
    % Ecuación del MOBIL:
    % acn-ac + p(ann - an + aon - ao) > ath

    %Primero anima los carriles por separado y anota las aceleraciones
    %Carril 1:
    Ao = IDM_engine(A1,dt);
    %Carril 2:
    An = IDM_engine(A2,dt);
    
    %con esto ya tenemos ac, an y ao de todos los coches, 
    %además de un paso completo de simulación si no hubiera cambios de carril. 

    Io = A1(:,7) < 0;
    In = A2(:,7) < 0; 
    %Cambio de carril 1 a carril 2:
    % En una primera versión (jajajajaja) lo voy a hacer con un "for"
    % En una segunda versión intentaré hacerlo completamente vectorial
    cambios = 0;
    if sum(~Io)>=2 && sum(~In)>=2 %Si hay al menos dos coches... 
        
        ids = unique(A1(:,7));
        ids(ids<0) = []; %elimina los dummies
        C = [ids, zeros(numel(ids),1)]; %Matriz de cambios
        
        %Va pasando por todos los coches menos los dummies. El coche
        %indexado se llama ACTUAL.
        
        for i=1:numel(ids) %INICIO DE LA EVALUACIÓN DE LOS CAMBIOS
            %Selecciona coche ACTUAL tanto antes como después de hacer IDM sin cambio de carril 
            actual_A1 = A1(:,7) == ids(i);
            actual_Ao = Ao(:,7) == ids(i);
            
            %OBTENCIÓN DE ao y CALCULO de aon
            %No existe seguidor y no influye en el decisor de MOBIL
            if actual_A1(1) == 1 
                aon = 0;
                ao  = 0;
            else %Existe seguidor y sí influye en el decisor de MOBIL 
                %identificadores de los coches implicados antes y después
                %del cambio
                follow_A1 = circshift(actual_A1,-1);
                follow_Ao = circshift(actual_Ao,-1);
                next_A1 = circshift(actual_A1,1);
                %next_Ao = circshift(actual_Ao,1);
                
                %Ejecuta una versión local del IDM con coche ACTUAL
                %retirado. Aquí podría usarse un modelo psicológico o
                %perceptual en coches normales, o parámetros alterados para
                %modelar la percepción humana. 
                temp = IDM_engine([A1(follow_A1,:); A1(next_A1,:)],dt); %ejecuta IDM con coche ACTUAL retirado 
                aon = temp(1,4);
                ao = Ao(follow_Ao,4);
            end
            
            %OBTENCIÓN DE ac Y an Y CÁLCULO DE acn Y ann
            ac = Ao(actual_Ao,4);
            %mete ACTUAL en A2
            temp = sortrows([A2; A1(actual_A1,:)],1);
            
            %localiza los coches implicados en nuevo vector
            actual_temp = temp(:,7) == ids(i);
            next_temp = circshift(actual_temp,1);

            %No hay nuevo seguidor y no influye en el decisor MOBIL
            if actual_temp(1) == 1
                ann = 0;
                an = 0;
                mini_vect = [temp(actual_temp,:);temp(next_temp,:)];
                gap1 = mini_vect(2,1) - mini_vect(1,1) - mini_vect(1,13); %gap delantero
                gap2 = gap_t+1; % mayor que gap trasero para que no influya;
                temp1 = IDM_engine(mini_vect,dt);
                acn = temp1(1,4);
                acn = 0;
            else %Existe nuevo seguidor y sí influye en el decisor de MOBIL
                follow_temp = circshift(actual_temp,-1);
                id_follow = temp(follow_temp,7);
                an = A2(A2(:,7) == id_follow,4);
                mini_vect = [temp(follow_temp,:);temp(actual_temp,:);temp(next_temp,:)];
                gap1 = mini_vect(3,1) - mini_vect(2,1) - mini_vect(2,13) ; %gap delantero
                gap2 = mini_vect(2,1) - mini_vect(1,1) - mini_vect(1,13);  %gap trasero
                temp = IDM_engine(mini_vect,dt);
                acn = temp(2,4);
                ann = temp(1,4);
            end
            
            % DECISOR DEL MOBIL:
            % INCENTIVO: acn-ac + p(ann - an + aon - ao) > ath  
            % SEGURIDAD: ann > b_mobil
            % ESPACIO: gap1 >= gap_d y gap2 >= gap_t
            
            if  ((acn - ac + p*(ann - an + aon - ao)) > ath) && (ann > b_mobil) && (gap1 >= gap_d) && (gap2 >= gap_t)
                C(i,2) = 1; 
            end
        end %FIN DE LA EVALUACIÓN DE LOS CAMBIOS
        
        %Ahora insertamos los coches marcados por C y ejecutamos IDM otra
        %vez (esto habrá que cambiarlo para hilar más fino, vectorial)
        cambios = sum(C(:,2));
        if  cambios > 0 %al menos uno cambia 
            I = ismember(A1(:,7),C(C(:,2) == 1,1)); %Matlab cinturon negro sexto dan
            A1(I,15) = 2; %Cambia la columna actual_lane             
            An = IDM_engine(sortrows([A2; A1(I,:)],1),dt);
            Ao = IDM_engine(A1(~I,:),dt);
        end 
    end   
end