function A = IDM_engine(A,dt)
    
    %A = [S T V a Ds Dv id s0 T a  ab v0 l  lane_x];
    %     1 2 3 4 5  6  7  8  9 10 11 12 13     14];
    
    I = A(:,7) < 0;
    dummies = A(I,:); 
    if sum(~I) 
        %gap y velocidad relativa
        A(:,5) = [diff(A(:,1)) - A(1:end-1,13); 0];
        A(:,6) = [-diff(A(:,3)); 0];
    
        %s*(v,Dv)
        temp_sa = A(:,8) + max(0, A(:,3).*A(:,9)+(A(:,3).*A(:,6))./A(:,11));
    
        %Aceleración
        A(:,4) = A(:,10).*(1-(A(:,3)./A(:,12)).^4-(temp_sa./A(:,5)).^2);
        A(I,4) = 0;
    
        %Integración balística:
        %velocidad
        A(:,3) = max(0, A(:,3) + A(:,4).*dt);
        A(I,3) = 0;

        %Posición
        A(:,1) = A(:,1) + max(0, A(:,3).*dt);
        A(I,1) = dummies(:,1);
    
        %Tiempo
        A(:,2) = A(:,2) + dt;
    
        %ordena todo (no hace falta sin adelantamientos y con dt pequeño)
        %A = sortrows(A,1);
    end
end