function [sensor, id, Nsv] = sensordata10(todo,prob,tipo,row)

switch tipo
    case 0 % ELEGIMOS una fraccion de sensores
        % signal length
        id = unique(todo(:,row));
        N = length(id);

        % numero de sensores
        Nsv = ceil(prob*N);

        % random 0-1 signal
        x = zeros(N,1);
        q = randperm(N)';
        x(q(1:Nsv)) = sign(rand(Nsv,1));
        x = logical(x);

        id = id(x);

        tf = ismember(todo(:,row),id);
        sensor = todo(tf,:);
        
    case 1 % Elegimos una fracción de puntos
        N = length(todo(:,row));
        
        % numero de puntos
        Nsv = ceil(prob*N);

        % random 0-1 signal
        x = zeros(N,1);
        q = randperm(N)';
        x(q(1:Nsv)) = sign(rand(Nsv,1));
        x = logical(x);
        
        sensor = todo(x,:);
        id = unique(sensor(:,row));
end