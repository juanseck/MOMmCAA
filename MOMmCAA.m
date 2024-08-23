%Multi-Objective MmCAA

function [visited_X,visited_F] = MOMmCAA(params,MultiObj)
    
% Parameters of MORECCA
smartcells_no  = params.ns;
neighbors_no  = params.neigh;
max_iteration  = params.maxgen;
prop               = params.proporcion;
num_dec_inf   = params.n_dec_inf;
num_dec_sup  = params.n_dec_sup;
rep_num         = params.nrep;
num_inter       = params.num_inter;

% Parameters of cost function
cost                = MultiObj.fun;
dim                = MultiObj.nVar;
lb                   = MultiObj.var_min(:);
ub                  = MultiObj.var_max(:);

% Initialization
smart_cells = repmat((ub-lb)',smartcells_no,1).*rand(smartcells_no,dim) + repmat(lb',smartcells_no,1);

% evaluate smart-cells vector by vector
%smart_cells_fit  = cost(smart_cells);

num_func = size(cost(smart_cells(1,:)),2);
smart_cells_fit = zeros(size(smart_cells,1),num_func);
for i = 1 : size(smart_cells,1)
    smart_cells_fit(i,:) = cost(smart_cells(i,:));
end



% Initilizate repository of non-dominated smart-cells
[repository,repository_fit]=initialize_rep(smart_cells,smart_cells_fit);
% Calculate hypercube features
[hypercube_limits, hyper_ind,hyper_density] = hypercube(repository_fit,num_inter,lb,ub);
%Plot current Pareto front and smart-cells
plot_solutions(smart_cells_fit,repository_fit,hypercube_limits,MultiObj)
%Plot only Pareto front 
plot_pareto(repository_fit,hypercube_limits,MultiObj)

% store visited points
visited_X = transpose(smart_cells);
visited_F = transpose(smart_cells_fit);


%Iteration loop
for i =1:max_iteration
    %Proporcion decreciente
    proporcion = prop-((prop-1)*(i/max_iteration));
    %Smart-cell loop
    for j=1:smartcells_no
        %Take each smart-cell
        smart=smart_cells(j,:);
        calif=smart_cells_fit(j,:);
        density = hypercube_den(calif,num_inter,hypercube_limits,hyper_ind,hyper_density);
        %Vecino 1, smart-cell aleatoria del repositorio
        indA=randi(size(repository,1));
        vecino1=repository(indA,:);
        califv1=repository_fit(indA,:);
        %Vecino 2, smart-cell aleatoria del repositorio
        indB=randi(size(repository,1));
        while indA==indB
            indB=randi(size(repository,1));
        end
        vecino2=repository(indB,:);
        califv2=repository_fit(indB,:);
        %Funcion para crear vecinos
        %neigh = MmCAA_n(neighbors_no,smart,calif,vecino1,califv1,vecino2,califv2,rand*proporcion,num_dec_inf,num_dec_sup,lb,ub);
        neigh = MmCAA_n(neighbors_no,smart,calif,vecino1,califv1,vecino2,califv2,proporcion-((proporcion-2)*(i/max_iteration)),num_dec_inf,num_dec_sup,lb,ub);
        %Select the best neighbor
        [best_n, best_nf, best_d] = best_neigh(neigh,cost,smartcells_no,hypercube_limits,hyper_ind,hyper_density);
        %Update the smart-cell if neighbor is better
        if sum(best_nf==inf)==0
            res = is_better(best_nf, best_d, calif, density);
            if res ==1
                smart_cells(j,:)=best_n;
                smart_cells_fit(j,:)=best_nf;
                %Update repository
                [repository,repository_fit,flag] = update_repository_old(repository,repository_fit,rep_num,best_n,best_nf,num_inter,hypercube_limits,hyper_ind,hyper_density);
                %If repositiory has been improved, update hypercube features
                if flag == 1
                    [hypercube_limits, hyper_ind,hyper_density] = hypercube(repository_fit,num_inter,lb,ub);
                end
            end
        end
    end
    %Plot current Pareto front and smart-cells
    plot_solutions(smart_cells_fit,repository_fit,hypercube_limits,MultiObj)
    %Plot only Pareto front
    plot_pareto(repository_fit,hypercube_limits,MultiObj)    
    % store visited points
    visited_X = [visited_X,transpose(smart_cells)];
    visited_F = [visited_F,transpose(smart_cells_fit)];
end
%Print repository and smart-cells fitness

end

% Function that returns 1 if x dominates y and 0 otherwise
function d = domination(x,y)
    d = all(x<=y,2) & any(x<y,2);
end

% Function for checking the domination between the smart-cells.
% Returns a vector indicating if smart_cells are dominated (1) or not (0)
function smart_cells_dom = dominatedSet(smart_cells_fit)
    ns = size(smart_cells_fit,1);
    smart_cells_dom = zeros(ns,1);
    %Set of permutations of two elements, from 1 to ns
    perm = nchoosek(1:ns,2); 
    %Rest of permutations swapping original elements
    perm = [perm; [perm(:,2) perm(:,1)]];
    %Check domination for every pair of elements in the permutation array
    d = domination(smart_cells_fit(perm(:,1),:),smart_cells_fit(perm(:,2),:));
    %Obtain set of dominated smart-cells
    dominated = unique(perm(d==1,2));
    smart_cells_dom(dominated) = 1;
end

%Function to initialize the repository with at least two non-dominated
%smart-cell. If there are least than two non-dominated smart-cells, the
%rest is taken at random 
function [repository,repository_fit]=initialize_rep(smart_cells,smart_cells_fit)
% Initilizate repository of non-dominated smart-cells
% Dominated smart-cells
smart_cells_dom = dominatedSet(smart_cells_fit);
% Index of non-dominated smart-cells
ind_scd = ~smart_cells_dom;
% If there are less than 2 non-dominated solutions, the rest are selected at random
aux1 = sum(ind_scd);
if aux1 < 2
    %Indices of dominated solutions
    ind1 = find(~ind_scd);
    %Select 2-aux1 indices at random
    ind2 = randperm(size(ind1,1),2-aux1);
    %Change the logical values of ind_scd in positions ind2
    ind_scd(ind2)=~ind_scd(ind2);
end
% Initilizate repositories of smart-cells and their fitness value
repository = smart_cells(ind_scd,:);
repository_fit = smart_cells_fit(ind_scd,:);
end


% Function that calculates density and list of smart-cells in each
% hypercube
function [hypercube_limits, hyper_ind,hyper_density] = hypercube(repository_fit,num_inter,lb,ub)
    % Limits of each hypercube
    dim_n = size(repository_fit,2);
    hypercube_limits = zeros(2,dim_n);
    for dim = 1:dim_n
        hypercube_limits(1,dim) = min(repository_fit(:,dim));
        hypercube_limits(2,dim) = max(repository_fit(:,dim));
        %Check valid limits
        if hypercube_limits(1,dim) == -inf
            hypercube_limits(1,dim)=lb(dim);
        end
        if hypercube_limits(2,dim) == inf
            hypercube_limits(2,dim) =ub(dim);
        end
    end
    
    
    % Calculate the hypercube indices of each smart-cell
    ns = size(repository_fit,1);
    %Hypercube indices of each smart-cell
    repository_hinds = zeros(ns,dim_n);
    %Hypercube density of each index for smart-cell
    repository_hden = zeros(ns,dim_n);
    %Hypercube global indices of each smart-cell
    repository_lind = zeros(ns,1);
    
    %Calculate hypercube segments and density for each smart-cell
    %Dimensions
    for i = 1:dim_n
        paso = (hypercube_limits(2,i) - hypercube_limits(1,i))/num_inter;
        %Smart-cells
        for j = 1:ns
            index = ceil((repository_fit(j,i)-hypercube_limits(1,i))/paso);
            if index == 0
                index = index+1;
            elseif index >  num_inter
                index = index-1;
            end
            repository_hinds(j,i)=index;
        end
        %Calculate density  of segments for each dimension
        %Segments
        for j = 1:num_inter
            %Number of solutions in segment i
            aux=repository_hinds(:,i)==j;
            seg_n=sum(aux);
            if seg_n>0
                repository_hden(aux,i)=seg_n;
            end
        end
    end
    
    %Calculate linear indices of hypercubes for clasification
    %Smart-cells
    for i = 1:ns
        %Convert indices to cell array to manage dimensions separately
        C = num2cell(repository_hinds(i,:));
        %Calculate linear index
        repository_lind(i,1)=sub2ind(num_inter*ones(1,dim_n),C{:});
    end

    %Array of unique hypercubes indices
    hyper_ind = unique(repository_lind);
    %Number of solutions in each hypercube 
    hyper_density = zeros(size(hyper_ind,1),1);
    %Number of solutions in each hypercube
    for i = 1:size(hyper_ind,1)
        hyper_density(i) = sum(repository_lind==hyper_ind(i));
    end
    
    %Array of number of dimensions by number of intervals
    hyper_ds=zeros(dim_n,num_inter);
    %Vector of dimensions
    vec=1:dim_n;
    %For every solution, take the segment of each dimension cost
    for i=1:ns
        %Segments in each dimensions of each solution
        segs=repository_hinds(i,:);
        %Absoulte indices
        indices=((segs-1)*dim_n)+vec;
        if ~isreal(indices)
            indices = real(indices);
        end
        %Increse the values in the indices by one in hyper_ds
        hyper_ds(indices)=hyper_ds(indices)+1;
    end
end

%Function to obtain density of the hypercube associated to a smart-cell
function density = hypercube_den(costo,num_inter,hypercube_limits,hyper_ind,hyper_density)
%Number of dimensions
dim_n=length(costo);
%Hypercube indices
hind = zeros(1,dim_n);
%For each value, check the hypercube it belongs
for i=1:dim_n
    value = costo(i);
    %If the value doesn't belong to an existent hypercube, density=0
    if value < hypercube_limits(1,i) || value > hypercube_limits(2,i)
        density = 0;
        return
    else
        %Calculate the index of value in the hypercube i
        paso = (hypercube_limits(2,i)-hypercube_limits(1,i))/num_inter;
        index = ceil((value-hypercube_limits(1,i))/paso);
        if index == 0
            index = index+1;
        elseif index >  num_inter
            index = index-1;
        end
        hind(i)=index;
    end
end
%Obtain the linear index of the hypercube indices
C = num2cell(hind);
lind=sub2ind(num_inter*ones(1,dim_n),C{:});
%Find the density
aux = hyper_ind==lind;
%If density is not listed, density=0
if sum(aux)==0
    density = 0;
else
    %Obtain the density of the hypercube 
    density = hyper_density(aux);
end
end

% Function to update repositories of smart-cells and fitness values
function [repository,repository_fit,flag] = update_repository_old(repository,repository_fit,rep_num,smart,smart_fit,num_inter,hypercube_limits,hyper_ind,hyper_density)
% Flag to indicate if repository has been updated
flag = 0;
% Number of solutions in the repository
rep_n = size(repository,1);
% If repository is not full and smart-cell is not dominated, add smart-cell
% to the repository and put flag=1
if rep_n < rep_num
    d = domination(repository_fit,smart_fit);
    dom = sum(d,1);
    if dom==0
        repository(rep_n+1,:) = smart;
        repository_fit(rep_n+1,:) = smart_fit;
        flag = 1;
    end
else
    % If repository is full, smart-cell is not dominated or there are
    % dominated solutions by the smart-cell in the repository:
    d1 = domination(smart_fit,repository_fit);
    dom1 = sum(d1,1);
    if dom1 > 0
        %  Obtain hypercube density of the smart-cell
        smart_den = hypercube_den(smart_fit,num_inter,hypercube_limits,hyper_ind,hyper_density);
        %  Take density of the dominated solutions
        index = find(d1);
        dom1_den = zeros(dom1,1);
        for i = 1:dom1
            dom1_den(i) = hypercube_den(repository_fit(index(i),:),num_inter,hypercube_limits,hyper_ind,hyper_density);
        end
        %  Take solution with maximum density and replace it with the
        %  smart-cell if it has greater density
        [max_den, max_index] = max(dom1_den);
        if max_den*2.1 >= smart_den
            repository(index(max_index),:) = smart;
            repository_fit(index(max_index),:) = smart_fit;
            %  Put flag=1
            flag=1;
        end
    end
    % Not dominated
    d2 = domination(repository_fit,smart_fit);
    dom2 = sum(d2,1);
    if dom2 == 0
        %  Obtain hypercube density of the smart-cell
        smart_den = hypercube_den(smart_fit,num_inter,hypercube_limits,hyper_ind,hyper_density);
        %  Take density of the other solutions
        dom2_den = zeros(rep_n,1);
        for i = 1:rep_n
            dom2_den(i) = hypercube_den(repository_fit(i,:),num_inter,hypercube_limits,hyper_ind,hyper_density);
        end
        %  Take solution with maximum density and replace it with the
        %  smart-cell if it has greater density
        [max_den, max_index] = max(dom2_den);
        if max_den*1.8 >= smart_den
            repository(max_index,:) = smart;
            repository_fit(max_index,:) = smart_fit;
            %  Put flag=1
            flag=1;
        end
    end
end
end

% Optimizer

%Function to calculate neighborhoods of a smart-cell, returns the set neighbors.
function [neigh] = MmCAA_n(neighbors_no,smart,calif,vecino1,califv1,vecino2,califv2,proporcion,num_dec_inf,num_dec_sup,lb,ub)
%Set of neighborhoods
neigh = zeros(neighbors_no, length(smart));
%Numero de reglas
nr=14;
%Neighborhoods loop
for i=1:neighbors_no
    %Elegir regla de forma aleatoria
    regla_ale=randi([1,nr]);
    %Reglas de mutación por mayoría o minoría
    if regla_ale == 1
        evolucion = regla_mayoria_solo(smart,proporcion); % Explotación
    elseif regla_ale == 2 %Explotación
        evolucion = regla_redondeo_solo(smart,calif,califv1,randi([num_dec_inf,num_dec_sup]));
    elseif regla_ale >= 3 && regla_ale <= 4 
        evolucion = regla_minoria_solo(smart,proporcion); % Exploración
    elseif regla_ale >= 5 && regla_ale <= 6
        evolucion = regla_mayoria_vecino(smart,calif,vecino1,califv1,proporcion); %Exploración
    elseif regla_ale >= 7 && regla_ale <= 8 
        evolucion = regla_minoria_vecino(smart,calif,vecino1,califv1,proporcion); %Exploración
    elseif regla_ale >= 9 && regla_ale <= 14
        % Seleccionar una solución elitista aleatoria
        evolucion = regla_minoria_tres_vecinos(smart,calif,vecino1,califv1,vecino2,califv2,proporcion); %Exploración
    end
    %Pasar valores de evolucion a intervalo valido
    ind_Neg=evolucion'<lb;
    if sum(ind_Neg)>0
        if length(lb)>1
            evolucion(ind_Neg)=lb(ind_Neg)+(((ub(ind_Neg)-lb(ind_Neg))/10)*rand);
        else
            evolucion(ind_Neg)=lb+((ub-lb)/10)*rand;
        end
    end
    ind_Pos=evolucion'>ub;
    if sum(ind_Pos)>0
        if length(ub)>1
            evolucion(ind_Pos)=ub(ind_Pos)-(((ub(ind_Pos)-lb(ind_Pos))/10)*rand);
        else
            evolucion(ind_Pos)=ub-((ub-lb)/10)*rand;
        end
    end
    %Add evolucion to list of neighborhoods
    neigh(i,:)=evolucion();
end
end

%Function to calculate the best neighbor of a set
function [best_n, best_nf, best_d] = best_neigh(neigh,fcosto,smartcells_no,hypercube_limits,hyper_ind,hyper_density)
%Number of neighbors
nn = size(neigh,1);
%Select the best neighbor
for i=1:nn
    sol = neigh(i,:);
    cost = fcosto(sol);
    %Obtain density of hypercube of the neighbor
    density = hypercube_den(cost,smartcells_no,hypercube_limits,hyper_ind,hyper_density);
    %Obtain best neighbor
    if i==1
        best_n = sol;
        best_nf = cost;
        best_d = density;
    else
        %Check if new solution is better than the current one
        res = is_better(cost, density, best_nf, best_d);
        if res==1 
            %Update the best neighbor
            best_n = sol;
            best_nf = cost;
            best_d = density;
        end
    end
end
end

%Function to check if a solution with cost1 and den1 is better than other one
%by domination and with cost2 and den2 hypercube density
function res = is_better(cost1, den1, cost2, den2)
%Initialize response as 0
res = 0;
if domination(cost1,cost2)==1 || domination(cost2,cost1)==0
    %Check if density is at most equal to the current one
    if den1<=den2
        %Update response as 1 meaning that solution 1 is better than
        %solution 2
        res = 1;
    end
end
end

%Function to plot the response of the smart-cells and the repository-fit
function plot_solutions(smart_fit,repository_fit,hypercube_limits,MultiObj)
%Number of dimensions
dim_n = size(repository_fit,2);
%Figure
figure(1);
%Case for a 2-dim Pareto front
if(dim_n==2)
    %Plot smart-cells
    plot(smart_fit(:,1),smart_fit(:,2),'ob');
    hold on;
    %Plot repository
    plot(repository_fit(:,1),repository_fit(:,2),'xr');
    %If exists, plot Pareto Front
    if(isfield(MultiObj,'truePF'))
        plot(MultiObj.truePF(:,1),MultiObj.truePF(:,2),'.','color',[0.9 0.9 0.9]);
    end
    %Set axis properties
    try
        set(gca,'xtick',hypercube_limits(:,1)','ytick',hypercube_limits(:,2)');
    end
    %axis([min(hypercube_limits(:,1)) max(hypercube_limits(:,1)) min(hypercube_limits(:,2)) max(hypercube_limits(:,2))]);
    grid on; 
    xlabel('f1'); 
    ylabel('f2');
    %Release figure
    hold off
end
if(dim_n==3)
    %Plot smart-cells
    plot3(smart_fit(:,1),smart_fit(:,2),smart_fit(:,3),'ob'); 
    hold on;
    %Plot repository
    plot3(repository_fit(:,1),repository_fit(:,2),repository_fit(:,3),'xr'); 
    %If exists, plot Pareto Front
    if(isfield(MultiObj,'truePF'))
        plot3(MultiObj.truePF(:,1),MultiObj.truePF(:,2),MultiObj.truePF(:,3),'.','color',[0.9 0.9 0.9]);
    end
    %Set axis properties
    set(gca,'xtick',hypercube_limits(:,1)','ytick',hypercube_limits(:,2)','ztick',hypercube_limits(:,3)');
    grid on;
    xlabel('f1');
    ylabel('f2');
    zlabel('f3');
    %axis square;
    %Release figure
    hold off
end
end

%Function to plot the response of the smart-cells and the repository-fit
function plot_pareto(repository_fit,hypercube_limits,MultiObj)
%Number of dimensions
dim_n = size(repository_fit,2);
%Figure
figure(2);
%Case for a 2-dim Pareto front
if(dim_n==2)
    %Plot repository
    plot(repository_fit(:,1),repository_fit(:,2),'xr');
    hold on;
    %If exists, plot Pareto Front
    if(isfield(MultiObj,'truePF'))
        plot(MultiObj.truePF(:,1),MultiObj.truePF(:,2),'.','color',[0.9 0.9 0.9]);
    end
    %Set axis properties
    try
        set(gca,'xtick',hypercube_limits(:,1)','ytick',hypercube_limits(:,2)');
    end
    %axis([min(hypercube_limits(:,1)) max(hypercube_limits(:,1)) min(hypercube_limits(:,2)) max(hypercube_limits(:,2))]);
    grid on; 
    xlabel('f1'); 
    ylabel('f2');
    %Release figure
    hold off
end
if(dim_n==3)
    %Plot repository
    plot3(repository_fit(:,1),repository_fit(:,2),repository_fit(:,3),'xr'); 
    hold on;
    %If exists, plot Pareto Front
    if(isfield(MultiObj,'truePF'))
        plot3(MultiObj.truePF(:,1),MultiObj.truePF(:,2),MultiObj.truePF(:,3),'.','color',[0.9 0.9 0.9]);
    end
    %Set axis properties
    set(gca,'xtick',hypercube_limits(:,1)','ytick',hypercube_limits(:,2)','ztick',hypercube_limits(:,3)');
    grid on;
    xlabel('f1');
    ylabel('f2');
    zlabel('f3');
    %axis square;
    %Release figure
    hold off
end
end

%Reglas

%Regla para acercar cada elemento del vector al valor mayoritario
function [evolucion] = regla_mayoria_solo(smart,proporcion)
evolucion=smart;
[x,e]=histcounts(smart);
[~,ind]=max(x);
valor=sum(e(ind:ind+1))/2;
%nueva distancia
dist=(smart-valor)*proporcion*rand;
evolucion=evolucion-dist;
end

%Regla para acercar cada elemento del vector al valor minoritario
function [evolucion] = regla_minoria_solo(smart,proporcion)
evolucion=smart;
[x,e]=histcounts(smart);
[~,ind]=min(x);
valor=sum(e(ind:ind+1))/2;
%nueva distancia
dist=(smart-valor)*proporcion*rand;
evolucion=evolucion-dist;
end

%Regla para redondeo en los valores de la smart_cell
function [evolucion] = regla_redondeo_solo(smart,calif,mejor_calif,num_dec)
evolucion=smart;
%suma de calificaciones
suma=calif+mejor_calif;
%ponderacion de calificación del vecino
pond=1-mean(calif./suma);
%Hacer redondeo
for i=1:length(smart)
    if rand<=pond
        evolucion(i)=round(evolucion(i),num_dec);
    end
end
end

%DOS VECINOS
%Regla para tomar mayoria del vecino para ponderar aleatoriamente dependiendo costos de los dos vecinos
function [evolucion] = regla_mayoria_vecino(smart,calif1,vecino,calif2,proporcion)
%La evolucion de la smart-cell 
evolucion=smart;
%suma de calificaciones
suma=calif1+calif2;
%ponderacion de calificación del vecino
pond=1-mean(calif2./suma);
%elemento mayoritaria del vecino
[x,e]=histcounts(vecino);
[~,ind]=max(x);
valor=sum(e(ind:ind+1))/2;
%Si se selecciona al vecino, se suma al elemento del vector multiplicado
%por un aletorio r entre -dist_mayor/2 y dist_mayor/2
r=(rand*proporcion)-(proporcion/2);
for i=1:length(smart)
    if rand<=pond
        evolucion(i)=evolucion(i)+(r*valor);
    end
end
end

%Regla para tomar minoría del vecino para ponderar aleatoriamente dependiendo costos de los dos vecinos
function [evolucion] = regla_minoria_vecino(smart,calif1,vecino,calif2,proporcion)
%La evolucion de la smart-cell 
evolucion=smart;
%suma de calificaciones
suma=calif1+calif2;
%ponderacion de calificación del vecino
pond=1-mean(calif2./suma);
%elemento mayoritaria del vecino
[x,e]=histcounts(vecino);
[~,ind]=min(x);
valor=sum(e(ind:ind+1))/2;
%Si se selecciona al vecino, se suma al elemento del vector multiplicado
%por un aletorio r entre -dist_mayor/2 y dist_mayor/2
r=(rand*proporcion)-(proporcion/2);
for i=1:length(smart)
    if rand<=pond
        evolucion(i)=evolucion(i)+(r*valor);
    end
end
end

%Regla para tomar minoría de tres vecinos ponderando aleatoriamente sus costos
function [evolucion] = regla_minoria_tres_vecinos(smart,calif,vecino1,calif1,vecino2,calif2,proporcion)
%La evolucion de la smart-cell 
evolucion=smart;
valor=smart;
%suma de calificaciones
suma=calif+calif1+calif2;
%ponderacion de calificación del vecino
pond=1-mean(calif2./suma);
%elemento mayoritario para cada posición en los tres vecinos
%Obtener distancias entre vecinos
dist1=abs(smart-vecino1);
dist2=abs(vecino1-vecino2);
dist3=abs(vecino2-smart);
dist=[dist1; dist2; dist3];
%Obtener indice de distancia menor por cada posicion
[~,pos]=min(dist);
%Obtener valor que no esté a mínima distancia
vecindad=[smart; vecino1; vecino2; smart];
for i=1:length(smart)
    ind=pos(i)-1;
    if ind==0
        ind=3;
    end
    valor(i)=vecindad(ind,i);
end
%Si se selecciona la posición para modificación, se suma al elemento del vector multiplicado
%por un aletorio r entre -dist_mayor/2 y dist_mayor/2
r=(rand*proporcion)-(proporcion/2);
for i=1:length(smart)
    if rand<=pond
        evolucion(i)=evolucion(i)+(r*valor(i));
    end
end
end
