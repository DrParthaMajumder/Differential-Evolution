clc
clear all
close all
format long g

%%
% [LB,UB,D,fobj] = Get_Functions_details('Ackley_F1');  % Fg=0
% [LB,UB,D,fobj] = Get_Functions_details('Beale_F2');   % Fg=0 
% [LB,UB,D,fobj] = Get_Functions_details('Bohachevsky_F3'); % Fg=0

% [LB,UB,D,fobj] = Get_Functions_details('Booth_F4');   % Fg=0;
% [LB,UB,D,fobj] = Get_Functions_details('BUKINN6_F5'); % Fg=0  
% [LB,UB,D,fobj] = Get_Functions_details('Colville_F6'); % Fg=0;  

% [LB,UB,D,fobj] = Get_Functions_details('Cross_In_Tray_F7');  % Fg=-2.06261;
% [LB,UB,D,fobj] = Get_Functions_details('DejongN5_1_F8');  % Fg=0
% [LB,UB,D,fobj] = Get_Functions_details('Dixonprice_F9'); % Fg=0;

% [LB,UB,D,fobj] = Get_Functions_details('Drop_Wave_F10');  % Fg=-1
% [LB,UB,D,fobj] = Get_Functions_details('EASOM1_F11');   % Fg=-1;
% [LB,UB,D,fobj] = Get_Functions_details('Eggholder_F12');  % Fg=959.6407;

% [LB,UB,D,fobj] = Get_Functions_details('GoldsteinPrice_F13');   % Fg=3;
% [LB,UB,D,fobj] = Get_Functions_details('GoldsteinPrice_Scaled_F14');  % Fg=3;
% [LB,UB,D,fobj] = Get_Functions_details('Griewank_F15');  % Fg=0

% [LB,UB,D,fobj] = Get_Functions_details('Hartmann_3D_F16');  % Fg=-3.67597355769227
% [LB,UB,D,fobj] = Get_Functions_details('Hartmann_4D_F17');  % F(g) = -3.135474
% [LB,UB,D,fobj] = Get_Functions_details('Hartmann_6D_F18');  % Fg=-3.04245773783059

% [LB,UB,D,fobj] = Get_Functions_details('Holder_Table_F19');  % Fg=-19.2085
% [LB,UB,D,fobj] = Get_Functions_details('Langermann_F20');  % Fg=-4.05404569816266;
% [LB,UB,D,fobj] = Get_Functions_details('Levy_F21');  %Fg=0;

% [LB,UB,D,fobj] = Get_Functions_details('LevyN13_F22');  % Fg=0;
% [LB,UB,D,fobj] = Get_Functions_details('Matyas_F23'); % Fg=0
% [LB,UB,D,fobj] = Get_Functions_details('Mccormick_F24'); % Fg=-1.9133; 

%  [LB,UB,D,fobj] = Get_Functions_details('Michalewicz1_F25'); % Fg=-9.66015;
% 26:Combinatorial Optimization
% 27: Combinatorial Optimization

% [LB,UB,D,fobj] = Get_Functions_details('Powell_F28'); % Fg=0 
% [LB,UB,D,fobj] = Get_Functions_details('Sum_Power_F29'); % Fg=0;
% [LB,UB,D,fobj] = Get_Functions_details('Rastrigin1_F30'); % Fg=0


% [LB,UB,D,fobj] = Get_Functions_details('Rosenbrock1_F31'); % Fg=0
% [LB,UB,D,fobj] = Get_Functions_details('Rotted_hyper_ellipsoid_F32');  %Fg=0
% [LB,UB,D,fobj] = Get_Functions_details('SchafferN2_F33'); % Fg=0

% [LB,UB,D,fobj] = Get_Functions_details('SchafferN4_F34'); % Fg=0.292579
[LB,UB,D,fobj] = Get_Functions_details('Schwef1_F35'); % Fg=0
% 36: Combinatorial Optimization



% [LB,UB,D,fobj] = Get_Functions_details('Shubert_F37'); % Fg=-186.7309
% [LB,UB,D,fobj] = Get_Functions_details('Six_Hump_Camel_F38');  % Fg=-1.031628453486
% [LB,UB,D,fobj] = Get_Functions_details('Sphere1_F39'); % Fg=0;

%[LB,UB,D,fobj] = Get_Functions_details('Styblinski_Tang_F40');  % Fg=-39.16599*D
% [LB,UB,D,fobj] = Get_Functions_details('Sum_Square_Function_F41'); % Fg=0
% [LB,UB,D,fobj] = Get_Functions_details('Three_Hump_Camel_F42');  % Fg=0;

% [LB,UB,D,fobj] = Get_Functions_details('Trid_F43'); % Fg=-50; D=6; % Fg=-200; D=10;


%%


%% DE Parameters
itmax=100;        % Maximum Number of Iterations
N=200;            % Population Size

beta_min=0.2;     % Lower Bound of Scaling Factor
beta_max=0.8;     % Upper Bound of Scaling Factor

pCR=0.2;          % Crossover Probability
VarSize=[1 D];

%% Bounds of Variable:

itmax=100; % Maximum numbef of iterations
N=100;

if length(LB)==1
for kk=1:1:D
    lb(1:N,kk)=LB; 
    ub(1:N,kk)=UB;
end
end

if length(LB)~=1
for kk=1:1:D
    lb(1:N,kk)=LB(kk); 
    ub(1:N,kk)=UB(kk);
    
end
end


x=lb+(ub-lb).*rand(N,D);
for ii=1:1:N
F(ii) = fobj(x(ii,:),D);
end
[F_g_best,pp]=min(F);
g_best=x(pp,:);


break_point=1;

%% Iteration

for it=1:1:itmax
    
    for ii=1:1:N
        x_vecT=x(ii,:);
        vect=randperm(N);
        vect(vect==ii)=[];
        
        ri1=vect(1);
        ri2=vect(2);
        ri3=vect(3);
        
        %% Mutation
        beta=unifrnd(beta_min,beta_max,VarSize);
        y=x(ri1,:)+beta.*(x(ri2,:)-x(ri3,:));    
        y = max(y, lb(1,:));
		y = min(y, ub(1,:));
        
        %% Crossover
        z_c=zeros(size(x_vecT));       %Crossover vector
        J0=randi([1 length(x_vecT)]);
        
        for jj=1:length(x_vecT)
            if jj==J0 || rand<=pCR
                z(jj)=y(jj);
            else
                z(jj)=x_vecT(jj);     
            end
        end
        
        Fz= fobj(z,D);      
        
        
        if  Fz<F(ii)
            x(ii,:)=z;
            F(ii)=Fz;
            
        if F(ii)<F_g_best
               F_g_best=F(ii);
               g_best=x(ii,:);
        end                      
        end   
    end
     for kk=1:1:D
         x(:,kk)=min(x(:,kk),ub(:,kk));
         x(:,kk)=max(x(:,kk),lb(:,kk));
     end
F_g_best_k(it)=F_g_best;
end
F_g_best=F_g_best
g_best=g_best

plot(F_g_best_k)
grid on;
break_point=2;









