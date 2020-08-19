% clear all;
 close all;


%% Description 26/07/2020 (Renan Liupekevicius Carnielli, ra 157139)

% This is an assignment of F604 Statistical Physics.
% University of Campinas.

 
% This script displays the behavior of an evolving 2D square lattice.
% Each cell has a spin up or a spin down. There is an interaction 
% between the spins using Ising Model.

%% Constants
kb   = 1.380649e-23;             % Boltzmann Constant
kbTc = 2/log(1+sqrt(2));      % Critical temperature (2D) times Kb



%% Parameters
J   = 1;      % ferromagnetic (>0) exchange interaction 
L   = 25;    % number of particles on an edge
N   = L^2;   % number of particles 
kbT = 2;     % Boltzmant contant Kb times temperature
Ne  = 10000;  % Number of flips           
BC  = "periodic"; % Boundary Condition "periodic" or "up"
IC  = "u";   % Inicial Condition "up", "down", "rand" or "previous"
    

 
%% Compute flips via metropolis algorithm  
switch IC
    case "previous"  
     
     % Parameters
     J   = parameters.J;     
     L   = parameters.L;      
     N   = L^2;               
     kbT = parameters.kbT;    
     BC  = parameters.BC ;    
     
     % Continue simulation
     [m,e,s,parameters] = isingmodel(J,L,kbT,Ne,BC,IC,1,m,e,s);
     
    otherwise
     % New simulation
     [m,e,s,parameters] = isingmodel(J,L,kbT,Ne,BC,IC,1);
end
 


%% Plot: magnetizaion, energy, suceptibilty and heat capacity
 
 % Set time-step Window to compute moving average
 fenetre = 10000;
 
 % Compute physical quantities from 'initial_index' until the end of the
 % vector
 initial_index =round(length(m)/2);
 
 % Compute magnetization (hopefully after equilibrium)
 m_mean  = mean(m(initial_index:2*end/3));
 
 % Compute energy per particle (hopefully after equilibrium)
 e_mean  = mean(e(initial_index:end));
 
 % Compute magnetization variance: mag. suceptibility 
 chi     = L^2 * 1/kbT     * var(m(initial_index:end)); 
 
 % Compute energy variance: Heat Capacity over kb:   C/kb^
 CokboN    = L^2 * 1/(kbT^2) * var(e(initial_index:2*end/3));
  
 % Magnetization plot
 figure(2)
 hold on;
 plot(m,'k');
 %plot(movmean(m,fenetre))
 box on;
 grid on;
 set(gca,'fontsize', 18);
 xlabel('Number of iterations');
 ylabel('Magnetization');
 title("$k_bT= $"+num2str(kbT)+...
       ";  $\chi =$ "+num2str(chi,'   %.2f')+ ...
       ";  $m    =$ "+num2str(m_mean,'%.2f'), ...
                            'Interpreter','latex');
 hold off;
 
 %Energy per particle plot
 figure(3)
 hold on;
 plot(e,'k');
 %plot(movmean(e,fenetre))
 box on;
 grid on;
 set(gca,'fontsize', 18);
 xlabel('Number of iterations');
 ylabel('Energy per particle');
 title("$k_bT= $"+num2str(kbT)+...
       ";   $C/N k_b =$ "+num2str(CokboN,'%.2f')+ ...
       ";   $e =$ "+num2str(e_mean,  '%.2f'),     ...
                                'Interpreter','latex');
 hold off;
 
 
 % Final state (snapshot) plot
 figure(7)
 titre="Lattice snapshot";
 imagesc(s(:,:,end));
 set(gca,'fontsize', 18);
 title(titre);
 colorbar('Ticks',[-1,1],'TickLabels',{'down','up'})
 axis equal off;
 drawnow;
 
 
 

 %% Store workspace file
 
 % set file name
 filename = "parameters_L"+ parameters.L +"_kbT"+ parameters.kbT +...
            "_Ne"+ parameters.Ne+ "_BC"+ parameters.BC;
 
 %save the workspace as matlab file
 save(filename+".mat", '-v7.3');
 
 
 
 
 %% Analytic expressions 
    
    % temperature vector and critical temperature
    kbT_vec   = 0:0.01:5;              % Temperature times Kb (vector)

    
    % Analytic magnectization of 2D ising model
    ma  = [(1-sinh(2*J./kbT_vec(1:227)).^(-4)).^(1/8) zeros(1,501-227)];
    
    
    % Analytic heat capacity
    C_okb_oN_analytic = - 2/pi * (2*J/kbTc)^2* log (abs(1-kbT_vec/kbTc));
    
    
    % Analytic energy per particle
    beta  = 1./kbT_vec;
    kappa = 2 * ( sinh(2 * beta *J) ) ./ ( cosh(2 * beta  * J) ).^2 ;
    
    K1    = zeros(1,length(beta));
    x     = 0:0.01:pi/2;
    
    for i= 1:length(beta)
    fun      = @(x)  1 ./ (    sqrt( 1- kappa(i)^2 * (sin(x)).^2 )   );
    K1(i)    = trapz(x,fun(x));
    end
    
    e_analytic = -2*J* tanh( 2 * beta * J) ...
                 - J * ( sinh(2 *beta * J).^2 -1 )...
                 ./( sinh(2 * beta * J) .* cosh(2 * beta * J) ) ...
                 .* ( 2 / pi .* K1 -1 );
             
 %% Plot analytic expressions    
             
    % Plot analytic magnetization
    figure(4)
    hold on;
    box on;
    plot(kbT_vec,ma,'k','LineWidth',2);
   % plot(kbT_vec,ones(size(kbT_vec)),'k--','LineWidth',2);
    grid on;
    set(gca,'fontsize', 25);
    xlabel(' K_bT');
    ylabel('$m$', 'interpreter','latex');
    hold off;
 

    %Plot analytic heat capacity
    figure(5)
    hold on;
    box on;
    plot(kbT_vec,C_okb_oN_analytic,'k','LineWidth',2);
    grid on;
    set(gca,'fontsize', 25);
    xlabel(' K_bT');
    ylabel('$C/N k_b$', 'interpreter','latex');
    hold off;
 
    
    %Plot analytic energy per particle
    figure(6)
    hold on;
    box on;
    plot(kbT_vec,e_analytic,'k','LineWidth',2);
    grid on;
    set(gca,'fontsize', 25);
    xlabel(' K_bT');
    ylabel('$E/N $', 'interpreter','latex');
    hold off;
   
    
  
 




