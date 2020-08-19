function [m,e,s,parameters] = isingmodel(J,L,kbT,Ne,BC,IC,Display,...
                                                      m_prev,e_prev,s_prev)
    %% Description 25/07/2020 (Renan Liupekevicius Carnielli, ra 157139)
   
    % This script computes the spin flips a certain amount of times for
    %  a 2D square with periodic boundary conditions or fixes spin-up 
    % boundary condition.
    % There is no external magnetic field applied to the system.

    % RETURN: 
    % m  - magnetization evolution
    % e  - energy per particle evolution
    % s -  spin matrix ( =lattice)
    % parameters  - Input

    % %% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J   ferromagnetic (>0) exchange interaction 
    % L   number of particles on an edge
    % kbT Boltzmant contant Kb times temperature
    % Ne  Number of iteractions (hopefully) sufficient to achive equil.          
    % BC  Boundary Condition "periodic" or "up"
    % IC  Inicial Condition "up", "down", "rand" or "previous
    % Display = "1" plot or "0" do not plot 
    
    %m_prev,e_prev,s_prev

    %% Fixed parameters
    % N       = L^2;                       %number of particles 
    % m(T) = | 0                          if T > Tc
    %        |(1-sinh(2J/kbT)^(-4))^(1/8) if T < Tc 

    %% Initialize
    tic
    
    % Magnetization and energy per particle
    m=[];
    e=[];
    
    % Initialize spin lattice
    switch IC
        case "up"   % Initialize spin-up lattice
            s   = ones(L,L, Ne);

        case "down" % Initialize spin down
            s   = -ones(L,L, Ne);

        case "rand" % Initialize random spin
            randspin = 2*randi([0 1],L) - ones(L,L);
            s        = repmat(randspin,1,1,Ne);

        case "previous" % Initialize with previous computed lattice
            % Number of input arguments
            switch nargin
                
                case 10
                m = m_prev;
                e = e_prev;
                s = repmat(s_prev(:,:,end),1,1,Ne);
                shape=size(s);
                Ntotal = shape(3) + Ne;
                
                otherwise
                error("Input arguments error."+...
                    " Try another initial condition IC.")
            end
        otherwise
                error("Choose the initial conditions 'up',"+...
                    " 'down', 'rand' or 'previous'")

     end

    
   % Indices randomly choosen
    randindex = randi(L,[2,Ne]);   % Random indices matrix


   % Random probability of flipping that is be compared to exp(-DE/kbT) 
    prob      = rand(Ne,1);


    
    % Initialize s_prev in case IC is different than "previous"
    if IC ~= "previous"
        s_prev=[];
    end
    
 


    %% Compute flips (hopefully) until equilibrium

switch Display
    
    case 1
            for k=2:Ne
                % Random index (i,j)
                i= randindex(1,k); j = randindex(2,k); 

                % Compute Delta E and flip/not flip
                DE = ComputeDeltaE(J,s,i,j,k,L,BC);   
                if DE < 0
                    s(i,j,k:end)=-s(i,j,k-1); % flip
                elseif prob(k) < exp(- DE / kbT)  
                    s(i,j,k:end)=-s(i,j,k-1); % flip
                else
                    % not flip
                end

                % Magnetization and Energy per particle evolution
                mk = Compute_m(s,k,L^2);
                ek = Compute_E(s,k,J,L,BC);

                % Magnetization and Energy per particle storage  
                m = [m mk];
                e = [e ek];

                % Selecting 3 digits for display
                mag               = num2str(mk,'%.3f');
                energyperparticle = num2str(ek,'%.3f');
                
                % Draw lattice
                titre="Computed "+ num2str(100*k/Ne,'%.f')+"%;"+...
                          " m = "+ mag+...
                         "; e = "+ energyperparticle + ...
                         "; kbT= "+ num2str(kbT,'%.2f');
                
                imagesc(s(:,:,k));
                set(gca,'fontsize', 18);
                title(titre);
                colorbar('Ticks',[-1,1],'TickLabels',{'down','up'})
                axis equal off;
                drawnow;
                
            end
        
    case 0
            for k=2:Ne

                % Random index (i,j)
                i= randindex(1,k); j = randindex(2,k); 

                % Compute Delta E and flip/not flip
                DE = ComputeDeltaE(J,s,i,j,k,L,"periodic");   
                if DE < 0
                    s(i,j,k:end)=-s(i,j,k-1); % flip
                elseif prob(k) < exp(- DE / kbT)  
                    s(i,j,k:end)=-s(i,j,k-1); % flip
                else
                    % not flip
                end

                % Magnetization and Energy per particle evolution
                mk = Compute_m(s,k,L^2);
                ek = Compute_E(s,k,J,L,"periodic");

                % Magnetization and Energy per particle storage  
                m = [m mk];
                e = [e ek];
            end
end


    %% return

    % Concatenate s
    if IC=="previous"
     s = cat(3,s_prev,s);
     shape=size(s);
     Ntotal = shape(3);
    else
        Ntotal=Ne;
    end

%     % Save variables
%     filename= w +'.mat';
%     save(filename,'m','e','s','parameters');
    
    % Save parameters of current execution
    parameters = struct( "L", L,"J",J,"kbT",kbT,...
                         "Ne",Ntotal,"BC",BC,"IC",IC);
   

    
    %Display message
    disp(" ising model has been succesfully computed");
    
    toc
end

function DE = ComputeDeltaE(J,s,i,j,k,l, BC)
%% Description
%
% Renan Liupekevicius Carnielli, ra 157139, 25/07/2020
%
% This function computes de energy variation from spin s to -s in a
% position (i,j) of a square lattice with periodic bounday condition.
%
% J - Exchange iteraction;
% S - Spin matrix that represents the lattice;
% i,j - Coordinates of the potential flipping spin;
% k - A snapshot of the system that will evolve to state
% k+1 by flipping a spin or not.
%
% shape = size(s);
% l     = shape(1);
% 
%% Equation 
% Initial first neighbors energy (before flip):
% Ei       = -J s    ( s_i+1j + s_i-1j + s_ij+1 +  s_ij-1 )
%
% Final first neighbors energy (after flip):
% Ef       = -J (-s) ( s_i+1j + s_i-1j + s_ij+1 +  s_ij-1 )
%
% Energy variation of the whole lattice is given by flipping s:
%
% DeltaE   = Ef - Ei
%          = 2Js ( s_i+1j + s_i-1j + s_ij+1 +  s_ij-1 )
%
%% Compute Delta E

switch BC
    
    case "periodic"
        % first lattice row (upper horizontal edge)
        if i == 1 
            if j == 1   % first element (upper left corner)
                DE =  s(i+1,j,k)+ s( l ,j,k) + s(i,j+1,k) + s(i, l ,k); 
            elseif j==l % last element  (upper right corner)
                DE =  s(i+1,j,k)+ s( l ,j,k) + s(i, 1 ,k) + s(i,j-1,k); 
            else        % middle elements
                DE =  s(i+1,j,k)+ s( l ,j,k) + s(i,j+1,k) + s(i,j-1,k);  
            end



        % last lattice row (lower horizontal edge)
        elseif i == l 
            if j == 1   % first element (lower left corner)
                DE = s( 1 ,j,k)+ s(i-1,j,k) + s(i,j+1,k) + s(i, l ,k); 
            elseif j==l % last element (lower right corner)
                DE = s( 1 ,j,k)+ s(i-1,j,k) + s(i, 1 ,k) + s(i,j-1,k); 
            else        % middle elements 
                DE = s( 1 ,j,k)+ s(i-1,j,k) + s(i,j+1,k) + s(i,j-1,k);          
            end



        % vertical edges and middle
        elseif j == 1     % first lattice column (middle elements) 
                DE = s(i+1,j,k)+ s(i-1,j,k) + s(i,j+1,k) + s(i, l ,k);



        elseif j == l % last lattice column (middle elements)
                DE  = s(i+1,j,k)+ s(i-1,j,k) + s(i, 1 ,k) + s(i,j-1,k) ;



        else          % middle elements
                DE  = s(i+1,j,k)+ s(i-1,j,k) + s(i,j+1,k) + s(i,j-1,k) ;

        end

    case "up"
        
        % first lattice row (upper horizontal edge)
        if i == 1 
            if j == 1   % first element (upper left corner)
                DE =  s(i+1,j,k)+      1     + s(i,j+1,k) +     1     ;
            elseif j==l % last element  (upper right corner)
                DE =  s(i+1,j,k)+      1     +      1     + s(i,j-1,k); 
            else        % middle elements (upper horizontal edge)
                DE =  s(i+1,j,k)+      1     + s(i,j+1,k) + s(i,j-1,k);  
            end



        % last lattice row (lower horizontal edge)
        elseif i == l 
            if j == 1   % first element (lower left corner)
                DE =      1    + s(i-1,j,k) + s(i,j+1,k) +      1    ; 
            elseif j==l % last element (lower right corner)
                DE =      1    + s(i-1,j,k) +      1     + s(i,j-1,k); 
            else        % middle elements (lower horizontal edge)
                DE =      1     + s(i-1,j,k) + s(i,j+1,k) + s(i,j-1,k);          
            end



        % vertical edges and middle
        elseif j == 1     % first lattice column (middle elements) 
                DE  = s(i+1,j,k)+ s(i-1,j,k) + s(i,j+1,k) +    1    ;



        elseif j == l % last lattice column (middle elements)
                DE  = s(i+1,j,k)+ s(i-1,j,k) +    1      + s(i,j-1,k) ;



        else          % middle elements
                DE  = s(i+1,j,k)+ s(i-1,j,k) + s(i,j+1,k) + s(i,j-1,k) ;
        end

end
%% Return
DE = 2*J*s(i,j,k)*DE;

end 


function E = Compute_E(s,k,J,L,BC)
%% Description
% Compute the energy with no external magnetic field using Ising's model.

switch BC
    
    case "periodic"

    % Edges (periodic BC)
    E = dot(s(1,:,k),s(L,:,k)) + dot(s(:,1,k),s(:,L,k));

    % Middle
    for i=1:L-1
        E = E + dot(s(i,:,k),s(i+1,:,k)) + dot(s(:,i,k),s(:,i+1,k));
    end
    
    case "up"
        
    E = sum(s(1,:,k)) + sum(s(L,:,k)) + sum(s(:,1,k)) + sum(s(:,L,k));

    for i=2:L-2

         E = E + dot(s(i,:,k),s(i+1,:,k)) + dot(s(:,i,k),s(:,i+1,k));
    end

end

E = -J*E/L^2;
end


function m = Compute_m(s,k,N)
%% Description 
% s - is a 2D-lattice of spins 1 or -1
% L - is an integer that represents the number of particles on an edge

 m = sum(s(:,:,k),"all") / N;
    
end


