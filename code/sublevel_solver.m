function [level] = sublevel_solver(level,timestep,tspan,initial,lvl,position)
  % Funzione che, presa la struttura level, il timestep corrente (del passo finale),
  % la struttura initial (definisce il tipo di dato iniziale)
  % il livello di lavoro e la posizione dell'elemento di lavoro, restituisce
  % la struttura level aggiornata in level{lvl,position}.solution = soluzione
  % calcolata nel sottolivello al timestep corrente
  
  x_max = level{lvl,position}.domain_position(4);
  x_min = level{lvl,position}.domain_position(2);
  y_max = level{lvl,position}.domain_position(1);
  y_min = level{lvl,position}.domain_position(3);
  
  m_x = ((level{lvl,position}.values(4)-level{lvl,position}.values(2))+1);
  m_y = ((level{lvl,position}.values(1)-level{lvl,position}.values(3))+1);
  h_y = (y_max-y_min)/(m_y-1);
  h_x = (x_max-x_min)/(m_x-1);
  
  x = linspace(x_min,x_max,m_x);               
  y = linspace(y_min,y_max,m_y);             

  
  [X,Y] = meshgrid(x,y);
  X = flip(X);
  Y = flip(Y);
  
  % meshgrid del livello padre del livello di lavoro
  father_pos = level{lvl,position}.father;
  X_father = level{lvl-1,father_pos}.Xmesh;
  Y_father = level{lvl-1,father_pos}.Ymesh;
  %x_max_father = level{lvl-1,father_pos}.domain_position(4);
  %x_min_father = level{lvl-1,father_pos}.domain_position(2);
  %y_max_father = level{lvl-1,father_pos}.domain_position(1);
  %y_min_father = level{lvl-1,father_pos}.domain_position(3);
  %m_x_father = (level{lvl-1,father_pos}.values(4)...
  %              -level{lvl-1,father_pos}.values(2))+1;
  %m_y_father = (level{lvl-1,father_pos}.values(1)...
  %              -level{lvl-1,father_pos}.values(3))+1;
  %x_father = linspace(x_min_father,x_max_father,m_x_father);
  %y_father = linspace(y_min_father,y_max_father,m_y_father);    
  %[X_father,Y_father] = meshgrid(x_father,y_father);
  %X_father = flip(X_father);
  %Y_father = flip(Y_father);
  
  
% CONDIZIONE INIZIALE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  U0 = NaN(m_y,m_x);
  % Se siamo al primo passo temporale, per condizione iniziale calcoliamo 
  % direttamente nel sottoproblema il dato iniziale
  if (timestep == 2)
    %fprintf("\nDato iniziale\n")
    a(1) = 0.340107907001967147606775060166996;
    b(1) = (((((722731008*a(1)-326467584)*a(1)-13427712)*a(1)+11551104)*a(1)+...
           834006)*a(1)-12183)/(((((2972712960*a(1)-1362493440)*a(1)-8744960)*...
           a(1)+19299840)*a(1)+8788080)*a(1)-174580);
    a(2) = a(1)*b(1) - 1/4*a(1);
    b(2) = ((((737280*a(1)+209920)*a(1)-134720)*a(1)-8960)*b(1)+...
           ((-364544*a(1)+70144)*a(1)+18256)*a(1)+393)/((2457600*a(1)-537600)*...
           a(1)-105600);
    a(3) = a(1)*(192*b(2)-48*b(1)+16*a(1)+5)/192;
    b(3) = ((61440*a(1)-9600)*b(2)+((-30720*a(1)+640)*a(1)+560)*b(1)+(8448*a(1)-...
           1056)*a(1)-21)/(368640*a(1)-92160);
    a(4) = a(1)*(4608*b(3)-1152*b(2)+(384*a(1)+120)*b(1)-128*a(1)-7)/4608;
    % Determinare la condizione da imporre , prendi dato iniziale
    % initial.type = "single vortex"
    %                "rotating vortex" 
    %                "traslating vortex"
    if (initial.type == "sin")
        R2 = (X-initial.x0).^2+(Y-initial.y0).^2;
        rho = (((a(4)*R2+a(3)).*R2+a(2)).*R2+a(1)).*R2./...
              ((((a(4)*R2+b(3)).*R2+b(2)).*R2+b(1)).*R2+1);
        U0 = sqrt(rho).*exp(-1i*atan2(Y-initial.y0,X-initial.x0));  
    else
      if (initial.type == "rot")
          R2_1 = (X-initial.x1).^2+(Y-initial.y1).^2;
          rho_1 = (((a(4)*R2_1+a(3)).*R2_1+a(2)).*R2_1+a(1)).*R2_1./...
                  ((((a(4)*R2_1+b(3)).*R2_1+b(2)).*R2_1+b(1)).*R2_1+1);
          U0_1 = sqrt(rho_1).*exp(1i*atan2(Y-initial.y1,X-initial.x1));
          R2_2 = (X-initial.x2).^2+(Y-initial.y2).^2;
          rho_2 = (((a(4)*R2_2+a(3)).*R2_2+a(2)).*R2_2+a(1)).*R2_2./...
                        ((((a(4)*R2_2+b(3)).*R2_2+b(2)).*R2_2+b(1)).*R2_2+1);
          U0_2 = sqrt(rho_2).*exp(1i*atan2(Y-initial.y2,X-initial.x2));
          U0 = U0_1.*U0_2;
      elseif(initial.type == "tra")
          R2_1 = (X-initial.x1).^2+(Y-initial.y1).^2;
          rho_1 = (((a(4)*R2_1+a(3)).*R2_1+a(2)).*R2_1+a(1)).*R2_1./...
                  ((((a(4)*R2_1+b(3)).*R2_1+b(2)).*R2_1+b(1)).*R2_1+1);
          U0_1 = sqrt(rho_1).*exp(-1i*atan2(Y-initial.y1,X-initial.x1));
          R2_2 = (X-initial.x2).^2+(Y-initial.y2).^2;
          rho_2 = (((a(4)*R2_2+a(3)).*R2_2+a(2)).*R2_2+a(1)).*R2_2./...
                        ((((a(4)*R2_2+b(3)).*R2_2+b(2)).*R2_2+b(1)).*R2_2+1);
          U0_2 = sqrt(rho_2).*exp(1i*atan2(Y-initial.y2,X-initial.x2));
          U0 = U0_1.*U0_2;
      else
          warning("Dato iniziale incorretto: controllare l'input in GP.m")
      end  
    end    
  else
    %fprintf("\nInterpolazione\n")
    %                              riprendiamo il dato iniziale della matrice padre prima del timestep
    %%fprintf("Il valore iniziale è per la posizione %d e tempo %d",position, timestep)
    %%fprintf("\n")
    U0 = interp2(X_father,Y_father,level{lvl-1,father_pos}.sol_before_timestep,...
                 X,Y,"cubic");
                 %pause(1)
                 %clc
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Per la condizione al bordo, determiniamo le condizioni al bordo come condizione
% di Dirichlet, determinata come interpolazione di Hermite fra la "soluzione" 
% del sottoproblema al tempo tspan(timestep-1) (CONDIZONE INIZIALE, INTERPOLATA O 
% CALCOLATA) e al tempo tspan(timestep) (CONDIZIONE INTERPOLATA)
% PER CALCOLARE LE DERIVATE NECESSARIE ALL'INTERPOLAZIONE DI HERMITE:
% INTERPOLAZIONE PER TEMPO INIZIALE E FINALE O NUOVO SISTEMA DISCRETIZZATO ?
% PRIMO APPROCCIO

% Calcoliamo in primis l'interpolazione della "soluzione" del sottoproblema al tempo 
% tspan(timestep), nei punti al bordo
% FORSE È MEGLIO TAGLIARE PRIMA DELL'INTERPOLAZIONE!!! PER IL MOMENTO LASCIAMO COSÌ

U_final = interp2(X_father,Y_father,level{lvl-1,father_pos}.sol_after_timestep,...
                 X,Y,"cubic");
% Analogamente, per le derivate temporali
U0_der = interp2(X_father,Y_father,level{lvl-1,father_pos}.der_before_timestep,...
                 X,Y,"cubic");
U_final_der = interp2(X_father,Y_father,level{lvl-1,father_pos}.der_after_timestep,...
                 X,Y,"cubic");
%fprintf("I Valori PER Passo temporale = %d, Elemento %d \n", timestep,position)
%U0
%pause
%U_final
%pause
% ModU0 = real(U0).^2+imag(U0).^2;
% ModU_final = real(U_final).^2+imag(U_final).^2;
% ModDiff = ModU_final-ModU0;
% max(abs(ModDiff(:)));
%pause
%U0_der
%pause
%U_final_der
%pause                
% Generiamo ora i vettori per imporre le condizioni al bordo, tramite
% Interpolazione di Hermite
U0_border = sparse(U0);
U_final_border = sparse(U_final);                 
U0_border(2:m_y-1,2:m_x-1) = 0;
U_final_border(2:m_y-1,2:m_x-1) = 0;
U0_der_border = sparse(U0_der);
U_final_der_border = sparse(U_final_der);
U0_der_border(2:m_y-1,2:m_x-1) = 0;
U_final_der_border(2:m_y-1,2:m_x-1) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sappiamo che, considerati i tempi t0 (tempo iniziale) e t1 (tempo finale), 
% il polinomio di Hermite è H(t) = U(t0)+U'(t0)(t-t0)+U[t0,t0,t1](t-t0)^2+
%                                  U[t0,t0,t1,t1](t-t0)^2(t-t1)
% ====> La derivata sarà H'(t) = U'(t0)+2*U[t0,t0,t1](t-t0)+U[t0,t0,t1,t1](2*(t-t0)(t-t1)+(t-t0)^2)
%U[t0,t0,t1] = (U[t0,t1]-U[t0,t0])/(t1-t0), con U[t0,t0] = U'(t0) e U[t0,t1] = (U(t1)-U(t0))/(t1-t0)
k = tspan(timestep)-tspan(timestep-1);
U_t0_t1 = (U_final_border-U0_border)./k;
U_t0_t0_t1 = (U_t0_t1-U0_der_border)./k;
%U[t0,t0,t1,t1] = (U[t0,t1,t1]-U[t0,t0,t1])/(t1-t0)
% U[t0,t1,t1] = (U[t1,t1]-U[t0,t1])/(t1-t0), con U[t1,t1] = U'(t1)
U_t0_t1_t1 = (U_final_der_border-U_t0_t1)./k;
U_t0_t0_t1_t1 = (U_t0_t1_t1-U_t0_t0_t1)./k; 

%G(tspan(timestep-1)) == U0_border
%G(tspan(timestep)) == U_final_border
der_G =@(t) sparse(U0_der_border(:)+2*U_t0_t0_t1(:)*(t-tspan(timestep-1))+...
            U_t0_t0_t1_t1(:)*(2*(t-tspan(timestep-1))*(t-tspan(timestep))+(t-tspan(timestep-1)).^2));
% Per la derivata finale c'è un po' di errore
%der_G(tspan(timestep))-U_final_der_border
%der_G(tspan(timestep-1)) == U0_der_border  
             
%pause

% ORA POSSIAMO PASSARE AL PROBLEMA VERO E PROPRIO, UTILIZZEREMO UN METODO 
% RUNGE-KUTTA 2
% Il problema è sempre ∂tΨ = (i/2)∇^2Ψ+(i/2)(1-|Ψ|^2)Ψ
% Discretizzando con le differenze finite
h_y = (y_max-y_min)/(m_y-1);
h_x = (x_max-x_min)/(m_x-1);


D_y_kron = spdiags(ones(m_y,1)*[1,-2,1],[-1,0,1],m_y,m_y);
D_y_complete = kron(speye(m_x,m_x),D_y_kron);
%size(D_y_complete)
D_x_kron = spdiags(ones(m_x,1)*[1,-2,1],[-1,0,1],m_x,m_x);
D_x_complete = kron(D_x_kron,speye(m_y,m_y));
%size(D_x_complete)
% Cancelliamo il bordo sinistro
D_y_complete(2:m_y-1,:) = zeros(m_y-2,m_y*m_x);
D_x_complete(2:m_y-1,:) = zeros(m_y-2,m_y*m_x);
% Cancelliamo il bordo superiore
D_y_complete(1:m_y:m_y*(m_x-1)+1,:) = zeros(m_x,m_y*m_x);
D_x_complete(1:m_y:m_y*(m_x-1)+1,:) = zeros(m_x,m_y*m_x);
% Cancelliamo il bordo inferiore
D_y_complete(m_y:m_y:m_y*m_x,:) = zeros(m_x,m_y*m_x);
D_x_complete(m_y:m_y:m_y*m_x,:) = zeros(m_x,m_y*m_x);
% Cancelliamo il bordo destro
D_y_complete(m_y*(m_x-1)+2:m_y*m_x-1,:) = zeros(m_y-2,m_y*m_x);
D_x_complete(m_y*(m_x-1)+2:m_y*m_x-1,:) = zeros(m_y-2,m_y*m_x);

D_x_complete = (1/(h_x.^2))*D_x_complete;
D_y_complete = (1/(h_y.^2))*D_y_complete;

Lap = ((1i)/2)*(D_x_complete+D_y_complete);

DelMat = speye(m_y*m_x,m_y*m_x);
DelMat(2:m_y-1,:) = zeros(m_y-2,m_y*m_x);
DelMat(1:m_y:m_y*(m_x-1)+1,:) = zeros(m_x,m_y*m_x);
DelMat(m_y:m_y:m_y*m_x,:) = zeros(m_x,m_y*m_x);
DelMat(m_y*(m_x-1)+2:m_y*m_x-1,:) = zeros(m_y-2,m_y*m_x);

F =@(t,u) Lap*u+((1i)/2)*((1-(real(DelMat*u).^2+imag(DelMat*u).^2)).*(DelMat*u))+((der_G(t)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metodo RK2, dimezzando il passo temporale
passi = 100;

k_rk = (tspan(timestep)-tspan(timestep-1))/(passi-1);
tspan_rk = linspace(tspan(timestep-1),tspan(timestep),passi);
U_i = U0(:);

for t = 2:passi
psi1 = U_i;
psi2 = U_i+k_rk*F(tspan_rk(t-1),U_i);
U_i = U_i+(k_rk/2)*F(tspan_rk(t-1),psi1)+(k_rk/2)*F(tspan_rk(t),psi2);
end  
U_f = reshape(U_i,m_y,m_x);
% Ora salviamo il valore iniziale e finale nel sottolivello
level{lvl,position}.sol_before_timestep = U0;
level{lvl,position}.sol_after_timestep = U_f;
level{lvl,position}.Xmesh = X;
level{lvl,position}.Ymesh = Y;
level{lvl,position}.modules2 = (real(U_f).^2+imag(U_f).^2);
