% Risolviamo l'equazione di Gross Pitaveskii tramite il metodo di
% Strang Splitting,
% ∂tΨ = (i/2)∇^2Ψ+(i/2)(1-|Ψ|^2)Ψ
% Condizioni al bordo di Neumann omogenee Ψ' = 0
% Condizioni iniziali U0
% Consideriamo un dominio quadrato
clear all
close all

m_x = 1000;
m_y = 1000;
m = m_x;
x_max = 10;
x_min = -10;
y_max = 10;
y_min = -10;
x = linspace(x_min,x_max,m_x);
y = linspace(y_min,y_max,m_y);
[A,B] = meshgrid(x,y);
A = flip(A);
B = flip(B);
% Discretizzazione spaziale analoga per entrambe le dimensioni
h = (x_max-x_min)/(m-1);

% N di livelli massimi
max_n_levels = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Condizione iniziale di un vortice stazionario centrato in un punto (x0,y0)
 initial.type = "sin"
x0 = 0; y0 = 0;
initial.x0 = x0;
initial.y0 = y0;
R2 = (A-x0).^2+(B-y0).^2;
rho = (((a(4)*R2+a(3)).*R2+a(2)).*R2+a(1)).*R2./...
       ((((a(4)*R2+b(3)).*R2+b(2)).*R2+b(1)).*R2+1);
psi = sqrt(rho).*exp(-1i*atan2(B-y0,A-x0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plottiamo |Ψ0|^2
InitMod = (real(psi).^2+imag(psi).^2);
U0 = psi;
set(gcf, 'WindowState', 'fullscreen');
%shading flat
Q = pcolor(A,B,InitMod)
set(Q, 'EdgeColor', 'none')
colormap parula
colorbar
axis square
saveas(gcf,'1000Sin.png')

tin = 0;
tstar = 120;
tm = 4000;
k = (tstar-tin)/(tm-1);
tspan = linspace(tin,tstar,tm);

D_x_kron = spdiags(ones(m,1)*[1,-2,1],[-1,0,1],m,m);
% Condizioni al bordo (Neumann omogenee)
D_x_kron(1,2) = 2;
D_x_kron(m,m-1) = 2;
D_x_kron = ((1i)/2)*(1/(h^2))*D_x_kron;
D_y_kron = D_x_kron;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo esplicito del Laplaciano secondo le differenze finite
LAP = kron(speye(m,m),D_y_kron)+kron(D_x_kron,speye(m,m));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcolo prima del ciclo degli esponenziali di matrice per il 2° passo di Strang Splitting
E_y = expm(k*D_y_kron);
E_x = transpose(expm(k*D_x_kron));
% Condizione iniziale
U = U0;

for p = 2:tm
  p
  tic
  % Salviamoci il valore della soluzione al passo precedente
  U_before = U;
  % Metodo di Strang Splitting
  U_onethird = exp((k/2)*(1i/2)*(1-(real(U).^2+imag(U).^2))).*U;
  U_twothird = E_y*U_onethird*E_x;
  U = exp((k/2)*((1i)/2)*(1-(real(U_twothird).^2+imag(U_twothird).^2))).*U_twothird;
  modules = sqrt(real(U).^2+imag(U).^2);
  % Flagging dei punti per cui |Ψ| <= 1/2 e dei punti limitrofi (secondo LUMR)
  flagged_modules = (modules <= 1/2);
  [row,col,flagged_Mat] = flagger2(flagged_modules);
  
  % Generiamo in partenza il livello 1, quello "completo" del problema
  s.values = [m_y,1;1,m_x];
  % Salviamoci i meshgrid del problema originale
  s.Xmesh = A;
  s.Ymesh = B;
  % Definiamo inoltre la posizione nel dominio del livello (le coordinate identificano
  % la posizione fisica nel dominio discretizzato)
  % Per livello 1
  s.domain_position = [y_max,y_min;x_min,x_max];
  % Per il livello 1 la soluzione nel timestep è quella calcolata mediante Strang-Splitting
  % (Successivamente verranno iniettati i valori nuovi dal livello 2)
  s.sol_before_timestep = U_before;
  s.sol_after_timestep = U;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Salviamo anche le derivate temporali prima del timestep e dopo del timestep
  % DOBBIAMO TRASFORMARE IN FORMA VETTORIALE ED UTILIZZARE IL LAPLACIANO COMPLETO
  s.der_before_timestep = reshape(LAP*U_before(:)+...
                          ((1i)/2)*(1-(real(U_before(:)).^2+imag(U_before(:)).^2))...
                          .*U_before(:),m,m);
  s.der_after_timestep = reshape(LAP*U(:)+...
                          ((1i)/2)*(1-(real(U(:)).^2+imag(U(:)).^2)).*U(:),m,m);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Salviamo il tutto nella struttura level
  level{1,1} = s;
  % Ora generiamo un array che definisca la prima posizione libera nel livello i-esimo
  first_free_in_lvl = ones(1,max_n_levels);
  % Il primo livello è stato generato, perciò aggiorniamo la prima componente
  first_free_in_lvl(1) = 2;
  
  % Per livello successivo
  % structure_refinement_lvl3(level,Matrice padre,Livello di lavoro,Posizione del padre in level{Livello-1,Posizione},first_free_in_lvl)
  [level,first_free_in_lvl] = structure_refinement_lvl3(level,flagged_Mat,2,1,first_free_in_lvl,row,col);

  % Ora, andiamo a risolvere il problema nei vari sottolivelli 2
  for no = 1:first_free_in_lvl(2)-1
    %       sublevel_solver(level,timestep,tspan,lvl,posizione in level{lvl,} dell'elemento di lavoro)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    level = sublevel_solver(level,p,tspan,initial,2,no);
    % Iniettiamo ora la soluzione del sottolivello nella soluzione padre
    %level = sol_injection(level,2,no);
    % PROTOTIPO
    level = sol_injection(level,2,no);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AGGIORNIAMO ORA LA DERIVATA AL TEMPO "FINALE"
    U = level{1,1}.sol_after_timestep;
    level{1,1}.der_after_timestep = reshape(LAP*U(:)+((1i)/2)*(1-(real(U(:)).^2+imag(U(:)).^2))...
                                    .*U(:),m,m);
  end

  U = level{1,1}.sol_after_timestep;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plotting della soluzione aggiornata
  modules2 = (real(U).^2+imag(U).^2);
  
  %axis equal
  if(p == 1000 | p == 2000 | p == 3000 | p == 4000)
      pause(2)
      set(gcf, 'WindowState', 'fullscreen');
      shading flat
      pcolor(A,B,modules2)
      %axis equal
      colormap parula
      %axis equal
      for no = 1:first_free_in_lvl(2)-1
         hold on
         %axis equal
         pcolor(level{2,no}.Xmesh,level{2,no}.Ymesh,level{2,no}.modules2)
         colormap parula
      end
      colorbar
      pause(2)
  end
  axis square
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(p == 1000)  
    saveas(gcf,'SecondTra.png')
  end
  if(p == 2000)
      saveas(gcf,'ThirdTra.png')
  end
  if(p == 3000)
      saveas(gcf,'FourthTra.png')
  end  
  if(p == 4000)  
    saveas(gcf,'FifthTra.png')
  end  
  %pause(0.1)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Elimininamo il contenuto di level e il resto al fine di ogni iterazione
  clear level;
  clear s;
  clear U_before;
  clear row;
  clear col;
  toc
end

