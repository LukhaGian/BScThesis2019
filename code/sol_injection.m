function [level] = sol_injection(level,lvl,position)
% Funzione che, data la struttura level, il livello di lavoro e la posizione 
% del sottolivello di lavoro in level{lvl,position} inietta la soluzione calcolata
% in level{lvl,position} nella soluzione padre in 
% level{lvl-1,level{lvl,position}.father}
father_pos = level{lvl,position}.father;
F_sol = level{lvl-1,father_pos}.sol_after_timestep;
S_sol = level{lvl,position}.sol_after_timestep;
n_y = size(S_sol,1);
n_x = size(S_sol,2);
% Ritagliamo la matrice figlia, al fine di poter iniettare la soluzione
Temp = S_sol(1:2:n_y,:);
Temp2 = Temp(:,1:2:n_x);
% Dobbiamo determinare la posizione di S_sol nel dominio del padre
% Per le colonne è facile:
left_col = (level{lvl,position}.values(2)+1)/2;
right_col = (level{lvl,position}.values(4)+1)/2;
% Per le righe, dobbiamo ribaltare l'ordinamento
n_y_father = size(F_sol,1);
left_row = (level{lvl,position}.values(1)+1)/2;
right_row = (level{lvl,position}.values(3)+1)/2;
left_row = n_y_father-left_row+1;
right_row = n_y_father-right_row+1;
% Iniettiamo la nuova soluzione
F_sol(left_row:right_row,left_col:right_col) = Temp2;
U = F_sol;
level{lvl-1,father_pos}.sol_after_timestep = U;
%pause