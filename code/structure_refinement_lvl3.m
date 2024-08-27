function [level,first_free_in_lvl] = structure_refinement_lvl3(level,A,lvl,father_pos,first_free_in_lvl,row,col)
  % Funzione che, passata la struttura level, la matrice A, il livello di lavoro e la posizione del padre
  % in level{lvl-1,...} genera il sottolivello di level a partire da A (matrice di flagging)
  % first_free_in_lvl è l'array in cui, nell'i-esima componente, è salvato il primo 
  % posto libero in level al livello i-esimo
  
  [n_rows,n_col] = size(A);
  % NECESSARIE PER DETERMINARE LA POSIZIONE DEI PUNTI NEL DOMINIO
  xspan = linspace(level{lvl-1,father_pos}.domain_position(2),level{lvl-1,father_pos}.domain_position(4),2*n_col-1);
  yspan = linspace(level{lvl-1,father_pos}.domain_position(3),level{lvl-1,father_pos}.domain_position(1),2*n_rows-1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 0 = physical boundary
  % 1 = the elements are inside the grid
% Generate the rectangles saved in the l
for m = 1:min(row)-1
   l{n_rows-m+1} = [];
end 
for i = min(row):n_rows-1
  row_el = 0;
  inside_rect = false;
  l{n_rows-i+1} = [];
  for j = min(col):n_col
    if(inside_rect == false && A(i,j) == 1 && A(i+1,j) == 1)
      inside_rect = true;
      left_value = [n_rows-i+1;j];
    % Right side
    else
      if(inside_rect == true && A(i,j) == 1 && A(i+1,j) == 1)
        if(j == n_col)
          right_value = [n_rows-(i+1)+1;j];
          row_el = row_el+1;
          %a.values = [left_value(1),right_value(1);left_value(2),right_value(2)];
          a.father = father_pos;
          a.values = [2*left_value(1)-1,2*right_value(1)-1;2*left_value(2)-1,2*right_value(2)-1];
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %                       riga                                colonna
          a.domain_position = [yspan(a.values(1)),yspan(a.values(3));xspan(a.values(2)),xspan(a.values(4))];
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          l{n_rows-i+1}{row_el} = a;
          inside_rect = false;
          %%%%%%%%%
          clear a;
          %%%%%%%%%
        else
          if(A(i,j+1) ~= 1 || A(i+1,j+1) ~= 1)
            right_value = [n_rows-(i+1)+1;j];
            row_el = row_el+1;
            %a.values = [left_value(1),right_value(1);left_value(2),right_value(2)];
            a.father = father_pos;
            a.values = [2*left_value(1)-1,2*right_value(1)-1;2*left_value(2)-1,2*right_value(2)-1];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     riga                                  colonna
            a.domain_position = [yspan(a.values(1)),yspan(a.values(3));xspan(a.values(2)),xspan(a.values(4))];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            l{n_rows-i+1}{row_el} = a;
            clear a;
            inside_rect = false;
          end 
        end
      end   
    end
  end
end

% Merge rectangles in common structures (Non efficiente, da migliorare)
for i = 2:n_rows-1
  % n_rows-max(row)
  % Only contigouos rectangles are allowed
    for p = 1:size(l{i},2)
      j = 1;
      while(j <= size(l{i+1},2))
        al = l{i}{p}.values(2);
        b = l{i+1}{j}.values(2);
        c = l{i}{p}.values(4);
        d = l{i+1}{j}.values(4);
        if(al == b && c == d)
          % Merge the rectangles in the upper row
          l{i+1}{j}.values(3:4) = l{i}{p}.values(3:4);
          l{i+1}{j}.domain_position(3:4) = l{i}{p}.domain_position(3:4);
          % Erase the useless lower cell
          l{i}{p} = [];
          break;
        else
          j = j+1;
        end
      end
    end
end

% Now we should erase the empty rows (after the merge)
for i = 2:n_rows-1
  n_el = size(l{i},2);
  if(n_el ~= 0)
    all_empty = true;
    % Check if there's only empty elements
    for j = 1:n_el
      if(isempty(l{i}{j}) == 0)
        all_empty = false;
        break;
      end
    end
    if(all_empty == true)
      l{i} = [];
    end
  end
end 

% Rearrange the data structure in order to define the level and the number of that
% element in the level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%no_element = 1;
no_element = first_free_in_lvl(lvl);
for i = n_rows:-1:2
  % If the row is not empty
  if(isempty(l{i}) ~= 1)
    n_el_4_row = size(l{i},2);
    for j = 1:n_el_4_row
      if(isempty(l{i}{j}) ~= 1)
        level{lvl,no_element} = l{i}{j};
        no_element = no_element+1;
      end
    end
  end
end
% Aggiorniamo la prima posizione disponibile in level{lvl,...}
first_free_in_lvl(lvl) = no_element;