function [row,col,flagged_Mat] = flagger2(M)
% Funzione che flagga i vicini dei punti flaggati (secondo LUMR)  
[m,n] = size(M);
flagged_Mat = sparse(zeros(m,n));
[row,col] = find(M);
no_flagged = size(row,1);
for(i = 1:no_flagged)
  q = row(i);
  p = col(i);
  %
  left_bound = 2;
  if(p == 1)
    left_bound = 0;
  elseif(p == 2)
    left_bound = 1;
  end
  right_bound = 2;
  if(p == n)
    right_bound = 0;
  elseif(p == n-1)
    right_bound = 1;
  end
    upper_bound = 2;
  if(q == 1)
    upper_bound = 0;
  elseif(q == 2)
    upper_bound = 1;
  end
  lower_bound = 2;
  if(q == m)
    lower_bound = 0;
  elseif(q == m-1)
    lower_bound = 1;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  row(i) = row(i)-upper_bound;
  col(i) = col(i)-left_bound;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  flagged_Mat(q-upper_bound:q+lower_bound,p-left_bound:p+right_bound) = 1; 
end  
   