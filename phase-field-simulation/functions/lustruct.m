% Calculates permute LU decomposition
function M = lustruct(M)
  [L,U,P,Q] = lu(M);
  M = struct('L',L,'U',U,'P',P,'Q',Q);
end