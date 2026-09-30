% Permutted LU solver for grid-defined functions
function x = lusolve(M,b)
  x = reshape(M.Q*(M.U\(M.L\(M.P*b(:)))),size(b));
end
