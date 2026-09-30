% Matrix-vector multiplication for grid-defined functions
function f = mult(L,u)
  f = reshape(L*u(:),size(u));
end