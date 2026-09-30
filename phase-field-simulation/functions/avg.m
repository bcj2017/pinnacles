% Face-cell/cell-face averaging operator
function uavg = avg(u)
  uavg = 0.5*(u(2:end,:) + u(1:end-1,:));
end