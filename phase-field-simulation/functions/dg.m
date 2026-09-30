% Phase-field potential
function val = dg(phi)
    val = 2*phi - 6*phi.^2 + 4*phi.^3;
end