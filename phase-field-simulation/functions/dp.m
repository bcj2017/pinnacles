% Phase-field potential
function val = dp(phi)
    val = 30*phi.^2 - 60*phi.^3 + 30*phi.^4;
end
