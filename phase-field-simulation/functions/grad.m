% Face-cell/cell-face gradient operator
function [ur,uz] = grad(u)

    global hz hp hc Nr %#ok<*GVMIS> 
    
    shp = size(u,1);
    du = diff(u,[],1);
    uz = diff(u,[],2)/hz; 

    if shp == Nr-1; ur = du./avg(hc);
    elseif shp == Nr; ur = du./hc;
    elseif shp == Nr+1; ur = du./hp; 
    elseif shp == Nr+2; ur = du./[hp(1);hp];
    end
    
end