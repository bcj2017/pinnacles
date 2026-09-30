
% Construct second order time-stepping matrices
function [Lu,Lv,Ltm,Ltp,Lphim,Lphip,Lumh,Luph,Lvmh,Lvph,Lp,Lq,K] = ...
  calcops(Nr,Nz,dt,Pr,m)

    global hz rc hp hm %#ok<*GVMIS> 
    
    r = 0.5*(rc(2:end) + rc(1:end-1));

    % Poisson matrix for pressure
    D = spdiags([hm./hp (hp./hm - hm./hp) -hp./hm]./(hp+hm)./rc,-1:1,Nr,Nr)';
    D2 = spdiags( 2*[hm, -(hp + hm), hp]./(hp.*hm.*(hp + hm)),-1:1,Nr,Nr)';

    Lrr = D2 + D;
    Lrr(end,end) = Lrr(end,end)/2; %correct Neumann condition at r = 0
    Lrr(end,end) = Lrr(end,end) + 1;
    
    Lzz = Kxx(Nz,hz,1,1);
    Lp = -(kron(speye(Nz),Lrr) + kron(Lzz,speye(Nr)));

    % Horizontal velocity
    hx = r(2:end) - r(1:end-1);
    hx = [hx(1);hx;hx(end)];
    hxk = hx(2:end);
    hxkm1 = hx(1:end-1);

    D2 = spdiags( 2*[hxkm1, -(hxk + hxkm1), hxk]./(hxk.*hxkm1.*(hxk + hxkm1)),-1:1,Nr-1,Nr-1)';
    D = spdiags([hxkm1./hxk (hxk./hxkm1 - hxkm1./hxk) -hxk./hxkm1]./(hxk+hxkm1)./r,-1:1,Nr-1,Nr-1)';
    D([1 end],:) = 0;
    Lrr = D2 + D + K0(Nr-1)./(r.^2);
    Lzz = Kxx(Nz,hz,3,3);

    Lu = (kron(speye(Nz),Lrr) + kron(Lzz,speye(Nr-1))); 
    Lumh = speye((Nr-1)*Nz) - 0.5*dt*Pr*Lu; 
    Luph = speye((Nr-1)*Nz) + 0.5*dt*Pr*Lu; 
    
    % vertical velocity
    D2 = spdiags( 2*[hm, -(hp + hm), hp]./(hp.*hm.*(hp + hm)),-1:1,Nr,Nr)';
    D = spdiags([hm./hp (hp./hm - hm./hp) -hp./hm]./(hp+hm)./rc,-1:1,Nr,Nr)';
    D([1 end],:) = 0;
    Lrr = D2 + D;
    Lrr(1,1) = Lrr(1,1)/2;
    Lrr(end,end) = 3*Lrr(end,end)/2;
    Lzz = Kxx(Nz-1,hz,2,2);

    Lv = kron(speye(Nz-1),Lrr) + kron(Lzz,speye(Nr)); 
    Lvmh = speye(Nr*(Nz-1)) - 0.5*dt*Pr*Lv; 
    Lvph = speye(Nr*(Nz-1)) + 0.5*dt*Pr*Lv; 
    
    % Temperature and phase field
    D = spdiags([hm./hp,(hp./hm - hm./hp),-hp./hm]./(hp+hm)./rc,-1:1,Nr,Nr)';
    D2 = spdiags( 2*[hm, -(hp + hm), hp]./(hp.*hm.*(hp + hm)),-1:1,Nr,Nr)';
    D([1 end],:) = 0;
    
    Lrr = D2 + D;
    
    Lrr(1,2) = 2*Lrr(1,2);
    Lrr(end,end) = Lrr(end,end)/2; %correct Neumann condition at r = 0
    Lzz = Kxx(Nz,hz,1,1);
    LL = (kron(speye(Nz),Lrr) + kron(Lzz,speye(Nr)));

    % Temperature and phase field    
    Ltm = speye(Nr*Nz) - 0.5*dt*LL;
    Ltp = speye(Nr*Nz) + 0.5*dt*LL;
    Lphim = speye(Nr*Nz) - 0.5*dt*m*LL;
    Lphip = speye(Nr*Nz) + 0.5*dt*m*LL;
    Lzz = Kxx(Nz,hz,1,1);
    Lq = kron(speye(Nz),Lrr) + kron(Lzz,speye(Nr));
    
end

% Identity without end points
function A = K0(n)
    A = spdiags([0;ones(n-2,1);0],0,n,n)';    
end
% Second derivative matrix
function A = Kxx(n,h,a11,ann)
    % a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
    A = spdiags([1 -a11 0;ones(n-2,1)*[1 -2 1];0 -ann 1],-1:1,n,n)'/h^2;
end