%------ HRBF based MAPS for Elliptic PDES ----------%

% Example 4 %

warning off;
clearvars;      close all;      clc;    format long;
tic;
%%---Data Points--%%
w0 = 9e-3;
c0 = 50; % MHRBF, MGA, GA, HRBF, c0 = 0.3
mm = 3;

lam = 1; 
lam2=lam^2;
%-----------------------%
nodes        = 'n4';       % n8 = 800, n4 = 400, n2 = 200
method       = 'MHRBF';     % MHRBF = MAPS GA+C3, HRBF = GA+C3, MPS = MAPS PS, MGA = MAPS GA, PS = PS RBF, GA = GA RBF
cases        = 'case1';     % case1 = dirichlet BCs (can add other type of BCS cases)
%------------ Domain parameters ---------- %
switch (nodes)
    case('n2')
        [interior_nodes, boundary_nodes] = generateCircleNodes([0,0], 1, 150, 50); 
        x=[interior_nodes(:,1);boundary_nodes(:,1)];
        y=[interior_nodes(:,2);boundary_nodes(:,2)];
        Xn     = [x(:), y(:)]; xx = Xn(:,1); yy = Xn(:,2); 
        Dbpi   = boundary_nodes;    % index of Dirichlet bpts
        NI     = interior_nodes;    % index of Interior pts
    case('n4')
        [interior_nodes, boundary_nodes] = generateCircleNodes([0,0], 1, 300, 100); 
        x=[interior_nodes(:,1);boundary_nodes(:,1)];
        y=[interior_nodes(:,2);boundary_nodes(:,2)];
        Xn     = [x(:), y(:)]; xx = Xn(:,1); yy = Xn(:,2);
        Dbpi   = boundary_nodes;    % index of Dirichlet bpts
        NI     = interior_nodes;    % index of Interior pts
    case('n8')
        [interior_nodes, boundary_nodes] = generateCircleNodes([0,0], 1, 600, 200); 
        x=[interior_nodes(:,1);boundary_nodes(:,1)];
        y=[interior_nodes(:,2);boundary_nodes(:,2)];
        Xn     = [x(:), y(:)]; xx = Xn(:,1); yy = Xn(:,2);
        Dbpi   = boundary_nodes;    % index of Dirichlet bpts
        NI     = interior_nodes;    % index of Interior pts
end
NN = length(Xn(:,1));
%--------- Test nodes -------------%
[Xa,Yb, ~] = circleUniformCenters(200,1,1);
XX = [Xa(:), Yb(:)];
N = length(XX(:,1)); % number of test nodes (1024)
%%****** shape parameter *******%%
c  = (c0/sqrt(NN))*ones(NN,1)';   % constant shape paramter
wp = (w0/sqrt(NN))*ones(NN,1)'; % constant wiegth parameter
        
%%****** HRBF *******%%
switch (cases)
    case ('case1')
        %%--soruce terms and boundary function--%%
        F=@(X,Y)    lam2*(cos(pi*X).*sinh(Y) - sin(pi*X).*cosh(Y)) ...
                    - cos(pi*X).*sinh(Y) + sin(pi*X).*cosh(Y) ...
                    + pi^2*cos(pi*X).*sinh(Y) - pi^2*sin(pi*X).*cosh(Y);      % source
        G=@(X,Y)   sin(pi*X).*cosh(Y) - cos(pi*X).*sinh(Y);                   % boundary 
        Uex=@(X,Y)  sin(pi*X).*cosh(Y) - cos(pi*X).*sinh(Y);                  % exact
end 

%--------------------------------%
switch (method)
    case ('MHRBF')
        %% GA+TPS-MAPS %%
        d0phai =@(X,c)      (-psi(1) + expint((c.*(X)).^2) + log((c.*(X+eps)).^2))./(4*c.^2) ...
                            +  wp.*((X.^(2*mm+2).*((mm+1).*log(abs(X)+eps)-1))./4/(mm+1).^3);
        d1phai =@(X,Y,c)    (Y./(X+eps)).*((X.^(2*mm + 2).*wp)./(4*(X+eps).*(mm + 1)^2) ...
                            - ((2*exp(-X.^2.*c.^2))./(X+eps) - 2./(X+eps))./(4*c.^2) ...
                            + (X.^(2*mm + 1).*wp.*(2*mm + 2).*(log(abs(X)+eps).*(mm + 1) - 1))./(4*(mm + 1).^3));
        d2phai =@(X,c)      exp(-(c.*X).^2) + wp.*(X+eps).^(2*mm).*log(abs(X)+eps);
    case ('HRBF')
        %% GA+C3 %%
        d0phai =@(X,c)      wp.*((X+eps).^(2*mm-1)) + exp(-(c.*X).^2);
        d1phai =@(X,Y,c)    Y.*(((((X+eps).^(2*mm - 2).*(2.*mm - 1))./(X+eps)).*wp - 2.*c.^2.*X.*exp(-c.^2.*X.^2))./(X+eps));
        d2phai =@(X,c)      wp.*(((X+eps).^(2.*mm - 2).*(2.*mm - 1))./(X+eps) + (X+eps).^(2.*mm - 3).*(2.*mm - 1).*(2.*mm - 2))...
                            + (- 2.*c.^2.*X.*exp(-c.^2.*X.^2))./(X+eps) - 2.*c.^2.*exp(-c.^2.*X.^2) + 4.*c.^4.*X.^2.*exp(-c.^2.*X.^2);
    case ('MPS')
        %% SpN-MAPS %%
        d0phai  =@(X,c)     ((X+eps).^(2*mm+1))./(2*mm+1).^2;
        d1phai  =@(X,Y,c)   ((2.*Y).*X.^(2*mm))./(2.*(X+eps));
        d2phai  =@(X,c)     (X+eps).^(2*mm-1);
    case ('MGA')
        %% GA-MAPS %%
        d0phai  =@(X,c)     (-psi(1) + expint((c.*(X+eps)).^2) + log((c.*(X+eps)).^2))./(4*c.^2);
        d1phai  =@(X,Y,c)   (c.^2.*(2.*Y).*sqrt((c.*X).^2 + 1) + (c.^2.*(2.*Y).*((c.*X).^2 + 4))./(2*sqrt((c.*X).^2 + 1)) ...
                            - (3*c.^2.*(2.*Y))./(2*(sqrt((c.*X).^2 + 1) + 1).*sqrt((c.*X).^2 + 1)))./(9*c.^2);
        d2phai  =@(X,c)     exp(-(c.*X).^2);
    case ('PS')
        %% SpN-RBF %%
        d0phai  =@(X,c)     (X+eps).^(2*mm-1);
        d1phai  =@(X,Y,c)   Y.*(((X+eps).^(2*mm - 2).*(2.*mm - 1))./(X+eps));
        d2phai  =@(X,c)     ((X+eps).^(2.*mm - 2).*(2.*mm - 1))./(X+eps) + (X+eps).^(2.*mm - 3).*(2.*mm - 1).*(2.*mm - 2);
    case ('GA')
        %% GA-RBF %%
        d0phai  =@(X,c)     exp(-(c.*X).^2);
        d1phai  =@(X,Y,c)   -2.*Y.*c.^2.*exp(-c.^2.*X.^2);
        d2phai  =@(X,c)     4.*c.^4.*X.^2.*exp(-c.^2.*X.^2) - 4.*c.^2.*exp(-c.^2.*X.^2);
end

%%------------System Matrices----------------------------------%%

rx = xx - xx';      ry = yy - yy';      r = sqrt(rx.^2 + ry.^2);
Rx = XX(:,1) - xx'; Ry = XX(:,2) - yy'; R = sqrt(Rx.^2 + Ry.^2);

A0 = d0phai(r,c); A0(isnan(A0) | isinf(A0)) = 0;
A1x = d1phai(r,rx,c); A1y = d1phai(r,ry,c); 
A2 = d2phai(r,c); 
H  = d0phai(R,c); H(isnan(H) | isinf(H)) = 0;

bb = zeros(NN,1);
bb(1:length(NI))    = F(NI(:,1),NI(:,2));
bb(length(NI)+1:NN) = G(Dbpi(:,1),Dbpi(:,2));

%% LSE solution %%
Ai=A2(1:length(NI),:) - lam2*A0(1:length(NI),:);
Ad=A0(length(NI)+1:NN,:);

A=[Ai;Ad];
AA= A + 2e-15*eye(NN); % [A, ones(1,NN)';  ones(1,NN), 0]
kA = cond(AA);
if kA > 2e+16
    lambda = pinv(AA)*bb;  % lsqminnorm(AA,bb); %
else
    lambda = mldivide(AA,bb);
end
uD=H*lambda; 
runT=toc;
uex=Uex(XX(:,1),XX(:,2));
format short e
PMAE = norm(uD-uex,inf);
RMSE = norm(uD-uex,2)/sqrt(NN);
[PMAE, RMSE, kA, runT]

figure;
td = delaunay(XX(:,1),XX(:,2));
trisurf(td,XX(:,1),XX(:,2),uD); grid off, xlabel x, ylabel y, zlabel wA, view([50 50]) % axis equal %  

%----------------- ThE EnD --------%
