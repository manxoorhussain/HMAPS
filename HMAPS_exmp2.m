%------ HRBF based MAPS for Elliptic PDES ----------%

% Example 2 %

warning off;
clearvars;      close all;      clc;    format long;
tic;

%%---Data Points--%%
w0 = 9e-3;
c0 = 0.3; % MHRBF, MGA, GA, HRBF, c0 = 0.3
mm = 2;

%-----------------------%
nodes        = 'n16';       % n16 = 1600, n9 = 900, n4 = 400, n1 = 100
method       = 'HMAPS';     % HMAPS = MAPS-GA+C3, HRBF = GA+C3, PMAPS = MAPS-PS, GMAPS = MAPS-GA, PSRBF = PS-RBF, GARBF = GA-RBF
BCS          = 'DtN';       % dirichlet-to-NeuMann BCs "other cases can be added" 
%------------ Domain parameters ---------- %
a = 0; b = 1;
switch (nodes)
    case('n16')
        [x,y]  = meshgrid(linspace(a,b,40),linspace(a,b,40));                  % n=1600
        Xn     = [x(:), y(:)];
        Dbpi   = find(Xn(:,1)==a | Xn(:,1)==b);                                % index of Dirichlet bpts
        Nbpi   = find((Xn(:,2)==a | Xn(:,2)==b) & ~(Xn(:,1)==a | Xn(:,1)==b)); % index of Neumann bpts
        NI     = find((Xn(:,1)~=a & Xn(:,1)~=b & Xn(:,2)~=a & Xn(:,2)~=b));    % index of Interior pts
    case('n9')
        [x,y]  = meshgrid(linspace(a,b,30),linspace(a,b,30));                  % n=900
        Xn     = [x(:), y(:)];
        Dbpi   = find(Xn(:,1)==a | Xn(:,1)==b);                                % index of Dirichlet bpts
        Nbpi   = find((Xn(:,2)==a | Xn(:,2)==b) & ~(Xn(:,1)==a | Xn(:,1)==b)); % index of Neumann bpts
        NI     = find((Xn(:,1)~=a & Xn(:,1)~=b & Xn(:,2)~=a & Xn(:,2)~=b));    % index of Interior pts
    case('n4')
        [x,y]  = meshgrid(linspace(a,b,20),linspace(a,b,20));                  % n=400
        Xn     = [x(:), y(:)];
        Dbpi   = find(Xn(:,1)==a | Xn(:,1)==b);                                % index of Dirichlet bpts
        Nbpi   = find((Xn(:,2)==a | Xn(:,2)==b) & ~(Xn(:,1)==a | Xn(:,1)==b)); % index of Neumann bpts
        NI     = find((Xn(:,1)~=a & Xn(:,1)~=b & Xn(:,2)~=a & Xn(:,2)~=b));    % index of Interior pts
    case('n1')
        [x,y]  = meshgrid(linspace(a,b,10),linspace(a,b,10));                  % n=200
        Xn     = [x(:), y(:)];
        Dbpi   = find(Xn(:,1)==a | Xn(:,1)==b);                                % index of Dirichlet bpts
        Nbpi   = find((Xn(:,2)==a | Xn(:,2)==b) & ~(Xn(:,1)==a | Xn(:,1)==b)); % index of Neumann bpts
        NI     = find((Xn(:,1)~=a & Xn(:,1)~=b & Xn(:,2)~=a & Xn(:,2)~=b));    % index of Interior pts
end
NN = length(Xn(:,1));
%----------------------%
[Xa,Yb] = meshgrid(linspace(a,b,32),linspace(a,b,32));
XX = [Xa(:), Yb(:)];
N = length(XX(:,1)); % number of test nodes (1024)
%%****** shape & wieght parameter *******%%
c  = (c0/sqrt(NN))*ones(NN,1)';   % constant shape paramter
wp = (w0/sqrt(NN))*ones(NN,1)'; % constant wiegth parameter
%---- functions BCS and source -----%
switch (BCS)
    case ('DtN')
        %%--soruce terms and boundary function--%%
        F=@(X,Y)   (130*(2*X - 2/5).^2)./((X - 1/5).^2 + (Y + 1/10).^2 + 65).^3 ...
                   + (130*(2*Y + 1/5).^2)./((X - 1/5).^2 + (Y + 1/10).^2 + 65).^3 ...
                   - 260./((X - 1/5).^2 + (Y + 1/10).^2 + 65).^2;                % source
        G1=@(X,Y)  65./((X - 1/5).^2 + (Y + 1/10).^2 + 65);                      % Dirichlet boundary data 
        G2=@(X,Y)  -(65.*(2.*Y + 1/5))./((X - 1/5).^2 + (Y + 1/10).^2 + 65).^2;  % Neumann boundary data  
        Uex=@(X,Y) 65./((X - 1/5).^2 + (Y + 1/10).^2 + 65);                      % exact
end 

%---- Square domain ----%
% for Halton nodes the next two-lines be un-commented (ref. Fig. 5) %
% XY0 = haltonset(2,'Skip',1e3,'Leap',1e2); XY = net(XY0,NN);
% Xn(NI,:) = XY(NI,:);
%%%%%%%%%%%%%%%%%%%%%%%%
xx=Xn(:,1); yy=Xn(:,2);
figure;
plot(xx(NI),yy(NI),'b.','linewidth',1.5); grid off; axis equal;
hold on
plot(xx(Dbpi),yy(Dbpi),'r.','linewidth',1.5); grid off; axis square;
hold on
plot(xx(Nbpi),yy(Nbpi),'k.','linewidth',1.5); grid off; axis square;
hold off
axis equal

%--------------------------------%
switch (method)
    case ('HMAPS')
        %% GA+C3-MAPS %%
        d0phai =@(X,c)      (-psi(1) + expint((c.*(X+eps)).^2) + log((c.*(X+eps)).^2))./(4*c.^2) +  wp.*(((X+eps).^(2*mm+1))./(2*mm+1).^2);
        d1phai =@(X,Y,c)    (c.^2.*(2.*Y).*sqrt((c.*X).^2 + 1) + (c.^2.*(2.*Y).*((c.*X).^2 + 4))./(2*sqrt((c.*X).^2 + 1)) ...
                            - (3*c.^2.*(2.*Y))./(2*(sqrt((c.*X).^2 + 1) + 1).*sqrt((c.*X).^2 + 1)))./(9*c.^2) ...
                            + wp.*(((2.*Y).*X.^(2*mm))./(2.*(X+eps)));
        d2phai =@(X,c)      exp(-(c.*X).^2) + wp.*(X+eps).^(2*mm-1);
    case ('HRBF')
        %% GA+C3 %%
        d0phai =@(X,c)      wp.*((X+eps).^(2*mm-1)) + exp(-(c.*X).^2);
        d1phai =@(X,Y,c)    Y.*(((((X+eps).^(2*mm - 2).*(2.*mm - 1))./(X+eps)).*wp - 2.*c.^2.*X.*exp(-c.^2.*X.^2))./(X+eps));
        d2phai =@(X,c)      wp.*(((X+eps).^(2.*mm - 2).*(2.*mm - 1))./(X+eps) + (X+eps).^(2.*mm - 3).*(2.*mm - 1).*(2.*mm - 2))...
                            + (- 2.*c.^2.*X.*exp(-c.^2.*X.^2))./(X+eps) - 2.*c.^2.*exp(-c.^2.*X.^2) + 4.*c.^4.*X.^2.*exp(-c.^2.*X.^2);
    case ('PMAPS')
        %% SpN-MAPS %%
        d0phai  =@(X,c)     ((X+eps).^(2*mm+1))./(2*mm+1).^2;
        d1phai  =@(X,Y,c)   ((2.*Y).*X.^(2*mm))./(2.*(X+eps));
        d2phai  =@(X,c)     (X+eps).^(2*mm-1);
    case ('GMAPS')
        %% GA-MAPS %%
        d0phai  =@(X,c)     (-psi(1) + expint((c.*(X+eps)).^2) + log((c.*(X+eps)).^2))./(4*c.^2);
        d1phai  =@(X,Y,c)   (c.^2.*(2.*Y).*sqrt((c.*X).^2 + 1) + (c.^2.*(2.*Y).*((c.*X).^2 + 4))./(2*sqrt((c.*X).^2 + 1)) ...
                            - (3*c.^2.*(2.*Y))./(2*(sqrt((c.*X).^2 + 1) + 1).*sqrt((c.*X).^2 + 1)))./(9*c.^2);
        d2phai  =@(X,c)     exp(-(c.*X).^2);
    case ('PSRBF')
        %% SpN-RBF %%
        d0phai  =@(X,c)     (X+eps).^(2*mm-1);
        d1phai  =@(X,Y,c)   Y.*(((X+eps).^(2*mm - 2).*(2.*mm - 1))./(X+eps));
        d2phai  =@(X,c)     ((X+eps).^(2.*mm - 2).*(2.*mm - 1))./(X+eps) + (X+eps).^(2.*mm - 3).*(2.*mm - 1).*(2.*mm - 2);
    case ('GARBF')
        %% GA-RBF %%
        d0phai  =@(X,c)     exp(-(c.*X).^2);
        d1phai  =@(X,Y,c)   -2.*Y.*c.^2.*exp(-c.^2.*X.^2);
        d2phai  =@(X,c)     4.*c.^4.*X.^2.*exp(-c.^2.*X.^2) - 4.*c.^2.*exp(-c.^2.*X.^2);
end

%%------------System Matrices Assembly ----------------------------------%%
b=zeros(NN,1); 
rx = xx - xx';      ry = yy - yy';      r = sqrt(rx.^2 + ry.^2);
Rx = XX(:,1) - xx'; Ry = XX(:,2) - yy'; R = sqrt(Rx.^2 + Ry.^2);

A0 = d0phai(r,c); 
A1 = d1phai(r,ry,c); 
A2 = d2phai(r,c); 
H  = d0phai(R,c);
Bi = F(xx(NI),yy(NI));      % interior points conditions
Bd = G1(xx(Dbpi),yy(Dbpi)); % boundary conditions
Bn = G2(xx(Nbpi),yy(Nbpi)); % boundary conditions
bb = [Bi;Bd;Bn];            % RHS vector

%------ LSE solution -----%
Ai=A2(NI,:);
Ad=A0(Dbpi,:);
An=A1(Nbpi,:);
A=[Ai;Ad;An];
AA= A + 2e-15*eye(size(A)); % MDI-scaling
kA = cond(AA);
if kA > 2e+16
    lambda = pinv(AA)*bb;  % can be replaced with "mldivide(AA,bb)"
else
    lambda = AA\bb;
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
subplot(1,2,1)
trisurf(td,XX(:,1),XX(:,2),uD); grid off, xlabel x, ylabel y, zlabel uA, view([50 50]) % axis equal %  
subplot(1,2,2)
trisurf(td,XX(:,1),XX(:,2),abs(uD-uex)); grid off, xlabel x, ylabel y, zlabel |error|, view([50 50]) %axis equal % 


%----------------- ThE EnD --------%
