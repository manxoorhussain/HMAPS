function [x,y, Nb] = circleUniformCenters(N,R,plt)
        if ~exist('plt','var'), plt = false;  end
        if ~exist('R','var'), R = 1;  end
        x(1) = 0; y(1) = 0;
        Ns = round( (sqrt(pi+4*(N-1)) - sqrt(pi)) /(2*sqrt(pi)) );       % number of circles
        K = pi*(Ns+1)/(N-1);            % constant used in each loop

        for i = 1:Ns-1
            ri = i/Ns;
            ni = round(2*pi*ri/K);      % number of points on circle i
            t = linspace(0, 2*pi, ni+1)';   t = t(1:ni);

            if mod(i,2)==0   % stagger the start of every other circle
                dt = t(2) - t(1);
                t = t + 0.5*dt;
            end

            x = [x; ri*cos(t)];   y = [y; ri*sin(t)];  
        end

        Nb = N - length(x);   % remaining points to be placed on the outter circle
        t = linspace(0, 2*pi, Nb + 1)';  t = t(1:Nb);
        x = [cos(t); x];  y = [sin(t); y];
        x = R*x;  y = R*y;                 % adjust to have radius R

        if plt, scatter( x, y,'b.'); end
    end