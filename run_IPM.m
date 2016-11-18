function [z,finalx] = run_IPM( H, c, A, b )
%IPM Solves a quadratic programming problem in the form
%of 0.5*H*x'*x + c'*x such that,
%Ax <= b and x >= 0
%IPM uses the prediction and central path method to find
%its way to the most optimal solution. It runs until a 200 iteration
%threshhold or a value threshhold

    [m,n] = size(A);
    x = ones(n,1);
    lamda = ones(m,1);
    s = ones(m,1);
    
    %Calculate residual values and mu
    rd0 = H*x + c + A'*lamda;
    rp0 = s + A*x - b;
    mu = sum(lamda.*s)/m;
    e = ones(m,1);
    rc = diag(s) * diag(lamda)*e;
    i = 0;
    eps = 1e-16;
    maxiter = 200;
    %Begin Algorithm
    while i <= maxiter && norm(rd0) >= eps && norm(rp0) >= eps && abs(mu) >= eps
        %Solve for xaff, lamdaaff, and saff to center the point
        big = [H A' zeros(n,m); A zeros(m,m) eye(m); zeros(m,n) diag(s) diag(lamda)];
        r = [-rd0; -rp0; -rc];
        aff = big\r;

        xaff = aff(1:n);
        lamdaaff = aff(n+1:n+m);
        saff = aff(n+m+1:n+m+m);
        
        %Compute alphaaffine
        alphaaff = 1;
        %Record indexes where lamdaaff < 0
        indexfindz = find ( lamdaaff <0) ;
        if (isempty ( indexfindz )==0)
            alphaaff = min (alphaaff , min(-lamda (indexfindz ) ./ lamdaaff ( indexfindz ) ) ) ;
            %steplength=steplength+(length(lamda(idxz)));
        end
         %Record indexes where saff < 0
        indexfinds = find ( saff  <0) ;
        if(isempty ( indexfinds )==0)
            alphaaff = min ( alphaaff , min(-s(indexfinds ) ./ saff(indexfinds) ) ) ;
            %steplength=steplength+(length(s(idxs)));
        end

        muaff = (s + alphaaff .* saff)'*(lamda + alphaaff .* lamdaaff)/m;
       % steplength=steplength+(length(saff));
        cent = (muaff/mu)^3;

        %Solve for the step direction (deltax, deltalamda, deltas)
        imprc = rc + saff.*lamdaaff - cent*mu*e;
        deltars = [-rd0; -rp0; -imprc];
        deltas = big\deltars;
        %steplength=2*steplength+(length(big))^3;

        deltax = deltas(1:n);
        deltalamda = deltas(n+1:n+m);
        deltas = deltas(n+m+1:n+m+m);
        
        %Computer stepsize
        stepsizealpha = 1 ;
        %Record all indexes where deltalamda < 0
        indexfindz = find( deltalamda<0) ;
        if(isempty ( indexfindz)==0)
            stepsizealpha = min ( stepsizealpha , min(-lamda ( indexfindz ) ./ deltalamda(indexfindz ) ) ) ;
           % steplength=steplength+(length(lamda(idxz)));
        end
        
         %Record all indexes where deltas < 0
        indexfinds = find(deltas<0) ;
        if ( isempty( indexfinds )==0)
            stepsizealpha = min ( stepsizealpha , min(-s ( indexfinds ) ./ deltas( indexfinds ) ) ) ;
            %steplength=steplength+(length(s(idxs)));
        end

        %Compute new points
        x = x + stepsizealpha*deltax;
        lamda = lamda + stepsizealpha*deltalamda;
        s = s + stepsizealpha*deltas;
        
        %Update residual values and mu
        rd0 = H*x + c + A'*lamda;
        rp0 = s + A*x -b;        
        rc = diag(s) * diag(lamda)*e; 
        mu = sum(s.*lamda)/m;
        i = i + 1;
    end
    finalx = x;
    z=((0.5*finalx')*(H*finalx))+(c'*finalx);%[xx, fval] = quadprog(H,c,A,b);
end

