function [W] = NucMT(mu,muu,nu,rsex,r,pl)
    % single nuclear locus

    M=2; % M is the number of discrete mitochondria per cell
    n=2;  % n is the number of possible nuclear alleles. n=3 here
          % P second index j = 1:a, 2:A
    x=2;  % x is the exponent in the fitness function. x>1 makes it concave which we want.
    s=1;  % selection strength
    pl=round(pl,3);
    
    % Initialize population P with random distributions
    % Females
    P1=zeros(M+1,n);
    p = rand([M+1 1]);
    P1(:,1)=1/2*p./sum(p); %a
    P1(:,2)=1/2*p./sum(p); %A    
    % Males
    P2=zeros(M+1,n);
    p = rand([M+1 1]);
    P2(:,1)=1/2*p./sum(p); %a
    P2(:,2)=1/2*p./sum(p); %A

    % Life cycle operators
    % Mitochondrial mutation
    U=zeros(M+1,M+1);
    for i=0:M
        for j=0:M
            k=0:M;
            U(i+1,j+1)=sum(binopdf(k,M-j,mu).*binopdf(k+j-i,j,muu));
        end
    end
    % Nuclear mutation 
    NU(1,1)=1-nu;
    NU(1,2)=nu;
    NU(2,1)=nu;
    NU(2,2)=1-nu;
    % Recombination rates between different types of zygotes
    % Re-index the zygote types. Gametes of type i (F) and j (M) form a zygote of
    % type ZI(i,j):
    ZI = [[1,2];[3,4]];
    R = diag(ones(4,1));
    R(ZI(1,2),ZI(1,2)) = 1-r;
    R(ZI(2,1),ZI(1,2)) = r;
    R(ZI(1,2),ZI(2,1)) = r;
    R(ZI(2,1),ZI(2,1)) = 1-r;

    w=zeros(M+1,n);
    for i=0:M
        % just count the number of mismatches
        % w = 1 - (total mismatches/2M)^x;
        w(i+1,1)=1 - s*( i / M     )^x; %a
        w(i+1,2)=1 - s*( (M-i) / M )^x; %A
    end
    
    % Resampling at Meiosis 1
    L1=zeros(2*M+1,2*M+1);
    for i=0:2*M
        for j=0:2*M
            L1(i+1,j+1)=hygepdf(i,4*M,2*j,2*M);
        end
    end
    
    % Resampling at Meiosis 2
    L2=zeros(M+1,2*M+1);
    for i=0:M
        for j=0:2*M
            L2(i+1,j+1)=hygepdf(i,2*M,j,M);
        end
    end
    
    % Mitochondrial resampling after cell fusion
    % With paternal leakage pl, the zygote has M+pl*M total mitochondria
    % PSI resamples that up to 2M, a proper number for the zygote
    PSI=zeros(0,1);
    for i=0:2*M
        for j=0:round(pl*M)+M
            PSI(i+1,j+1)=binopdf(i,2*M,j/(round(pl*M)+M));
        end
    end
    
    % Mitochondrial resampling at cell fusion for the paternal cell
    % With paternal leakage it should only contribute pl*M mitochondria
    % PSIP reamples the number of mito mutants without replacement
    PSIP=zeros(0,1);
    for i=0:pl*M
        for j=0:M
            PSIP(i+1,j+1)=hygepdf(i,M,j,round(pl*M));
        end
    end
    
    % Mitochondrial resampling for asexual doubling and clonal division
    A=zeros(0,1);
    for i=0:M
        for j=0:M
            A(i+1,j+1)=hygepdf(i,2*M,2*j,M);
        end
    end
    
    AD=zeros(0,1);
    for i=0:2*M
        for j=0:M
            AD(i+1,j+1)=binopdf(i,2*M,j/M);
        end
    end
    
    for k=1:5000
        
        P1 = U*P1;
        P2 = U*P2;
        P1 = (NU.'*P1.').';
        P2 = (NU.'*P2.').';
        % Selection
        P1 = (w.*P1) / sum(w.*P1,'All');
        P2 = (w.*P2) / sum(w.*P2,'All');

        % Sexual distributions S
        for i=1:2 % from female a A
            for j=1:2 % from male a A
                Z(:,ZI(i,j))=PSI*(conv(P1(:,i),PSIP*P2(:,j))*2);
            end
        end
        Z = Z/sum(Z,'All');
        % Recombine
        Z=(R*Z.').';

        S1(:,1)=L2*L1*( Z(:,ZI(1,1))+Z(:,ZI(1,2)) );
        S1(:,2)=L2*L1*( Z(:,ZI(2,2))+Z(:,ZI(2,1)) );
        S2(:,1)=L2*L1*( Z(:,ZI(1,1))+Z(:,ZI(2,1)) );
        S2(:,2)=L2*L1*( Z(:,ZI(2,2))+Z(:,ZI(1,2)) );
        
        % asexual and sexual combined
        P1=L2*L1*AD*P1*(1-rsex)+S1*rsex;
        P2=L2*L1*AD*P2*(1-rsex)+S2*rsex;

        % record mean fitness
        W1(k) = sum(w.*P1,'All'); 
        W2(k) = sum(w.*P2,'All'); 
       
    end
    W = [W1(end); W2(end)]; 
end
    