function [W] = NucMTPL(mu,muu,nu,rsex,r,pl,plr,fr,sex)

    pl=round(pl,2);
    plr=round(plr,2);
    M=20; % M is the number of discrete mitochondria per cell
    n=2;  % number of possible nuclear valures: aa aA AA
    x=2;  % x is the exponent in the fitness function. x>1 makes it concave which is what we want.
    s=1;  % selection strength
    
    % Initialize population P with random distributions
    % Maternal population, females
    P1=zeros(M+1,n);
    p = rand([M+1 2]);
    P1(:,[1:2])=1/2 * p./sum(p);
        
    % Paternal population, males
    P2=zeros(M+1,n);
    p = rand([M+1 2]);
    P2(:,[1:2])=1/2*p./sum(p); %Maa %MaA %MAA
    
    % Paternal population, males mutants
    P2R=zeros(M+1,n);
    P2R(:,[1:2])=0*1/2*p./sum(p);
    
    % Life cycle operators
    % Mitochondrial mutation
    % Both ways, rates mu and muu
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
    
    % w is a matrix version of the fitness function
    % w(i+1,j) is the fitness of a cell with i mito mutations and j nuclear
    % genotype
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
    % PSIP resamples the number of mito mutants without replacement
    PSIP=zeros(0,1);
    for i=0:pl*M
        for j=0:M
            PSIP(i+1,j+1)=hygepdf(i,M,j,round(pl*M));
        end
    end

    PSIr=zeros(0,1);
    for i=0:2*M
        for j=0:round(plr*M)+M
            PSIr(i+1,j+1)=binopdf(i,2*M,j/(round(plr*M)+M));
        end
    end

    PSIPr=zeros(0,1);
    for i=0:plr*M
        for j=0:M
            PSIPr(i+1,j+1)=hygepdf(i,M,j,round(plr*M));
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
    
    % Now iterate the life cycle for many many generations k
    for k=1:2000

        if k==1
            P2R=fr*P2;
            P2=P2-fr*P2;
        end

        % Mutate mitochondria from a to A, from A to a
        P1 = U*P1;
        P2 = U*P2;
        P2R = U*P2R;
        
        % Nuclear mutation between aa AA aA 
        P1 = (NU.'*P1.').';
        P2 = (NU.'*P2.').';
        P2R = (NU.'*P2R.').';
        
        % Selection
        P1 = (w.*P1) / sum(w.*P1,'All');
        P2P2R = ([w,w].*[P2,P2R]) / sum([w,w].*[P2,P2R],'All');
        P2 = P2P2R(:,1:2);
        P2R = P2P2R(:,3:4);

        % Sexual distributions S
        for i=1:2 % from female a A
            for j=1:2 % from male a A
                if sex==0
                 Z(:,ZI(i,j))=PSI*(conv(P1(:,i),PSIP*P2(:,j))*2);
                 ZR(:,ZI(i,j))=PSIr*(conv(P1(:,i),PSIPr*P2R(:,j))*2);
                end
                if sex==1
                 Z(:,ZI(i,j))=PSI*(conv(PSIP*P1(:,i),P2(:,j))*2);
                 ZR(:,ZI(i,j))=PSIr*(conv(PSIPr*P1(:,i),P2R(:,j))*2);
                end
                
            end
        end
        Zsum = sum(Z+ZR,'All');
        Z = Z/Zsum;
        ZR = ZR/Zsum;
        
        % Recombine
        Z=(R*Z.').';
        ZR=(R*ZR.').';

        S1(:,1)=L2*L1*( Z(:,ZI(1,1))+Z(:,ZI(1,2)) + ZR(:,ZI(1,1))+ZR(:,ZI(1,2)) );
        S1(:,2)=L2*L1*( Z(:,ZI(2,2))+Z(:,ZI(2,1)) + ZR(:,ZI(2,2))+ZR(:,ZI(2,1)) );
        
        S2(:,1)=L2*L1*( Z(:,ZI(1,1))+Z(:,ZI(2,1)) );
        S2(:,2)=L2*L1*( Z(:,ZI(2,2))+Z(:,ZI(1,2)) );
        
        S2R(:,1)=L2*L1*( ZR(:,ZI(1,1))+ZR(:,ZI(2,1)) );
        S2R(:,2)=L2*L1*( ZR(:,ZI(2,2))+ZR(:,ZI(1,2)) );
        
        % asexual and sexual combined
        P1=L2*L1*AD*P1*(1-rsex)+S1*rsex;
        P2=L2*L1*AD*P2*(1-rsex)+S2*rsex;
        P2R=L2*L1*AD*P2R*(1-rsex)+S2R*rsex;
        
        % if we need to keep mutant frequency constant to measure fitness
        % advantage, set it to fr here:
        P2 = P2/sum(P2,'All')*(1-fr);
        P2R = P2R/sum(P2R,'All')*(fr);
        
        % record mean fitness or frequencies
        W1(k) = sum(w.*P1,'All') ;
        W2(k) = sum([w,w].*[P2,P2R],'All');
        % fitness advantage of the mutant to return
        W(k) = sum(w.*P2R,'All')/sum(P2R,'All')-sum(w.*P2,'All')/sum(P2,'All') ;
        % or frequency of the mutant to return
        %W(k) = sum(P2R,'All');

    end
    W=W(end);
end
    