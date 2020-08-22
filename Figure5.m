R=[0.05:0.05:0.95,1.0];
l=[0:0.05:1];
n=size(R,2);
nn=size(l,2);
Wl=ones(1,nn);
Wf=zeros(1,n);
Wm=zeros(1,n);
dr=0.05;
for j=1:nn
    
    for i=1:n
        i
        if NucMTSexRate(0.01,0.01,0.01,R(i),R(i)-dr, 0.1,l(j), 0.01,0)<0
            Wl(j)=R(i)
            break
        end
    end    
end

plot(Wl,l)
