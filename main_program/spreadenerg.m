function energ_vec=spreadenerg(N,loc)
    energ_vec=zeros(1,N);
    energ_vec(loc)=1/3;
    if(loc==1 || loc==N)
        if(loc==1)
            energ_vec(loc+1)=1/3;
            energ_vec(N)=1/3;
        else
            energ_vec(loc-1)=1/3;
            energ_vec(1)=1/3;
        end
        
    else
        energ_vec(loc-1)=1/3;
        energ_vec(loc+1)=1/3;
    end
end