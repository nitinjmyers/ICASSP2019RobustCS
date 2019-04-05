function w=designawv(N,energ,qbit, Nrefine)
    w=exp(1i*11*pi*[0:1:N-1].^2/N)/sqrt(N);
    des_mag=sqrt(energ);
    for k=1:1:Nrefine
        DFTw=fft(w)/sqrt(N);
        DFTw=exp(1i*phase(DFTw)).*des_mag;
        w=exp(1i*phase(ifft(DFTw)*sqrt(N))); 
    end
 %   w=quantz(w,qbit); Quantize later - vary this and study
    w=w/norm(w);
end