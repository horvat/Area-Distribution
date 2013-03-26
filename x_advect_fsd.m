function DA = x_advect_fsd(A,D,V1,V2,A1,A2)
    
    L = max(D);
    Ftop = V1*(A1);
    Fbottom = V2*(A);
    DA = (Fbottom-Ftop)/L; 
    
end