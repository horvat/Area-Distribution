function DA = horvsmooth(DA,A)

DA0 = DA; 

normsmooth = 1; 

if normsmooth == 0
if max(DA(2:end)) ~= 0
s1 = sum(DA(2:end));
DA(2:end) = smooth(DA(2:end))';
s2 = sum(DA(2:end));
DA(2:end) = (s1/s2)*DA(2:end);
end
else
    DA1 = smooth(DA)'; 
    DA = norm(DA)/norm(DA1) * DA1; 
end


[C,I] = min(A + DA);

if C < 0 
    DA = DA0;
end

end