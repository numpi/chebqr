function s = format_number(f)
%FORMAT_NUMBER 

if abs(f) <= 10 && abs(f) >= 1 || f == 0
    s = sprintf('%1.2f', f);
else
    exponent = floor(log10(abs(f)));
    
    f = f * 10^(-exponent);
    
    s = sprintf('%2.2f \\cdot 10^{%d}', f, exponent);
end

 