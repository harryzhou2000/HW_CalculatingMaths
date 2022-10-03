
I0b = 1- 0.3679;
I0 =  1- exp(-1);
Is = nan(1,10);
Isb = Is;
Is(1) = I0;
Isb(1) = I0b;


for i = 1:9
   Is(i + 1) = 1 - i * Is(i);
   Isb(i + 1) = round(1 - round(i * Isb(i), 4, 'significant'),4,'significant');
end
d = [Is', Isb'];
