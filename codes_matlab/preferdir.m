function ang = preferdir(vect)

[r1, a1] = resultant(vect,1);
[r2, a2] = resultant(vect,2);

if r1 >= r2
    ang = a1;
else
    if degreediff(a1,a2,1) < degreediff(a1,a2+180,1)
        ang = a2;
    else
        ang = a2 + 180;
    end
end
