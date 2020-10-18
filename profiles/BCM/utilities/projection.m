function proj = projection(v,u)
%PROJECTION Projects the vector v into the line spanned by u.
proj = (v'*u)/(u'*u)*u;
end

