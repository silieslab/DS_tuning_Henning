function mm = mm_Collapse(mm, ix)

% if all classes are to be purged, randomly save one class
if(all(ix))
    ix(1) = false;
end

mm.nClasses = mm.nClasses - sum(ix);

mm.a(ix)    = [];
mm.a        = mm.a ./ sum(mm.a);
mm.class    = removecells(mm.class, ix);

end

