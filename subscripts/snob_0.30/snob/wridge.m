function b = wridge(X, y, w, lambda)

%% weighted ridge regression
%% see https://www.nber.org/papers/w0011
n = sum(w);
xw = bsxfun(@times, X, w);
ym = sum(w .* y) / n;   % weighted mean of y
xm = sum(xw, 1) ./ n;   % weighted mean of each column of X

ystd = y - ym;                      % subtract weighted mean from X and y
xstd = bsxfun(@minus, X, xm);

xstdw = bsxfun(@times, xstd, sqrt(w));
ystdw = w .* ystd;
s2 = sum(bsxfun(@times, xstd.^2, w), 1);    % weighted length of x

XtX = xstdw'*xstdw + lambda*diag(s2);
Xty = (xstd'*ystdw);
b = XtX \ Xty;
b = [ym - xm*b; b];

end