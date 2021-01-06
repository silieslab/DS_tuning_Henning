function f = logdet(A)

f = 2 * sum(log(diag(chol(A))));

end