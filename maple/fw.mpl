w  := t -> 1.0*(K_FACTOR_C - t)/(K_FACTOR_C + t):
fw := (t, n) -> sum(coeff_a[i]*w(t)^(i-1), i=1..n):