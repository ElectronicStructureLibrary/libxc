(* type: work_gga_x *)

theta0 := 1.0008:
theta1 := 0.1926:
theta2 := 1.8962:

f0 := s -> s^2/(1 + s)^2:
f  := x -> theta0 + f0(X2S*x)* (theta1 + f0(X2S*x) * theta2):
