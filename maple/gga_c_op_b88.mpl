(* type: work_gga_c *)

# Not that the files have to be included in this specific order
$define gga_x_b88_params
$include "gga_x_b88.mpl"

$include "op.mpl"

op_qab := 2.3670:
op_f   := xs -> f_b88(xs):
