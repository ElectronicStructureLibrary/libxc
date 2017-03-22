(* type: work_gga_c *)

# Not that the files have to be included in this specific order
$define gga_x_pbe_params
$include "gga_x_pbe.mpl"

$include "op.mpl"

op_qab := 2.3789:
op_f   := xs -> f_pbe(xs):
