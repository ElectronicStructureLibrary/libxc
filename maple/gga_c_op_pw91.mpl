(* type: work_gga_c *)

# Not that the files have to be included in this specific order
$define gga_x_pw91_params
$include "gga_x_pw91.mpl"

$include "op.mpl"

op_qab := 2.3706:
op_f   := xs -> f_pw91(xs):
