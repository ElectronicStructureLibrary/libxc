my $info = [
  [
    ['r_a_f','f(r_a_rs, r_a_z)']
  ], [
    ['r_a_dfdrs',    'diff(f(r_a_rs, r_a_z), r_a_rs)'],
    ['r_a_dfdz',     'diff(f(r_a_rs, r_a_z), r_a_z)']
  ], [
    ['r_a_d2fdrs2',  'diff(f(r_a_rs, r_a_z), r_a_rs$2)'],
    ['r_a_d2fdrsz',  'diff(f(r_a_rs, r_a_z), r_a_rs, r_a_z)'],
    ['r_a_d2fdz2',   'diff(f(r_a_rs, r_a_z), r_a_z$2)']
  ], [
    ['r_a_d3fdrs3',  'diff(f(r_a_rs, r_a_z), r_a_rs$3)'],
    ['r_a_d3fdrs2z', 'diff(f(r_a_rs, r_a_z), r_a_rs$2, r_a_z)'],
    ['r_a_d3fdrsz2', 'diff(f(r_a_rs, r_a_z), r_a_rs, r_a_z$2)'],
    ['r_a_d3fdz3',   'diff(f(r_a_rs, r_a_z), r_a_z$3)']
  ]
];
