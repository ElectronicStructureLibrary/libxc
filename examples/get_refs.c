#include <stdio.h>

#include <xc.h>

int main()
{
  xc_func_type func;

  xc_func_init(&func, XC_GGA_X_B88, XC_UNPOLARIZED);

  printf("The functional '%s' is defined in the reference(s):\n%s\n", func.info->name, func.info->refs);

  xc_func_end(&func);
}
