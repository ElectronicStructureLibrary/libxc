## Copyright (C) 2017 M. Marques
##
## This Source Code Form is subject to the terms of the Mozilla Public
## License, v. 2.0. If a copy of the MPL was not distributed with this
## file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

hgga_order_4 = 

hgga_order_3 = 

dist_noinst_HEADERS = $(hgga_order_4) $(hgga_order_3)

maple2c:
	@for funct in $(hgga_order_4); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.py --srcdir $(top_srcdir) --functional $${funct%.*} --maxorder 4; \
	done
	@for funct in $(hgga_order_3); do \
	  echo Processing: $${funct%.*}; \
	  $(top_srcdir)/scripts/maple2c.py --srcdir $(top_srcdir) --functional $${funct%.*} --maxorder 3; \
	done
