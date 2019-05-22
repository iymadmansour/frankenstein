#*****************************************************************************
#
#  $Id: Makefile,v 1.1 2005/04/29 20:50:20 iwan Exp $
#
#*****************************************************************************

include $(EGS_CONFIG)
include $(SPEC_DIR)egspp.spec
include $(SPEC_DIR)egspp_$(my_machine).conf

DEFS = $(DEF1) -DBUILD_ELLIPSOIDS_DLL

library = egs_ellipsoids
lib_files = egs_ellipsoids

include $(SPEC_DIR)egspp_libs.spec

$(make_depend)
