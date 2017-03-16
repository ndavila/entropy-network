#///////////////////////////////////////////////////////////////////////////////
#  Copyright (c) 2017 Clemson University.
# 
#  This file was originally written by Bradley S. Meyer.
# 
#  This is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#  USA
# 
#///////////////////////////////////////////////////////////////////////////////

#///////////////////////////////////////////////////////////////////////////////
#//!
#//! \file Makefile
#//! \brief A makefile to generate network examples.
#//!
#///////////////////////////////////////////////////////////////////////////////

ifndef NUCNET_TARGET
NUCNET_TARGET = ../..
endif

NNT_DIR = $(NUCNET_TARGET)/nnt
USER_DIR = $(NUCNET_TARGET)/user
BUILD_DIR = $(NUCNET_TARGET)/build

ifndef DATA_DIR
DATA_DIR = $(NUCNET_TARGET)/data_pub
endif

ifndef OBJDIR
OBJDIR = $(NUCNET_TARGET)/obj
endif

ifndef VENDORDIR
VENDORDIR = $(NUCNET_TARGET)/vendor
endif

#///////////////////////////////////////////////////////////////////////////////
# End of lines to be edited.
#///////////////////////////////////////////////////////////////////////////////

include $(BUILD_DIR)/Makefile

include $(BUILD_DIR)/Makefile.sparse

include $(USER_DIR)/Makefile.inc

VPATH = $(BUILD_DIR):$(NNT_DIR):$(USER_DIR)

#===============================================================================
# Initial objects.
#===============================================================================

MY_HYDRO_OBJ = $(OBJDIR)/my_hydro_helper.o                 \

$(MY_HYDRO_OBJ): $(OBJDIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

NETWORK_OBJS = $(WN_OBJ)        \
               $(NNT_OBJ)	\
               $(HYDRO_OBJ)	\
               $(MY_HYDRO_OBJ)	\
               $(SOLVE_OBJ)	\
               $(USER_OBJ)      \

#===============================================================================
# Use Sparskit2, if desired.  NNT_USE_SPARSKIT2 is an environment variable.
# In a bash shell, set this by typing at the command line, for example,
# 'export NNT_USE_SPARSKIT2=1'.  If you use Sparskit2, you must set zone
# properties "iterative solver method" (for example, gmres) and
# "t9 for iterative solver" (for example, 2.0).  This latter property is
# the t9 below which to user Sparskit2.
#===============================================================================

ifdef NNT_USE_SPARSKIT2
  CFLAGS += -DSPARSKIT2
  NETWORK_OBJS += $(SP_OBJ) $(ILU_OBJ)
  NET_DEP = sparse
  FLIBS= -L$(SPARSKITDIR) -lskit -lstdc++
  MC = $(FF)
  ifeq ($(GC), clang++)
    CC+= -stdlib=libstdc++
  endif
else
  MC = $(CC)
endif

#===============================================================================
# Final network dependencies.
#===============================================================================

NET_DEP += $(NETWORK_OBJS)

#===============================================================================
# Executables.
#===============================================================================

NETWORK_EXEC = run_entropy

$(NETWORK_EXEC): $(NET_DEP)
	$(CC) -c -o $(OBJDIR)/$@.o $@.cpp
	$(MC) $(NETWORK_OBJS) $(OBJDIR)/$@.o -o $(BINDIR)/$@ $(CLIBS) $(FLIBS)

.PHONY all_entropy : $(NETWORK_EXEC)

#===============================================================================
# Clean up.
#===============================================================================

.PHONY: clean_entropy cleanall_entropy

clean_entropy: 
	rm -f $(NETWORK_OBJS)

cleanall_entropy: clean_entropy
	rm -f $(BINDIR)/$(NETWORK_EXEC) $(BINDIR)/$(NETWORK_EXEC).exe

#===============================================================================
# Define.
#===============================================================================

NETWORK_DEF = yes
