CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-g

#################### DEPENDENCIES ####################

DIPS_INCLUDE = $(shell root-config --cflags)
DIPS_LIBS = $(shell root-config --ldflags) $(shell root-config --libs)

CPPFLAGS+= -I $(DIPS_INCLUDE)
LDFLAGS += $(DIPS_LIBS) 

######################################################

SRCS=seedsProduction.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: tool

tool: $(OBJS)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) -o PSeeds $(OBJS) $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) *~ .depend

include .depend
