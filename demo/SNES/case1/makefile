
export PETSC_DIR=/home/wangzewei/PETSc/petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


CFLAGS += -pedantic -std=c99


SRCS := $(wildcard *.c)


OBJS := $(SRCS:.c=.o)


e: $(OBJS)
	-${CLINKER} -o e $(OBJS) ${PETSC_LIB}
	${RM} $(OBJS)


rune_1:
	-@../testit.sh e "" 1 1

test_e: rune_1

test: test_e


clean::
	@rm -f *~ e $(OBJS) *tmp

distclean: clean

.PHONY: clean distclean rune_1 test_e test

