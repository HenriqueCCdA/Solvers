#!/bin/sh
PATH_INCLUDE="./include"    
PATH_LIB="./lib"    
PRENAME=solver
CC=gcc
FC=ifort   
OPENMP=yes
OT=O3
DEBUG=no
#BLAS=/usr/local/atlas/lib

#------------------gerando o nome do excutavel-------------
ifeq ($(CC),icc)
  ifeq ($(OPENMP),yes)
    COMPILER_NAME=intel_omp_$(OT)
  else
    COMPILER_NAME=intel_$(OT)
  endif  
endif

ifeq ($(CC),gcc)
  ifeq ($(OPENMP),yes)
    COMPILER_NAME=gnu_omp_$(OT)
  else
    COMPILER_NAME=gnu_$(OT)
  endif  
endif
  
NAME+=$(PRENAME)_$(COMPILER_NAME)
#-------------------Fontes--------------------------------
fontes = \
src/Iters.c\
src/Time.c\
src/HccaBlas.c\
src/WriteMtx.c\
src/Graph.c\
src/Gauss.c\
src/Main.c
#-------------------Flags necessarios--------------------------------
NFLAGS=-I$(PATH_INCLUDE) -L$(PATH_LIB) 
LDFLAGS=-lmmio

#--------------------compiladores------------------------------------
# intel ifort
ifeq ($(CC),icc)
  LDFLAGS += 
  OFLAGS  +=  -ansi -std=c99 -mkl:parallel 
  ifeq ($(OPENMP),yes)
    OFLAGS  += -openmp
  endif
endif
# gnu gcc
ifeq ($(CC),gcc)
  LDFLAGS +=  -lm 
  OFLAGS  +=  -Wall -ansi -std=c99 -pedantic-errors
  ifeq ($(OPENMP),yes)
    OFLAGS  += -fopenmp
  endif
endif
#--------------------------------------------------------------------
#---------------------------Debug------------------------------------
ifeq ($(DEBUG),yes)
  OFLAGS += -g -$(OT)	
else
  OFLAGS += -$(OT) 
endif
#--------------------------------------------------------------------
FCFLAGS= $(NFLAGS) $(OFLAGS) -L$(BLAS) -L$(LAPACK)
CFLAGS = $(NFLAGS) $(OFLAGS) 


.SUFFIXES: 
.SUFFIXES: .c .h .o
OBJS  = $(patsubst %.c,%.o ,$(fontes))

build:	$(OBJS) 
	ls bin || mkdir -p bin
	$(CC) $(CFLAGS) $(OBJS) -o bin/$(NAME) $(LDFLAGS)

$(OBJS):%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $@

tags:
	ctags -R src/*.c include/*.h

.PHONY: cleantags
cleantags:
	@rm -fv tags
	
.PHONY: clean
clean:  
	@rm -fv src/*.o
	@rm -fv ../src/*.o
	@rm -fv bin/$(NAME)

.PHONY: cleanall
cleanall:  
	@rm -fv tags
	@rm -fv src/*.o
	@rm -fv ../src/*.o
	@rm -fv bin/$(NAME)


.PHONY: help
help:
	@echo "Autor :$(AUTHOR)                              "
	@echo "Makefile para prepar para sitemas linux.      "
	@echo -e "\E[7;32mOpcoes:\E[1;0m                      "
	@echo "build         - compila o prepar              "
	@echo "build_modules - gera os modulos               "
	@echo "tags          - gera os tags                  "
	@echo "cleantags     - limpa os tags                 "
	@echo "clean         - limpa os obj, bin e mod       "
	@echo "cleaall       - limpa tudo obj,bin,mod e tags "

# DO NOT DELETE

# DO NOT DELETE
