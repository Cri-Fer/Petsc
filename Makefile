PETSC_DIR = $(CURDIR)/petsc
PETSC_ARCH = arch-linux-opt

CXX = mpicxx
SRC = main.cpp
OUT = main.exe

INCLUDES = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LIBS     = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -Wl,-rpath=$(PETSC_DIR)/$(PETSC_ARCH)/lib


.PHONY: deps build run clean

build: deps $(OUT)
	@echo ">> Compilazione completata ✅"

$(OUT): $(SRC)
	$(CXX) $(SRC) -o $(OUT) $(INCLUDES) $(LIBS)

run: $(OUT)
	mpirun -n 4 ./$(OUT) 
clean:
	rm -f $(OUT)

deps:
	@echo ">> Verifica installazione PETSc..."
	@if [ ! -d "$(PETSC_DIR)/.git" ]; then \
		echo ">> Clonazione di PETSc..."; \
		git clone -b release https://gitlab.com/petsc/petsc.git "$(PETSC_DIR)"; \
	fi
	@echo ">> Configurazione e compilazione di PETSc..."
	cd $(PETSC_DIR) && ./configure \
	--with-cc=mpicc \
	--with-cxx=mpicxx \
	--with-fc=mpif90 \
	--download-f2cblaslapack=1 \
	--download-hypre=1 \
	--with-debugging=0 \
	PETSC_ARCH=$(PETSC_ARCH)
	cd $(PETSC_DIR) && make all PETSC_ARCH=$(PETSC_ARCH) -j$(shell nproc)
	cd $(PETSC_DIR) && make check PETSC_ARCH=$(PETSC_ARCH)
	@echo ">> PETSc installato correttamente! ✅"
