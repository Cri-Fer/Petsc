#include <petscksp.h>
#include <iostream>
#include <vector> 

/*
Questo codice usa le matrici estrapolate da Matlab e fa il calcolo. Controlla che il calcolo sia corretto
In fondo ci sono le varie norme di quanto i vari valori si discostano dalla soluzione di matlab e la exact solution
*/

int main(int argc, char **argv) {
  PetscInitialize(&argc, &argv, NULL, NULL);
  /*
    =======================================================================
                    START OF MATRIX AND VECTOR BUILDING
    =======================================================================
  */
  Mat A;
  Vec b, uhm; // uhm soluzione di matlab. uex exact solution
  // Open and load the matrix
  PetscViewer viewerA;
  int N = 1000;
  int p = 1;

  std::string pathA = "FilesLap/A_" + std::to_string(N) + "p" + std::to_string(p) + ".dat";
  std::string pathF = "FilesLap/F_" + std::to_string(N) + "p" + std::to_string(p) + ".dat";
  std::string pathUm = "FilesLap/Um_" + std::to_string(N) + "p" + std::to_string(p) + ".dat";
  std::string patherr = "err/lap/err_L2_" + std::to_string(N) + "p" + std::to_string(p) + ".dat";


  PetscViewerBinaryOpen(PETSC_COMM_WORLD, pathA.c_str() , FILE_MODE_READ, &viewerA);
 
  MatCreate(PETSC_COMM_WORLD, &A);       // <--- CREA IL CONTENITORE PRIMA
  MatLoad(A, viewerA);                   // <--- POI CARICA I DATI
  PetscViewerDestroy(&viewerA);

  // Open and load the vector
  PetscViewer viewerB;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, pathF.c_str(), FILE_MODE_READ, &viewerB);

  VecCreate(PETSC_COMM_WORLD, &b);       // <--- CREA IL CONTENITORE PRIMA
  VecLoad(b, viewerB);
  PetscViewerDestroy(&viewerB);

  PetscViewer viewerC;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, pathUm.c_str(), FILE_MODE_READ, &viewerC);

  VecCreate(PETSC_COMM_WORLD, &uhm);       // <--- CREA IL CONTENITORE PRIMA
  VecLoad(uhm, viewerC);
  PetscViewerDestroy(&viewerC);

  /*
    =======================================================================
                      EDN OF MATRIX AND VECTOR BUILDING
    ======================================================================= 
  */

  Vec uh, err; // uh soluzione approssimata da petsc
  KSP ksp;
  std::vector<PetscReal> residual = {1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-12, 1e-14, 1e-16, 1e-20};
  Vec err_L2 ;
  PetscViewer vectorViewer;

  VecDuplicate(b, &uh);
  VecDuplicate(b, &err);
  VecCreate(PETSC_COMM_WORLD, &err_L2);
  VecSetSizes(err_L2, PETSC_DECIDE, residual.size());
  VecSetFromOptions(err_L2);

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, A, A);
  KSPSetType(ksp, KSPCG); //gradiente coniugato come solver

  PC pc;
  KSPGetPC(ksp, &pc); // Dice "lega il precond di ksp a pc"
  //PCSetType(pc, PCGAMG); // Setta il precondizionatore a default AMG
  //KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // dettagli del Krilov method
  PCSetType(pc, PCHYPRE);
  PCHYPRESetType(pc, "boomeramg");  // Precond: AMG Hypre
  PetscInt its;
  PetscReal rnorm, e;

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  for (PetscInt i = 0; i < residual.size(); ++i) {
    VecSet(uh, 0.0);

    PetscPrintf(PETSC_COMM_WORLD, "\n======== Test %.d ========\n", i + 1);

    KSPSetTolerances(ksp, residual[i], 1e-30, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSolve(ksp, b, uh);

    KSPGetResidualNorm(ksp, &rnorm);
    PetscPrintf(PETSC_COMM_WORLD, "- Final residual norm = %.6e\n", (double)rnorm);
    PetscPrintf(PETSC_COMM_WORLD, "- Residual \t = %e\n", (double)residual[i]);
    VecWAXPY(err, -1.0, uh, uhm);
    VecNorm(err, NORM_2, &e);
    VecSetValues(err_L2, 1, &i, &e, INSERT_VALUES);


        
    KSPGetIterationNumber(ksp, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Converged in %d iterations\n", (int)its);
    
    PetscPrintf(PETSC_COMM_WORLD, "Convergence reason: %s\n", KSPConvergedReasons[reason]);
    PetscPrintf(PETSC_COMM_WORLD, "|| uhm - uh||_2 = %.5e \n", double(e));
  }

  VecAssemblyBegin(err_L2);
  VecAssemblyEnd(err_L2);
  // The solution conververges in one step because LU factorization is not an iterative method, but a direct one
  // so we are directly solving the system, and not finding a possible one, and at each iteration we change it
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, patherr.c_str(), FILE_MODE_WRITE, &vectorViewer);
  VecView(err_L2, vectorViewer);
  VecDestroy(&uh);
  VecDestroy(&b);
  MatDestroy(&A);
  KSPDestroy(&ksp);


  PetscFinalize();
  return 0;
}
