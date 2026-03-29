CC = gcc

CFLAGS = -std=c11 -Wall -pedantic -O3

OBJ =	AndersonDarlingTest.o Brent_fmin.o bessel_ratio.o dlrice.o eigens.o gofsignif.o  inverse.o \
	logBessel.o loglikelihood.o mat_vec.o matvecops.o nelder_mead_min.o order.o \
	quantile.o rice_emcluster.o rice_em_rndinit.o rice_em_variance.o run_rice_local_estimates.o \
	sorted.o run_rice_em_sample.o 

testrice_em: testrice_em.c $(OBJ)
	$(CC) $(CFLAGS) -o testrice_em testrice_em.c $(OBJ) -lm -lRmath -llapack

testevaluate_rice_em: testevaluate_rice_em.c $(OBJ)
	$(CC) $(CFLAGS) -o testevaluate_rice_em testevaluate_rice_em.c $(OBJ) -lm -lRmath -llapack

testrice_em_3d: testrice_em_3d.c $(OBJ)
	$(CC) $(CFLAGS) -o testrice_em_3d testrice_em_3d.c $(OBJ) -lm -lRmath -llapack

testlocal_sigma: testlocal_sigma.c $(OBJ)
	$(CC) $(CFLAGS) -o testlocal_sigma testlocal_sigma.c $(OBJ) -lm -lRmath -llapack

testevaluate_local_sigma: testevaluate_local_sigma.c $(OBJ)
	$(CC) $(CFLAGS) -o testevaluate_local_sigma testevaluate_local_sigma.c $(OBJ) -lm -lRmath -llapack

test_local_sigma_3d: test_local_sigma_3d.c $(OBJ)
	$(CC) $(CFLAGS) -o test_local_sigma_3d test_local_sigma_3d.c $(OBJ) -lm -lRmath -llapack


clean:	
	rm -rf testrice_em testevaluate_rice_em testrice_em_3d testevaluate_local_sigma testevaluate_local_sigma_3d testlocal_sigma test_local_sigma_3d $(OBJ) 

all: testrice_em testevaluate_rice_em testrice_em_3d testlocal_sigma testevaluate_local_sigma test_local_sigma_3d

