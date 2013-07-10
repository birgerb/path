
#include <math.h>
#include <apop.h>
/* #include <./trace.h> */

apop_data   *DataGl = NULL;
int *useFinishedGl = NULL;
int NoOfCombGl = 0;
int NoOfEstGl = 0;
int NoOfUseGl = 0;

int BBfindVariable(apop_data *Data, char *Name, int *Pos)
{
	int i = 0;

	*Pos = -1;
	while ((i < Data->matrix->size2) && (*Pos == -1)) {
	       	if (strcmp(Name, Data->names->column[i]) == 0) {
			*Pos = i;
		}
		i++;
        }
}

int BBmergeUse(int *ThisUseMax, apop_data *ThisDataLevel, int *UseCombined)
{
	int i = 0, pos = -1;

	if ((ThisUseMax != (int *)NULL) && (UseCombined != (int *)NULL)) {
		for (i = 1; i < ThisDataLevel->matrix->size2; i++) {
			if (ThisUseMax[i] == 0) { 
				BBfindVariable(DataGl, ThisDataLevel->names->column[i], &pos);
				if (pos != -1) {
					UseCombined[pos] = 0;
					/*
						printf("BBmergeUse: UseCombined: ");
        					fflush(stdout);
						BBprintUse(UseCombined, DataGl->matrix->size2, -1, 1);
					*/
				}
				else {
					/*
						printf("ERROR in BBmergeUse: Cannot find %s\n", ThisDataLevel->names->column[i]);
						fflush(stdout);
					*/
				}
			}
		}
	}
	else if (ThisUseMax == (int *)NULL) {
		/*
			printf("ERROR in BBmergeUse: ThisUseMax is NULL.\n");
			fflush(stdout);
		*/
	}
	else if (UseCombined == (int *)NULL) {
		/*
			printf("ERROR in BBmergeUse: UseCombined is NULL.\n");
			fflush(stdout);
		*/
	}
}


int BBcopyUse(int *Use, int *UseHere, int Size)
{
	int i;

	for (i = 0; i < Size; i++) {
		UseHere[i] = Use[i];
	}

}

int BBinitUse(int *Use, int Size)
{
	int i;

	for (i = 0; i < Size; i++) {
		Use[i] = 1;
	}

}

int BBcheckUse(int *Use, int Size, int *InUse)
{
	int i = 0;

	*InUse = 1;
	while ((i < Size) && (*InUse == 1)) {
		if (Use[i] == 0) {
			*InUse = 0;
		}
		i++;
	}

}

int BBgetUseIndep(int *Use, int Size, int *IndepSize, int *FirstIndepPos)
{
	int i = 1;

	*IndepSize = 0;
	*FirstIndepPos = 0;
	while (i < Size) {
		if (Use[i] == 0) {
			if (*FirstIndepPos == 0) {
				*FirstIndepPos = i;
			}
			*IndepSize = *IndepSize + 1;
		}
		i++;
	}

}

int BBprintIndent(int No)
{
	int i;
	int retVal = 1;

	for (i = 0; i < No; i++) {
		printf("+ ");
		fflush(stdout);
	}
error:
return(retVal);
}

int BBprintUse( int *Use, int Size, int NoOfComb, int Always )
{
	int i;
	int retVal = 1;

	if (Use == (int *)NULL) {
		printf( "ERROR: Use is NULL\n" );
		fflush(stdout);
		retVal = 0;
		goto error;
	}
	if ((NoOfCombGl % 10000 == 0) || (Always == 1)) {
		fprintf(stderr, "NoOfComb: %d: ", NoOfComb);
		for (i = 0; i < Size; i++) {
			if (Use[i] == 1) {
				fprintf(stderr, " ");
			}
			else if (Use[i] == 0) {
				fprintf(stderr, "x");
			}
			else {
				fprintf(stderr, "?");
				fprintf(stderr, "%d", Use[i]);
			}
		}
		fprintf(stderr, ".\n");
		fflush(stderr);
	}
error:
return(retVal);
}


int BBcalcFirstPart(apop_data *OrigData, int *Use, int OrigIndep, int Control)
{
	apop_data       *dataDep = NULL;
	apop_data       *dataIndep = NULL;
	apop_data       *newData = NULL;
	apop_data       *dataClean = NULL;
	apop_model   *est = NULL;
	int i, j, indep, control;
	double tVal, coeff;
	int sign = 1;
	int *useNext = NULL;
	int matrixLarge = 3;

	printf("---------------------------------------------\n");
	fflush(stdout);
	printf("NR. %d-%d\n", OrigIndep, Control);
	fflush(stdout);
	printf("NORSK: Inkluderer testvariabel. Regresjon av første part fra original uavhengig nr.: %d til kontrollvar nr.: %d: ", OrigIndep, Control);
	fflush(stdout);
       	dataDep   = apop_data_copy(OrigData);
       	dataIndep   = apop_data_copy(OrigData);
	BBinitUse(Use, OrigData->matrix->size2);
	Use[Control] = 0;
	apop_data_rm_columns(dataDep, Use);
	memset(Use, 1, DataGl->matrix->size2 * sizeof(int));
	Use[OrigIndep] = 0;
	apop_data_rm_columns(dataIndep, Use);
	newData = apop_data_stack(dataDep, dataIndep, 'c', 'n');
       	dataClean   = apop_data_listwise_delete(newData); 
	if ((dataClean != NULL) && (dataClean->matrix->size1 > 30)) {
		est  = apop_estimate(dataClean, apop_ols); 
		for (i=1; i < est->parameters->matrix->size1; i++) { 
			tVal = apop_data_get(est->parameters, i, 2); 		
			coeff = gsl_vector_get(est->parameters->vector, i); 		
			printf("NORSK:         %.2f (%.2f)      \\ \n", coeff, tVal); 
		}
	}
	fflush(stdout);
	apop_data_free(newData);
	apop_data_free(dataDep);
	apop_data_free(dataIndep);
	apop_data_free(dataClean);
	fflush(stdout);
}

int BBcontrol(double OrigCoeff, double OrigT, char *OrigDepName, char *OrigIndepName)
{
	apop_data       *dataDepend = NULL, *dataIndepend = NULL, *dataCopy = NULL;
	apop_data       *dataClean = NULL, *dataControl = NULL;
	apop_model   *est = NULL;
	int i, j, indep, control;
	int *useCopy;
	double tVal, coeff;
	int sign = 1;
	int *dataControlUse = NULL;
	int matrixLarge = 3;
	char *indepLevelStr = NULL, *levelStr = NULL;
	char nameCopy[200];
	int indepLevel, noOfIndep = 0, pos, indepPos, origIndepPos;

	dataDepend = apop_data_copy(DataGl);
	useCopy = calloc(DataGl->matrix->size2 , sizeof(int));
	memset(useCopy, 0, DataGl->matrix->size2 * sizeof(int));
	BBinitUse(useCopy, DataGl->matrix->size2);
	BBfindVariable(DataGl, OrigDepName, &pos);
	if (pos > -1) {
		useCopy[pos] = 0;
	}
	else {
		useCopy[0] = 0;
	}
       	apop_data_rm_columns(dataDepend, useCopy);
	dataIndepend = apop_data_copy(DataGl);
	BBfindVariable(DataGl, OrigIndepName, &indepPos);
	if (indepPos > -1) {
		useCopy[indepPos] = 0;
		indepLevelStr = strtok(OrigIndepName, "_");
		if (indepLevelStr != (char *)NULL) {
			indepLevel = atoi(indepLevelStr);
			for (j = 0; j < DataGl->matrix->size2; j++) {
	       			strncpy(nameCopy, DataGl->names->column[j], 198);
	       			levelStr = strtok(nameCopy, "_");
				if ((levelStr != (char *)NULL) && (atoi(levelStr) > indepLevel)) {
				/**/
					printf("Found indep: %s\n", DataGl->names->column[j]);
					fflush(stdout);
				/**/
					useCopy[j] = 0;
					noOfIndep++;
				}
				else {
					useCopy[j] = 1;
				}
        		}
			if (noOfIndep > 0) {
				printf("Control use: ");
				fflush(stdout);
				BBprintUse(useCopy, DataGl->matrix->size2, -1, 1);
        			apop_data_rm_columns(dataIndepend, useCopy);
        			dataControl = apop_data_stack(dataDepend, dataIndepend, 'c', 'n');
				dataControlUse = calloc(dataControl->matrix->size2 , sizeof(int));
				BBfindVariable(dataControl, OrigIndepName, &origIndepPos);
				if (origIndepPos < 0) {
					printf("ERROR. Cannot find OrigIndepName in dataControl!\n");
					fflush(stdout);
					goto error;
				}
				for (indep = 1; indep < dataControl->matrix->size2; indep++) {
					if (indep != origIndepPos) {
						printf("---------------------------------------------\n");
						fflush(stdout);
						printf("NR. %d-%d\n", origIndepPos, indep);
						fflush(stdout);
						BBinitUse(dataControlUse, dataControl->matrix->size2);
						dataControlUse[0] = 0;
						dataControlUse[origIndepPos] = 0;
						dataControlUse[indep] = 0;
						printf("Control use: ");
						fflush(stdout);
						BBprintUse(dataControlUse, dataControl->matrix->size2, -1, 1);
						printf("NORSK: Inkluderer testvariabel. Regresjon mellom avhengig variabel: %s og uavhengige variabel nr. %d: %s med testvariabel nr. %d: ", dataControl->names->column[0], origIndepPos, OrigIndepName, indep);
						fflush(stdout);
       						dataCopy   = apop_data_copy(dataControl);
						apop_data_rm_columns(dataCopy, dataControlUse); 
						for (j = 1; j < dataCopy->matrix->size2; j++) {
							if (strcmp(OrigIndepName, dataCopy->names->column[j]) != 0) {
								printf("%s, ", dataCopy->names->column[j]);
							}
							printf("\n");
							fflush(stdout);
						}
        					dataClean   = apop_data_listwise_delete(dataCopy); 
						if ((dataClean == NULL) || (dataClean->matrix->size2 > 2)) {
						matrixLarge = 3;
					}		
					else {
						matrixLarge = 2;
					}
					if ((dataClean != NULL) && (dataClean->matrix->size1 > 30) && (matrixLarge == 3)) {
						est  = apop_estimate(dataClean, apop_ols); 
						sign = 1;
						control = 0;
						for (i=1; i < est->parameters->matrix->size1; i++) { 
							tVal = apop_data_get(est->parameters, i, 2); 		
							coeff = gsl_vector_get(est->parameters->vector, i); 		
							if ((tVal < -1000.0) || (tVal > 1000.0)) {
								sign = 2;
								printf("NORSK: finner en t-verdi som er ekstremt høy for: %s. Mest sannsynlig singular matrise slik at OLS har brutt sammen: %f\n", est->parameters->names->row[i], tVal); 
								fflush(stdout);
							}
							else if ((tVal < 2.0) && (tVal > -2.0)) {
								sign = 0;
								printf("NORSK: finner en ikke ok t-verdi for: %s: %f\n", est->parameters->names->row[i], tVal); 
								fflush(stdout);
							}
							else {
								printf("NORSK: finner ok t-verdi for: %s: %f\n", est->parameters->names->row[i], tVal);  
								fflush(stdout);
							}
							if ((strcmp(OrigIndepName, est->parameters->names->row[i]) == 0) && (OrigT >= 2.0) && (tVal < 0.5)) { 
								control = 1;
							}
							if ((strcmp(OrigIndepName, est->parameters->names->row[i]) == 0) && (OrigT <= -2.0) && (tVal > -0.5)) { 
								control = 1;
							}
							if ((strcmp(OrigIndepName, est->parameters->names->row[i]) == 0) && (OrigT <= 0.5) && (tVal >= 2.0)) { 
								control = 2;
							}
							if ((strcmp(OrigIndepName, est->parameters->names->row[i]) == 0) && (OrigT >= -0.5) && (tVal <= -2.0)) { 
								control = 2;
							}
						}
						if ((sign == 0) && (control == 1))  {
							printf("+++++++++ Forsterker ++++++++++\n"); 
							printf("NR. %d-%d\n", origIndepPos, indep);
							printf("!!!!!: Opprinnelig modell: %s ---> %.2f (%.2f) ---> %s\n", OrigIndepName, OrigCoeff, OrigT, OrigDepName); 
							printf("!!!!!: testvariabelen sender t for: %s under 0,5 og forklarer dermed opprinnelig sammenheng.\n", OrigDepName); 
							for (i=1; i < est->parameters->matrix->size1; i++) {	 
								if (strcmp(OrigDepName, est->parameters->names->row[i]) != 0) {
									tVal = apop_data_get(est->parameters, i, 2); 		
									coeff = gsl_vector_get(est->parameters->vector, i); 		
									printf("!!!!!:                 %s \n", est->parameters->names->row[i]); 
									printf("!!!!!:                /        \\\n"); 
									BBcalcFirstPart(dataControl, dataControlUse, origIndepPos, indep);
									printf("!!!!!:         /                %.2f (%.2f) \n", coeff, tVal); 
									printf("!!!!!:        /                     \\\n"); 
								}
								fflush(stdout);
							}
							for (i=1; i < est->parameters->matrix->size1; i++) { 
								if (strcmp(OrigIndepName, est->parameters->names->row[i]) == 0) {
									tVal = apop_data_get(est->parameters, i, 2); 		
									coeff = gsl_vector_get(est->parameters->vector, i); 		
									printf("!!!!!: %s ----> %.2f (%.2f) ----> %s\n", est->parameters->names->row[i], coeff, tVal, OrigDepName); 
								}
							}
							fflush(stdout);
						}	
						if ((control == 2))  {
							printf("--------- Bremser ----------\n"); 
							printf("NR. %d-%d\n", origIndepPos, indep);
							printf("!!!!!: Opprinnelig modell: %s ---> %.2f (%.2f) ---> %s\n", OrigIndepName, OrigCoeff, OrigT, OrigDepName); 
							printf("!!!!!: testvariabelen øker: %s over 0,5 og forklarer hvorfor det ikke var opprinnelig sammenheng.\n", OrigDepName); 
							for (i=1; i < est->parameters->matrix->size1; i++) { 
								if (strcmp(OrigDepName, est->parameters->names->row[i]) != 0) {
									tVal = apop_data_get(est->parameters, i, 2); 		
									coeff = gsl_vector_get(est->parameters->vector, i); 		
									printf("NORSK:                 %s \n", est->parameters->names->row[i]); 
									printf("NORSK:                /        \\\n"); 
									BBcalcFirstPart(dataControl, dataControlUse, origIndepPos, indep);
									printf("NORSK:         /                %.2f (%.2f) \n", coeff, tVal); 
									printf("NORSK:        /                     \\\n"); 
								}
								fflush(stdout);
							}
							for (i=1; i < est->parameters->matrix->size1; i++) { 
								if (strcmp(OrigIndepName, est->parameters->names->row[i]) == 0) {
									tVal = apop_data_get(est->parameters, i, 2); 		
									coeff = gsl_vector_get(est->parameters->vector, i); 		
									printf("NORSK: %s ----> %.2f (%.2f) ----> %s\n", est->parameters->names->row[i], coeff, tVal, OrigDepName); 
								}
							}
							fflush(stdout);
						}
						if ((sign == 1) && (control != 0))  {
							printf("....................................\n"); 
							printf("NORSK: Signifikant modell\n");
							apop_model_show(est);
							fflush(stdout);
						}
						if ((sign == 1) && (control == 0))  {
							printf("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY\n"); 
							printf("NR. %d-%d\n", origIndepPos, indep);
							fflush(stdout);
							printf("NORSK: testvariabelen sender ikke: %s under 0,5 men den nye modellen er signifikant\n", OrigDepName); 
							apop_model_show(est);
							fflush(stdout);
						}
						apop_model_free(est);
						est = NULL;
					}
					else {
						printf("NORSK: Uavhengig varabel nr.: %d gir liten eller tom matrise\n", indep);
					}	
					if (dataClean != NULL) {
						apop_data_free(dataClean);
						dataClean = NULL;
					}
					if (dataCopy != NULL) {
						apop_data_free(dataCopy);
						dataCopy = NULL;
					}
				}
			}
		}
	}
}
	printf("NORSK: slutten av control-regresjon for avhengig variabel: %s\n", dataControl->names->column[0]);
error:
	fflush(stdout);
	
}


int BBregress(int *Use, int J, int CombLevel, int LevelDep, apop_data *DataLevel, int **UseMax, double *MaxLL, int *MaxNoOfComb)
{
	apop_data       *dataCopy = NULL;
	apop_data       *dataClean = NULL;
	apop_model   *est = NULL;
	int i, j;
	double tVal, origT, coeff;
	char origName[200];
	int sign = 1;
	int *useNext = NULL;
	int matrixLarge = 3;
	int indepSize = 0, firstIndepPos = 0;
	/* double logLikelihood = 0.0; */

	/* printf("BBregress start\n");
	fflush(stdout); 
	*/
	CombLevel++;
	NoOfCombGl ++;
	BBprintUse(Use, DataLevel->matrix->size2, NoOfCombGl, 0);
					/*
					printf("................................................................BBregress start, level: %d, dependent var: %s, indep vars: ", CombLevel, DataLevel->names->column[0]);
					fflush(stdout);
					*/
       	dataCopy   = apop_data_copy(DataLevel);
	if (CombLevel != 1) {
		apop_data_rm_columns(dataCopy, Use); //leaks memory, fixed
	}
				/*
				for (j = 1; j < dataCopy->matrix->size2; j++) {
					printf("%s, ", dataCopy->names->column[j]);
					fflush(stdout);
				}
				printf("................................................\n");
				fflush(stdout);
				*/
        dataClean   = apop_data_listwise_delete(dataCopy); //leaks memory, fixed
	if ((dataClean == NULL) || (dataClean->matrix->size2 > 2)) {
		matrixLarge = 3;
	}
	else {
		matrixLarge = 2;
	}
	if ((dataClean != NULL) && (dataClean->matrix->size1 > 15)) {
		est  = apop_estimate(dataClean, apop_ols); // Leaks memory
		NoOfEstGl ++;
		sign = 1;
		for (i=1; i < est->parameters->matrix->size1; i++) { 
			tVal = apop_data_get(est->parameters, i, 2); 		
			coeff = gsl_vector_get(est->parameters->vector, i);
    			origT = tVal;
                        strcpy(origName, est->parameters->names->row[1]);
                        printf("origName: %s\n", origName);
                        fflush(stdout);
			if (isnan(tVal)) {
				sign = 2;
			}
			else if ((tVal < -1000.0) || (tVal > 1000.0)) {
				sign = 2;
				BBprintIndent(LevelDep);
				printf("NORSK: finner en t-verdi som er ekstremt høy for: %s. Mest sannsynlig singular matrise (kollineæritet) slik at OLS har brutt sammen: %f\n", est->parameters->names->row[i], tVal); 
				fflush(stdout);
			}
			else if ((tVal < 1.6) && (tVal > -1.6)) {
				sign = 0;
							/*
							 printf("t-verdi under: %f\n", tVal); 
							 fflush(stdout);
							*/
			}
			else {  
				// OK
			}
		}
		if (sign == 1) {
			/*
			printf("before logLikelihood\n");
			fflush(stdout);
			logLikelihood = apop_log_likelihood(dataClean, est);
			printf("logLikelihood: %.2f\n", logLikelihood);
			printf("after logLikelihood\n");
			fflush(stdout);
			*/
			if ((est->llikelihood <= 0.0) && (est->llikelihood > -170.0)) {
				printf("\n");
				BBprintIndent(LevelDep);
				printf("Dependent var: %s -----------------------  ", dataCopy->names->column[0]);
				if (est->llikelihood > *MaxLL) {
					if (*UseMax != (int *)NULL) {
						free(*UseMax);
						*UseMax = NULL;
					}
					*UseMax = calloc(DataLevel->matrix->size2 , sizeof(int));
					memset(*UseMax, 0, DataLevel->matrix->size2 * sizeof(int));
					*MaxLL = est->llikelihood;
					*MaxNoOfComb = NoOfCombGl;
					BBprintIndent(LevelDep);
					printf("*** MAX LL så langt *** ");
					fflush(stdout);
					BBcopyUse(Use, *UseMax, DataLevel->matrix->size2);
					BBprintUse(*UseMax, DataLevel->matrix->size2, *MaxNoOfComb, 1);
				}
				BBprintIndent(LevelDep);
				printf("Regresjonsnr: %d", NoOfEstGl);
				printf("   Kombinasjonsnr: %d\n", NoOfCombGl);
				for (i=1; i < est->parameters->matrix->size1; i++) { 
					tVal = apop_data_get(est->parameters, i, 2); 		
					coeff = gsl_vector_get(est->parameters->vector, i);
					/*
					printf(": finner en ok koeffisient for: %s: %.2f (t-verdi: %.2f)\n", est->parameters->names->row[i], coeff, tVal);
					fflush(stdout);
					*/
				}
				BBprintIndent(LevelDep);
				apop_model_show(est);
				fflush(stdout);
				BBprintIndent(LevelDep);
				apop_data_show(apop_estimate_r_squared(est));
				fflush(stdout);
			}
			else {
				BBprintIndent(LevelDep);
				printf("logLikelihood over 0: %.2f\n", est->llikelihood);
				fflush(stdout);
			}
	                if (sign != 2) {
printf("AAA\n"); fflush(stdout);
				BBgetUseIndep(Use, DataLevel->matrix->size2, &indepSize, &firstIndepPos);
printf("BBB\n"); fflush(stdout);
				if (indepSize == 1) {
printf("CCC\n"); fflush(stdout);
       		                 	BBcontrol(coeff, origT, DataLevel->names->column[0], origName);
printf("DDD\n"); fflush(stdout);
				}
                	}
		}
		apop_model_free(est);
		est = NULL;
	}
	if (dataClean != NULL) {
		apop_data_free(dataClean);
		dataClean = NULL;
	}
	if (dataCopy != NULL) {
		apop_data_free(dataCopy);
		dataCopy = NULL;
	}
	if (matrixLarge > 2) {
		useNext = calloc(DataLevel->matrix->size2 , sizeof(int));
		NoOfUseGl ++;
		/* 
			printf("NoOfUseGl: %d\n", NoOfUseGl);
			fflush(stdout); 
		*/
		for (j = J - 1; j > 0; j--) { 
			memset(useNext, 0, DataLevel->matrix->size2 * sizeof(int));
			BBcopyUse(Use, useNext, DataLevel->matrix->size2);
			useNext[j] = 1;
			// printf("v");
			BBregress(useNext, j, CombLevel, LevelDep, DataLevel, UseMax, MaxLL, MaxNoOfComb);
		}
		if (useNext != NULL) {
			free(useNext);
			useNext = NULL;
			NoOfUseGl --;
			/* 
				printf("NoOfUseGl: %d\n", NoOfUseGl);
				fflush(stdout); 
			*/
		}
	}
	CombLevel--;
	/*
		printf("................................................................BBregress ..end, level: %d, dependent var: %s\n", Level, DataLevel->names->column[0]);
		fflush(stdout);
	*/
	 // printf("^");
	/*
		printf("BBregress end\n");
		fflush(stdout); 
	*/
}


int BBlevel(int *Use, int J, apop_data *DataDep, int LevelDep, int LevelIndep, int **ReturnUseMax, apop_data **ReturnDataLevel)
{
	apop_data   *dataLevel = NULL, *thisDataLevel = NULL, *dataLevelCombined = NULL, *dataIndep = NULL, *dataIndepCombined = NULL, *newDataDep = NULL;
	int i, j, *levelUse = NULL, *dataDepUse = NULL, *useMax = NULL, *useMaxCombined = NULL, *thisUseMax = NULL, *useCombined = NULL, *levelUseCombined = NULL;
	char *levelStr = NULL;
	char nameCopy[200];
	double maxLogLikelihood = -1000.0;
	int maxNoOfComb = 0;
	int noOfIndep = 0;
	int pos = -1, inUse = 1;

	/*
		printf("BBlevel start\n");
		fflush(stdout);
	*/
	*ReturnUseMax = NULL;
	*ReturnDataLevel = NULL;
	useCombined = calloc(DataGl->matrix->size2 , sizeof(int));
	memset(useCombined, 0, DataGl->matrix->size2 * sizeof(int));
	BBinitUse(useCombined, DataGl->matrix->size2);
	/*
		printf("useCombined 0: ");
       		fflush(stdout);
		BBprintUse(useCombined, DataGl->matrix->size2, maxNoOfComb, 1);
	*/
	for (j = 0; j < J; j++) {
	       	strncpy(nameCopy, DataGl->names->column[j], 198);
	       	levelStr = strtok(nameCopy, "_");
		if ((levelStr != (char *)NULL) && (atoi(levelStr) == (LevelIndep + 1))) {
			/*
				printf("Found indep: %s\n", DataGl->names->column[j]);
				fflush(stdout);
			*/
			Use[j] = 0;
			noOfIndep++;
		}
		else {
			Use[j] = 1;
		}
        }
	if (noOfIndep > 0) {
		dataIndep = apop_data_copy(DataGl);
        	apop_data_rm_columns(dataIndep, Use);
        	dataLevel = apop_data_stack(DataDep, dataIndep, 'c', 'n');
		levelUse = calloc(dataLevel->matrix->size2 , sizeof(int));
		memset(levelUse, 0, dataLevel->matrix->size2 * sizeof(int));
		printf("levelUse: ");
       		fflush(stdout);
		BBprintUse(levelUse, dataLevel->matrix->size2, maxNoOfComb, 1);
		maxLogLikelihood = -1000.0;
		BBprintIndent(LevelDep);
               	printf("****************************************\n");
		BBprintIndent(LevelDep);
               	printf("Steg %d i denne stien som er avhengig variabel %s  mot kombinasjon av uavhengige variabeler på nivå %d\n", LevelDep, dataLevel->names->column[0], LevelIndep + 1);
		fflush(stdout);
		/* 
			BBregress(levelUse, dataLevel->matrix->size2, Level, dataLevel, &useMax, &maxLogLikelihood, &maxNoOfComb);
		*/
		BBregress(levelUse, dataLevel->matrix->size2, 0, LevelDep, dataLevel, &useMax, &maxLogLikelihood, &maxNoOfComb);
		*ReturnUseMax = useMax;
		*ReturnDataLevel = dataLevel;
		printf("useMax: ");
        	fflush(stdout);
		BBprintUse(useMax, dataLevel->matrix->size2, maxNoOfComb, 1);
		printf("ReturnUseMax: ");
        	fflush(stdout);
		BBprintUse(*ReturnUseMax, dataLevel->matrix->size2, maxNoOfComb, 1);
		if (useMax != (int *)NULL) {
			for (i = 1; i < dataLevel->matrix->size2; i++) {
				if (useMax[i] == 0) { 
					BBfindVariable(DataGl, dataLevel->names->column[i], &pos);
					if ((pos != -1) && (useFinishedGl[pos] != 0)) {
						useFinishedGl[pos] = 0;
        					newDataDep = apop_data_copy(dataLevel);
						dataDepUse = calloc(dataLevel->matrix->size2 , sizeof(int));
        					memset(dataDepUse, 1, dataLevel->matrix->size2 * sizeof(int));
        					dataDepUse[i] = 0;
        					apop_data_rm_columns(newDataDep, dataDepUse);
        					memset(Use, 1, DataGl->matrix->size2 * sizeof(int));
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 1, &thisUseMax, &thisDataLevel);
						/*
							printf("thisUseMax: ");
        						fflush(stdout);
							if (thisDataLevel != NULL) {
								BBprintUse(thisUseMax, thisDataLevel->matrix->size2, maxNoOfComb, 1);
							}
							else {
								printf("thisDataLevel is NULL");
							}
						*/
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 2, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 3, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 4, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 5, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 6, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 7, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 8, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 9, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBlevel(Use, J, newDataDep, LevelDep + 1, LevelIndep + 10, &thisUseMax, &thisDataLevel);
						BBmergeUse(thisUseMax, thisDataLevel, useCombined);
						BBcheckUse(useCombined, DataGl->matrix->size2, &inUse);
						/*
							printf("inUse: %d", inUse); 
							fflush(stdout);
						*/
						if (inUse == 0) {
							dataIndepCombined = apop_data_copy(DataGl);
        						apop_data_rm_columns(dataIndepCombined, useCombined);
        						dataLevelCombined = apop_data_stack(newDataDep, dataIndepCombined, 'c', 'n');
							levelUseCombined = calloc(dataLevelCombined->matrix->size2, sizeof(int));
							memset(levelUseCombined, 0, dataLevelCombined->matrix->size2 * sizeof(int));
							maxLogLikelihood = -1000.0;
							BBprintIndent(LevelDep);
               						printf("****************************************\n");
							BBprintIndent(LevelDep);
               						printf("Steg %d i denne stien som er avhengig variabel %s  mot kombinasjon av uavhengige variabeler med MAX LL på alle nivåer\n", LevelDep, dataLevelCombined->names->column[0]);
							fflush(stdout);
							fprintf(stderr, "useCombined:");
							fflush(stderr);
							BBprintUse(useCombined, DataGl->matrix->size2, maxNoOfComb, 1);
							/* 
								printf("dataIndepCombined size: %d", dataIndepCombined->matrix->size2); 
								fflush(stdout);
							*/
							BBregress(levelUseCombined, dataLevelCombined->matrix->size2, 0, LevelDep, dataLevelCombined, &useMaxCombined, &maxLogLikelihood, &maxNoOfComb);
							/*
								printf("useMax: ");
        							fflush(stdout);
								BBprintUse(useMax, dataLevel->matrix->size2, maxNoOfComb, 1);
							*/
						}
					}
					else if (pos == -1) {
						printf("*********************************\n");
						printf("ERROR: Finner ikke posisjon til avhengig variabel!!!!\n", DataGl->names->column[i]);
						printf("*********************************\n");
						fflush(stdout);
					}
					else if (useFinishedGl[pos] == 0) {
						printf("*********************************\n");
						printf("Sti fra %s er allerede undersøkt.\n", DataGl->names->column[pos]);
						printf("*********************************\n");
						fflush(stdout);
					}
				}
			}
		}
	}
	/*
		printf("BBlevel end\n");
		fflush(stdout); 
	*/
}

int BBfindDep(int J, apop_data *DataDep)
{
	int i, j;
	char nameCopy[200];
	int *dataDepUse = NULL;
	int foundDep = 0;

        printf(": ****************************************\n");
	printf("BBfindDep start\n");
	fflush(stdout); 
	dataDepUse = calloc(DataGl->matrix->size2 , sizeof(int));
        memset(dataDepUse, 1, DataGl->matrix->size2 * sizeof(int));
	for (j = 0; j < J; j++) {
	       	strncpy(nameCopy, DataGl->names->column[j], 198);
		if ((DataGl->names->column[j] != (char *)NULL) && (DataGl->names->column[j][0] == '0') && (foundDep == 0)) {
			dataDepUse[j] = 0;
			useFinishedGl[j] = 0;
			foundDep = 1;
               		printf(": Fant avhengig variabel: %s \n", DataGl->names->column[j]);
		}
        }
	if (foundDep == 0) {
        	dataDepUse[0] = 0;
		useFinishedGl[0] = 0;
               	printf(": Fant ikke avhengig variabel. Bruker den første: %s \n", DataGl->names->column[0]);
	}
        apop_data_rm_columns(DataDep, dataDepUse);
	printf("BBfindDep end\n");
	fflush(stdout); 
}



int main(void)
{
	apop_data       *data = NULL, *dataCopy = NULL, *dataDep = NULL, *dataIndep = NULL, *newData = NULL;
	apop_model   *est = NULL;
	apop_data       *test = NULL;
	FILE        *f;
	char        outfile[]   = "scatter.gplot";
	int i, j, *use = NULL;
	double tVal;
	apop_data   *thisDataLevel = NULL, *dataLevelCombined = NULL, *dataIndepCombined = NULL;
	int *thisUseMax = NULL, *useCombined = NULL, *levelUseCombined = NULL, *useMax = NULL, maxNoOfComb = 0;
	double maxLogLikelihood = -1000.0;
	int inUse = 1;

	gsl_set_error_handler_off();

	DataGl = apop_text_to_data("pathdata.csv");
       	apop_data_print(DataGl); 
	use = calloc(DataGl->matrix->size2 , sizeof(int));
	memset(use, 1, DataGl->matrix->size2 * sizeof(int));
	useFinishedGl = calloc(DataGl->matrix->size2 , sizeof(int));
	memset(useFinishedGl, 1, DataGl->matrix->size2 * sizeof(int));
	useCombined = calloc(DataGl->matrix->size2 , sizeof(int));
	memset(useCombined, 0, DataGl->matrix->size2 * sizeof(int));
	BBinitUse(useCombined, DataGl->matrix->size2);
	printf("useCombined main: ");
       	fflush(stdout);
	BBprintUse(useCombined, DataGl->matrix->size2, maxNoOfComb, 1);
        dataDep   = apop_data_copy(DataGl);
	BBfindDep(DataGl->matrix->size2, dataDep);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 0, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 1, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 2, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 3, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 4, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 5, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 6, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBlevel(use, DataGl->matrix->size2, dataDep, 0, 7, &thisUseMax, &thisDataLevel);
	BBmergeUse(thisUseMax, thisDataLevel, useCombined);
	BBcheckUse(useCombined, DataGl->matrix->size2, &inUse);
	printf("inUse: %d", inUse); 
	fflush(stdout);
	if (inUse == 0) {
		dataIndepCombined = apop_data_copy(DataGl);
       		apop_data_rm_columns(dataIndepCombined, useCombined);
		printf("useCombined: ");
       		fflush(stdout);
		BBprintUse(useCombined, DataGl->matrix->size2, maxNoOfComb, 1);
		/*
			printf("dataIndepCombined->matrix->size2: %d", dataIndepCombined->matrix->size2); 
			fflush(stdout);
		*/
       		dataLevelCombined = apop_data_stack(dataDep, dataIndepCombined, 'c', 'n');
		levelUseCombined = calloc(dataLevelCombined->matrix->size2, sizeof(int));
		memset(levelUseCombined, 0, dataLevelCombined->matrix->size2 * sizeof(int));
		maxLogLikelihood = -1000.0;
		BBprintIndent(0);
       		printf("****************************************\n");
		BBprintIndent(0);
        	printf("Steg %d i denne stien som er avhengig variabel %s  mot kombinasjon av uavhengige variabeler på alle nivåer\n", 0, dataDep->names->column[0]);
		fflush(stdout);
		BBregress(levelUseCombined, dataLevelCombined->matrix->size2, 0, 0, dataLevelCombined, &useMax, &maxLogLikelihood, &maxNoOfComb);
		/*
			printf("useMax: ");
        		fflush(stdout);
			BBprintUse(useMax, dataLevel->matrix->size2, maxNoOfComb, 1);
		*/
	}
	/*
       		dataCopy = apop_data_copy(data);
		for (i = 1; i < data->matrix->size2; i++) {
       			dataDep   = apop_data_copy(data);
       			dataIndep   = apop_data_copy(data);
			memset(use, 1, data->matrix->size2 * sizeof(int));
			use[i] = 0;
			apop_data_rm_columns(dataDep, use);
			memset(use, 0, data->matrix->size2 * sizeof(int));
			use[i] = 1;
			apop_data_rm_columns(dataIndep, use);
			newData = apop_data_stack(dataDep, dataIndep, 'c', 'n');
			BBregress(newData, newData->matrix->size2, 1);
			apop_data_free(newData);
			apop_data_free(dataDep);
			apop_data_free(dataIndep);
		}
	*/
	printf("FERDIG -----------------------  ");
	printf("Antall regresjoner: %d", NoOfEstGl);
	printf("   Antall kombinasjoner: %d (inkludert for lite datagrunnlag på grunn av tomme verdier.)\n", NoOfCombGl);
	fflush(stdout); 
	fflush(stdout); 
	fflush(stdout); 
	fflush(stdout); 


return 0;
}

