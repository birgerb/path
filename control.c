#include <apop.h>

static apop_data   *DataGl = NULL;
int NoOfEstGl = 0;
int NoOfUseGl = 0;

int BBcopyUse(int *Use, int *UseHere, int Size)
{
	int i;

	for (i = 0; i < Size; i++) {
		UseHere[i] = Use[i];
	}

}

int BBprintUse(int *Use, int Size)
{
	int i;

	for (i = 0; i < Size; i++) {
		if (Use[i] == 0) {
			printf("x");
		}
		else {
			printf(" ");
		}
	}
	printf(".\n");
}


int BBcalcFirstPart(int *Use, int OrigIndep, int Control)
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

	/*
	printf("---------------------------------------------\n");
	fflush(stdout);
	printf("NR. %d-%d\n", OrigIndep, Control);
	fflush(stdout);
	printf("NORSK: Inkluderer testvariabel. Regresjon av første part fra original uavhengig nr.: %d til kontrollvar nr.: %d: ", OrigIndep, Control);
	fflush(stdout);
	*/
       	dataDep   = apop_data_copy(DataGl);
       	dataIndep   = apop_data_copy(DataGl);
	memset(Use, 1, DataGl->matrix->size2 * sizeof(int));
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

int BBcontrol(int *Use, int OrigIndep, double OrigCoeff, double OrigT, char *OrigName)
{
	apop_data       *dataCopy = NULL;
	apop_data       *dataClean = NULL;
	apop_model   *est = NULL;
	int i, j, indep, control;
	double tVal, coeff;
	int sign = 1;
	int *useNext = NULL;
	int matrixLarge = 3;

	for (indep = 1; indep < DataGl->matrix->size2; indep++) {
		if (indep != OrigIndep) {
			/*
			printf("---------------------------------------------\n");
			fflush(stdout);
			printf("NR. %d-%d\n", OrigIndep, indep);
			fflush(stdout);
			*/
			memset(Use, 1, DataGl->matrix->size2 * sizeof(int));
			Use[0] = 0;
			Use[OrigIndep] = 0;
			Use[indep] = 0;
			/*
			BBprintUse(Use, DataGl->matrix->size2);
			printf("NORSK: Inkluderer testvariabel. Regresjon mellom avhengig variabel: %s og uavhengige variabel nr. %d: %s med testvariabel nr. %d: ", DataGl->names->column[0], OrigIndep, OrigName, indep);
			fflush(stdout);
			*/
       			dataCopy   = apop_data_copy(DataGl);
			apop_data_rm_columns(dataCopy, Use); 
			/*
			for (j = 1; j < dataCopy->matrix->size2; j++) {
				if (strcmp(OrigName, dataCopy->names->column[j]) != 0) {
					printf("%s, ", dataCopy->names->column[j]);
				}
			}
			printf("\n");
			fflush(stdout);
			*/
        		dataClean   = apop_data_listwise_delete(dataCopy); 
			if ((dataClean == NULL) || (dataClean->matrix->size2 > 2)) {
				matrixLarge = 3;
			}
			else {
				matrixLarge = 2;
			}
			if ((dataClean != NULL) && (dataClean->matrix->size1 > 30)) {
				est  = apop_estimate(dataClean, apop_ols); 
				sign = 1;
				control = 0;
				for (i=1; i < est->parameters->matrix->size1; i++) { 
					tVal = apop_data_get(est->parameters, i, 2); 		
					coeff = gsl_vector_get(est->parameters->vector, i); 		
					if ((tVal < -1000.0) || (tVal > 1000.0)) {
						sign = 2;
						/*
						printf("NORSK: finner en t-verdi som er ekstremt høy for: %s. Mest sannsynlig singular matrise slik at OLS har brutt sammen: %f\n", est->parameters->names->row[i], tVal); 
						fflush(stdout);
						*/
					}
					else if ((tVal < 2.0) && (tVal > -2.0)) {
						sign = 0;
						/*
						printf("NORSK: finner en ikke ok t-verdi for: %s: %f\n", est->parameters->names->row[i], tVal); 
						fflush(stdout);
						*/
					}
					else {
						/*
						printf("NORSK: finner ok t-verdi for: %s: %f\n", est->parameters->names->row[i], tVal);  
						fflush(stdout);
						*/
					}
					if ((strcmp(OrigName, est->parameters->names->row[i]) == 0) && (OrigT >= 2.0) && (tVal < 0.5)) { 
						control = 1;
					}
					if ((strcmp(OrigName, est->parameters->names->row[i]) == 0) && (OrigT <= -2.0) && (tVal > -0.5)) { 
						control = 1;
					}
					if ((strcmp(OrigName, est->parameters->names->row[i]) == 0) && (OrigT <= 0.5) && (tVal >= 2.0)) { 
						control = 2;
					}
					if ((strcmp(OrigName, est->parameters->names->row[i]) == 0) && (OrigT >= -0.5) && (tVal <= -2.0)) { 
						control = 2;
					}
				}
				if ((sign == 0) && (control == 1))  {
					printf("+++++++++ Forsterker ++++++++++\n"); 
					printf("NR. %d-%d\n", OrigIndep, indep);
					printf("NORSK: Opprinnelig modell: %s ---> %.2f (%.2f) ---> %s\n", OrigName, OrigCoeff, OrigT, DataGl->names->column[0]); 
					printf("NORSK: testvariabelen sender t for: %s under 0,5 og forklarer dermed opprinnelig sammenheng.\n", OrigName); 
					for (i=1; i < est->parameters->matrix->size1; i++) { 
						if (strcmp(OrigName, est->parameters->names->row[i]) != 0) {
							tVal = apop_data_get(est->parameters, i, 2); 		
							coeff = gsl_vector_get(est->parameters->vector, i); 		
							printf("NORSK:                 %s \n", est->parameters->names->row[i]); 
							printf("NORSK:                /        \\\n"); 
							BBcalcFirstPart(Use, OrigIndep, indep);
							printf("NORSK:         /                %.2f (%.2f) \n", coeff, tVal); 
							printf("NORSK:        /                     \\\n"); 
						}
						fflush(stdout);
					}
					for (i=1; i < est->parameters->matrix->size1; i++) { 
						if (strcmp(OrigName, est->parameters->names->row[i]) == 0) {
							tVal = apop_data_get(est->parameters, i, 2); 		
							coeff = gsl_vector_get(est->parameters->vector, i); 		
							printf("NORSK: %s ----> %.2f (%.2f) ----> %s\n", est->parameters->names->row[i], coeff, tVal, DataGl->names->column[0]); 
						}
					}
					fflush(stdout);
				}
				if ((control == 2))  {
					printf("--------- Bremser ----------\n"); 
					printf("NR. %d-%d\n", OrigIndep, indep);
					printf("NORSK: Opprinnelig modell: %s ---> %.2f (%.2f) ---> %s\n", OrigName, OrigCoeff, OrigT, DataGl->names->column[0]); 
					printf("NORSK: testvariabelen øker: %s over 0,5 og forklarer hvorfor det ikke var opprinnelig sammenheng.\n", OrigName); 
					for (i=1; i < est->parameters->matrix->size1; i++) { 
						if (strcmp(OrigName, est->parameters->names->row[i]) != 0) {
							tVal = apop_data_get(est->parameters, i, 2); 		
							coeff = gsl_vector_get(est->parameters->vector, i); 		
							printf("NORSK:                 %s \n", est->parameters->names->row[i]); 
							printf("NORSK:                /        \\\n"); 
							BBcalcFirstPart(Use, OrigIndep, indep);
							printf("NORSK:         /                %.2f (%.2f) \n", coeff, tVal); 
							printf("NORSK:        /                     \\\n"); 
						}
						fflush(stdout);
					}
					for (i=1; i < est->parameters->matrix->size1; i++) { 
						if (strcmp(OrigName, est->parameters->names->row[i]) == 0) {
							tVal = apop_data_get(est->parameters, i, 2); 		
							coeff = gsl_vector_get(est->parameters->vector, i); 		
							printf("NORSK: %s ----> %.2f (%.2f) ----> %s\n", est->parameters->names->row[i], coeff, tVal, DataGl->names->column[0]); 
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
					printf("NR. %d-%d\n", OrigIndep, indep);
					fflush(stdout);
					printf("NORSK: testvariabelen sender ikke: %s under 0,5 men den nye modellen er signifikant\n", OrigName); 
					apop_model_show(est);
					fflush(stdout);
				}
				apop_model_free(est);
				est = NULL;
			}
			else {
				/*
				printf("NORSK: Uavhengig varabel nr.: %d gir liten eller tom matrise\n", indep);
				*/
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
	printf("NORSK: slutten av regresjon for avhengig variabel: %s\n", DataGl->names->column[0]);
	fflush(stdout);
	
}

int BBtest(int *Use)
{
	apop_data       *dataCopy = NULL;
	apop_data       *dataClean = NULL;
	apop_model   *est = NULL;
	int i, j, indep;
	double tVal, origT, coeff;
	char origName[200];
	int sign = 1;
	int *useNext = NULL;
	int matrixLarge = 3;

	printf("A\n");
	fflush(stdout);
	for (indep = 1; indep < DataGl->matrix->size2; indep++) {
		origName[0] = '\0';
		printf("---------------------------------------------\n");
		fflush(stdout);
		printf("NR. %d\n", indep);
		fflush(stdout);
		memset(Use, 1, DataGl->matrix->size2 * sizeof(int));
		Use[0] = 0;
		Use[indep] = 0;
		BBprintUse(Use, DataGl->matrix->size2);
		printf("NORSK: Regresjon mellom avhengig variabel: %s og uavhengig variabel nr. %d: ", DataGl->names->column[0], indep);
		fflush(stdout);
       		dataCopy   = apop_data_copy(DataGl);
		apop_data_rm_columns(dataCopy, Use); 
		for (j = 1; j < dataCopy->matrix->size2; j++) {
			printf("%s, ", dataCopy->names->column[j]);
		}
		printf("\n");
		fflush(stdout);
        	dataClean   = apop_data_listwise_delete(dataCopy); 
		if ((dataClean == NULL) || (dataClean->matrix->size2 > 2)) {
			matrixLarge = 3;
		}
		else {
			matrixLarge = 2;
		}
		if ((dataClean != NULL) && (dataClean->matrix->size1 > 30)) {
			est  = apop_estimate(dataClean, apop_ols); 
			sign = 1;
			for (i=1; i < est->parameters->matrix->size1; i++) { 
				tVal = apop_data_get(est->parameters, i, 2); 		
				coeff = gsl_vector_get(est->parameters->vector, i); 		
				origT = tVal;
				strcpy(origName, est->parameters->names->row[1]); 
				printf("origName: %s\n", origName);
				fflush(stdout);
				if ((tVal < -1000.0) || (tVal > 1000.0)) {
					sign = 2;
					printf("NORSK: finner en t-verdi som er ekstremt høy for: %s. Mest sannsynlig singular matrise slik at OLS har brutt sammen: %f\n", est->parameters->names->row[i], tVal); 
					fflush(stdout);
				}
				else if ((tVal < 2.0) && (tVal > -2.0)) {
					sign = 0;
					printf("NORSK: finner en ikke ok koeffisient for: %s: %.2f (t-verdi: %.2f)\n", est->parameters->names->row[i], coeff, tVal); 
					fflush(stdout);
				}
				else {
					printf("NORSK: finner en ok koeffisient for: %s: %.2f (t-verdi: %.2f)\n", est->parameters->names->row[i], coeff, tVal); 
					fflush(stdout);
				}
			}
			if (sign == 1) {
				printf("\n");
				printf("NORSK: Signifikant modell\n");
				apop_model_show(est);
				fflush(stdout);
			}
			apop_model_free(est);
			est = NULL;
		}
		else {
			printf("NORSK: Uavhengig varabel nr.: %d gir liten eller tom matrise\n", indep);
		}
		if (sign != 2) {
			BBcontrol(Use, indep, coeff, origT, origName);
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
	printf("NORSK: slutten av regresjon for avhengig variabel: %s\n", DataGl->names->column[0]);
	fflush(stdout);
}


int main(void){
	apop_data       *data = NULL, *dataCopy = NULL, *dataDep = NULL, *dataIndep = NULL, *newData = NULL;
	apop_model   *est = NULL;
	apop_data       *test = NULL;
	FILE        *f;
	char        outfile[]   = "scatter.gplot";
	int i, j, *use = NULL, *use2 = NULL;
	double tVal;

	gsl_set_error_handler_off();

	DataGl = apop_text_to_data("data");
       	apop_data_print(DataGl); 

	use = calloc(DataGl->matrix->size2 , sizeof(int));
	memset(use, 1, DataGl->matrix->size2 * sizeof(int));
	printf("O\n");
	fflush(stdout);
	BBtest(use);

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


return 0;
}

