#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//DEFINES
#define SUCCESS 1
#define FAILURE -1
#define PROGRESS "===================="
#define PRG_WIDTH 20

//FUNCTIONS DECLARATION
int rows_number(char *);
double mean(double *, int);
void cumsum(double *, double *, int);
void slice_vec(double *, double *, int, int);
int lin_fit(int, const double *, const double *, double *, double *);
int get_num_elements(char *, char);
int get_scales_vec(char *, char *, int *, int);

//MAIN
int main(int argc, char **argv)
{
    //inputs
    int input_args = 4;
    if(argc < input_args){
        printf("Not enough input arguments!\n");
        return FAILURE;
    }else if(argc > input_args){
        printf("Too many input arguments!\n");
        return FAILURE;
    }
    //file with data (single-column only)
    char file_name[255];
    memset(file_name, 0x00, sizeof(file_name));
    sprintf(file_name, "%s", argv[1]);
    // polynomial order fixed to 1
    int pol_ord = 1;
    //scales are inputed as a string of ints separated by ","
    int num_elem = get_num_elements(argv[2], ',');
    if(num_elem == FAILURE){
        printf("ERROR: Scales separator found at the beginning or end of the string!\n");
        return FAILURE;
    }
    int *scales;
    scales = calloc(num_elem, sizeof(int));
    if(!scales){
        printf("MALLOC ERROR (scales)\n");
        return FAILURE;
    }
    int res = get_scales_vec(argv[2], ",", scales, pol_ord);
    if(res == FAILURE){
        printf("ERROR: Individual scales must be greater than %d!\n", pol_ord+2);
        free(scales);
        return FAILURE;
    }else{
        // path where to save results file
        char path_tot[255];
        memset(path_tot, 0x00, sizeof(path_tot));
        sprintf(path_tot, "%s", argv[3]);
        //file length
        int N = rows_number(file_name);
        if(N == FAILURE){
            printf("Cannot open input file.\n");
            return FAILURE;
        }
        //time series
        double *pn;
        pn = calloc(N, sizeof(double));
        if(!pn){
            printf("MALLOC ERROR (pn)\n");
            return FAILURE;
        }
        FILE *f;
        f = fopen(file_name, "r");
        if(!f){
            printf("Cannot open input file.\n");
            return FAILURE;
        }
        for(int i = 0; i < N; i++)
            fscanf(f, "%lf", pn+i);
        fclose(f);
        printf("Data length: <%d>\n", N);
        //time vector
        double *t;
        t = calloc(N, sizeof(double));
        if(!t){
            printf("MALLOC ERROR (t)\n");
            return FAILURE;
        }
        for(int i = 0; i < N; i++)
            t[i] = (double)(i+1);
        //time series minus its mean
        double a_ave = mean(pn, N);
        double *pn_nomean;
        pn_nomean = calloc(N, sizeof(double));
        if(!pn_nomean){
            printf("MALLOC ERROR (pn_nomean)\n");
            return FAILURE;
        }
        for(int i = 0; i < N; i++)
            pn_nomean[i] = pn[i] - a_ave;
        //cumulative sum
        double *X;
        X = calloc(N, sizeof(double));
        if(!X){
            printf("MALLOC ERROR (X)\n");
            return FAILURE;
        }
        printf("Data integration...");
        cumsum(pn_nomean, X, N);
        printf("done!\n");
        //defining parameters
        int scmin = 4;
        int max_win = 5;
        int scmax = N / max_win;
        int scres = 20;
        double step = (scmax - scmin) /(double) scres;
        int range0 = scres + 1;
        int *scale0;
        scale0 = calloc(range0, sizeof(int));
        if(!scale0){
            printf("MALLOC ERROR (scale0)");
            return FAILURE;
        }
        for(int i = 0; i < range0; i++)
            scale0[i] = (int)(scmin+i*step);
        //Fluctuations at q = 0
        printf("Computing fluctuations at q=0...");
        int num_seg, start_lim, end_lim;
	double ang_coeff, intercept;
        double *t_fit, *X_fit, *diff_vec, *RMS0, *Fq0, *fit_coeffs0;
        int flctLen = N / scmin;
        t_fit = calloc(scmax, sizeof(double));
        if(!t_fit){
            printf("MALLOC ERROR (t_fit)\n");
            return FAILURE;
        }
        X_fit = calloc(scmax, sizeof(double));
        if(!X_fit){
            printf("MALLOC ERROR (X_fit)\n");
            return FAILURE;
        }
        diff_vec = calloc(scmax, sizeof(double));
        if(!diff_vec){
            printf("MALLOC ERROR (diff_vec)\n");
            return FAILURE;
        }
        RMS0 = calloc(flctLen, sizeof(double));
        if(!RMS0){
            printf("MALLOC ERROR (RMS0)\n");
            return FAILURE;
        }
        Fq0 = calloc(range0, sizeof(double));
        if(!Fq0){
            printf("MALLOC ERROR (Fq0)\n");
            return FAILURE;
        }
        fit_coeffs0 = calloc(pol_ord+1, sizeof(double));
        if(!fit_coeffs0){
            printf("MALLOC ERROR (fit_coeffs0)\n");
            return FAILURE;
        }
        for(int i = 0; i < range0; i++){
            num_seg = N / scale0[i];
            for(int v = 0; v < num_seg; v++){
                start_lim = v * scale0[i];
                end_lim = (v + 1) * scale0[i];
                slice_vec(t, t_fit, start_lim, end_lim);
                slice_vec(X, X_fit, start_lim, end_lim);
                lin_fit(scale0[i], t_fit, X_fit, &ang_coeff, &intercept);
                for(int j = 0; j < scale0[i]; j++)
                    diff_vec[j] = pow(X_fit[j]-(intercept+ang_coeff*t_fit[j]), 2.0);
                RMS0[v] = log(mean(diff_vec, scale0[i]));
            }
            Fq0[i] = exp(0.5 * mean(RMS0, num_seg));
        }
        printf("done\n");
        //ht for every scale
        for(int i = 0; i < num_elem; i++){
            int scale = scales[i];
            int HtLen = N - scale + 1;
	    double ang_coeff, intercept;
            double *X_fitH, *t_fitH, *diff_vecH, *RMS, *fit_coeffs;
            X_fitH = calloc(scale, sizeof(double));
            if(!X_fitH){
                printf("MALLOC ERROR (X_fitH)\n");
                return FAILURE;
            }
            t_fitH = calloc(scale, sizeof(double));
            if(!t_fitH){
                printf("MALLOC ERROR (t_fitH)\n");
                return FAILURE;
            }
            diff_vecH = calloc(scale, sizeof(double));
            if(!diff_vecH){
                printf("MALLOC ERROR (diff_vecH)\n");
                return FAILURE;
            }
            RMS = calloc(HtLen, sizeof(double));
            if(!RMS){
                printf("MALLOC ERROR (RMS)\n");
                return FAILURE;
            }
            fit_coeffs = calloc(pol_ord+1, sizeof(double));
            if(!fit_coeffs){
                printf("MALLOC ERROR (fit_coeffs)\n");
                return FAILURE;
            }
            for(int v = 0; v <= N-scale; v++){
                double perc = v * 100 / (double)HtLen;
                int prg = (v * PRG_WIDTH) / HtLen;
                printf("Computing fluctuations at scale <%d> => [%.*s%*s] %.2lf%%\r", scale, prg, PROGRESS, PRG_WIDTH-prg, "", perc);
                fflush(stdout);
                start_lim = v;
                end_lim = v + scale - 1;
                for(int j = 0; j < scale; j++)
                    t_fitH[j] = (double)(start_lim+j);
                slice_vec(X, X_fitH, start_lim, end_lim);
                lin_fit(scale, t_fitH, X_fitH, &ang_coeff, &intercept);
                for(int j = 0; j < scale; j++)
                    diff_vecH[j] = pow(X_fitH[j]-(intercept+ang_coeff*t_fitH[j]), 2.0);
                RMS[v] = sqrt(mean(diff_vecH, scale));
            }
            printf("Computing fluctuations at scale <%d> => [%s] 100.00%%\r", scale, PROGRESS);
            fflush(stdout);
            printf("\n");
            //local hurst exponent
            printf("Computing local hurst exponent at scale <%d>...", scale);
            double *log_scale0, *log_Fq0, *Hq0_fit;
            log_scale0 = calloc(range0, sizeof(double));
            if(!log_scale0){
                printf("MALLOC ERROR (log_scale0)\n");
                return FAILURE;
            }
            log_Fq0 = calloc(range0, sizeof(double));
            if(!log_Fq0){
                printf("MALLOC ERROR (log_Fq0)\n");
                return FAILURE;
            }
            Hq0_fit = calloc(2, sizeof(double));
            if(!Hq0_fit){
                printf("MALLOC ERROR (Hq0_fit)\n");
                return FAILURE;
            }
            for(int i = 0; i < range0; i++){
                log_scale0[i] = log(scale0[i]);
                log_Fq0[i] = log(Fq0[i]);
            }
            polynomialFit(range0, 2, log_scale0, log_Fq0, Hq0_fit);
            double Regfit, logscale;
            double *resRMS, *Ht;
            double Hq0_intercept = Hq0_fit[0];
            double Hq0 = Hq0_fit[1];
            Regfit = Hq0_intercept + Hq0 * log(scale);
            logscale = log(HtLen) - log(scale);
            resRMS = calloc(HtLen, sizeof(double));
            if(!resRMS){
                printf("MALLOC ERROR (resRMS)\n");
                return FAILURE;
            }
            Ht = calloc(HtLen, sizeof(double));
            if(!Ht){
                printf("MALLOC ERROR (Ht)\n");
                return FAILURE;
            }
            for(int i = 0; i < HtLen; i++){
                resRMS[i] = Regfit - log(RMS[i]);
                Ht[i] = resRMS[i] / (double)logscale + Hq0;
            }
            printf("done\n");
            //local Hurst exponent file
            printf("Output file for scale <%d>...", scale);
            char path_file[300];
            memset(path_file, 0x00, sizeof(path_file));
            sprintf(path_file, "%s/Ht_scale=%d.txt", path_tot, scale);
            f = fopen(path_file, "w");
            if(!f){
                printf("Cannot open output file.\n");
                return FAILURE;
            }
            for(int i = 0; i < HtLen; i++)
                fprintf(f, "%lf\n", Ht[i]);
            fclose(f);
            printf("done\n");
            //free memory
            free(X_fitH); free(t_fitH); free(diff_vecH);
            free(RMS); free(fit_coeffs); free(log_scale0);
            free(log_Fq0); free(resRMS); free(Ht); free(Hq0_fit);
        }
        //free memory
        free(scales); free(pn); free(t);
        free(pn_nomean); free(X); free(scale0);
        free(t_fit); free(X_fit); free(diff_vec);
        free(RMS0); free(Fq0); free(fit_coeffs0);
    }
    return 0;
}

//FUNCTIONS
int rows_number(char *file_name)
{
    FILE *f;
    int stop;
    int lines = 0;
    f = fopen(file_name,"r");
    if(!f){
        printf("Cannot open file %s\n", file_name);
        return FAILURE;
    }
    while(!feof(f)){
        stop = fgetc(f);
        if(stop == '\n')
            lines++;
    }
    fclose(f);
    return lines;
}

double mean(double *vec, int L)
{
    double avg = 0.0;
    for(int i = 0; i < L; i++)
        avg += vec[i];
    avg /= (double)L;
    return avg;
}

void cumsum(double *vec, double *sum_vec, int L)
{
    sum_vec[0] = vec[0];
    for(int i = 1; i < L; i++)
        sum_vec[i] = sum_vec[i-1] + vec[i];
}

void slice_vec(double *all_vec, double *sliced_vec, int start, int end)
{
    for(int i = 0; i <= (end-start); i++)
        sliced_vec[i] = all_vec[start+i];
}

int lin_fit(int L, const double *x, const double *y, double *m, double *q)
{
    double sumx = 0.0;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    double sumy = 0.0;
    double sumy2 = 0.0;
    for(int i = 0; i < L; i++){
        sumx += x[i];
        sumx2 += x[i] * x[i];
        sumxy += x[i] * y[i];
        sumy += y[i];
        sumy2 += y[i] * y[i];
    }
    double denom = (L * sumx2 - sumx * sumx);
    if(denom == 0){
        *m = 0;
        *q = 0;
        return 1;
    }
    *m = (L * sumxy - sumx * sumy) / denom;
    *q = (sumy * sumx2 - sumx * sumxy) / denom;
    return 0;
}

int get_num_elements(char *str, char del)
{
	char c;
	int len, i;
	len = 0;
	i = 0;
	c = *str;
    if(c == del)
        return FAILURE;
	while(c != '\0'){
		if(c == del){
			len += 1;
            do{
                i++;
            }while(*(str+i) == del);
        }
        else
            i++;
		c = *(str+i);
	}
    if(*(str+i-1) == del)
        return FAILURE;
	return len+1;
}

int get_scales_vec(char *scales_str, char *del, int *scales_vec, int pol_deg)
{
    char *str_part;
	str_part = strtok(scales_str, del);
    int idx = 0;
	while(str_part != NULL){
        if(strlen(str_part) == 0)
            return FAILURE;
        int int_scale = atoi(str_part);
        if(int_scale < pol_deg+2)
            return FAILURE;
        else{
            scales_vec[idx] = int_scale;
    		idx++;
    		str_part = strtok(NULL, del);
        }
	}
    return SUCCESS;
}
