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
int lin_fit(int, const double *, const double *, double *, double *, double *);

//MAIN
int main(int argc, char **argv)
{
    //inputs
    int input_args = 6;
    if(argc < input_args){
        printf("Not enough input arguments!\n");
        return FAILURE;
    }else if(argc > input_args){
        printf("Too many input arguments!\n");
        return FAILURE;
    }
    char file_name[255];
    memset(file_name, 0x00, sizeof(file_name));
    sprintf(file_name, "%s", argv[1]);
    int min_win = atoi(argv[2]);
    int ord = atoi(argv[3]);
    int rev_seg = atoi(argv[4]);
    char path_tot[255];
    memset(path_tot, 0x00, sizeof(path_tot));
    sprintf(path_tot, "%s/dfa_mw=%d_pol=%d_bck=%d.txt", argv[5], min_win, ord, rev_seg);
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
    double *y;
    y = calloc(N, sizeof(double));
    if(!y){
        printf("MALLOC ERROR (y)\n");
        return FAILURE;
    }
    printf("Data integration...");
    cumsum(pn_nomean, y, N);
    printf("done!\n");
    //defining parameters
    int max_win = 5;
    int end_dfa = N / max_win;
    int *s;
    int range_dfa = end_dfa - min_win + 1;
    s = calloc(range_dfa, sizeof(int));
    if(!s){
        printf("MALLOC ERROR (s)\n");
        return FAILURE;
    }
    for(int i = 0; i < range_dfa; i++)
        s[i] = i + min_win;
    //fluctuations vector and other arrays
    double *F, *F_nu1, *F_nu2, *t_fit, *y_fit, *diff_vec;
    int F_len = N / min_win;
    F = calloc(range_dfa, sizeof(double));
    if(!F){
        printf("MALLOC ERROR (F)\n");
        return FAILURE;
    }
    F_nu1 = calloc(F_len, sizeof(double));
    if(!F_nu1){
        printf("MALLOC ERROR (F_nu1)\n");
        return FAILURE;
    }
    F_nu2 = calloc(F_len, sizeof(double));
    if(!F_nu2){
        printf("MALLOC ERROR (F_nu2)\n");
        return FAILURE;
    }
    t_fit = calloc(end_dfa, sizeof(double));
    if(!t_fit){
        printf("MALLOC ERROR (t_fit)\n");
        return FAILURE;
    }
    y_fit = calloc(end_dfa, sizeof(double));
    if(!y_fit){
        printf("MALLOC ERROR (y_fit)\n");
        return FAILURE;
    }
    diff_vec = calloc(end_dfa, sizeof(double));
    if(!diff_vec){
        printf("MALLOC ERROR (diff_vec)\n");
        return FAILURE;
    }
    //computation
    int start_lim, end_lim;
    double ang_coeff, intercept, r_coeff;
    for(int i = 0; i < range_dfa; i++){
        int N_s = N / s[i];
        double perc = i * 100 / (double)range_dfa;
        int prg = (i * PRG_WIDTH) / range_dfa;
        printf("Computing fluctuations => [%.*s%*s] %.2lf%%\r", prg, PROGRESS, PRG_WIDTH-prg, "", perc);
        fflush(stdout);
        for(int v = 0; v < N_s; v++){
            start_lim = v * s[i];
            end_lim = (v + 1) * s[i] - 1;
            slice_vec(t, t_fit, start_lim, end_lim);
            slice_vec(y, y_fit, start_lim, end_lim);
            lin_fit(s[i], t_fit, y_fit, &ang_coeff, &intercept, &r_coeff);
            for(int j = 0; j < s[i]; j++)
                diff_vec[j] = pow((y_fit[j] - (intercept + ang_coeff * t_fit[j])), 2.0);
            F_nu1[v] = mean(diff_vec, s[i]);
        }
        if(rev_seg == 1){
            for(int v = 0; v < N_s; v++){
                start_lim = v * s[i] + (N - N_s * s[i]);
                end_lim = (v + 1) * s[i] + (N - N_s * s[i]);
                slice_vec(t, t_fit, start_lim, end_lim);
                slice_vec(y, y_fit, start_lim, end_lim);
                lin_fit(s[i], t_fit, y_fit, &ang_coeff, &intercept, &r_coeff);
                for(int j = 0; j < s[i]; j++)
                    diff_vec[j] = pow((y_fit[j] - (intercept + ang_coeff * t_fit[j])), 2.0);
                F_nu2[v] = mean(diff_vec, s[i]);
            }
            F[i] = sqrt((mean(F_nu1, N_s) + mean(F_nu2, N_s)) / (double)2);
        }else{
            F[i] = sqrt(mean(F_nu1, N_s));
        }
    }
    printf("Computing fluctuations => [%s] 100.00%%\r", PROGRESS);
    fflush(stdout);
    printf("\n");
    free(F_nu1); free(F_nu2);
    free(t_fit); free(y_fit);
    free(diff_vec);
    //HURST EXPONENT
    printf("Computing hurst exponent...");
    double *log_s, *log_F, H, H_intercept, H_rcoeff;
    log_s = calloc(range_dfa, sizeof(double));
	if(!log_s){
        printf("MALLOC ERROR (log_s)\n");
        return FAILURE;
    }
    log_F = calloc(range_dfa, sizeof(double));
    if(!log_F){
        printf("MALLOC ERROR (log_F)\n");
        return FAILURE;
    }
    for(int i = 0; i < range_dfa; i++){
        log_s[i] = log(s[i]);
        log_F[i] = log(F[i]);
    }
    lin_fit(range_dfa, log_s, log_F, &H, &H_intercept, &H_rcoeff);
    printf("done\n");
    printf("Output file...");
    f = fopen(path_tot, "w");
    if(!f){
        printf("Cannot open output file.\n");
        return FAILURE;
    }
    for(int i = 0; i < range_dfa; i++)
        fprintf(f, "%d %lf %lf %lf\n", s[i], F[i], H_intercept+log_s[i]*H, H);
    fclose(f);
    printf("done\n");
    free(F); free(log_s); free(log_F);
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

int lin_fit(int L, const double *x, const double *y, double *m, double *q, double *r)
{
    double sumx = 0.0;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    double sumy = 0.0;
    double sumy2 = 0.0;
    for(int i = 0;i < L; i++){
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
        if(r)
            *r = 0;
        return 1;
    }
    *m = (L * sumxy - sumx * sumy) / denom;
    *q = (sumy * sumx2 - sumx * sumxy) / denom;
    if(r != NULL){
        *r = (sumxy - sumx * sumy / (double)L) /
              sqrt((sumx2 - (sumx * sumx) / L) *
              (sumy2 - (sumy * sumy) / L));
    }
    return 0;
}
