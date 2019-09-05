#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_multifit.h>

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
void polynomialFit(int obs, int degree, double *dx, double *dy, double *store);

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
    double *fit_coeffs;
    fit_coeffs = calloc(ord+1, sizeof(double));
    if(!fit_coeffs){
        printf("MALLOC ERROR (fit_coeffs)\n");
        return FAILURE;
    }
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
            polynomialFit(s[i], ord+1, t_fit, y_fit, fit_coeffs);
            for(int j = 0; j < s[i]; j++){
                for(int k = 0; k < ord+1; k++)
                    y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
                diff_vec[j] = pow(y_fit[j], 2.0);
            }
            F_nu1[v] = mean(diff_vec, s[i]);
        }
        if(rev_seg == 1){
            for(int v = 0; v < N_s; v++){
                start_lim = v * s[i] + (N - N_s * s[i]);
                end_lim = (v + 1) * s[i] + (N - N_s * s[i]);
                slice_vec(t, t_fit, start_lim, end_lim);
                slice_vec(y, y_fit, start_lim, end_lim);
                polynomialFit(s[i], ord+1, t_fit, y_fit, fit_coeffs);
                for(int j = 0; j < s[i]; j++){
                    for(int k = 0; k < ord+1; k++)
                        y_fit[j] -= fit_coeffs[k] * pow(t_fit[j], k);
                    diff_vec[j] = pow(y_fit[j], 2.0);
                }
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
    free(diff_vec); free(fit_coeffs);
    //HURST EXPONENT
    printf("Computing hurst exponent...");
    double *log_s, *log_F, *H_fit;
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
    H_fit = calloc(2, sizeof(double));
    if(!H_fit){
        printf("MALLOC ERROR (H_fit)\n");
        return FAILURE;
    }
    for(int i = 0; i < range_dfa; i++){
        log_s[i] = log(s[i]);
        log_F[i] = log(F[i]);
    }
    polynomialFit(range_dfa, 2, log_s, log_F, H_fit);
    printf("done\n");
    printf("Output file...");
    f = fopen(path_tot, "w");
    if(!f){
        printf("Cannot open output file.\n");
        return FAILURE;
    }
    for(int i = 0; i < range_dfa; i++)
        fprintf(f, "%d %lf %lf %lf\n", s[i], F[i], H_fit[0]+log_s[i]*H_fit[1], H_fit[1]);
    fclose(f);
    printf("done\n");
    free(F); free(log_s); free(log_F); free(H_fit);
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

void polynomialFit(int obs, int degree, double *dx, double *dy, double *store)
{
    gsl_multifit_linear_workspace *ws;
    gsl_matrix *cov, *X;
    gsl_vector *y, *c;
    double chisq;
    int i, j;
    //alloc
    X = gsl_matrix_alloc(obs, degree);
    y = gsl_vector_alloc(obs);
    c = gsl_vector_alloc(degree);
    cov = gsl_matrix_alloc(degree, degree);
    //computation
    for(i = 0; i < obs; i++){
        for(j = 0; j < degree; j++)
            gsl_matrix_set(X, i, j, pow(dx[i], j));
        gsl_vector_set(y, i, dy[i]);
    }
    ws = gsl_multifit_linear_alloc(obs, degree);
    gsl_multifit_linear(X, y, c, cov, &chisq, ws);
    for(i = 0; i < degree; i++)
        store[i] = gsl_vector_get(c, i);
    //free memory
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);
}
