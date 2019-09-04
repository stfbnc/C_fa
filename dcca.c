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
int dcca_comp(double *, double *, double *, int, int, int, int, double *, int);
int rows_number(char *);
double mean(double *, int, int);
void cumsum(double *, double *, int);
void slice_vec(double *, double *, int, int);
int lin_fit(int, const double *, const double *, double *, double *, double *);

//MAIN
int main(int argc, char **argv)
{
    //inputs
    int input_args = 7;
    if(argc < input_args){
        printf("Not enough input arguments!\n");
        return FAILURE;
    }else if(argc > input_args){
        printf("Too many input arguments!\n");
        return FAILURE;
    }
    char file_name1[255]; //path+name for the first file
    memset(file_name1, 0x00, sizeof(file_name1));
    sprintf(file_name1, "%s", argv[1]);
    char file_name2[255]; //path+name for the second file
    memset(file_name2, 0x00, sizeof(file_name2));
    sprintf(file_name2, "%s", argv[2]);
    int smallest_win_size = atoi(argv[3]); //minimum size for the sliding window
    int biggest_win_num = atoi(argv[4]); //number of windows in the last step
    int pol_ord = atoi(argv[5]); //polynomial order for fit
    char path_out[255]; //path+name for the output file
    memset(path_out, 0x00, sizeof(path_out));
    sprintf(path_out, "%s/rho_mw=%d_Mw=%d_pol=%d.txt", argv[6], smallest_win_size, biggest_win_num, pol_ord);
    //files length
    int N1 = rows_number(file_name1);
    if(N1 == FAILURE){
        printf("Cannot open input file (1).\n");
        return FAILURE;
    }
    int N2 = rows_number(file_name2);
    if(N2 == FAILURE){
        printf("Cannot open input file (2).\n");
        return FAILURE;
    }
    //time series
    double *ts1, *ts2;
    ts1 = calloc(N1, sizeof(double));
    if(!ts1){
        printf("MALLOC ERROR (ts1)\n");
        return FAILURE;
    }
    ts2 = calloc(N2, sizeof(double));
    if(!ts2){
        printf("MALLOC ERROR (ts2)\n");
        return FAILURE;
    }
    FILE *f1, *f2;
    f1 = fopen(file_name1, "r");
    if(!f1){
        printf("Cannot open input file (1).\n");
        return FAILURE;
    }
    for(int i = 0; i < N1; i++)
        fscanf(f1, "%lf", ts1+i);
    fclose(f1);
    f2 = fopen(file_name2, "r");
    if(!f2){
        printf("Cannot open input file (2).\n");
        return FAILURE;
    }
    for(int i = 0; i < N2; i++)
        fscanf(f2, "%lf", ts2+i);
    fclose(f2);
    //time series must have the same length
    int L = N1;
    if(N1 < N2){
        printf("Series 2 is longer than series 1. Data length reduced to: <%d>\n", L);
        slice_vec(ts2, ts2, 0, N1-1);
    }else if(N1 > N2){
        slice_vec(ts1, ts1, 0, N2-1);
        L = N2;
        printf("Series 1 is longer than series 2. Data length reduced to: <%d>\n", L);
    }else{
        printf("Data length: <%d>\n", L);
    }
    //time vector
    double *time;
    time = calloc(L, sizeof(double));
    if(!time){
        printf("MALLOC ERROR (time)\n");
        return FAILURE;
    }
    for(int i = 0; i < L; i++)
        time[i] = (double)(i+1);
    //time series minus its mean
    double ts1_ave = mean(ts1, L, L);
    double ts2_ave = mean(ts2, L, L);
    double *ts1_nomean, *ts2_nomean;
    ts1_nomean = calloc(L, sizeof(double));
    if(!ts1_nomean){
        printf("MALLOC ERROR (ts1_nomean)\n");
        return FAILURE;
    }
    ts2_nomean = calloc(L, sizeof(double));
    if(!ts2_nomean){
        printf("MALLOC ERROR (ts2_nomean)\n");
        return FAILURE;
    }
    for(int i = 0; i < L; i++){
        ts1_nomean[i] = ts1[i] - ts1_ave;
        ts2_nomean[i] = ts2[i] - ts2_ave;
    }
    //cumulative sum
    double *ts1_cum, *ts2_cum;
    ts1_cum = calloc(L, sizeof(double));
    if(!ts1_cum){
        printf("MALLOC ERROR (ts1_cum)\n");
        return FAILURE;
    }
    ts2_cum = calloc(L, sizeof(double));
    if(!ts2_cum){
        printf("MALLOC ERROR (ts2_cum)\n");
        return FAILURE;
    }
    printf("Data integration...");
    cumsum(ts1_nomean, ts1_cum, L);
    cumsum(ts2_nomean, ts2_cum, L);
    printf("done!\n");
    //dcca computation
    int biggest_win_size = L / biggest_win_num;
    int win_range = biggest_win_size - smallest_win_size + 1;
    double *F_DCCA, *F_DFA_x, *F_DFA_y, *rho;
    F_DCCA = calloc(win_range, sizeof(double));
    F_DFA_x = calloc(win_range, sizeof(double));
    F_DFA_y = calloc(win_range, sizeof(double));
    rho = calloc(win_range, sizeof(double));
    printf("Computing DCCA between series 1 and series 2...\n");
    int ret = dcca_comp(time, ts1_cum, ts2_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_DCCA, win_range);
    if(ret == FAILURE)
        return FAILURE;
    printf("Computing DCCA between series 1 and series 1...\n");
    ret = dcca_comp(time, ts1_cum, ts1_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_DFA_x, win_range);
    if(ret == FAILURE)
        return FAILURE;
    printf("Computing DCCA between series 2 and series 2...\n");
    ret = dcca_comp(time, ts2_cum, ts2_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_DFA_y, win_range);
    if(ret == FAILURE)
        return FAILURE;
    FILE *fOut;
    fOut = fopen(path_out, "w");
    if(!fOut){
        printf("Cannot open output file.\n");
        return FAILURE;
    }
    printf("Computing rho...");
    for(int i = 0; i < win_range; i++){
        F_DFA_x[i] = sqrt(F_DFA_x[i]);
        F_DFA_y[i] = sqrt(F_DFA_y[i]);
        rho[i] = F_DCCA[i] / (double)(F_DFA_x[i] * F_DFA_y[i]);
        fprintf(fOut, "%lf %lf %lf %lf\n", F_DCCA[i], F_DFA_x[i], F_DFA_y[i], rho[i]);
    }
    printf("done\n");
    printf("Output file...");
    fclose(fOut);
    printf("done\n");
    //free memory
    free(ts1); free(ts2); free(time);
    free(ts1_nomean); free(ts2_nomean);
    free(ts1_cum); free(ts2_cum);
    free(F_DCCA); free(F_DFA_x); free(F_DFA_y); free(rho);
    return 0;
}

//FUNCTIONS
int dcca_comp(double *t, double *y1, double *y2, int N, int min_win_size, int last_win_len, int ord, double *F, int range_dcca)
{
    //fluctuations vector and other arrays
    int F_len = N - min_win_size;
    double *F_nu, *t_fit, *y_fit1, *y_fit2, *diff_vec;
    F_nu = calloc(F_len, sizeof(double));
    if(!F_nu){
        printf("MALLOC ERROR (F_nu)\n");
        return FAILURE;
    }
    t_fit = calloc((last_win_len+1), sizeof(double));
    if(!t_fit){
        printf("MALLOC ERROR (t_fit)\n");
        return FAILURE;
    }
    y_fit1 = calloc((last_win_len+1), sizeof(double));
    if(!y_fit1){
        printf("MALLOC ERROR (y_fit1)\n");
        return FAILURE;
    }
    y_fit2 = calloc((last_win_len+1), sizeof(double));
    if(!y_fit2){
        printf("MALLOC ERROR (y_fit2)\n");
        return FAILURE;
    }
    diff_vec = calloc((last_win_len+1), sizeof(double));
    if(!diff_vec){
        printf("MALLOC ERROR (diff_vec)\n");
        return FAILURE;
    }
    //computation
    int s, N_s, start_lim, end_lim;
    double ang_coeff1, intercept1, r_coeff1, ang_coeff2, intercept2, r_coeff2;
    for(int i = 0; i < range_dcca; i++){
        double perc = i * 100 / (double)range_dcca;
        int prg = (i * PRG_WIDTH) / range_dcca;
        printf(" [%.*s%*s] %.2lf%%\r", prg, PROGRESS, PRG_WIDTH-prg, "", perc);
        fflush(stdout);
        s = i + min_win_size;
        N_s = N - s;
        for(int v = 0; v < N_s; v++){
            start_lim = v;
            end_lim = v + s;
            slice_vec(t, t_fit, start_lim, end_lim);
            slice_vec(y1, y_fit1, start_lim, end_lim);
            slice_vec(y2, y_fit2, start_lim, end_lim);
            lin_fit(s+1, t_fit, y_fit1, &ang_coeff1, &intercept1, &r_coeff1);
            lin_fit(s+1, t_fit, y_fit2, &ang_coeff2, &intercept2, &r_coeff2);
            for(int j = 0; j <= s; j++)
                diff_vec[j] = (y_fit1[j] - (intercept1 + ang_coeff1 * t_fit[j])) * (y_fit2[j] - (intercept2 + ang_coeff2 * t_fit[j]));
            F_nu[v] = mean(diff_vec, s+1, s-1);
        }
        F[i] = mean(F_nu, N_s, N_s);
    }
    printf(" [%s] 100.00%%\r", PROGRESS);
    fflush(stdout);
    printf("\n");
    free(F_nu); free(t_fit); free(y_fit1); free(y_fit2); free(diff_vec);
    return SUCCESS;
}

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

double mean(double *vec, int vecL, int L)
{
    double avg = 0.0;
    for(int i = 0; i < vecL; i++)
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
    for(int i = 0; i < L; i++){
        sumx += *(x + i);
        sumx2 += *(x + i) * *(x + i);
        sumxy += *(x + i) * *(y + i);
        sumy += *(y + i);
        sumy2 += *(y + i) * *(y + i);
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
