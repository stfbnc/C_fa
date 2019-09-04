#include "includes.h"

void dcca(double *ts1, int N1, double *ts2, int N2, int smallest_win_size, int biggest_win_num, int pol_ord, char *path_out, int rank, int size)
{
    //time series must have the same length
    int L = N1;
    if(N1 < N2){
        slice_vec(ts2, ts2, 0, N1-1);
    }else if(N1 > N2){
        slice_vec(ts1, ts1, 0, N2-1);
        L = N2;
    }
    //time vector
    double *time;
    time = calloc(L, sizeof(double));
    if(!time){
        printf("MALLOC ERROR (time)\n");
        return 99;
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
        return 99;
    }
    ts2_nomean = calloc(L, sizeof(double));
    if(!ts2_nomean){
        printf("MALLOC ERROR (ts2_nomean)\n");
        return 99;
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
        return 99;
    }
    ts2_cum = calloc(L, sizeof(double));
    if(!ts2_cum){
        printf("MALLOC ERROR (ts2_cum)\n");
        return 99;
    }
    cumsum(ts1_nomean, ts1_cum, L);
    cumsum(ts2_nomean, ts2_cum, L);
    //dcca computation
    int biggest_win_size = L / biggest_win_num;
    int win_range = biggest_win_size - smallest_win_size + 1;
    int F_proc_size = win_range / size;
    int remd_size = win_range % size;
    double *F_proc;
    F_proc = calloc(F_proc_size, sizeof(double));
    if(!F_proc){
        printf("MALLOC ERROR (F_proc)\n");
        return 99;
    }
    double *F_DCCA, *F_DFA_x, *F_DFA_y, *rho;
    if(rank == 0){
        F_DCCA = calloc(win_range, sizeof(double));
        if(!F_DCCA){
            printf("MALLOC ERROR (F_DCCA)\n");
            return 99;
        }
        F_DFA_x = calloc(win_range, sizeof(double));
        if(!F_DFA_x){
            printf("MALLOC ERROR (F_DFA_x)\n");
            return 99;
        }
        F_DFA_y = calloc(win_range, sizeof(double));
        if(!F_DFA_y){
            printf("MALLOC ERROR (F_DFA_y)\n");
            return 99;
        }
        rho = calloc(win_range, sizeof(double));
        if(!rho){
            printf("MALLOC ERROR (rho)\n");
            return 99;
        }
    }
    dcca_comp(time, ts1_cum, ts2_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_proc, F_proc_size, rank);
    MPI_Gather(F_proc, F_proc_size, MPI_DOUBLE, F_DCCA, F_proc_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    dcca_comp(time, ts1_cum, ts1_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_proc, F_proc_size, rank);
    MPI_Gather(F_proc, F_proc_size, MPI_DOUBLE, F_DFA_x, F_proc_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    dcca_comp(time, ts2_cum, ts2_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_proc, F_proc_size, rank);
    MPI_Gather(F_proc, F_proc_size, MPI_DOUBLE, F_DFA_y, F_proc_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(rank == 0){
        double *F_rmnd_DCCA, *F_rmnd_DFA_x, *F_rmnd_DFA_y;
        F_rmnd_DCCA = calloc(F_proc_size, sizeof(double));
        if(!F_rmnd_DCCA){
            printf("MALLOC ERROR (F_rmnd_DCCA)\n");
            return 99;
        }
        F_rmnd_DFA_x = calloc(F_proc_size, sizeof(double));
        if(!F_rmnd_DFA_x){
            printf("MALLOC ERROR (F_rmnd_DFA_x)\n");
            return 99;
        }
        F_rmnd_DFA_y = calloc(F_proc_size, sizeof(double));
        if(!F_rmnd_DFA_y){
            printf("MALLOC ERROR (F_rmnd_DFA_y)\n");
            return 99;
        }
        dcca_comp(time, ts1_cum, ts2_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_rmnd_DCCA, F_proc_size, size);
        dcca_comp(time, ts1_cum, ts1_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_rmnd_DFA_x, F_proc_size, size);
        dcca_comp(time, ts2_cum, ts2_cum, L, smallest_win_size, biggest_win_size, pol_ord, F_rmnd_DFA_y, F_proc_size, size);
        FILE *fOut;
        fOut = fopen(path_out, "w");
        for(int i = 0; i < (size*F_proc_size); i++){
            F_DFA_x[i] = sqrt(F_DFA_x[i]);
            F_DFA_y[i] = sqrt(F_DFA_y[i]);
            rho[i] = F_DCCA[i] / (double)(F_DFA_x[i] * F_DFA_y[i]);
            fprintf(fOut, "%lf\n", rho[i]);
        }
        for(int i = 0; i < remd_size; i++){
            F_rmnd_DFA_x[i] = sqrt(F_rmnd_DFA_x[i]);
            F_rmnd_DFA_y[i] = sqrt(F_rmnd_DFA_y[i]);
            rho[i] = F_rmnd_DCCA[i] / (double)(F_rmnd_DFA_x[i] * F_rmnd_DFA_y[i]);
            fprintf(fOut, "%lf\n", rho[i]);
        }
        fclose(fOut);
        free(F_rmnd_DCCA); free(F_rmnd_DFA_x); free(F_rmnd_DFA_y);
        free(F_DCCA); free(F_DFA_x); free(F_DFA_y); free(rho);
    }
    free(time);
    free(ts1_nomean); free(ts2_nomean);
    free(ts1_cum); free(ts2_cum);
    free(F_proc);
}

void dcca_comp(double *t, double *y1, double *y2, int N, int min_win_size, int last_win_len, int ord, double *F, int proc_arr_size, int proc_num)
{
    //fluctuations vector and other arrays
    int F_len = N - min_win_size;
    double *F_nu, *t_fit, *y_fit1, *y_fit2, *diff_vec;
    F_nu = calloc(F_len, sizeof(double));
    if(!F_nu){
        printf("MALLOC ERROR (F_nu)\n");
        return 99;
    }
    t_fit = calloc(last_win_len+1, sizeof(double));
    if(!t_fit){
        printf("MALLOC ERROR (t_fit)\n");
        return 99;
    }
    y_fit1 = calloc(last_win_len+1, sizeof(double));
    if(!y_fit1){
        printf("MALLOC ERROR (y_fit1)\n");
        return 99;
    }
    y_fit2 = calloc(last_win_len+1, sizeof(double));
    if(!y_fit2){
        printf("MALLOC ERROR (y_fit2)\n");
        return 99;
    }
    diff_vec = calloc(last_win_len+1, sizeof(double));
    if(!diff_vec){
        printf("MALLOC ERROR (diff_vec)\n");
        return 99;
    }
    //computation
    int s, N_s, start_lim, end_lim;
    double ang_coeff1, intercept1, r_coeff1, ang_coeff2, intercept2, r_coeff2;
    for(int i = 0; i < proc_arr_size; i++){
        s = (proc_num * proc_arr_size) + i + min_win_size;
        if(s < last_win_len+1){
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
        }else{
            break;
        }
    }
    free(F_nu); free(t_fit); free(y_fit1); free(y_fit2); free(diff_vec);
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
