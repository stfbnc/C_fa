#include "includes.h"

void up_down_lims(char * file_name, int N, double alpha)
{
    int numln = rows_number(file_name) / N; //length of every of the N simulations
    double *rho_mtx;
    rho_mtx = calloc(numln*N, sizeof(double));
    if(!rho_mtx){
        printf("MALLOC ERROR (rho_mtx)\n");
        return 99;
    }
    FILE *f;
    f = fopen(file_name, "r");
    for(int i = 0; i < numln*N; i++)
        fscanf(f, "%lf\n", rho_mtx + i);
    fclose(f);
    double *down_lim, *up_lim, *rho_fixed;
    down_lim = calloc(numln, sizeof(double));
    if(!down_lim){
        printf("MALLOC ERROR (down_lim)\n");
        return 99;
    }
    up_lim = calloc(numln, sizeof(double));
    if(!up_lim){
        printf("MALLOC ERROR (up_lim)\n");
        return 99;
    }
    rho_fixed = calloc(N, sizeof(double));
    if(!rho_fixed){
        printf("MALLOC ERROR (rho_fixed)\n");
        return 99;
    }
    int k;
	f = fopen("down_up_lims.txt", "w");
    for(int i = 0; i < numln; i++){
        k = -1;
        for(int j = i; j < numln*N; j += numln){
            k++;
            rho_fixed[k] = rho_mtx[j];
        }
		histogram_conf(rho_fixed, N, alpha, down_lim+i, up_lim+i);
		fprintf(f, "%lf %lf\n", down_lim[i], up_lim[i]);
    }
	fclose(f);
    free(rho_mtx); free(down_lim); free(up_lim); free(rho_fixed);
}

int rows_number(char *file_name)
{
    FILE *f;
    int stop;
    int lines = 0;
    f = fopen(file_name,"r");
    while(!feof(f)){
        stop = fgetc(f);
        if(stop == '\n')
            lines++;
    }
    fclose(f);
    return lines;
}

double max(double *vec, int len)
{
    double M = vec[0];
    for(int i = 1; i < len; i++)
        if(vec[i] > M)
            M = vec[i];
    return M;
}

double min(double *vec, int len)
{
    double m = vec[0];
    for(int i = 1; i < len; i++)
        if(vec[i] < m)
            m = vec[i];
    return m;
}

void histogram_conf(double *vec, int len, double alpha, double *left_lim, double *right_lim)
{
	double min_val, max_val, bin_width;
	double *norm_hist;
    int j, num_bins;
	min_val = min(vec, len);
	max_val = max(vec, len);
	num_bins = (int)(2.0 * pow(len, 1.0/3.0));
	bin_width = (max_val - min_val) /(double) num_bins;
	norm_hist = calloc(num_bins, sizeof(double));
    if(!norm_hist){
        printf("MALLOC ERROR (norm_hist)\n");
        return 99;
    }
	for(int i = 0; i < len; i++){
        j = floor((vec[i] - min_val) / bin_width);
		norm_hist[j] += 1;
    }
    for(int i = 0; i < num_bins; i++)
        norm_hist[i] /= (double)(bin_width * len);
	double sum_area = 0.0, left_val, right_val;
	int idx = 0;
	while(sum_area <= alpha && idx < (num_bins / 2)){
        left_val = norm_hist[idx];
        right_val = norm_hist[num_bins-1-idx];
		sum_area += ((left_val + right_val) * bin_width);
        idx++;
	}
    idx++;
	*left_lim = min_val + (idx * bin_width);
	*right_lim = max_val - (idx * bin_width);
	free(norm_hist);
}
