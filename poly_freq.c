/*
 Copyright (C) 2023 Tuomas Hamala

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 For any other inquiries, send an email to tuomas.hamala@gmail.com

 ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

 Program for estimating allele frequencies from mixed ploidy VCF files.
 Output will be either population-specific allele frequencies or allele counts in the format required by BayPass.

 Compiling: gcc poly_freq.c -o poly_freq -lm

 Usage:
 -vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.
 -pops [file] Tab delimited file listing individuals to use and their populations (format: individual id, population id).
 -sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.
 -mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default > 0.
 -maf [double] Minimum minor allele frequency allowed. Default 0.
 -r2 [int] [int] [double] Excludes sites based on squared genotypic correlation. Requires a window size in number of SNPs, a step size in number of SNPs, and a maximum r2 value. Optional.
 -out [int] Whether to output allele frequencies (0) or allele counts in the BayPass format (1). Default 0.
 -info [string] If -out is 1, records populations and locations of used SNPs into this file. Default 'info.txt'.

 Example:
 ./poly_freq -vcf in.vcf -pops pops.txt -sites 4fold.sites -mis 0.8 -maf 0.05 -r2 100 50 0.1 -out 1 -info 4fold_ld_pruned.info > 4fold_ld_pruned.baypass
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n\n"

typedef struct {
    int idx;
    char ind[200], pop[200];
} Pop_s;

typedef struct {
    int pos;
    char chr[100];
} Site_s;

typedef struct {
    int pos, ok;
    double *geno, **counts;
    char chr[100];
} SNP_s;

void openFiles(int argc, char *argv[]);
Pop_s *readPops(FILE *pop_file, FILE *out_file, int out, int *n, int *m);
Site_s *readSites(FILE *site_file, int *n);
void readVcf(FILE *vcf_file, FILE *out_file, Pop_s *pops, Site_s *sites, int win, int step, int out, int ind_n, int pop_n, int site_n, double mis, double maf, double r2);
void estLD(SNP_s *snps, int win, int ind_n, double r2);
double estR2(double geno1[], double geno2[], int n);
void printOut(FILE *out_file, double **counts, char chr[], int pos, int out, int n);
int isNumeric(const char *s);
void stringTerminator(char *string);
void printHelp(void);

int main(int argc, char *argv[]) {
    int second = 0, minute = 0, hour = 0;
    time_t timer = 0;

    timer = time(NULL);
    openFiles(argc, argv);
    second = time(NULL) - timer;
    minute = second / 60;
    hour = second / 3600;

    fprintf(stderr, "Done!");
    if(hour > 0)
        fprintf(stderr, "\nElapsed time: %i h, %i min & %i sec\n\n", hour, minute - hour * 60, second - minute * 60);
    else if(minute > 0)
        fprintf(stderr, "\nElapset time: %i min & %i sec\n\n", minute, second - minute * 60);
    else if(second > 5)
        fprintf(stderr, "\nElapsed time: %i sec\n\n", second);
    else
        fprintf(stderr, "\n\n");

    return 0;
}

void openFiles(int argc, char *argv[]) {
    int i, win = 0, step = 0, out = 0, ind_n = 0, pop_n = 0, site_n = 0;
    double mis = 0, maf = 0, r2 = 1;
    char info[200] = "info.txt";
    Pop_s *pops = NULL;
    Site_s *sites = NULL;
    FILE *vcf_file = NULL, *pop_file = NULL, *site_file = NULL, *out_file = NULL;

    if(argc == 1) {
        printHelp();
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "\nParameters:\n");

    for(i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-vcf") == 0) {
            if((vcf_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-vcf %s\n", argv[i]);
        } else if(strcmp(argv[i], "-pops") == 0) {
            if((pop_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-pops %s\n", argv[i]);
        } else if(strcmp(argv[i], "-sites") == 0) {
            if((site_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sites %s\n", argv[i]);
        } else if(strcmp(argv[i], "-mis") == 0) {
            if(isNumeric(argv[++i])) {
                mis = atof(argv[i]);
                if(mis < 0 || mis > 1) {
                    fprintf(stderr, "\nERROR: Invalid value for -mis [double]!\n\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "\nERROR: Invalid value for -mis [double]!\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-mis %s\n", argv[i]);
        } else if(strcmp(argv[i], "-maf") == 0) {
            if(isNumeric(argv[++i])) {
                maf = atof(argv[i]);
                if(maf < 0 || maf > 1) {
                    fprintf(stderr, "\nERROR: Invalid value for -maf [double]!\n\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "\nERROR: Invalid value for -maf [double]!\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-maf %s\n", argv[i]);
        } else if(strcmp(argv[i], "-r2") == 0) {
            if(isNumeric(argv[++i])) {
                win = atoi(argv[i]);
                if(win < 1) {
                    fprintf(stderr, "\nERROR: Invalid value for the -r2 window size [int]!\n\n");
                    exit(EXIT_FAILURE);
                }
            }
            if(isNumeric(argv[++i])) {
                step = atoi(argv[i]);
                if(step > win || step < 1) {
                    fprintf(stderr, "\nERROR: Invalid value for the -r2 step size [int]!\n\n");
                    exit(EXIT_FAILURE);
                }
            }
            if(isNumeric(argv[++i])) {
                r2 = atof(argv[i]);
                if(r2 < 0 || r2 > 1) {
                    fprintf(stderr, "\nERROR: Invalid value for -r2 [double]!\n\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "\nERROR: Invalid value for -r2 [int] [int] [double]!\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-r2 %i %i %s\n", win, step, argv[i]);
        } else if(strcmp(argv[i], "-out") == 0) {
            if(isNumeric(argv[++i]))
                out = atoi(argv[i]);
            if(out != 0 && out != 1) {
                fprintf(stderr, "\nERROR: Invalid value for -out [int]! Allowed are 0 (allele frequencies) and 1 (allele counts).\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-out %s\n", argv[i]);
        } else if(strcmp(argv[i], "-info") == 0) {
            strncpy(info, argv[++i], 199);
            fprintf(stderr, "\t-info %s\n", argv[i]);
        } else if(strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            fprintf(stderr, "\t%s\n", argv[i]);
            printHelp();
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    fprintf(stderr, "\n");

    if(vcf_file == NULL || pop_file == NULL) {
        fprintf(stderr, "\nERROR: -vcf [file] and -pops [file] are required!\n\n");
        exit(EXIT_FAILURE);
    }
    if(r2 < 1 && maf == 0) {
        fprintf(stderr, "Warning: Doing LD-pruning, setting -maf to 0.05\n\n");
        maf = 0.05;
    }
    if(out == 1) {
        if((out_file = fopen(info, "w")) == NULL) {
            fprintf(stderr, "\n\nERROR: Cannot create file '%s'\n\n", info);
            exit(EXIT_FAILURE);
        }
    }
    if(site_file != NULL)
        sites = readSites(site_file, &site_n);
    pops = readPops(pop_file, out_file, out, &ind_n, &pop_n);
    readVcf(vcf_file, out_file, pops, sites, win, step, out, ind_n, pop_n, site_n, mis, maf, r2);

    if(out == 1)
        fclose(out_file);
}

Pop_s *readPops(FILE *pop_file, FILE *out_file, int out, int *n, int *m) {
    int i;
    double list_i = 200, pops_i = 50;
    char *line = NULL, **pops = NULL;
    Pop_s *list = NULL;
    size_t len = 0;
    ssize_t read;

    if((list = malloc(list_i * sizeof(Pop_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    while((read = getline(&line, &len, pop_file)) != -1) {
        if(line[0] == '\n' || line[0] == '#')
            continue;
        stringTerminator(line);
        strncpy(list[*n].ind, strtok(line, "\t"), 199);
        strncpy(list[*n].pop, strtok(NULL, "\t"), 199);
        if(pops == NULL) {
            if((pops = malloc(pops_i * sizeof(char *))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
            for(i = 0; i < pops_i; i++) {
                if((pops[i] = calloc(200, sizeof(char))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
            }
            strcpy(pops[*m], list[*n].pop);
            *m = *m + 1;
        }
        for(i = 0; i <= *m; i++) {
            if(strcmp(pops[i], list[*n].pop) == 0) {
                list[*n].idx = i;
                break;
            }
        }
        if(i > *m) {
            list[*n].idx = *m;
            strcpy(pops[*m], list[*n].pop);
            *m = *m + 1;
            if(*m >= pops_i) {
                pops_i += 20;
                if((pops = realloc(pops, pops_i * sizeof(char *))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
                for(i = *m; i < pops_i; i++) {
                    if((pops[i] = calloc(200, sizeof(char))) == NULL) {
                        fprintf(stderr, merror);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        *n = *n + 1;
        if(*n >= list_i) {
            list_i += 100;
            if((list = realloc(list, list_i * sizeof(Pop_s))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
        }
    }
    if(out == 0) {
        printf("\t");
        for(i = 0; i < *m; i++) {
            if(i < *m - 1)
                printf("%s\t", pops[i]);
            else
                printf("%s\n", pops[i]);
        }
    } else {
        fprintf(out_file, "#");
        for(i = 0; i < *m; i++) {
            if(i < *m - 1)
                fprintf(out_file, "%s\t", pops[i]);
            else
                fprintf(out_file, "%s\n", pops[i]);
        }
    }
    free(line);
    for(i = 0; i < pops_i; i++)
        free(pops[i]);
    free(pops);
    fclose(pop_file);

    return list;
}

Site_s *readSites(FILE *site_file, int *n) {
    double list_i = 1e6;
    char *line = NULL;
    Site_s *list = NULL;
    size_t len = 0;
    ssize_t read;

    if((list = malloc(list_i * sizeof(Site_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    while((read = getline(&line, &len, site_file)) != -1) {
        if(line[0] == '\n' || line[0] == '#')
            continue;
        strncpy(list[*n].chr, strtok(line, "\t"), 99);
        list[*n].pos = atoi(strtok(NULL, "\t"));
        if(*n > 0) {
            if(strcmp(list[*n].chr, list[*n - 1].chr) < 0) {
                fprintf(stderr, "ERROR: Site file is not sorted. Use: sort -k1,1 -k2,2n list.sites > sorted.sites\n\n");
                exit(EXIT_FAILURE);
            } else if(strcmp(list[*n].chr, list[*n - 1].chr) == 0) {
                if(list[*n].pos < list[*n - 1].pos) {
                    fprintf(stderr, "ERROR: Site file is not sorted. Use: sort -k1,1 -k2,2n list.sites > sorted.sites\n\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        *n = *n + 1;
        if(*n >= list_i) {
            list_i += 1e5;
            if((list = realloc(list, list_i * sizeof(Site_s))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
        }
    }

    free(line);
    fclose(site_file);

    return list;
}

void readVcf(FILE *vcf_file, FILE *out_file, Pop_s *pops, Site_s *sites, int win, int step, int out, int ind_n, int pop_n, int site_n, double mis, double maf, double r2) {
    int i, j = 0, k = 0, pos = 0, ok = 0, stop = 0, geno_n = 0, ind_i = 0, site_i = 0, win_n = 0, win_i = 0, step_i = 0, snp_i = 0, *pop_l = NULL;
    double ploidy = 0, mis_i = 0, alt_i = 0, hap_i = 0, **counts = NULL;
    char chr[100], hap[50], *line = NULL, *temp = NULL, *end = NULL;
    SNP_s *snps = NULL;
    size_t len = 0;
    ssize_t read;

    while((read = getline(&line, &len, vcf_file)) != -1) {
        if(line[0] == '\n' || (line[0] == '#' && line[1] == '#'))
            continue;
        temp = strtok_r(line, "\t", &end);
        k = 1;
        if(strcmp(temp, "#CHROM") == 0) {
            if(r2 < 1) {
                if((snps = malloc(win * sizeof(SNP_s))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
                geno_n = (int)read;
                for(i = 0; i < win; i++) {
                    if((snps[i].geno = malloc(geno_n * sizeof(double))) == NULL) {
                        fprintf(stderr, merror);
                        exit(EXIT_FAILURE);
                    }
                    if((snps[i].counts = malloc(pop_n * sizeof(double *))) == NULL) {
                        fprintf(stderr, merror);
                        exit(EXIT_FAILURE);
                    }
                    for(j = 0; j < pop_n; j++) {
                        if((snps[i].counts[j] = malloc(2 * sizeof(double))) == NULL) {
                            fprintf(stderr, merror);
                            exit(EXIT_FAILURE);
                        }
                    }
                }
            } else {
                if((counts = malloc(pop_n * sizeof(double *))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
                for(i = 0; i < pop_n; i++) {
                    if((counts[i] = malloc(2 * sizeof(double))) == NULL) {
                        fprintf(stderr, merror);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            if((pop_l = malloc(read * sizeof(int))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
            memset(pop_l, -1, read * sizeof(int));
            while(temp != NULL) {
                if(k > 9) {
                    stringTerminator(temp);
                    for(i = 0; i < ind_n; i++) {
                        if(strcmp(temp, pops[i].ind) == 0) {
                            pop_l[k] = pops[i].idx;
                            ind_i++;
                        }
                    }
                }
                temp = strtok_r(NULL, "\t", &end);
                k++;
            }
            if(ind_i == 0) {
                fprintf(stderr, "\nERROR: Individuals in pops file were not found in the VCF file!\n\n");
                exit(EXIT_FAILURE);
            }
            if(ind_i < ind_n) {
                fprintf(stderr, "Warning: pops file contains individuals that are not in the VCF file\n\n");
                ind_n = ind_i;
            }
            continue;
        }
        ind_i = 0;
        mis_i = 0;
        alt_i = 0;
        hap_i = 0;
        if(r2 < 1) {
            memset(snps[win_i].geno, 0, geno_n * sizeof(double));
            for(i = 0; i < pop_n; i++)
                memset(snps[win_i].counts[i], 0, 2 * sizeof(double));
        } else {
            for(i = 0; i < pop_n; i++)
                memset(counts[i], 0, 2 * sizeof(double));
        }
        while(temp != NULL) {
            if(k == 1)
                strncpy(chr, temp, 99);
            else if(k == 2) {
                pos = atoi(temp);
                if(site_n > 0) {
                    ok = 0;
                    while(site_i < site_n) {
                        if(strcmp(chr, sites[site_i].chr) == 0) {
                            if(pos == sites[site_i].pos) {
                                ok = 1;
                                break;
                            } else if(pos < sites[site_i].pos)
                                break;
                        } else if(strcmp(chr, sites[site_i].chr) < 0)
                            break;
                        site_i++;
                    }
                    if(ok == 0)
                        break;
                }
                if(r2 < 1) {
                    if(win_n > 0 && strcmp(snps[0].chr, chr) != 0) {
                        estLD(snps, win_n + 1, ind_n, r2);
                        for(i = 0; i < win; i++) {
                            if(snps[win_i].ok == 1) {
                                printOut(out_file, snps[win_i].counts, snps[win_i].chr, snps[win_i].pos, out, pop_n);
                                snps[win_i].ok = 0;
                                snp_i++;
                            }
                            win_i++;
                            if(win_i == win)
                                win_i = 0;
                        }
                        win_n = 0;
                        win_i = 0;
                        step_i = 0;
                        memset(snps[win_i].geno, 0, geno_n * sizeof(double));
                        for(i = 0; i < pop_n; i++)
                            memset(snps[win_i].counts[i], 0, 2 * sizeof(double));
                    }
                    strcpy(snps[win_i].chr, chr);
                    snps[win_i].pos = pos;
                    snps[win_i].ok = -1;
                }
            } else if(k > 9) {
                if(pop_l[k] == -1) {
                    temp = strtok_r(NULL, "\t", &end);
                    k++;
                    continue;
                }
                if(strchr(temp, ':') != NULL)
                    strncpy(hap, strtok(temp, ":"), 49);
                else {
                    stringTerminator(temp);
                    strncpy(hap, temp, 49);
                }
                if(hap[0] == '.') {
                    if(r2 < 1)
                        snps[win_i].geno[ind_i++] = -1;
                    mis_i++;
                    temp = strtok_r(NULL, "\t", &end);
                    k++;
                    continue;
                }
                switch(strlen(hap)) {
                    case 3:
                        ploidy = 2;
                        stop = 2;
                        break;
                    case 7:
                        ploidy = 4;
                        stop = 6;
                        break;
                    case 11:
                        ploidy = 6;
                        stop = 10;
                        break;
                    case 15:
                        ploidy = 8;
                        stop = 14;
                        break;
                    default:
                        fprintf(stderr, "\nERROR: Allowed ploidy-levels are 2, 4, 6, and 8!\n\n");
                        exit(EXIT_FAILURE);
                }
                for(i = 0; i <= stop; i += 2) {
                    if(hap[i] == '0' || hap[i] == '1') {
                        if(r2 < 1) {
                            snps[win_i].geno[ind_i] += hap[i] - '0';
                            snps[win_i].counts[pop_l[k]][1] += hap[i] - '0';
                        } else
                            counts[pop_l[k]][1] += hap[i] - '0';
                        alt_i += hap[i] - '0';
                    } else {
                        fprintf(stderr, "\nERROR: Unknown alleles found at site %s:%i! Only 0 and 1 are allowed.\n\n", chr, pos);
                        exit(EXIT_FAILURE);
                    }
                }
                if(r2 < 1)
                    snps[win_i].counts[pop_l[k]][0] += ploidy;
                else
                    counts[pop_l[k]][0] += ploidy;
                hap_i += ploidy;
                ind_i++;
            }
            temp = strtok_r(NULL, "\t", &end);
            k++;
        }
        if(temp == NULL) {
            if(mis_i / ind_n > 1 - mis || mis_i == ind_n) {
                if(r2 < 1)
                    snps[win_i].chr[0] = '\0';
                continue;
            }
            if(alt_i / hap_i < maf || alt_i / hap_i > 1 - maf) {
                if(r2 < 1)
                    snps[win_i].chr[0] = '\0';
                continue;
            }
            if(r2 < 1) {
                if((win_n == win - 1 && step_i >= step) || (win == step && win_i == win - 1)) {
                    estLD(snps, win_n + 1, ind_n, r2);
                    step_i = 0;
                }
                if(win_n < win - 1)
                    win_n++;
                step_i++;
                win_i++;
                if(win_i == win)
                    win_i = 0;
                if(win_n == win - 1) {
                    if(snps[win_i].ok == 1) {
                        printOut(out_file, snps[win_i].counts, snps[win_i].chr, snps[win_i].pos, out, pop_n);
                        snps[win_i].ok = 0;
                        snp_i++;
                    }
                }
            } else {
                printOut(out_file, counts, chr, pos, out, pop_n);
                snp_i++;
            }
        }
    }
    if(r2 < 1) {
        estLD(snps, win_n + 1, ind_n, r2);
        for(i = 0; i < win; i++) {
            if(snps[win_i].ok == 1) {
                printOut(out_file, snps[win_i].counts, snps[win_i].chr, snps[win_i].pos, out, pop_n);
                snps[win_i].ok = 0;
                snp_i++;
            }
            win_i++;
            if(win_i == win)
                win_i = 0;
        }
    }

    if(isatty(1))
        fprintf(stderr, "\n");
    fprintf(stderr, "Kept %i variants\n\n", snp_i);

    if(r2 < 1) {
        for(i = 0; i < win; i++) {
            free(snps[i].geno);
            for(j = 0; j < pop_n; j++)
                free(snps[i].counts[j]);
            free(snps[i].counts);
        }
        free(snps);
    } else {
        for(i = 0; i < pop_n; i++)
            free(counts[i]);
        free(counts);
    }
    free(pop_l);
    free(pops);
    if(site_n > 0)
        free(sites);
    free(line);
    fclose(vcf_file);
}

void estLD(SNP_s *snps, int win, int ind_n, double r2) {
    int i, j;
    for(i = 0; i < win; i++) {
        if(snps[i].chr[0] == '\0')
            continue;
        for(j = i + 1; j < win; j++) {
            if(snps[j].chr[0] == '\0')
                continue;
            if(strcmp(snps[i].chr, snps[j].chr) != 0)
                continue;
            /*
            fprintf(stderr, "%s %i: ", snps[i].chr, snps[i].pos);
            for(k = 0; k < ind_n; k++)
                fprintf(stderr, "%.0f ", snps[i].geno[k]);
            fprintf(stderr, "\n");
            fprintf(stderr, "%s %i: ", snps[j].chr, snps[j].pos);
            for(k = 0; k < ind_n; k++)
                fprintf(stderr, "%.0f ", snps[j].geno[k]);
            fprintf(stderr, "\n");
            fprintf(stderr, "%f\n\n", estR2(snps[i].geno, snps[j].geno, ind_n));
            */
            if(estR2(snps[i].geno, snps[j].geno, ind_n) > r2)
                break;
        }
        if(j == win && (snps[i].ok == -1 || snps[i].ok == 1))
            snps[i].ok = 1;
        else
            snps[i].ok = 0;
    }
}

double estR2(double geno1[], double geno2[], int n) {
    int i, m = 0;
    double sum1 = 0, sum2 = 0, sum12 = 0, sqsum1 = 0, sqsum2 = 0, r = 0;
    m = n;
    for(i = 0; i < m; i++) {
        if(geno1[i] == -1 || geno2[i] == -1) {
            n--;
            continue;
        }
        sum1 += geno1[i];
        sum2 += geno2[i];
        sum12 += geno1[i] * geno2[i];
        sqsum1 += geno1[i] * geno1[i];
        sqsum2 += geno2[i] * geno2[i];
    }
    r = (n * sum12 - sum1 * sum2) / sqrt((n * sqsum1 - sum1 * sum1) * (n * sqsum2 - sum2 * sum2));
    return r * r;
}

void printOut(FILE *out_file, double **counts, char chr[], int pos, int out, int n) {
    int i;
    if(out == 0)
        printf("%s:%i\t", chr, pos);
    else
        fprintf(out_file, "%s\t%i\n", chr, pos);
    for(i = 0; i < n; i++) {
        if(out == 0) {
            if(i < n - 1)
                printf("%f\t", counts[i][1] / (counts[i][0]));
            else
                printf("%f\n", counts[i][1] / (counts[i][0]));
        } else {
            if(i < n - 1)
                printf("%.0f %.0f ", counts[i][0] - counts[i][1], counts[i][1]);
            else
                printf("%.0f %.0f\n", counts[i][0] - counts[i][1], counts[i][1]);
        }
    }
}

int isNumeric(const char *s) {
    char *p;
    if(s == NULL || *s == '\0' || isspace(*s))
        return 0;
    strtod(s, &p);
    return *p == '\0';
}

void stringTerminator(char *string) {
    string[strcspn(string, "\n")] = 0;
}

void printHelp(void) {
    fprintf(stderr, "\nProgram for estimating allele frequencies from mixed ploidy VCF files.\nOutput will be either population-specific allele frequencies or allele counts in the format required by BayPass.\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "-vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.\n");
    fprintf(stderr, "-pops [file] Tab delimited file listing individuals to use and their populations (format: individual id, population id).\n");
    fprintf(stderr, "-sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.\n");
    fprintf(stderr, "-mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default > 0.\n");
    fprintf(stderr, "-maf [double] Minimum minor allele frequency allowed. Default 0.\n");
    fprintf(stderr, "-r2 [int] [int] [double] Excludes sites based on squared genotypic correlation. Requires a window size in number of SNPs, a step size in number of SNPs, and a maximum r2 value. Optional.\n");
    fprintf(stderr, "-out [int] Whether to output allele frequencies (0) or allele counts in the BayPass format (1). Default 0.\n");
    fprintf(stderr, "-info [string] If -out is 1, records populations and locations of used SNPs into this file. Default 'info.txt'.\n\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "./poly_freq -vcf in.vcf -pops pops.txt -sites 4fold.sites -mis 0.8 -maf 0.05 -r2 100 50 0.1 -out 1 -info 4fold_ld_pruned.info > 4fold_ld_pruned.baypass\n\n");
}