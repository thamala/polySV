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

 Program for conducting LD-pruning on mixed ploidy VCF files.

 Compiling: gcc prune_ld.c -o prune_ld -lm

 Usage:
 -vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.
 -sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.
 -r2 [int] [int] [double] Excludes sites based on squared genotypic correlation. Requires a window size in number of SNPs, a step size in number of SNPs, and a maximum r2 value.
 -mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default 0.6.
 -maf [double] Minimum minor allele frequency allowed. Default 0.05.

 Example:
 ./prune_ld -vcf in.vcf -sites 4fold.sites -mis 0.8 -maf 0.05 -r2 100 50 0.1 > 4fold_ld_pruned.vcf
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
    int pos;
    char chr[100];
} Site_s;

typedef struct {
    int pos, ok;
    double *geno;
    char ref, alt, id[100], chr[100], **hap;
} SNP_s;

void openFiles(int argc, char *argv[]);
Site_s *readSites(FILE *site_file, int *n);
void readVcf(FILE *vcf_file, Site_s *sites, int win, int step, int site_n, double mis, double maf, double r2);
void estLD(SNP_s *snps, int win, int ind_n, double r2);
double estR2(double geno1[], double geno2[], int n);
void printOut(SNP_s snp, int n);
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
    int i, win = 0, step = 0, out = 0, site_n = 0;
    double mis = 0.6, maf = 0.05, r2 = -1;
    Site_s *sites = NULL;
    FILE *vcf_file = NULL, *site_file = NULL;

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

    if(vcf_file == NULL || win == 0 || step == 0 || r2 == -1) {
        fprintf(stderr, "\nERROR: -vcf [file] and -r2 [int] [int] [double] are required!\n\n");
        exit(EXIT_FAILURE);
    }
    if(maf == 0) {
        fprintf(stderr, "Warning: Doing LD-pruning, setting -maf to 0.05\n\n");
        maf = 0.05;
    }
    if(mis == 0) {
        fprintf(stderr, "Warning: Doing LD-pruning, setting -mis to 0.6\n\n");
        mis = 0.6;
    }
    if(site_file != NULL)
        sites = readSites(site_file, &site_n);
    readVcf(vcf_file, sites, win, step, site_n, mis, maf, r2);
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

void readVcf(FILE *vcf_file, Site_s *sites, int win, int step, int site_n, double mis, double maf, double r2) {
    int i, j = 0, k = 0, pos = 0, ok = 0, stop = 0, geno_n = 0, ind_i = 0, ind_n = 0, site_i = 0, win_n = 0, win_i = 0, step_i = 0, snp_i = 0;
    double ploidy = 0, mis_i = 0, alt_i = 0, hap_i = 0;
    char chr[100], *line = NULL, *temp = NULL, *end = NULL;
    SNP_s *snps = NULL;
    size_t len = 0;
    ssize_t read;

    while((read = getline(&line, &len, vcf_file)) != -1) {
        if(line[0] == '\n')
            continue;
        if(line[0] == '#') {
            printf("%s", line);
            continue;
        }
        if(snps == NULL) {
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
                if((snps[i].hap = malloc(geno_n * sizeof(char *))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
                for(j = 0; j < geno_n; j++) {
                    if((snps[i].hap[j] = malloc(50 * sizeof(char))) == NULL) {
                        fprintf(stderr, merror);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        temp = strtok_r(line, "\t", &end);
        k = 1;
        ind_i = 0;
        mis_i = 0;
        alt_i = 0;
        hap_i = 0;
        memset(snps[win_i].geno, 0, geno_n * sizeof(double));
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
                if(win_n > 0 && strcmp(snps[0].chr, chr) != 0) {
                    estLD(snps, win_n + 1, ind_n, r2);
                    for(i = 0; i < win; i++) {
                        if(snps[win_i].ok == 1) {
                            printOut(snps[win_i], ind_n);
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
                }
                strcpy(snps[win_i].chr, chr);
                snps[win_i].pos = pos;
                snps[win_i].ok = -1;
            } else if(k == 3)
                strncpy(snps[win_i].id, temp, 99);
            else if(k == 4)
                snps[win_i].ref = temp[0];
            else if(k == 5)
                snps[win_i].alt = temp[0];
            else if(k > 9) {
                if(strchr(temp, ':') != NULL)
                    strncpy(snps[win_i].hap[ind_i], strtok(temp, ":"), 49);
                else {
                    stringTerminator(temp);
                    strncpy(snps[win_i].hap[ind_i], temp, 49);
                }
                if(snps[win_i].hap[ind_i][0] == '.') {
                    snps[win_i].geno[ind_i++] = -1;
                    mis_i++;
                    temp = strtok_r(NULL, "\t", &end);
                    k++;
                    continue;
                }
                switch(strlen(snps[win_i].hap[ind_i])) {
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
                    if(snps[win_i].hap[ind_i][i] == '0' || snps[win_i].hap[ind_i][i] == '1') {
                        snps[win_i].geno[ind_i] += snps[win_i].hap[ind_i][i] - '0';
                        alt_i += snps[win_i].hap[ind_i][i] - '0';
                    } else {
                        fprintf(stderr, "\nERROR: Unknown alleles found at site %s:%i! Only 0 and 1 are allowed.\n\n", chr, pos);
                        exit(EXIT_FAILURE);
                    }
                }
                hap_i += ploidy;
                ind_i++;
            }
            temp = strtok_r(NULL, "\t", &end);
            k++;
        }
        if(temp == NULL) {
            if(ind_n == 0)
                ind_n = ind_i;
            if(mis_i / ind_n > 1 - mis || mis_i == ind_n) {
                snps[win_i].chr[0] = '\0';
                continue;
            }
            if(alt_i / hap_i < maf || alt_i / hap_i > 1 - maf) {
                snps[win_i].chr[0] = '\0';
                continue;
            }
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
                    printOut(snps[win_i], ind_n);
                    snps[win_i].ok = 0;
                    snp_i++;
                }
            }
        }
    }
    estLD(snps, win_n + 1, ind_n, r2);
    for(i = 0; i < win; i++) {
        if(snps[win_i].ok == 1) {
            printOut(snps[win_i], ind_n);
            snps[win_i].ok = 0;
            snp_i++;
        }
        win_i++;
        if(win_i == win)
            win_i = 0;
    }

    if(isatty(1))
        fprintf(stderr, "\n");
    fprintf(stderr, "After pruning, kept %i variants\n\n", snp_i);

    for(i = 0; i < win; i++) {
        free(snps[i].geno);
        for(j = 0; j < geno_n; j++)
            free(snps[i].hap[j]);
        free(snps[i].hap);
    }
    free(snps);
    if(site_n > 0)
        free(sites);
    free(line);
    fclose(vcf_file);
}

void estLD(SNP_s *snps, int win, int ind_n, double r2) {
    int i, j, k;
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

void printOut(SNP_s snp, int n) {
    int i;
    printf("%s\t%i\t%s\t%c\t%c\t.\tPASS\t.\tGT:FT\t", snp.chr, snp.pos, snp.id, snp.ref, snp.alt);
    for(i = 0; i < n; i++) {
        if(i < n - 1)
            printf("%s:PASS\t", snp.hap[i]);
        else
            printf("%s:PASS\n", snp.hap[i]);
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
    fprintf(stderr, "\nProgram for conducting LD-pruning on mixed ploidy VCF files.\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "-vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.\n");
    fprintf(stderr, "-sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.\n");
    fprintf(stderr, "-r2 [int] [int] [double] Excludes sites based on squared genotypic correlation. Requires a window size in number of SNPs, a step size in number of SNPs, and a maximum r2 value.\n");
    fprintf(stderr, "-mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default 0.6.\n");
    fprintf(stderr, "-maf [double] Minimum minor allele frequency allowed. Default 0.05.\n\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "./prune_ld -vcf in.vcf -sites 4fold.sites -mis 0.8 -maf 0.05 -r2 100 50 0.1 > 4fold_ld_pruned.vcf\n\n");
}