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

 Program for estimating SFS from mixed ploidy VCF files. Missing alleles are imputed by drawing them from a Bernoulli distribution.

 Compiling: gcc poly_sfs.c -o poly_sfs -lm

 Usage:
 -vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.
 -inds [file] File listing individuals to use. Optional.
 -sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.
 -mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default 0.6.
 -seed [int] Seed number used for imputation. Default is a random seed.

 Example:
 ./poly_sfs -vcf in.vcf -inds inds.txt -sites 4fold.sites -mis 0.8 -seed 1524796 > out.sfs
*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define merror "ERROR: System out of memory\n\n"

typedef struct {
    int pos;
    char chr[100];
} Site_s;

void openFiles(int argc, char *argv[]);
char **readInds(FILE *ind_file, int *n);
Site_s *readSites(FILE *site_file, int *n);
void readVcf(FILE *vcf_file, char **inds, Site_s *sites, int ind_n, int site_n, long int seed, double mis);
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
    int i, ind_n = 0, site_n = 0;
    long int seed = 0;
    double mis = 0.6;
    char **inds = NULL;
    Site_s *sites = NULL;
    FILE *vcf_file = NULL, *ind_file = NULL, *site_file = NULL;

    if(argc == 1) {
        printHelp();
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "\nParameters:\n");

    for(i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-vcf") == 0) {
            if((vcf_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-vcf %s\n", argv[i]);
        } else if(strcmp(argv[i], "-inds") == 0) {
            if((ind_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-inds %s\n", argv[i]);
        } else if(strcmp(argv[i], "-sites") == 0) {
            if((site_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sites %s\n", argv[i]);
        } else if(strcmp(argv[i], "-mis") == 0) {
            if(isNumeric(argv[++i])) {
                mis = atof(argv[i]);
                if(mis < 0 || mis > 1) {
                    fprintf(stderr, "ERROR: Invalid value for -mis [double]!\n\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "ERROR: Invalid value for -mis [double]!\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-mis %s\n", argv[i]);
        } else if(strcmp(argv[i], "-seed") == 0) {
            if(isNumeric(argv[++i]))
                seed = atoi(argv[i]);
            else {
                fprintf(stderr, "ERROR: Invalid value for -seed [int]!\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-seed %s\n", argv[i]);
        } else if(strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            fprintf(stderr, "\t%s\n", argv[i]);
            printHelp();
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "ERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    fprintf(stderr, "\n");

    if(vcf_file == NULL) {
        fprintf(stderr, "ERROR: -vcf [file] is required!\n\n");
        exit(EXIT_FAILURE);
    }
    if(mis < 0.6)
        fprintf(stderr, "Warning: When over 40%% missing data is allowed, imputation is unreliable\n\n");
    if(ind_file != NULL)
        inds = readInds(ind_file, &ind_n);
    if(site_file != NULL)
        sites = readSites(site_file, &site_n);
    readVcf(vcf_file, inds, sites, ind_n, site_n, seed, mis);
}

char **readInds(FILE *ind_file, int *n) {
    double list_i = 100;
    char *line = NULL, **list = NULL;
    size_t len = 0;
    ssize_t read;

    if((list = malloc(list_i * sizeof(char *))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    while((read = getline(&line, &len, ind_file)) != -1) {
        if(line[0] == '\n' || line[0] == '#')
            continue;
        if((list[*n] = malloc(read * sizeof(char))) == NULL) {
            fprintf(stderr, merror);
            exit(EXIT_FAILURE);
        }
        stringTerminator(line);
        strcpy(list[*n], line);
        *n = *n + 1;
        if(*n >= list_i) {
            list_i += 50;
            if((list = realloc(list, list_i * sizeof(char *))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
        }
    }

    free(line);
    fclose(ind_file);

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

void readVcf(FILE *vcf_file, char **inds, Site_s *sites, int ind_n, int site_n, long int seed, double mis) {
    int i, j = 0, pos = 0, ok = 0, stop = 0, ind_i = 0, site_i = 0, mis_i = 0, *ind_l = NULL;
    double alt_i = 0, hap_i = 0, hap_n = 0, p = 0, *sfs = NULL;
    char chr[100], hap[50], *line = NULL, *temp = NULL, *end = NULL;
    size_t len = 0;
    ssize_t read;

    if(seed == 0) {
        seed = (long int)time(NULL);
        fprintf(stderr, "Seed number used for imputation: %ld\n\n", seed);
    }
    srand(seed);

    while((read = getline(&line, &len, vcf_file)) != -1) {
        if(line[0] == '\n' || (line[0] == '#' && line[1] == '#'))
            continue;
        temp = strtok_r(line, "\t", &end);
        j = 1;
        if(strcmp(temp, "#CHROM") == 0) {
            if(ind_n == 0)
                continue;
            if((ind_l = calloc(read, sizeof(int))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
            while(temp != NULL) {
                if(j > 9) {
                    stringTerminator(temp);
                    for(i = 0; i < ind_n; i++) {
                        if(strcmp(temp, inds[i]) == 0) {
                            ind_l[j] = 1;
                            ind_i++;
                        }
                    }
                }
                temp = strtok_r(NULL, "\t", &end);
                j++;
            }
            if(ind_i == 0) {
                fprintf(stderr, "ERROR: Individuals in -ind file were not found in the VCF file!\n\n");
                exit(EXIT_FAILURE);
            }
            if(ind_i < ind_n)
                fprintf(stderr, "Warning: -ind file contain individuals that are not in the VCF file\n\n");
            continue;
        }
        alt_i = 0;
        hap_i = 0;
        while(temp != NULL) {
            if(j == 1)
                strncpy(chr, temp, 99);
            else if(j == 2) {
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
            } else if(j > 9) {
                if(ind_n > 0 && ind_l[j] == 0) {
                    temp = strtok_r(NULL, "\t", &end);
                    j++;
                    continue;
                }
                if(strchr(temp, ':') != NULL)
                    strncpy(hap, strtok(temp, ":"), 49);
                else {
                    stringTerminator(temp);
                    strncpy(hap, temp, 49);
                }
                if(sfs == NULL) {
                    switch(strlen(hap)) {
                        case 3:
                            hap_n += 2;
                            break;
                        case 7:
                            hap_n += 4;
                            break;
                        case 11:
                            hap_n += 6;
                            break;
                        case 15:
                            hap_n += 8;
                            break;
                        default:
                            fprintf(stderr, "ERROR: Allowed ploidy-levels are 2, 4, 6, and 8!\n\n");
                            exit(EXIT_FAILURE);
                    }
                }
                if(hap[0] == '.') {
                    temp = strtok_r(NULL, "\t", &end);
                    j++;
                    continue;
                }
                switch(strlen(hap)) {
                    case 3:
                        stop = 2;
                        hap_i += 2;
                        break;
                    case 7:
                        stop = 6;
                        hap_i += 4;
                        break;
                    case 11:
                        stop = 10;
                        hap_i += 6;
                        break;
                    case 15:
                        stop = 14;
                        hap_i += 8;
                        break;
                    default:
                        fprintf(stderr, "ERROR: Allowed ploidy-levels are 2, 4, 6, and 8!\n\n");
                        exit(EXIT_FAILURE);
                }
                for(i = 0; i <= stop; i += 2) {
                    if(hap[i] == '0' || hap[i] == '1')
                        alt_i += hap[i] - '0';
                    else {
                        fprintf(stderr, "ERROR: Unknown alleles found at site %s:%i! Only 0 and 1 are allowed.\n\n", chr, pos);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            temp = strtok_r(NULL, "\t", &end);
            j++;
        }
        if(temp == NULL) {
            if(sfs == NULL) {
                if((sfs = calloc(hap_n + 1, sizeof(double))) == NULL) {
                    fprintf(stderr, merror);
                    exit(EXIT_FAILURE);
                }
            }
            if(hap_i / hap_n < mis)
                continue;
            if(hap_i < hap_n) {
                p = alt_i / hap_i;
                mis_i = hap_n - hap_i;
                if(p == 1)
                    alt_i += mis_i;
                else if(p > 0) {
                    for(i = 0; i < mis_i; i++) {
                        if((double)rand() / RAND_MAX < p)
                            alt_i++;
                    }
                }
            }
            sfs[(int)alt_i]++;
        }
    }
    if(sfs == NULL)
        fprintf(stderr, "Warning: SFS is empty. Please check your input files!\n\n");
    else {
        for(i = 0; i <= hap_n; i++) {
            if(i < hap_n)
                printf("%.0f,", sfs[i]);
            else
                printf("%.0f\n", sfs[i]);
        }
        if(isatty(1))
            fprintf(stderr, "\n");
        free(sfs);
    }
    if(ind_n > 0) {
        for(i = 0; i < ind_n; i++)
            free(inds[i]);
        free(inds);
        free(ind_l);
    }
    if(site_n > 0)
        free(sites);
    free(line);
    fclose(vcf_file);
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
    fprintf(stderr, "\nProgram for estimating SFS from mixed ploidy VCF files.\nMissing alleles are imputed by drawing them from a Bernoulli distribution.\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "-vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.\n");
    fprintf(stderr, "-inds [file] File listing individuals to use. Optional.\n");
    fprintf(stderr, "-sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.\n");
    fprintf(stderr, "-mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default 0.6.\n");
    fprintf(stderr, "-seed [int] Seed number used for imputation. Default is a random seed.\n\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "./poly_sfs -vcf in.vcf -inds inds.txt -sites 4fold.sites -mis 0.8 -seed 1524796 > out.sfs\n\n");
}