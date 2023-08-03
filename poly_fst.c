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

 Program for estimating pairwise Fst and Dxy from mixed ploidy VCF files.

 Compiling: gcc poly_fst.c -o poly_fst -lm

 Usage:
 -vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.
 -pop1 [file] File listing individuals from population 1.
 -pop2 [file] File listing individuals from population 2.
 -sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.
 -genes [file] Tab delimited file listing genes to use (format: chr, start, end, id). Output will be Fst/Dxy calculated for each gene. Optional.
 -mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default > 0.
 -maf [double] Minimum minor allele frequency allowed. Default 0.
 -stat [string] Whether to calculate 'fst' or 'dxy'. Default 'fst'. Note that dxy requires invariant sites to be included in the VCF file.
 -out [int] Whether to print full output (0) or genome-wide estimate only (1). Default 0.

 Example:
 ./poly_fst -vcf in.vcf -pop1 pop1.txt -pop2 pop2.txt -sites 4fold.sites -genes genes.txt -mis 0.8 -stat dxy > out_gene.dxy
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

typedef struct {
    int start, end;
    double hw, hb, n;
    char chr[100], id[200];
} Gene_s;

void openFiles(int argc, char *argv[]);
char **readInds(FILE *ind_file, int *n);
Site_s *readSites(FILE *site_file, int *n);
Gene_s *readGenes(FILE *gene_file, int *n);
void readVcf(FILE *vcf_file, char **pop1, char **pop2, Site_s *sites, Gene_s *genes, int stat, int out, int pop1_n, int pop2_n, int site_n, int gene_n, double mis, double maf);
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
    int i, stat = 0, pop1_n = 0, pop2_n = 0, site_n = 0, gene_n = 0, out = 0;
    double mis = 0, maf = 0;
    char temp[10], **pop1 = NULL, **pop2 = NULL;
    Site_s *sites = NULL;
    Gene_s *genes = NULL;
    FILE *vcf_file = NULL, *pop1_file = NULL, *pop2_file = NULL, *site_file = NULL, *gene_file = NULL;

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
        } else if(strcmp(argv[i], "-pop1") == 0) {
            if((pop1_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-pop1 %s\n", argv[i]);
        } else if(strcmp(argv[i], "-pop2") == 0) {
            if((pop2_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-pop2 %s\n", argv[i]);
        } else if(strcmp(argv[i], "-sites") == 0) {
            if((site_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sites %s\n", argv[i]);
        } else if(strcmp(argv[i], "-genes") == 0) {
            if((gene_file = fopen(argv[++i], "r")) == NULL) {
                fprintf(stderr, "ERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-genes %s\n", argv[i]);
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
        } else if(strcmp(argv[i], "-maf") == 0) {
            if(isNumeric(argv[++i])) {
                maf = atof(argv[i]);
                if(maf < 0 || maf > 1) {
                    fprintf(stderr, "ERROR: Invalid value for -maf [double]!\n\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "ERROR: Invalid value for -maf [double]!\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-maf %s\n", argv[i]);
        } else if(strcmp(argv[i], "-stat") == 0) {
            strncpy(temp, argv[++i], 9);
            if(strcmp(temp, "fst") == 0)
                stat = 0;
            else if(strcmp(temp, "dxy") == 0)
                stat = 1;
            else {
                fprintf(stderr, "ERROR: Invalid input for -stat [string]! Allowed are 'fst' and 'dxy'\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-stat %s\n", argv[i]);
        } else if(strcmp(argv[i], "-out") == 0) {
            if(isNumeric(argv[++i])) {
                out = atoi(argv[i]);
                if(out != 0 && out != 1) {
                    fprintf(stderr, "ERROR: Invalid value for -out [int]! Only 0 (full output) and 1 (genome-wide output) are allowed\n\n");
                    exit(EXIT_FAILURE);
                }
            } else {
                fprintf(stderr, "ERROR: Invalid value for -out [int]! Only 0 (full output) and 1 (genome-wide output) are allowed\n\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-out %s\n", argv[i]);
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

    if(vcf_file == NULL || pop1_file == NULL || pop2_file == NULL) {
        fprintf(stderr, "ERROR: -vcf [file] -pop1 [file] -pop2 [file] are required!\n\n");
        exit(EXIT_FAILURE);
    }
    pop1 = readInds(pop1_file, &pop1_n);
    pop2 = readInds(pop2_file, &pop2_n);
    if(site_file != NULL)
        sites = readSites(site_file, &site_n);
    if(gene_file != NULL)
        genes = readGenes(gene_file, &gene_n);
    readVcf(vcf_file, pop1, pop2, sites, genes, stat, out, pop1_n, pop2_n, site_n, gene_n, mis, maf);
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

Gene_s *readGenes(FILE *gene_file, int *n) {
    double list_i = 1e4;
    char *line = NULL;
    Gene_s *list = NULL;
    size_t len = 0;
    ssize_t read;

    if((list = malloc(list_i * sizeof(Gene_s))) == NULL) {
        fprintf(stderr, merror);
        exit(EXIT_FAILURE);
    }
    while((read = getline(&line, &len, gene_file)) != -1) {
        if(line[0] == '\n')
            continue;
        stringTerminator(line);
        strncpy(list[*n].chr, strtok(line, "\t"), 99);
        list[*n].start = atoi(strtok(NULL, "\t"));
        list[*n].end = atoi(strtok(NULL, "\t"));
        strncpy(list[*n].id, strtok(NULL, "\t"), 199);
        list[*n].hw = 0;
        list[*n].hb = 0;
        list[*n].n = 0;
        *n = *n + 1;
        if(*n >= list_i) {
            list_i += 1000;
            if((list = realloc(list, list_i * sizeof(Gene_s))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
        }
    }

    free(line);
    fclose(gene_file);

    return list;
}

void readVcf(FILE *vcf_file, char **pop1, char **pop2, Site_s *sites, Gene_s *genes, int stat, int out, int pop1_n, int pop2_n, int site_n, int gene_n, double mis, double maf) {
    int i, j = 0, k = 0, pos = 0, ok = 0, stop = 0, pop_i = 0, site_i = 0, gene_i = 0, *pop_l = NULL;
    double ploidy = 0, mis1_i = 0, mis2_i = 0, pop1_i = 0, pop2_i = 0, p1 = 0, p2 = 0, n1 = 0, n2 = 0, hw = 0, hb = 0, tot_hw = 0, tot_hb = 0, tot_n = 0;
    char chr[100], hap[50], *line = NULL, *temp = NULL, *end = NULL;
    size_t len = 0;
    ssize_t read;

    while((read = getline(&line, &len, vcf_file)) != -1) {
        if(line[0] == '\n' || (line[0] == '#' && line[1] == '#'))
            continue;
        temp = strtok_r(line, "\t", &end);
        j = 1;
        if(strcmp(temp, "#CHROM") == 0) {
            if((pop_l = calloc(read, sizeof(int))) == NULL) {
                fprintf(stderr, merror);
                exit(EXIT_FAILURE);
            }
            while(temp != NULL) {
                if(j > 9) {
                    stringTerminator(temp);
                    for(i = 0; i < pop1_n; i++) {
                        if(strcmp(temp, pop1[i]) == 0) {
                            pop_l[j] = 1;
                            pop_i++;
                        }
                    }
                    for(i = 0; i < pop2_n; i++) {
                        if(strcmp(temp, pop2[i]) == 0) {
                            pop_l[j] = 2;
                            pop_i++;
                        }
                    }
                }
                temp = strtok_r(NULL, "\t", &end);
                j++;
            }
            if(pop_i == 0) {
                fprintf(stderr, "ERROR: Individuals in -pop1 and -pop2 files were not found in the VCF file!\n\n");
                exit(EXIT_FAILURE);
            }
            if(pop_i < pop1_n + pop2_n)
                fprintf(stderr, "Warning: -pop1 and -pop2 files contain individuals that are not in the VCF file\n\n");
            continue;
        }
        mis1_i = 0;
        mis2_i = 0;
        pop1_i = 0;
        pop2_i = 0;
        p1 = 0;
        p2 = 0;
        n1 = 0;
        n2 = 0;
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
                if(gene_n > 0) {
                    ok = 0;
                    while(gene_i < gene_n) {
                        if(strcmp(chr, genes[gene_i].chr) == 0) {
                            if(pos <= genes[gene_i].end && pos >= genes[gene_i].start) {
                                ok = 1;
                                break;
                            } else if(pos < genes[gene_i].start)
                                break;
                        } else if(strcmp(chr, genes[gene_i].chr) < 0)
                            break;
                        gene_i++;
                    }
                    if(ok == 0)
                        break;
                }
            } else if(j > 9) {
                if(pop_l[j] == 0) {
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
                if(hap[0] == '.') {
                    if(pop_l[j] == 1)
                        mis1_i++;
                    else if(pop_l[j] == 2)
                        mis2_i++;
                    else {
                        fprintf(stderr, "ERROR: Unknown population labels!\n\n");
                        exit(EXIT_FAILURE);
                    }
                    temp = strtok_r(NULL, "\t", &end);
                    j++;
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
                        fprintf(stderr, "ERROR: Allowed ploidy-levels are 2, 4, 6, and 8!\n\n");
                        exit(EXIT_FAILURE);
                }
                if(pop_l[j] == 1) {
                    n1 += ploidy;
                    pop1_i++;
                } else if(pop_l[j] == 2) {
                    n2 += ploidy;
                    pop2_i++;
                } else {
                    fprintf(stderr, "ERROR: Unknown population labels!\n\n");
                    exit(EXIT_FAILURE);
                }
                for(i = 0; i <= stop; i += 2) {
                    if(hap[i] == '0' || hap[i] == '1') {
                        if(pop_l[j] == 1)
                            p1 += hap[i] - '0';
                        else if(pop_l[j] == 2)
                            p2 += hap[i] - '0';
                        else {
                            fprintf(stderr, "ERROR: Unknown population labels!\n\n");
                            exit(EXIT_FAILURE);
                        }
                    } else {
                        fprintf(stderr, "ERROR: Unknown alleles found at site %s:%i! Only 0 and 1 are allowed.\n\n", chr, pos);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            temp = strtok_r(NULL, "\t", &end);
            j++;
        }
        if(temp == NULL) {
            if(pop1_i == 0 || pop2_i == 0)
                continue;
            if(pop1_i / (pop1_i + mis1_i) < mis || pop2_i / (pop2_i + mis2_i) < mis)
                continue;
            p1 /= n1;
            p2 /= n2;
            if(p1 < maf || p1 > 1 - maf || p2 < maf || p2 > 1 - maf)
                continue;
            if(stat == 0 && p1 == 0 && p2 == 0)
                continue;
            hw = (p1 - p2) * (p1 - p2) - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1);
            hb = p1 * (1 - p2) + p2 * (1 - p1);
            tot_hw += hw;
            tot_hb += hb;
            tot_n++;
            if(gene_n == 0 && out == 0) {
                if(stat == 1)
                    printf("%s\t%i\t%f\n", chr, pos, hb);
                else if(isnan(hw / hb) == 0)
                    printf("%s\t%i\t%f\n", chr, pos, hw / hb);
            } else if(out == 0) {
                for(i = k; i < gene_n; i++) {
                    if(strcmp(chr, genes[i].chr) == 0) {
                        if(pos <= genes[i].end && pos >= genes[i].start) {
                            genes[i].hw += hw;
                            genes[i].hb += hb;
                            genes[i].n++;
                        } else if(pos < genes[i].start) {
                            k = i;
                            for(j = 1; j <= i; j++) {
                                if(pos <= genes[i - j].end && pos >= genes[i - j].start)
                                    k = i - j;
                                else if(pos > genes[i - j].end) {
                                    if(i - j - 1 >= 0) {
                                        if(pos > genes[i - j - 1].end)
                                            break;
                                    } else
                                        break;
                                }
                            }
                            break;
                        }
                    } else if(strcmp(chr, genes[i].chr) < 0) {
                        k = i;
                        for(j = 1; j <= i; j++) {
                            if(strcmp(chr, genes[i - j].chr) == 0) {
                                if(pos <= genes[i - j].end && pos >= genes[i - j].start)
                                    k = i - j;
                                else if(pos > genes[i - j].end) {
                                    if(i - j - 1 >= 0) {
                                        if(pos > genes[i - j - 1].end)
                                            break;
                                    } else
                                        break;
                                }
                            } else if(strcmp(chr, genes[i - j].chr) > 0)
                                break;
                        }
                        break;
                    }
                }
            }
        }
    }
    if(gene_n > 0 && out == 0) {
        for(i = 0; i < gene_n; i++) {
            printf("%s\t", genes[i].id);
            if(stat == 1)
                printf("%f\t%0.f\n", genes[i].hb / genes[i].n, genes[i].n);
            else
                printf("%f\t%.0f\n", genes[i].hw / genes[i].hb, genes[i].n);
        }
    } else if(out == 1) {
        if(stat == 1)
            printf("%f\n", tot_hb / tot_n);
        else
            printf("%f\n", tot_hw / tot_hb);
    }

    if(isatty(1))
        fprintf(stderr, "\n");
    if(stat == 1)
        fprintf(stderr, "Average Dxy = %f\nTotal sites = %.0f\n", tot_hb / tot_n, tot_n);
    else
        fprintf(stderr, "Average weighted Fst = %f\nTotal sites = %.0f\n\n", tot_hw / tot_hb, tot_n);

    for(i = 0; i < pop1_n; i++)
        free(pop1[i]);
    free(pop1);
    for(i = 0; i < pop2_n; i++)
        free(pop2[i]);
    free(pop2);
    free(pop_l);
    if(site_n > 0)
        free(sites);
    if(gene_n > 0)
        free(genes);
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
    fprintf(stderr, "\nProgram for estimating pairwise Fst and Dxy from mixed ploidy VCF files.\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "-vcf [file] VCF file containing biallelic sites. Allowed ploidies are 2, 4, 6, and 8.\n");
    fprintf(stderr, "-pop1 [file] File listing individuals from population 1.\n");
    fprintf(stderr, "-pop2 [file] File listing individuals from population 2.\n");
    fprintf(stderr, "-sites [file] Tab delimited file listing sites to use (format: chr, pos). Optional.\n");
    fprintf(stderr, "-genes [file] Tab delimited file listing genes to use (format: chr, start, end, id). Output will be Fst/Dxy calculated for each gene. Optional.\n");
    fprintf(stderr, "-mis [double] Excludes sites based of the proportion of missing data (0 = all missing allowed, 1 = no missing data allowed). Default > 0.\n");
    fprintf(stderr, "-maf [double] Minimum minor allele frequency allowed. Default 0.\n");
    fprintf(stderr, "-stat [string] Whether to calculate 'fst' or 'dxy'. Default 'fst'. Note that dxy requires invariant sites to be included in the VCF file.\n");
    fprintf(stderr, "-out [int] Whether to print full output (0) or genome-wide estimate only (1). Default 0.\n\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "./poly_fst -vcf in.vcf -pop1 pop1.txt -pop2 pop2.txt -sites 4fold.sites -genes genes.txt -mis 0.8 -stat dxy > out_gene.dxy\n\n");
}