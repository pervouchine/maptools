//	Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
//	This file is a part of the IRBIS package.
//	IRBIS package is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//	
//	IRBIS package is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/dir.h>
#include "genutils.h"

#define INPUT_BED 0
#define INPUT_TSV 1
#define INPUT_USV 2

#define OUTPUT_MAF 0
#define OUTPUT_TSV 1

struct record {
    int beg;
    int end;
    int str;
    char name[100];
    struct record *next;
};

int main(int argc, char* argv[]) {
    char idx_file_name[MAXBUFFLENGTH];
    char dbx_file_name[MAXBUFFLENGTH];

    char spacer[MAXBUFFLENGTH]=".....";
    int spacer_length = 5;

    FILE *idx_file;
    FILE *dbx_file;


    long offset;
    long seqlen;

    long length_limit = MAXLONGBUFFLENGTH;

    long half_length_limit;
    char name[MAXBUFFLENGTH];

    char chr_name[MAXBUFFLENGTH];
    char buff[MAXBUFFLENGTH];
    char longbuff[MAXLONGBUFFLENGTH + 100];

    int out_type = OUTPUT_MAF;
    int inp_type = INPUT_BED;

    int i,j,q,k,n,s;
    char c;
    long a, b;

    long** beg;
    long** end;
    int** str;
    char*** ids;
    long pos,p,l;
    long margin1=0;
    long margin2=0;
    int warnings=0;
    int cut=0;

    int record_count[MAXCHR];
    int record_idx[MAXCHR];

    if(argc==1) {
	fprintf(stderr, "This utility does sequence retrieval from the 4-bit repository (see transf) given the input file of intervals\n");
	fprintf(stderr, "Last update by Dmitri Pervouchine (dp@crg.eu) on Apr 2, 2014\n");
	fprintf(stderr, "Usage: %s -dbx <file> -idx <file> [-limit <length>] [-margins <margin1> <margin2>] [-quiet]\n",argv[0]);
        fprintf(stderr," -dbx <database_file>\n -idx <database_index_file>\n -limit <length_limit> [default=%i]\n", length_limit);
        fprintf(stderr," -margins <margin1> <margin2> [default=%i %i]\n -spacer <spacer> [default=%s]\n -quiet suppress verbose output [default=NO]\n", margin1, margin2, spacer);
	fprintf(stderr," -inp_type 0=BED, 1=TSV [default=%i]\n -out_type 0=MAF, 1=TSV [default=%i]\n", inp_type, out_type);
	exit(1);
    }

    timestamp_set();
    for(i=1;i<argc;i++) {

        if(strcmp(argv[i],"-dbx")==0) {
            sscanf(argv[++i], "%s", &dbx_file_name[0]);
        }

        if(strcmp(argv[i],"-idx")==0) {
            sscanf(argv[++i], "%s", &idx_file_name[0]);
        }

        if(strcmp(argv[i],"-limit")==0) {
            sscanf(argv[++i], "%li", &length_limit);
        }

        if(strcmp(argv[i],"-margin")==0) {
            sscanf(argv[++i], "%li", &margin1);
	    margin2 = margin1;
        }

        if(strcmp(argv[i],"-margins")==0) {
            sscanf(argv[++i], "%li", &margin1);
            sscanf(argv[++i], "%li", &margin2);
        }

	if(strcmp(argv[i],"-spacer")==0) {
	    sscanf(argv[++i], "%s", &spacer[0]);
	}

        if(strcmp(argv[i],"-inp_type")==0) {
            sscanf(argv[++i], "%i", &inp_type);
        }

	if(strcmp(argv[i],"-out_type")==0) {
	    sscanf(argv[++i], "%i", &out_type);
	}

        if(strcmp(argv[i],"-quiet")==0) {
            verbose = 0;
        }

    }

    half_length_limit = length_limit/2;

    if(spacer[0]=='0') spacer[0]=0;
    spacer_length = strlen(spacer);

    struct record** curr[MAXCHR];
    struct record*  data[MAXCHR];
    struct record* ptr;

    for(i=0;i<MAXCHR;i++) {
	data[i] = NULL;
        curr[i] = &(data[i]);
    }

    while(fgets(buff,MAXBUFFLENGTH,stdin)) {
        if(strlen(buff)<2) break;
        if(inp_type==INPUT_USV) {
            for(k=0;k<strlen(buff);k++) {if(buff[k]=='_') buff[k]=' ';}
        }
      	sscanf(buff,"%s" , &chr_name[0]);
      	i = assign_code(chr_name);

	ptr = (struct record*)malloc(sizeof(struct record)+1);
        switch(inp_type) {
            case INPUT_BED :
                sscanf(buff,"%*s %li %li %s %*s %c" ,&ptr->beg, &ptr->end, &ptr->name[0], &c);
                break;
            case INPUT_USV :
            case INPUT_TSV :
                sscanf(buff,"%*s %li %li %c" ,&ptr->beg, &ptr->end, &c);
                strcpy(&ptr->name[0], (char*)".");
                break;
        }
	ptr->str = strand_c2i(c);
	if(ptr->beg > ptr->end) {
	    pos = ptr->beg;
	    ptr->beg = ptr->end;
	    ptr->end = pos;
	}
	ptr->next = NULL;

	*(curr[i]) = ptr;
	curr[i] = &(ptr->next);
    }

    if(verbose) fprintf(stderr,"[<%s,%s",idx_file_name,dbx_file_name);
    idx_file = fopen(idx_file_name,"r");
    dbx_file = fopen(dbx_file_name,"r");
    if(idx_file == NULL || dbx_file == NULL) {
        fprintf(stderr,"[ERROR: cannot access %s or %s, exiting]\n", idx_file_name, dbx_file_name);
        exit(1);
    }

    offset = 0;
    while(fgets(buff,MAXBUFFLENGTH,idx_file)) {
        if(strlen(buff)<2) break;
        sscanf(buff,"%s" , &name[0]);
	while(fgets(buff,MAXBUFFLENGTH,idx_file)) {
            if(strlen(buff)<2) break;
	    sscanf(buff,"%s %li" , &chr_name[0], &seqlen);
	    i = get_chr_code(chr_name);
	    while(data[i]!=NULL) {
	    	for(s=-1; s<=1; s+=2) {
		    if(s == data[i]->str || data[i]->str == 0) {

			a = data[i]->beg - margin1;
			b = data[i]->end + margin2;

			l = b - a;

			if(l <= 0 || b > seqlen) {
		    	    warnings++;
		    	    continue;
			}

			p = a - 1;
			if(l<length_limit) {
		    	    fget_segment(longbuff, dbx_file, offset, p, l);
			}
			else {
		    	    fget_segment(longbuff, dbx_file, offset, p, half_length_limit);
		    	    strcat(longbuff + half_length_limit, spacer);
                    	    fget_segment(longbuff + half_length_limit + spacer_length, dbx_file, offset, p + l - half_length_limit, half_length_limit);
		    	    cut++;
			}

			if(s<0) {
		    	    rev1(longbuff);
			}

			switch(out_type) {
		    	    case OUTPUT_TSV : printf("%s_%i_%i_%c\t%s\n",chr_name, data[i]->beg, data[i]->end, strand_i2c(s), longbuff);
					      break;
		    	    case OUTPUT_MAF : printf("%s\t%s\t%li\t%li\t%c\t%li\t%s\n", data[i]->name, chr_name, (s>0 ? data[i]->beg : seqlen - data[i]->end), l, 
								strand_i2c(s), seqlen, longbuff);
					      break;
			}
		    }
		}
		data[i] = data[i]->next;
	    }
            offset+= (seqlen % 8 == 0) ? seqlen/8 : (seqlen/8 + 1);
	}
    }
    fclose(idx_file);
    fclose(dbx_file);
    if(verbose) fprintf(stderr,"]\n");
    if(verbose && warnings) fprintf(stderr,"[WARNING: %i segments out of range]\n", warnings);
    if(verbose && warnings) fprintf(stderr,"[WARNING: %i segments shortened due to length limit]\n", cut);
    timestamp_report();
    exit(0);
}
